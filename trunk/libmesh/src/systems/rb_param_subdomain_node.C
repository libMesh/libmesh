// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "rb_param_subdomain_node.h"
#include "rb_param_subdomain_tree.h"
#include "rb_system.h"
#include "parallel.h"
#include "libmesh_logging.h"
#include "utility.h"

#include <ctime>

namespace libMesh
{

RBParamSubdomainNode::RBParamSubdomainNode(RBParamSubdomainTree& tree_in,
                                           const std::vector<Real>& anchor_in)
        : left_child(NULL),
        right_child(NULL),
        _tree(tree_in),
        _rb_system(tree_in._rb_system),
        anchor(anchor_in),
        model_number(-1),
        training_set_initialized(false)
{
    libmesh_assert( get_n_params() != 0);

    // Clear the training set to begin with
    training_set.resize(get_n_params());
    for (unsigned int i=0; i<get_n_params(); i++)
        training_set[i].clear();
}

RBParamSubdomainNode::~RBParamSubdomainNode()
{
    if (left_child)
    {
        delete left_child;
        left_child = NULL;
    }
    if (right_child)
    {
        delete right_child;
        right_child = NULL;
    }
}

bool RBParamSubdomainNode::operator== (const RBParamSubdomainNode& node) const
{
    START_LOG("operator==", "RBParamSubdomainNode");

    // We consider two RBParamSubdomainNodes to be equal if they have the
    // same model number, anchor and training set

    if (this->n_local_training_parameters() == node.n_local_training_parameters())
    {
        bool equal = true;
        for (unsigned int i=0; i<get_n_params(); i++)
            for (unsigned int j=0; j<this->n_local_training_parameters(); j++)
                equal = equal && (this->training_set[i][j] == node.training_set[i][j]);

        STOP_LOG("operator==", "RBParamSubdomainNode");
        return equal && (this->model_number == node.model_number) && (this->anchor == node.anchor);
    }

    STOP_LOG("operator==", "RBParamSubdomainNode");

    return false;
}

unsigned int RBParamSubdomainNode::n_local_training_parameters() const
{
    return training_set[0].size();
}

unsigned int RBParamSubdomainNode::n_global_training_parameters() const
{
    unsigned int global_training_set_size = this->n_local_training_parameters();
    Parallel::sum(global_training_set_size);
    return global_training_set_size;
}

void RBParamSubdomainNode::hp_greedy()
{
    _rb_system.clear_basis_function_dependent_data();

    // Load the (full or subsampled) training set
    if(_tree.n_subsampled_training_points >= n_global_training_parameters())
    {
      _rb_system.load_training_set( training_set );
    }
    else
    {
      std::vector< std::vector<Number> > subsampled_training_set = get_subsampled_training_set();
      _rb_system.load_training_set( subsampled_training_set );
    }

    _rb_system.set_current_parameters( this->anchor );
    _rb_system.set_training_tolerance(_tree.h_tol);

    // Save _rb_system's Nmax and set Nmax to N_bar
    const unsigned int initial_Nmax = _rb_system.get_Nmax();
    _rb_system.set_Nmax(_tree.N_bar);

    Real greedy_bound = _rb_system.train_reduced_basis();

    // Restore _rb_system's Nmax to initial_Nmax after h-stage greedy
    _rb_system.set_Nmax(initial_Nmax);


    if ( greedy_bound > _tree.h_tol) // recursive call to hp_greedy
    {
        libMesh::out << "h tolerance not satisfied, splitting subdomain..." << std::endl;
        split_this_subdomain(true);

        left_child->hp_greedy();
        right_child->hp_greedy();
    }
    else // terminate branch, populate the model with standard p-type, write out subelement data
    {
        libMesh::out << "h tolerance satisfied, performing p-refinement..." << std::endl;
        greedy_bound = perform_p_stage(greedy_bound);
        if (greedy_bound > _tree.p_tol)
        {
            libMesh::out << "p tolerance not satisfied, splitting subdomain..." << std::endl;
            split_this_subdomain(false);

            left_child->hp_greedy();
            right_child->hp_greedy();
        }
        else
        {
            libMesh::out << "p tolerance satisfied, subdomain " << _tree.leaf_node_index << " is a leaf node..." << std::endl;

            // Finally, write out the data for this subdomain
            this->model_number = _tree.leaf_node_index;
            write_subdomain_data_to_files();
            _tree.leaf_node_index++;
        }
    }
}

void RBParamSubdomainNode::split_this_subdomain(bool )
{
    START_LOG("split_this_subdomain()", "RBParamSubdomainNode");

    this->add_child( _rb_system.get_greedy_parameter(0), RBParamSubdomainNode::LEFT);
    this->add_child( _rb_system.get_greedy_parameter(1), RBParamSubdomainNode::RIGHT);

    // Set the distance_between_anchors member of the children
    left_child->distance_between_anchors = right_child->distance_between_anchors = 
      distance(left_child->anchor, right_child->anchor);

    // Make sure that the anchor points of the children are different
    if (left_child->distance_between_anchors == 0.)
    {
        libMesh::out << "Error: Anchor points for children are equal!"
                     << std::endl;
        libmesh_error();
    }

    this->initialize_child_training_sets();

    STOP_LOG("split_this_subdomain()", "RBParamSubdomainNode");
}

Real RBParamSubdomainNode::perform_p_stage(Real greedy_bound)
{
    START_LOG("perform_p_stage()", "RBParamSubdomainNode");

    // Continue the greedy process on this subdomain, i.e.
    // we do not discard the basis functions generated for
    // this subdomain in the h-refinement phase
    _rb_system.set_training_tolerance(_tree.p_tol);

    // Checking if p-tol is already satisfied or Nmax has been reached
    // if not do another (standard) greedy
    if ( (greedy_bound > _tree.p_tol) || (_rb_system.get_n_basis_functions() < _rb_system.get_Nmax()) )
    {
        greedy_bound = _rb_system.train_reduced_basis();
    }

    STOP_LOG("perform_p_stage()", "RBParamSubdomainNode");

    return greedy_bound;
}

void RBParamSubdomainNode::write_subdomain_data_to_files()
{
    START_LOG("write_subdomain_data_to_files()", "RBParamSubdomainNode");

    libMesh::out << "Writing out RB data for leaf subdomain "
                 << _tree.leaf_node_index << std::endl;

    std::stringstream dir_name_stream;
    dir_name_stream << "offline_data_hp" << _tree.leaf_node_index;
    const std::string& directory_name = dir_name_stream.str();
    _rb_system.write_offline_data_to_files(directory_name);

    STOP_LOG("write_subdomain_data_to_files()", "RBParamSubdomainNode");
}

void RBParamSubdomainNode::add_child(const std::vector<Real>& new_anchor, Child c)
{
    START_LOG("add_child()", "RBParamSubdomainNode");

    if (c == LEFT)
    {
        if (left_child != NULL)
        {
            libMesh::out << "Error: Child already exists!"
                         << std::endl;
            libmesh_error();
        }
        else
        {
            // Allocate a new child node
            left_child = new RBParamSubdomainNode(_tree, new_anchor);
        }
    }

    if (c == RIGHT)
    {
        if (right_child != NULL)
        {
            libMesh::out << "Error: Child already exists!"
                         << std::endl;
            libmesh_error();
        }
        else
        {
            // Allocate a new child node
            right_child = new RBParamSubdomainNode(_tree, new_anchor);
        }
    }

    STOP_LOG("add_child()", "RBParamSubdomainNode");
}


void RBParamSubdomainNode::refine_training_set(const unsigned int new_local_training_set_size)
{
    START_LOG("refine_training_set()", "RBParamSubdomainNode");

    // seed the random number generator with the system time
    // and the processor ID so that the seed is different
    // on different processors
    std::srand((unsigned)( std::time(0)*(1+libMesh::processor_id()) ));

    // Use a deterministic seed, useful for testing purposes so that
    // we obtain reproducible training sets
    //std::srand((unsigned)( (1+libMesh::processor_id()) ));

    std::vector<Real> random_point(training_set.size());
    std::vector< std::vector<Real> > corners = get_training_bbox();

    while (n_local_training_parameters() < new_local_training_set_size)
    {
        Real ran;

        for (unsigned int i=0; i<anchor.size(); i++)
        {
            ran = ((double)std::rand())/RAND_MAX; // in range [0,1]
            random_point[i] = corners[0][i] + (corners[1][i]-corners[0][i])*ran;
        }

        if ( (_tree.determine_subdomain(random_point)) == this )
        {
            for (unsigned int i=0; i<anchor.size(); i++)
                training_set[i].push_back(random_point[i]);
        }
    }

    STOP_LOG("refine_training_set()", "RBParamSubdomainNode");
}

void RBParamSubdomainNode::initialize_child_training_sets()
{
    START_LOG("initialize_child_training_sets()", "RBParamSubdomainNode");
    
    if ( (left_child == NULL) || (right_child == NULL) )
    {
        libMesh::out << "ERROR: Children cannot be NULL in initialize_child_training_sets()."
                     << std::endl;
        libmesh_error();
    }

    libmesh_assert(!left_child->training_set_initialized &&
                   !right_child->training_set_initialized);

    left_child->training_set.resize(anchor.size());
    right_child->training_set.resize(anchor.size());

    for (unsigned int i=0; i<anchor.size(); i++)
    {
        left_child->training_set[i].clear();
        right_child->training_set[i].clear();
    }

    // Split this training set among the children
    for (unsigned int i=0; i<n_local_training_parameters(); i++)
    {
        std::vector<Real> next_param = get_local_training_parameter(i);
        Real left_dist  = left_child->dist_from_anchor(next_param);
        Real right_dist = right_child->dist_from_anchor(next_param);
        if ( left_dist <= right_dist )
        {
            for (unsigned int i=0; i<next_param.size(); i++)
                left_child->training_set[i].push_back(next_param[i]);
        }
        else
        {
            for (unsigned int i=0; i<next_param.size(); i++)
                right_child->training_set[i].push_back(next_param[i]);
        }
    }

    // Make sure that each child has at least one training point
    if ( (left_child->n_global_training_parameters()  == 0) ||
         (right_child->n_global_training_parameters() == 0) )
    {
        libMesh::out << "Error: Child training set is empty after initialization." << std::endl
                     << "At least the training set should contain the anchor point!" << std::endl;
        libmesh_error();
    }
    
    // Calculate the target number of local training samples. The global target
    // is min(_tree.n_subsampled_training_points, n_global_training_parameters())
    // since we just want to ignore n_subsampled_training_points if it is larger
    // than the original training set
    unsigned int target_n_local_training_samples;
    if( _tree.n_subsampled_training_points >= n_global_training_parameters() )
    {
        target_n_local_training_samples = this->n_local_training_parameters();
    }
    else
    {
        // Calculate the minimum number of training points that this processor must
        // contain in order for this subdomain to have at least _tree.n_subsampled_training_points
        // in total
        unsigned int quotient  = _tree.n_subsampled_training_points/libMesh::n_processors();
        unsigned int remainder = _tree.n_subsampled_training_points%libMesh::n_processors();
        target_n_local_training_samples =
          (libMesh::processor_id() < remainder) ? (quotient + 1) : quotient;
    }

    // Finally, perform the actual training set refinement
    if (left_child->n_local_training_parameters() < target_n_local_training_samples)
    {
        left_child->refine_training_set(target_n_local_training_samples);
    }

    if (right_child->n_local_training_parameters() < target_n_local_training_samples)
    {
        right_child->refine_training_set(target_n_local_training_samples);
    }
    
    // possibly clear the parent node's training set once the children
    // sets have been initialized
    if(_tree.clear_training_sets_during_hp)
    {
        clear_training_set();
    }

    left_child->training_set_initialized  = true;
    right_child->training_set_initialized = true;

    STOP_LOG("initialize_child_training_sets()", "RBParamSubdomainNode");
}



void RBParamSubdomainNode::copy_training_set_from_system()
{
    START_LOG("copy_training_set_from_system()", "RBParamSubdomainNode");
    
    training_set.resize(_rb_system.get_n_params());

    for (unsigned int i=0; i<training_set.size(); i++)
    {
      training_set[i].clear();
      for (unsigned int j=_rb_system.training_parameters[i]->first_local_index();
                        j<_rb_system.training_parameters[i]->last_local_index();
                        j++)
        training_set[i].push_back( (*_rb_system.training_parameters[i])(j) );
    }

    training_set_initialized = true;

    STOP_LOG("copy_training_set_from_system()", "RBParamSubdomainNode");
}

Real RBParamSubdomainNode::distance(const std::vector<Real>& p1, const std::vector<Real>& p2) const
{
    START_LOG("distance()", "RBParamSubdomainNode");
    
    if (p1.size() != p2.size())
    {
        libMesh::out << "Error: Input vectors must have the same size in RBParamSubdomain::distance."
                     << std::endl;
        libmesh_error();
    }

    Real sum = 0.;
    for (unsigned int i=0; i<p1.size(); i++)
    {
        // Map anchor[i] and new_param[i] to (0,1) before computing distance
        Real jacobian = _rb_system.get_parameter_max(i) - _rb_system.get_parameter_min(i);
        Real p1_mapped = (p1[i] - _rb_system.get_parameter_min(i))/jacobian;
        Real p2_mapped = (p2[i] - _rb_system.get_parameter_min(i))/jacobian;
        sum += libmesh_norm(p1_mapped - p2_mapped);
    }
    
    STOP_LOG("distance()", "RBParamSubdomainNode");

    return std::sqrt(sum);
}

Real RBParamSubdomainNode::dist_from_anchor(const std::vector<Real>& new_param) const
{
    return distance(anchor, new_param);
}

std::vector<Real> RBParamSubdomainNode::get_local_training_parameter(unsigned int i)
{
    START_LOG("get_local_training_parameter()", "RBParamSubdomainNode");
    
    if (i >= n_local_training_parameters())
    {
        libMesh::out << "Error: Argument is too large in get_training_parameter."
                     << std::endl;
    }

    std::vector<Real> param(training_set.size());
    for (unsigned int j=0; j<param.size(); j++)
        param[j] = libmesh_real(training_set[j][i]);

    STOP_LOG("get_local_training_parameter()", "RBParamSubdomainNode");

    return param;
}

std::vector< std::vector<Real> > RBParamSubdomainNode::get_training_bbox()
{
    START_LOG("get_training_bbox()", "RBParamSubdomainNode");

    std::vector <std::vector<Real> > corner_values(2);
    corner_values[0].resize(training_set.size());
    corner_values[1].resize(training_set.size());

    for (unsigned int i=0; i<training_set.size(); i++)
    {
        corner_values[0][i] = anchor[i] - _tree.bbox_margin*this->distance_between_anchors;
        corner_values[1][i] = anchor[i] + _tree.bbox_margin*this->distance_between_anchors;

        // Might have to set corner_values to the parameter domain boundaries
        if (corner_values[0][i] < _rb_system.get_parameter_min(i))
            corner_values[0][i] = _rb_system.get_parameter_min(i);
        if (corner_values[1][i] > _rb_system.get_parameter_max(i))
            corner_values[1][i] = _rb_system.get_parameter_max(i);
    }


#ifdef DEBUG
    // Check whether the bbox is only one point; in this case print an error message
    bool trivial_bbox = true;
    for (unsigned int i=0; i<training_set.size(); i++)
    {
        if (corner_values[0][i] != corner_values[1][i])
            trivial_bbox = false;
        break;
    }
    libmesh_assert (!trivial_bbox);
#endif

    STOP_LOG("get_training_bbox()", "RBParamSubdomainNode");

    return corner_values;
}

std::vector< std::vector<Number> > RBParamSubdomainNode::get_subsampled_training_set()
{
  START_LOG("get_subsampled_training_set()", "RBParamSubdomainNode");

  // Return the full training set if the subsampled set size is big enough
  if(_tree.n_subsampled_training_points >= n_global_training_parameters())
  {
    STOP_LOG("get_subsampled_training_set()", "RBParamSubdomainNode");
    return training_set;
  }
  
  // Otherwise, build a new set of size n_subsampled_training_points
  unsigned int quotient  = _tree.n_subsampled_training_points/libMesh::n_processors();
  unsigned int remainder = _tree.n_subsampled_training_points%libMesh::n_processors();
  unsigned int target_n_local_training_samples =
    (libMesh::processor_id() < remainder) ? (quotient + 1) : quotient;

  // Seed the random number generator
  std::srand((unsigned)( std::time(0)*(1+libMesh::processor_id()) ));
  std::vector< std::vector<Number> > subsampled_training_set(anchor.size());
  for(unsigned int i=0; i<subsampled_training_set.size(); i++)
    subsampled_training_set[i].clear();
  
  // Create a vector of the local indices
  std::vector<unsigned int> local_index_set( n_local_training_parameters() );
  Utility::iota(local_index_set.begin(), local_index_set.end(), 0);

  while(subsampled_training_set[0].size() < target_n_local_training_samples)
  {
    unsigned int position = std::rand() % local_index_set.size();
    unsigned int random_index = local_index_set[position];
    
    // Make sure that we don't add repeated parameters to subsampled_training_set
    std::vector<Real> sample = get_local_training_parameter(random_index);
    local_index_set.erase(local_index_set.begin() + position);

    for (unsigned int j=0; j<subsampled_training_set.size(); j++)
    {
      subsampled_training_set[j].push_back( sample[j] );
    }
  }

  STOP_LOG("get_subsampled_training_set()", "RBParamSubdomainNode");

  return subsampled_training_set;
}

void RBParamSubdomainNode::clear_training_set()
{
  START_LOG("clear_training_set()", "RBParamSubdomainNode");

  // Use the "swapping idiom" to clear the capacity of each
  // vector in the training set.
  for(unsigned int i=0; i<training_set.size(); i++)
    std::vector<Number>().swap(training_set[i]);

  STOP_LOG("clear_training_set()", "RBParamSubdomainNode");
}

} // namespace libMesh


