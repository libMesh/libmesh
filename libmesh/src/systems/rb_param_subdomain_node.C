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
#include "transient_rb_system.h"
#include "parallel.h"

#include <ctime>

RBParamSubdomainNode::RBParamSubdomainNode(RBParamSubdomainTree& tree_in, const std::vector<Real>& anchor_in)
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
  for(unsigned int i=0; i<get_n_params(); i++)
    training_set[i].resize(0);
}

RBParamSubdomainNode::~RBParamSubdomainNode()
{
  if(left_child)
  {
    delete left_child;
    left_child = NULL;
  }
  if(right_child)
  {
    delete right_child;
    right_child = NULL;
  }
}

bool RBParamSubdomainNode::operator== (const RBParamSubdomainNode& node) const
{
  // We consider two RBParamSubdomainNodes to be equal if they have the
  // same anchor and training set

  if(this->n_local_training_parameters() == node.n_local_training_parameters())
  {
    bool equal = true;
    for(unsigned int i=0; i<get_n_params(); i++)
      for(unsigned int j=0; j<this->n_local_training_parameters(); j++)
        equal = equal && (this->training_set[i][j] == node.training_set[i][j]);

    return equal && (this->anchor == node.anchor);
  }

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

void RBParamSubdomainNode::hp_greedy(Real h_tol, Real p_tol, unsigned int Nbar, Real conserv_factor, bool use_delta_N_in_h_stage)
{
  _rb_system.clear_basis_function_dependent_data();
  _rb_system.load_training_set( this->training_set );
  _rb_system.set_current_parameters( this->anchor );
  _rb_system.set_training_tolerance(h_tol);

  // Check if conservative factor is being used (for the transient case)
  Real greedy_bound;
  if (conserv_factor > 0.)
    {
      // This cast has to succeed since we checked that _rb_system is a TransientRBSystem
      // in the RBParamSubdomainTree constructor
      TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

      // Set the maximum number of truth solve to Nbar in the time-dependent case
      trans_rb.set_max_truth_solves(Nbar);

      if(!use_delta_N_in_h_stage)
      {
        trans_rb.set_POD_tol(h_tol/conserv_factor);
      }
      else
      {
        trans_rb.set_POD_tol(-1.);
      }

      // delta_N might be changed in basis training in the transient case
      // (i.e., if we're using POD_tol, or if we hit Nmax)
      // hence we save and reload delta_N
      unsigned int saved_delta_N = trans_rb.get_delta_N();
      greedy_bound = trans_rb.train_reduced_basis();
      // Reload delta_N
      trans_rb.set_delta_N(saved_delta_N);

      trans_rb.set_current_parameters(this->anchor);
      Real RB_error = trans_rb.RB_solve(trans_rb.get_n_basis_functions());
      if(RB_error > h_tol/conserv_factor)
      {
        std::cout << "Error: The h-tolerance was not satisfied at the "
                  << "anchor point hence h-type refinement may not converge."
                  << std::endl;
        libmesh_error();
      }
    }
  else
    {
      const unsigned int initial_Nmax = _rb_system.get_Nmax();
      // Set Nmax to Nbar in the time-dependent case
      _rb_system.set_Nmax(Nbar);

      greedy_bound = _rb_system.train_reduced_basis();

      // Restore _rb_system's Nmax to initial_Nmax after h-stage greedy
      _rb_system.set_Nmax(initial_Nmax);
    }



  if( greedy_bound > h_tol) // recursive call to hp_greedy
  {
    std::cout << "h tolerance not satisfied, splitting subdomain..." << std::endl;
    split_this_subdomain(true);

    left_child->hp_greedy(h_tol,p_tol,Nbar,conserv_factor,use_delta_N_in_h_stage);
    right_child->hp_greedy(h_tol,p_tol,Nbar,conserv_factor,use_delta_N_in_h_stage);
  }
  else // terminate branch, populate the model with standard p-type,write out subelement data
  {
    std::cout << "h tolerance satisfied, performing p-refinement..." << std::endl;
    greedy_bound = perform_p_stage(greedy_bound, p_tol, conserv_factor);
    if (greedy_bound > p_tol)
    {
      std::cout << "p tolerance not satisfied, splitting subdomain..." << std::endl;
      split_this_subdomain(false);

      left_child->hp_greedy(h_tol,p_tol,Nbar,conserv_factor,use_delta_N_in_h_stage);
      right_child->hp_greedy(h_tol,p_tol,Nbar,conserv_factor,use_delta_N_in_h_stage);
    }
    else
    {
      std::cout << "p tolerance satisfied, subdomain " << _tree.leaf_node_index << " is a leaf node..." << std::endl;
      
      // Finally, write out the data for this subdomain
      write_subdomain_data_to_files();
      _tree.leaf_node_index++;
    }
  }
}

void RBParamSubdomainNode::split_this_subdomain(bool h_stage_split)
{
  this->refine_training_set(_tree.refined_training_set_size);

  this->add_child( _rb_system.get_greedy_parameter(0), RBParamSubdomainNode::LEFT);
  this->add_child( _rb_system.get_greedy_parameter(1), RBParamSubdomainNode::RIGHT);

  bool anchors_are_equal = true;
  for (unsigned int i = 0; i < left_child->anchor.size(); i++)
    anchors_are_equal = ( anchors_are_equal && (left_child->anchor[i] == right_child->anchor[i]) );

  if (h_stage_split)
  {
    if (anchors_are_equal)
    {
      std::cout << "Error: Anchor points for children are equal!"
                << std::endl;
        libmesh_error();
      }
    }
    else
    {
      if (anchors_are_equal)
      {
        for(unsigned int i=2; i< _rb_system.greedy_param_list.size() ; i++)
        {
          bool parameters_are_equal = true;
          for (unsigned int  j = 0; j < left_child->anchor.size(); j++)
          {
            parameters_are_equal = ( parameters_are_equal && (left_child->anchor[j] == _rb_system.get_greedy_parameter(i)[j]));
          }
          if (!parameters_are_equal)
          {
            right_child->anchor = _rb_system.get_greedy_parameter(i);
            anchors_are_equal = false;
            break;
          }
        }
        if (anchors_are_equal)
        {
          std::cout << "Error: Unable to find distinct anchors in additional splitting step." << std::endl;
          libmesh_error();
        }
      }
    }

  this->initialize_child_training_sets();

  if( (left_child->n_global_training_parameters()  == 0) ||
      (right_child->n_global_training_parameters() == 0) )
    {
      std::cout << "Error: Child training set is empty after initialization. At least the training set should contain the anchor point!" << std::endl;
      libmesh_error();
    }

  if ( (left_child->n_global_training_parameters() < _tree.min_training_set_size) ||
       (right_child->n_global_training_parameters() < _tree.min_training_set_size) )
  {
    // Determine a splitting that gives min_training_set_size global training parameters
    unsigned int n_local_training_samples;
    unsigned int quotient  = _tree.min_training_set_size/libMesh::n_processors();
    unsigned int remainder = _tree.min_training_set_size%libMesh::n_processors();
    if(libMesh::processor_id() < remainder)
      n_local_training_samples = (quotient + 1);
    else
      n_local_training_samples = quotient;

    if (left_child->n_global_training_parameters() < _tree.min_training_set_size)
    {
      left_child->refine_training_set(n_local_training_samples);
    }

    if (right_child->n_global_training_parameters() < _tree.min_training_set_size)
    {
      right_child->refine_training_set(n_local_training_samples);
    }
  }
}

Real RBParamSubdomainNode::perform_p_stage(Real greedy_bound, Real p_tol, Real conserv_factor)
{
  // Continue the greedy process on this subdomain, i.e.
  // we do not discard the basis functions generated for
  // this subdomain in the h-refinement phase
  _rb_system.set_training_tolerance(p_tol);

  // Clear POD_tol for p-type greedy
  if (conserv_factor > 0.)
	{
	  TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

    // Ignore max_truth_solves and POD-tol in the p-stage
    trans_rb.set_POD_tol(-1.);
    trans_rb.set_max_truth_solves(-1);

    // Clear the reduced basis and reinitialize the greedy to the anchor point
    trans_rb.clear_basis_function_dependent_data();
    trans_rb.set_current_parameters( this->anchor );
	}

  // Checking if p-tol is already satisfied or Nmax has been reached
  // if not do another (standard) greedy
  if ( (greedy_bound > p_tol) || (_rb_system.get_n_basis_functions() < _rb_system.get_Nmax()) )
	{
	  greedy_bound = _rb_system.train_reduced_basis();
	}

  return greedy_bound;
}

void RBParamSubdomainNode::write_subdomain_data_to_files()
{
  std::cout << "Writing out RB data for leaf subdomain "
            << _tree.leaf_node_index << std::endl;

  std::stringstream dir_name_stream;
  dir_name_stream << "offline_data_hp" << _tree.leaf_node_index;
  const std::string& directory_name = dir_name_stream.str();
  _rb_system.write_offline_data_to_files(directory_name);
}

void RBParamSubdomainNode::add_child(const std::vector<Real>& new_anchor, Child c)
{
  if(c == LEFT)
  {
    if(left_child != NULL)
    {
      std::cout << "Error: Child already exists!"
                << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      left_child = new RBParamSubdomainNode(_tree, new_anchor);
    }
  }

  if(c == RIGHT)
  {
    if(right_child != NULL)
    {
      std::cout << "Error: Child already exists!"
                << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      right_child = new RBParamSubdomainNode(_tree, new_anchor);
    }
  }
}

void RBParamSubdomainNode::refine_training_set(const unsigned int new_local_training_set_size)
{

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

      for(unsigned int i=0; i<anchor.size(); i++)
	{
	  ran = ((double)std::rand())/RAND_MAX; // in range [0,1]
	  random_point[i] = corners[0][i] + (corners[1][i]-corners[0][i])*ran;
	}

      if( (_tree.determine_subdomain(random_point)) == this )
      {
	for(unsigned int i=0; i<anchor.size(); i++)
	  training_set[i].push_back(random_point[i]);
      }
    }

//   // Print out the training sample
//   for (unsigned int j=0; j <n_training_parameters(); j++)
//     std :: cout << training_set[0][j] << " " << training_set[1][j] << std :: endl;
}

void RBParamSubdomainNode::initialize_child_training_sets()
{
  if( (left_child == NULL) || (right_child == NULL) )
  {
    std::cout << "ERROR: Children cannot be NULL in initialize_child_training_sets()."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(!left_child->training_set_initialized &&
                 !right_child->training_set_initialized);

  left_child->training_set.resize(anchor.size());
  right_child->training_set.resize(anchor.size());

  for(unsigned int i=0; i<anchor.size(); i++)
  {
    left_child->training_set[i].clear();
    right_child->training_set[i].clear();
  }

  for(unsigned int i=0; i<n_local_training_parameters(); i++)
  {
    std::vector<Real> next_param = get_local_training_parameter(i);
    Real left_dist  = left_child->dist_from_anchor(next_param);
    Real right_dist = right_child->dist_from_anchor(next_param);
    if( left_dist <= right_dist )
    {
      for(unsigned int i=0; i<next_param.size(); i++)
        left_child->training_set[i].push_back(next_param[i]);
    }
    else
    {
      for(unsigned int i=0; i<next_param.size(); i++)
        right_child->training_set[i].push_back(next_param[i]);
    }
  }

  left_child->training_set_initialized  = true;
  right_child->training_set_initialized = true;
}



void RBParamSubdomainNode::copy_training_set_from_system()
{
  training_set.resize(_rb_system.get_n_params());

  for(unsigned int i=0; i<training_set.size(); i++)
    _rb_system.training_parameters[i]->localize(training_set[i]);

  training_set_initialized = true;
}

Real RBParamSubdomainNode::dist_from_anchor(const std::vector<Real>& new_param) const
{
  if(new_param.size() != anchor.size())
  {
    std::cout << "Error: Input vector is of incorrect size in dist_from_anchor."
              << std::endl;
    libmesh_error();
  }

  Real sum = 0.;
  for(unsigned int i=0; i<anchor.size(); i++)
  {
    // Map anchor[i] and new_param[i] to (0,1) before computing distance
    Real jacobian = _rb_system.get_parameter_max(i) - _rb_system.get_parameter_min(i);
    Real anchor_mapped = (anchor[i] - _rb_system.get_parameter_min(i))/jacobian;
    Real new_param_mapped = (new_param[i] - _rb_system.get_parameter_min(i))/jacobian;
    sum += std::pow((anchor_mapped - new_param_mapped), 2.);
    //sum += std::pow((anchor[i] - new_param[i]), 2.);
  }

  return std::sqrt(sum);
}

std::vector<Real> RBParamSubdomainNode::get_local_training_parameter(unsigned int i)
{
  if(i >= n_local_training_parameters())
  {
    std::cout << "Error: Argument is too large in get_training_parameter."
              << std::endl;
  }

  std::vector<Real> param(training_set.size());
  for(unsigned int j=0; j<param.size(); j++)
    param[j] = training_set[j][i];

  return param;
}

std::vector< std::vector<Real> > RBParamSubdomainNode::get_training_bbox()
{
  std::vector <Real> max_values(training_set.size());
  std::vector <Real> min_values(training_set.size());
  std::vector <Real> range(training_set.size());

  std::vector <std::vector<Real> > corner_values(2);
  corner_values[0].resize(training_set.size());
  corner_values[1].resize(training_set.size());

  for(unsigned int i=0; i<training_set.size(); i++)
    {
      min_values[i] = training_set[i][0];
      max_values[i] = training_set[i][0];
    }

  for(unsigned int i=0; i<training_set.size(); i++)
    {
      for(unsigned int j=1; j<training_set[0].size(); j++)
	{
	  if (training_set[i][j] < min_values[i])
	    min_values[i] = training_set[i][j];

	  if (training_set[i][j] > max_values[i])
	    max_values[i] = training_set[i][j];
	}
    }

  for(unsigned int i=0; i<training_set.size(); i++)
    {
      range[i] = max_values[i] - min_values[i];
      corner_values[0][i] = min_values[i] - _tree.bbox_margin*range[i];
      corner_values[1][i] = max_values[i] + _tree.bbox_margin*range[i];

      if (corner_values[0][i] < _rb_system.get_parameter_min(i))
	corner_values[0][i] = _rb_system.get_parameter_min(i);
      if (corner_values[1][i] > _rb_system.get_parameter_max(i))
	corner_values[1][i] = _rb_system.get_parameter_max(i);
    }

  return corner_values;
}


