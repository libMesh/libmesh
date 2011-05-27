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

#include "transient_rb_param_subdomain_node.h"
#include "transient_rb_param_subdomain_tree.h"
#include "transient_rb_system.h"
#include "libmesh_logging.h"

namespace libMesh
{

TransientRBParamSubdomainNode::TransientRBParamSubdomainNode(TransientRBParamSubdomainTree& tree_in,
        const std::vector<Real>& anchor_in)
        : Parent(tree_in, anchor_in)
{ }

void TransientRBParamSubdomainNode::add_child(const std::vector<Real>& new_anchor, Child c)
{
  START_LOG("add_child()", "TransientRBParamSubdomainNode");

  // cast the tree reference to a TransientRBParamSubdomainTree
  TransientRBParamSubdomainTree& trans_tree =
    libmesh_cast_ref<TransientRBParamSubdomainTree&>(_tree);

  if(c == LEFT)
  {
    if(left_child != NULL)
    {
      libMesh::out << "Error: Child already exists!"
                   << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      left_child = new TransientRBParamSubdomainNode(trans_tree, new_anchor);
    }
  }

  if(c == RIGHT)
  {
    if(right_child != NULL)
    {
      libMesh::out << "Error: Child already exists!"
                   << std::endl;
      libmesh_error();
    }
    else
    {
      // Allocate a new child node
      right_child = new TransientRBParamSubdomainNode(trans_tree, new_anchor);
    }
  }

  STOP_LOG("add_child()", "TransientRBParamSubdomainNode");
}

void TransientRBParamSubdomainNode::hp_greedy(bool store_basis_functions)
{
    _rb_system.rb_eval->clear();

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

    Real greedy_bound;

    // These casts have to succeed
    TransientRBParamSubdomainTree& trans_tree =
      libmesh_cast_ref<TransientRBParamSubdomainTree&>(_tree);

    TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

    // Set the maximum number of truth solves to N_bar in the time-dependent case
    trans_rb.set_max_truth_solves(_tree.N_bar);

    if (!trans_tree.use_delta_N_in_h_stage)
    {
        trans_rb.set_POD_tol(_tree.h_tol/trans_tree.conserv_factor);
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
    Real RB_error = trans_rb.rb_eval->RB_solve(trans_rb.rb_eval->get_n_basis_functions());
    if (RB_error > _tree.h_tol/trans_tree.conserv_factor)
    {
        libMesh::out << "Error: The h-tolerance was not satisfied at the "
                     << "anchor point hence h-type refinement may not converge."
                     << std::endl;
        libmesh_error();
    }




    if ( greedy_bound > _tree.h_tol) // recursive call to hp_greedy
    {
        libMesh::out << "h tolerance not satisfied, splitting subdomain..." << std::endl;
        split_this_subdomain(true);

        left_child->hp_greedy(store_basis_functions);
        right_child->hp_greedy(store_basis_functions);
    }
    else // terminate branch, populate the model with standard p-type,write out subelement data
    {
        libMesh::out << "h tolerance satisfied, performing p-refinement..." << std::endl;
        greedy_bound = perform_p_stage(greedy_bound);
        if (greedy_bound > _tree.p_tol)
        {
            libMesh::out << "p tolerance not satisfied, splitting subdomain..." << std::endl;
            split_this_subdomain(false);

            left_child->hp_greedy(store_basis_functions);
            right_child->hp_greedy(store_basis_functions);
        }
        else
        {
            libMesh::out << "p tolerance satisfied, subdomain "
                         << _tree.leaf_node_index << " is a leaf node..."
                         << std::endl;

            // Finally, write out the data for this subdomain
            write_subdomain_data_to_files(store_basis_functions);
            _tree.leaf_node_index++;
        }
    }
}

Real TransientRBParamSubdomainNode::perform_p_stage(Real greedy_bound)
{
    START_LOG("perform_p_stage()", "TransientRBParamSubdomainNode");
    
    // Continue the greedy process on this subdomain, i.e.
    // we do not discard the basis functions generated for
    // this subdomain in the h-refinement phase
    _rb_system.set_training_tolerance(_tree.p_tol);

    TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

    // Ignore max_truth_solves and POD-tol in the p-stage
    trans_rb.set_POD_tol(-1.);
    trans_rb.set_max_truth_solves(-1);

    // Clear the reduced basis and reinitialize the greedy to the anchor point
    trans_rb.clear();
    trans_rb.set_current_parameters( this->anchor );

    // Checking if p-tol is already satisfied or Nmax has been reached
    // if not do another (standard) greedy
    if ( (greedy_bound > _tree.p_tol) ||
         (_rb_system.rb_eval->get_n_basis_functions() < _rb_system.get_Nmax()) )
    {
        greedy_bound = _rb_system.train_reduced_basis();
    }

    STOP_LOG("perform_p_stage()", "TransientRBParamSubdomainNode");

    return greedy_bound;
}

void TransientRBParamSubdomainNode::split_this_subdomain(bool h_stage_split)
{
    START_LOG("split_this_subdomain()", "TransientRBParamSubdomainNode");

    // These first few lines are the same as RBParamSubdomainNode::split_this_subdomain
    this->add_child( _rb_system.get_greedy_parameter(0), RBParamSubdomainNode::LEFT);
    this->add_child( _rb_system.get_greedy_parameter(1), RBParamSubdomainNode::RIGHT);

    // Compute distance between the children anchor points, and pass to children (JLE 2010-09-16)
    Real distance_between_children_anchors = 0.;
    for (unsigned int i = 0; i < left_child->anchor.size(); i++)
    {
        distance_between_children_anchors += std::pow((left_child->anchor[i] - right_child->anchor[i]),Real(2.));
    }
    distance_between_children_anchors = std::sqrt(distance_between_children_anchors);
    left_child->distance_between_anchors = distance_between_children_anchors;
    right_child->distance_between_anchors = distance_between_children_anchors;


    // We now need some code specific to the transient case because it
    // is possible that we have repeated selection of training points
    // in the transient case (due to the POD-Greedy) and hence we
    // may have distance_between_children_anchors == 0.
    bool anchors_are_equal = (distance_between_children_anchors == 0.);

    if (h_stage_split)
    {
        if (anchors_are_equal)
        {
            libMesh::out << "Error: Anchor points for children are equal!"
                         << std::endl;
            libmesh_error();
        }
    }
    else
    {
        if (anchors_are_equal)
        {
            for (unsigned int i=2; i< _rb_system.rb_eval->greedy_param_list.size() ; i++)
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
            
            // anchors_are_equal has been updated, check if we have found different point.
            if(anchors_are_equal) 
            {
                libMesh::out << "Error: Unable to find distinct anchors in additional splitting step." << std::endl;
                libmesh_error();
            }
        }
    }

    this->initialize_child_training_sets();

    STOP_LOG("split_this_subdomain()", "TransientRBParamSubdomainNode");
}

} // namespace libMesh
