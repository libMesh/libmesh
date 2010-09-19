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

namespace libMesh
{

TransientRBParamSubdomainNode::TransientRBParamSubdomainNode(TransientRBParamSubdomainTree& tree_in,
        const std::vector<Real>& anchor_in)
        : Parent(tree_in, anchor_in)
{ }

void TransientRBParamSubdomainNode::add_child(const std::vector<Real>& new_anchor, Child c)
{
  // cast the tree reference to a TransientRBParamSubdomainTree
  TransientRBParamSubdomainTree& trans_tree =
    libmesh_cast_ref<TransientRBParamSubdomainTree&>(_tree);

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
      left_child = new TransientRBParamSubdomainNode(trans_tree, new_anchor);
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
      right_child = new TransientRBParamSubdomainNode(trans_tree, new_anchor);
    }
  }
}

void TransientRBParamSubdomainNode::hp_greedy(Real h_tol, Real p_tol, unsigned int Nbar)
{
    _rb_system.clear_basis_function_dependent_data();
    _rb_system.load_training_set( this->training_set );
    _rb_system.set_current_parameters( this->anchor );
    _rb_system.set_training_tolerance(h_tol);

    Real greedy_bound;

    // These casts have to succeed
    TransientRBParamSubdomainTree& trans_tree =
      libmesh_cast_ref<TransientRBParamSubdomainTree&>(_tree);

    TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

    // Set the maximum number of truth solve to Nbar in the time-dependent case
    trans_rb.set_max_truth_solves(Nbar);

    if (!trans_tree.use_delta_N_in_h_stage)
    {
        trans_rb.set_POD_tol(h_tol/trans_tree.conserv_factor);
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
    if (RB_error > h_tol/trans_tree.conserv_factor)
    {
        std::cout << "Error: The h-tolerance was not satisfied at the "
        << "anchor point hence h-type refinement may not converge."
        << std::endl;
        libmesh_error();
    }




    if ( greedy_bound > h_tol) // recursive call to hp_greedy
    {
        std::cout << "h tolerance not satisfied, splitting subdomain..." << std::endl;
        split_this_subdomain(true);

        left_child->hp_greedy(h_tol,p_tol,Nbar);
        right_child->hp_greedy(h_tol,p_tol,Nbar);
    }
    else // terminate branch, populate the model with standard p-type,write out subelement data
    {
        std::cout << "h tolerance satisfied, performing p-refinement..." << std::endl;
        greedy_bound = perform_p_stage(greedy_bound, p_tol);
        if (greedy_bound > p_tol)
        {
            std::cout << "p tolerance not satisfied, splitting subdomain..." << std::endl;
            split_this_subdomain(false);

            left_child->hp_greedy(h_tol,p_tol,Nbar);
            right_child->hp_greedy(h_tol,p_tol,Nbar);
        }
        else
        {
            std::cout << "p tolerance satisfied, subdomain "
                      << _tree.leaf_node_index << " is a leaf node..."
                      << std::endl;

            // Finally, write out the data for this subdomain
            write_subdomain_data_to_files();
            _tree.leaf_node_index++;
        }
    }
}

Real TransientRBParamSubdomainNode::perform_p_stage(Real greedy_bound, Real p_tol)
{
    // Continue the greedy process on this subdomain, i.e.
    // we do not discard the basis functions generated for
    // this subdomain in the h-refinement phase
    _rb_system.set_training_tolerance(p_tol);

    TransientRBSystem& trans_rb = libmesh_cast_ref<TransientRBSystem&>(_rb_system);

    // Ignore max_truth_solves and POD-tol in the p-stage
    trans_rb.set_POD_tol(-1.);
    trans_rb.set_max_truth_solves(-1);

    // Clear the reduced basis and reinitialize the greedy to the anchor point
    trans_rb.clear_basis_function_dependent_data();
    trans_rb.set_current_parameters( this->anchor );

    // Checking if p-tol is already satisfied or Nmax has been reached
    // if not do another (standard) greedy
    if ( (greedy_bound > p_tol) || (_rb_system.get_n_basis_functions() < _rb_system.get_Nmax()) )
    {
        greedy_bound = _rb_system.train_reduced_basis();
    }

    return greedy_bound;
}

} // namespace libMesh
