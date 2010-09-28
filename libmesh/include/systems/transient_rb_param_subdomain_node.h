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

#ifndef __transient_rb_param_subdomain_node_h__
#define __transient_rb_param_subdomain_node_h__

#include "rb_param_subdomain_node.h"

namespace libMesh
{

// Forward declaration
class TransientRBParamSubdomainTree;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBParamSubdomainNode implements a "node" of an "hp-tree" structure
 * which is used to generate a hierarchical partition of the parameter
 * domain. Each RBParamSubdomainNode is associated with a reduced basis
 * model that is developed on that particular subdomain. In this class,
 * the name "hp" refers to a reduced basis approach in which we first
 * subdivide the parameter domain ("h" refinement) and then enrich the
 * reduced basis in each subdomain ("p" refinement).
 *
 * @author David J. Knezevic and Jens L. Eftang, 2009
 */

// ------------------------------------------------------------
// RBParamSubdomainNode class definition

class TransientRBParamSubdomainNode : public RBParamSubdomainNode
{
public:

  /**
   * Constructor. Initializes required data structures.
   */
  TransientRBParamSubdomainNode (TransientRBParamSubdomainTree& tree, const std::vector<Real>& anchor);
  
  /**
   * Convenient typedef for the Parent class.
   */
  typedef RBParamSubdomainNode Parent;
  
  /**
   * Add a new (LEFT or RIGHT) child RBParamSubdomainNode. Overloaded
   * in order to add QNTransientRBParamSubdomainNode children.
   */
  virtual void add_child(const std::vector<Real>& child_anchor, Child c);

  /**
   * Recursive function that performs the "hp" greedy algorithm.
   */
  virtual void hp_greedy();

  /**
   * Split the current subdomain into two new subdomains.
   * The argument h_stage_split indicates whether or not
   * this is called during the h-stage or the p-stage.
   * Overloaded because we need specific behavior in case
   * the POD-Greedy selected parameters more than once;
   * need to be careful not to set the same anchor point
   * in both child nodes.
   */
  virtual void split_this_subdomain(bool h_stage_split);

  /**
   * This function performs the "p" stage of the "hp"
   * greedy algorithm.
   */
  virtual Real perform_p_stage(Real greedy_bound);

};

} // namespace libMesh

#endif
