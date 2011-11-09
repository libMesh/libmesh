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

#ifndef __transient_rb_param_subdomain_tree_h__
#define __transient_rb_param_subdomain_tree_h__

#include "rb_param_subdomain_tree.h"
#include "transient_rb_construction.h"

namespace libMesh
{

// Forward declaration
class RBParamSubdomainNode;

/**
 * This class is part of the rbOOmit framework.
 *
 * TransientRBParamSubdomainTree extends the adaptive partition of parameter
 * space functionality to the case of time-dependent problems.
 *
 * @author David J. Knezevic and Jens L. Eftang, 2009
 */

// ------------------------------------------------------------
// TransientRBParamSubdomainTree class definition

class TransientRBParamSubdomainTree : public RBParamSubdomainTree
{
public:

  /**
   * Constructor. Initializes required
   * data structures.
   */
  TransientRBParamSubdomainTree (TransientRBConstruction& rb_construction_in,
                                 const std::string& parameters_filename);
  
  /**
   * Build the root_node which, in this case, is a TransientRBParamNode.
   */
  virtual void build_root_node();
  
  /**
   * A convenient typedef for the parent class.
   */
  typedef RBParamSubdomainTree Parent;

  /**
   * Conservativity factor to determine POD_tol from h_tol
   */
  Real conserv_factor;

  /**
   * Flag that indicates whether we use delta_N or the POD-tolerance
   * to govern how many basis functions we add in a time-dependent
   * hp-greedy.
   */
  bool use_delta_N_in_h_stage;

};

} // namespace libMesh

#endif
