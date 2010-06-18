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

#ifndef __qn_transient_rb_param_subdomain_node_h__
#define __qn_transient_rb_param_subdomain_node_h__

// Configuration data
#include "libmesh_config.h"

// This class depends on QNTransientSCMSystem, which is only
// available if SLEPc is defined.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "rb_param_subdomain_node.h"

// Forward declaration
class QNTransientSCMSystem;
class QNTransientRBParamSubdomainTree;

/**
 * This class is part of the rbOOmit framework.
 *
 * This class extends RBParamSubdomainNode to add some extra
 * functionality for quadratically-nonlinear time-dependent
 * problems.
 *
 * @author David J. Knezevic, 2009
 */

// ------------------------------------------------------------
// RBParamSubdomainNode class definition

class QNTransientRBParamSubdomainNode : public RBParamSubdomainNode
{
public:

  /**
   * Constructor. Initializes required
   * data structures.
   */
  QNTransientRBParamSubdomainNode (QNTransientRBParamSubdomainTree& tree, const std::vector<Real>& anchor);

  /**
   * Destructor.
   */
  virtual ~QNTransientRBParamSubdomainNode () {}

  /**
   * The type of the parent.
   */
  typedef RBParamSubdomainNode Parent;

  /**
   * Add a new (LEFT or RIGHT) child RBParamSubdomainNode. Overloaded
   * in order to add QNTransientRBParamSubdomainNode children.
   */
  virtual void add_child(const std::vector<Real>& child_anchor, Child c);

  /**
   * Overload in order to perform a per-subdomain SCM
   * (which is not needed during the hp refinement) and then
   * write out the SCM data as well.
   */
  virtual void write_subdomain_data_to_files();

  /**
   * Reference to the QNTransientSCMSystem that is used
   * to perform the SCM on each subdomain.
   */
  QNTransientSCMSystem& _scm_system;

};

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
