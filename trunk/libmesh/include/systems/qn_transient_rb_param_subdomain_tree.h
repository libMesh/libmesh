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

#ifndef __qn_transient_rb_param_subdomain_tree_h__
#define __qn_transient_rb_param_subdomain_tree_h__

// Configuration data
#include "libmesh_config.h"

// This class depends on QNTransientSCMSystem, which is only
// available if SLEPc is defined.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "rb_param_subdomain_tree.h"

namespace libMesh
{

// Forward declaration
class RBSystem;
class QNTransientSCMSystem;

/**
 * This class is part of the rbOOmit framework.
 *
 * QNTransientRBParamSubdomainTree extends the "hp-tree"
 * data structure to also store a reference to a
 * QNTransientSCMSystem, which we need to use on each
 * leaf node of the tree in the time-dependent quadratically
 * nonlinear case.
 *
 * @author David J. Knezevic, 2009
 */

// ------------------------------------------------------------
// RBParamSubdomainTree class definition

class QNTransientRBParamSubdomainTree : public RBParamSubdomainTree
{
public:

  /**
   * Constructor. Initializes required
   * data structures.
   */
  QNTransientRBParamSubdomainTree (RBSystem& system,
                    const std::string& parameters_filename,
                    QNTransientSCMSystem& scm_system);

  /**
   * Destructor.
   */
  virtual ~QNTransientRBParamSubdomainTree() {}

  /**
   * The type of the parent.
   */
  typedef RBParamSubdomainTree Parent;

  /**
   * Build a QNTransientRBParamSubdomainNode. Overloaded virtual function.
   */
  virtual void build_root_node();

  /**
   * Reference to the RBSystem that is used
   * to actually perform the truth and RB solves.
   */
  QNTransientSCMSystem & _scm_system;

};

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
