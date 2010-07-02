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

// Configuration data
#include "libmesh_config.h"

#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "qn_transient_rb_param_subdomain_tree.h"
#include "qn_transient_rb_param_subdomain_node.h"
#include "qn_transient_scm_system.h"
#include "rb_system.h"

namespace libMesh
{


QNTransientRBParamSubdomainTree::QNTransientRBParamSubdomainTree(RBSystem& system,
                    const std::string& parameters_filename,
                    QNTransientSCMSystem& scm_system)
  : Parent(system, parameters_filename),
    _scm_system(scm_system)
{ }

void QNTransientRBParamSubdomainTree::build_root_node()
{
  libmesh_assert(root_node == NULL);
  root_node = new QNTransientRBParamSubdomainNode(*this, _rb_system.get_current_parameters());
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
