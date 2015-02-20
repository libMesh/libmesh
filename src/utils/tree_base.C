// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/tree_base.h"
#include "libmesh/mesh_base.h"

namespace libMesh
{
  const Elem* TreeBase::find_element(const Point& p,
                                     const std::set<subdomain_id_type> *allowed_subdomains,
                                     Real relative_tol) const
  {
    return this->find_element(p,this->mesh.mesh_dimension(),allowed_subdomains,relative_tol);
  }

} // namespace libMesh
