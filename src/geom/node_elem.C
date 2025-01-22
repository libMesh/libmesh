// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/node_elem.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// NodeElem class static member initializations
const int NodeElem::num_nodes;
const int NodeElem::num_sides;

// ------------------------------------------------------------
// NodeElem class member functions

// If we're computing on a node, we only need CONSTANT approximation
// to get the full function space
Order NodeElem::default_order() const
{
  return CONSTANT;
}

// But we can "support" any order in 0D, really...
Order NodeElem::supported_nodal_order() const
{
  return MAXIMUM;
}



void NodeElem::connectivity(const unsigned int,
                            const IOPackage,
                            std::vector<dof_id_type> &) const
{
  libmesh_not_implemented();
}

bool NodeElem::contains_point(const Point & p, Real tol) const
{
  return (this->point(0) - p).norm() < tol;
}

bool NodeElem::close_to_point(const Point & p, Real tol) const
{
  return this->contains_point(p, tol);
}

#ifdef LIBMESH_ENABLE_AMR

const Real NodeElem::_embedding_matrix[1][1][1] =
  {
    // embedding matrix for child 0
    {
      // 0
      {1.0}, // 0
    }
  };

#endif

} // namespace libMesh
