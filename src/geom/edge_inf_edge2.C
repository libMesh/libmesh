// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfEdge2 class member functions

bool InfEdge2::is_vertex(const unsigned int i) const
{
  if (i)
    return false;
  return true;
}

bool InfEdge2::is_edge(const unsigned int i) const
{
  if (i)
    return true;
  return false;
}

bool InfEdge2::is_face(const unsigned int) const
{
  return false;
}

bool InfEdge2::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, 1);
  return (s == n);
}

bool InfEdge2::is_node_on_edge(const unsigned int,
                               const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}

Order InfEdge2::default_order() const
{
  return FIRST;
}

void InfEdge2::connectivity(const unsigned int libmesh_dbg_var(se),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert_equal_to (se, 0);
  libmesh_assert_less (se, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}

} // namespace libMesh


#endif
