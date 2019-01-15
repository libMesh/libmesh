// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/edge_edge2.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// Edge2 class static member initializations
const int Edge2::num_nodes;
const int Edge2::num_children;

#ifdef LIBMESH_ENABLE_AMR

const float Edge2::_embedding_matrix[Edge2::num_children][Edge2::num_nodes][Edge2::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1
      {1.0, 0.0}, // 0
      {0.5, 0.5}  // 1
    },

    // embedding matrix for child 1
    {
      // 0    1
      {0.5, 0.5}, // 0
      {0.0, 1.0}  // 1
    }
  };

#endif

bool Edge2::is_vertex(const unsigned int) const
{
  return true;
}

bool Edge2::is_edge(const unsigned int) const
{
  return false;
}

bool Edge2::is_face(const unsigned int) const
{
  return false;
}

bool Edge2::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, Edge2::num_nodes);
  return (s == n);
}

bool Edge2::is_node_on_edge(const unsigned int,
                            const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



Order Edge2::default_order() const
{
  return FIRST;
}



void Edge2::connectivity(const unsigned int libmesh_dbg_var(sc),
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_equal_to (sc, 0);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
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


Real Edge2::volume () const
{
  // OK, so this is probably overkill, since it is equivalent to
  // Elem::hmax() for the Edge2, but here it is nonetheless...
  return (this->point(1) - this->point(0)).norm();
}



dof_id_type Edge2::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1));
}

} // namespace libMesh
