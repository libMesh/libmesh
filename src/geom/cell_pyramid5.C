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


// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid5 class static member initializations
const unsigned int Pyramid5::side_nodes_map[5][4] =
  {
    {0, 1, 4, 99}, // Side 0
    {1, 2, 4, 99}, // Side 1
    {2, 3, 4, 99}, // Side 2
    {3, 0, 4, 99}, // Side 3
    {0, 3, 2,  1}  // Side 4
  };

const unsigned int Pyramid5::edge_nodes_map[8][2] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 4}, // Edge 5
    {2, 4}, // Edge 6
    {3, 4}  // Edge 7
  };



// ------------------------------------------------------------
// Pyramid5 class member functions

bool Pyramid5::is_vertex(const unsigned int) const
{
  return true;
}

bool Pyramid5::is_edge(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_face(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Pyramid5::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Pyramid5::has_affine_map() const
{
  //  Point v = this->point(3) - this->point(0);
  //  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
  return false;
}



Order Pyramid5::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> Pyramid5::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:
        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<Tri3,Pyramid5>>(this,i);

        case 4:
          return libmesh_make_unique<Side<Quad4,Pyramid5>>(this,i);

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Return value
      std::unique_ptr<Elem> face;

      switch (i)
        {
        case 0: // triangular face 1
        case 1: // triangular face 2
        case 2: // triangular face 3
        case 3: // triangular face 4
          {
            face = libmesh_make_unique<Tri3>();
            break;
          }
        case 4: // the quad face at z=0
          {
            face = libmesh_make_unique<Quad4>();
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(Pyramid5::side_nodes_map[i][n]);

      return face;
    }
}



std::unique_ptr<Elem> Pyramid5::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2,Pyramid5>>(this,i);
}



void Pyramid5::connectivity(const unsigned int libmesh_dbg_var(sc),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        conn[4] = this->node_id(4)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(4)+1;
        conn[7] = this->node_id(4)+1;
        return;
      }

    case VTK:
      {
        conn.resize(5);
        conn[0] = this->node_id(3);
        conn[1] = this->node_id(2);
        conn[2] = this->node_id(1);
        conn[3] = this->node_id(0);
        conn[4] = this->node_id(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


Real Pyramid5::volume () const
{
  // The pyramid with a bilinear base has volume given by the
  // formula in: "Calculation of the Volume of a General Hexahedron
  // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
  Point
    x0 = point(0), x1 = point(1), x2 = point(2),
    x3 = point(3), x4 = point(4);

  // Construct various edge and diagonal vectors.
  Point v40 = x0 - x4;
  Point v13 = x3 - x1;
  Point v02 = x2 - x0;
  Point v03 = x3 - x0;
  Point v01 = x1 - x0;

  // Finally, ready to return the volume!
  return
    triple_product(v40, v13, v02) / 6. +
    triple_product(v02, v01, v03) / 12.;
}

} // namespace libMesh
