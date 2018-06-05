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
#include "libmesh/face_inf_quad6.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/side.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// InfQuad6 class static member initializations
const unsigned int InfQuad6::side_nodes_map[3][3] =
  {
    {0, 1, 4},  // Side 0
    {1, 3, 99}, // Side 1
    {0, 2, 99}  // Side 2
  };


// ------------------------------------------------------------
// InfQuad6 class member functions

bool InfQuad6::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool InfQuad6::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}

bool InfQuad6::is_face(const unsigned int) const
{
  return false;
}

bool InfQuad6::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

#ifdef LIBMESH_ENABLE_AMR

const float InfQuad6::_embedding_matrix[2][6][6] =
  {
    // embedding matrix for child 0
    {
      //     0       1       2       3       4       5th parent node
      {    1.0,    0.0,    0.0,    0.0,    0.0,    0.0 }, // 0th child node
      {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 1
      {    0.0,    0.0,    1.0,    0.0,    0.0,    0.0 }, // 2
      {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 3
      {  0.375, -0.125,    0.0,    0.0,   0.75,    0.0 }, // 4
      {    0.0,    0.0,  0.375, -0.125,    0.0,   0.75 }  // 5
    },

    // embedding matrix for child 1
    {
      //     0       1       2       3       4       5th parent node
      {    0.0,    0.0,    0.0,    0.0,    1.0,    0.0 }, // 0th child node
      {    0.0,    1.0,    0.0,    0.0,    0.0,    0.0 }, // 1
      {    0.0,    0.0,    0.0,    0.0,    0.0,    1.0 }, // 2
      {    0.0,    0.0,    0.0,    1.0,    0.0,    0.0 }, // 3
      { -0.125,  0.375,    0.0,    0.0,   0.75,    0.0 }, // 4
      {    0.0,    0.0, -0.125,  0.375,    0.0,   0.75 }  // 5
    }
  };

#endif




Order InfQuad6::default_order() const
{
  return SECOND;
}



dof_id_type InfQuad6::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
      // Edge3 side
    case 0:
      return this->compute_key (this->node_id(4));

      // InfEdge
    case 1:
    case 2:
      return InfQuad::key(s);

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }
}



unsigned int InfQuad6::which_node_am_i(unsigned int side,
                                       unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert ((side == 0 && side_node < 3) || (side_node < 2));

  return InfQuad6::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> InfQuad6::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  // libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:
          return libmesh_make_unique<Side<Edge3,InfQuad6>>(this,i);

        case 1:
        case 2:
          return libmesh_make_unique<Side<InfEdge2,InfQuad6>>(this,i);

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Return value
      std::unique_ptr<Elem> edge;

      switch (i)
        {
        case 0:
          {
            edge = libmesh_make_unique<Edge3>();
            break;
          }

          // adjacent to another infinite element
        case 1:
        case 2:
          {
            edge = libmesh_make_unique<InfEdge2>();
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<edge->n_nodes(); ++n)
        edge->set_node(n) = this->node_ptr(InfQuad6::side_nodes_map[i][n]);

      return edge;
    }
}




void InfQuad6::connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        switch(sf)
          {
          case 0:
            // linear sub-quad 0
            conn[0] = this->node_id(0)+1;
            conn[1] = this->node_id(4)+1;
            conn[2] = this->node_id(5)+1;
            conn[3] = this->node_id(2)+1;

            return;

          case 1:
            // linear sub-quad 1
            conn[0] = this->node_id(4)+1;
            conn[1] = this->node_id(1)+1;
            conn[2] = this->node_id(3)+1;
            conn[3] = this->node_id(5)+1;

            return;

          default:
            libmesh_error_msg("Invalid sf = " << sf);
          }
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




unsigned short int InfQuad6::second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



const unsigned short int InfQuad6::_second_order_adjacent_vertices[2][2] =
  {
    {0, 1}, // vertices adjacent to node 4
    {2, 3}  // vertices adjacent to node 5
  };



std::pair<unsigned short int, unsigned short int>
InfQuad6::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  return std::pair<unsigned short int, unsigned short int>
    (0, 2*n-7);
}

} // namespace libMesh




#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
