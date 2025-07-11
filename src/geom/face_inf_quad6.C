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



#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/face_inf_quad6.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{




// ------------------------------------------------------------
// InfQuad6 class static member initializations
const int InfQuad6::num_nodes;
const int InfQuad6::nodes_per_side;

const unsigned int InfQuad6::side_nodes_map[InfQuad6::num_sides][InfQuad6::nodes_per_side] =
  {
    {0, 1, 4},  // Side 0
    {1, 3, 99}, // Side 1
    {0, 2, 99}  // Side 2
  };


// ------------------------------------------------------------
// InfQuad6 class member functions

bool InfQuad6::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfQuad6::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s == 0) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
InfQuad6::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
}

#ifdef LIBMESH_ENABLE_AMR

const Real InfQuad6::_embedding_matrix[InfQuad6::num_children][InfQuad6::num_nodes][InfQuad6::num_nodes] =
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



unsigned int InfQuad6::local_side_node(unsigned int side,
                                       unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert ((side == 0 && side_node < InfQuad6::nodes_per_side) || (side_node < 2));

  return InfQuad6::side_nodes_map[side][side_node];
}



std::unique_ptr<Elem> InfQuad6::build_side_ptr (const unsigned int i)
{
  // libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> edge;

  switch (i)
    {
    case 0:
      {
        edge = std::make_unique<Edge3>();
        break;
      }

      // adjacent to another infinite element
    case 1:
    case 2:
      {
        edge = std::make_unique<InfEdge2>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(InfQuad6::side_nodes_map[i][n]));

  edge->set_interior_parent(this);
  edge->inherit_data_from(*this);

  return edge;
}



void InfQuad6::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  // Think of a unit cube: (-1,1) x (-1,1) x (1,1)
  switch (i)
    {
      // the base face
    case 0:
      {
        if (!side.get() || side->type() != EDGE3)
          {
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }

      // connecting to another infinite element
    case 1:
    case 2:
      {
        if (!side.get() || side->type() != INFEDGE2)
          {
            side = this->build_side_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  side->subdomain_id() = this->subdomain_id();
  side->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  side->set_p_level(this->p_level());
#endif
  // Set the nodes
  for (auto n : side->node_index_range())
    side->set_node(n, this->node_ptr(InfQuad6::side_nodes_map[i][n]));
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


ElemType
InfQuad6::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 3);
  if (s == 0)
    return EDGE3;
  return INFEDGE2;
}


void InfQuad6::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2neighbors(1,2);
  swap2boundarysides(1,2,boundary_info);
  swap2boundaryedges(1,2,boundary_info);
}


} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
