// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes cont'd
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/fe_interface.h"
#include "libmesh/side.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfPrism6 class static member initializations
const int InfPrism6::num_nodes;
const int InfPrism6::num_sides;
const int InfPrism6::num_edges;
const int InfPrism6::num_children;
const int InfPrism6::nodes_per_side;
const int InfPrism6::nodes_per_edge;

const unsigned int InfPrism6::side_nodes_map[InfPrism6::num_sides][InfPrism6::nodes_per_side] =
  {
    { 0, 1, 2, 99}, // Side 0
    { 0, 1, 3, 4},  // Side 1
    { 1, 2, 4, 5},  // Side 2
    { 2, 0, 5, 3}   // Side 3
  };

const unsigned int InfPrism6::edge_nodes_map[InfPrism6::num_edges][InfPrism6::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {0, 2}, // Edge 2
    {0, 3}, // Edge 3
    {1, 4}, // Edge 4
    {2, 5}  // Edge 5
  };

// ------------------------------------------------------------
// InfPrism6 class member functions

bool InfPrism6::is_vertex(const unsigned int i) const
{
  if (i < 3)
    return true;
  return false;
}

bool InfPrism6::is_edge(const unsigned int i) const
{
  if (i < 3)
    return false;
  return true;
}

bool InfPrism6::is_face(const unsigned int) const
{
  return false;
}

bool InfPrism6::is_node_on_side(const unsigned int n,
                                const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfPrism6::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  auto trim = (s > 0) ? 0 : 1;
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s]) - trim};
}

std::vector<unsigned>
InfPrism6::nodes_on_edge(const unsigned int e) const
{
  libmesh_assert_less(e, n_edges());
  return {std::begin(edge_nodes_map[e]), std::end(edge_nodes_map[e])};
}

bool InfPrism6::is_node_on_edge(const unsigned int n,
                                const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}



Order InfPrism6::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> InfPrism6::build_side_ptr (const unsigned int i,
                                                 bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face;
  if (proxy)
    {
#ifdef LIBMESH_ENABLE_DEPRECATED
      libmesh_deprecated();
      switch (i)
        {
          // base
        case 0:
          {
            face = libmesh_make_unique<Side<Tri3,InfPrism6>>(this,i);
            break;
          }

          // ifem sides
        case 1:
        case 2:
        case 3:
          {
            face = libmesh_make_unique<Side<InfQuad4,InfPrism6>>(this,i);
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
#else
      libmesh_error();
#endif // LIBMESH_ENABLE_DEPRECATED
    }
  else
    {
      switch (i)
        {
        case 0: // the triangular face at z=-1, base face
          {
            face = libmesh_make_unique<Tri3>();
            break;
          }

        case 1: // the quad face at y=0
        case 2: // the other quad face
        case 3: // the quad face at x=0
          {
            face = libmesh_make_unique<InfQuad4>();
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      // Set the nodes
      for (auto n : face->node_index_range())
        face->set_node(n) = this->node_ptr(InfPrism6::side_nodes_map[i][n]);
    }

#ifdef LIBMESH_ENABLE_DEPRECATED
  if (!proxy) // proxy sides used to leave parent() set
#endif
    face->set_parent(nullptr);
  face->set_interior_parent(this);

  face->subdomain_id() = this->subdomain_id();
  face->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  face->set_p_level(this->p_level());
#endif

  return face;
}


void InfPrism6::build_side_ptr (std::unique_ptr<Elem> & side,
                                const unsigned int i)
{
  this->side_ptr(side, i);
}


std::unique_ptr<Elem> InfPrism6::build_edge_ptr (const unsigned int i)
{
  if (i < 3)
    return this->simple_build_edge_ptr<Edge2,InfPrism6>(i);

  // infinite edges
  return this->simple_build_edge_ptr<InfEdge2,InfPrism6>(i);
}



void InfPrism6::build_edge_ptr (std::unique_ptr<Elem> & edge,
                                const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  switch (i)
    {
      // the base edges
    case 0:
    case 1:
    case 2:
      {
        if (!edge.get() || edge->type() != EDGE2)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

      // the infinite edges
    case 3:
    case 4:
    case 5:
      {
        if (!edge.get() || edge->type() != INFEDGE2)
          {
            edge = this->build_edge_ptr(i);
            return;
          }
        break;
      }

    default:
      libmesh_error_msg("Invalid edge i = " << i);
    }

  edge->subdomain_id() = this->subdomain_id();
  edge->set_mapping_type(this->mapping_type());
#ifdef LIBMESH_ENABLE_AMR
  edge->set_p_level(this->p_level());
#endif

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n) = this->node_ptr(InfPrism6::edge_nodes_map[i][n]);
}



void InfPrism6::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[3] = this->node_id(2)+1;
        conn[4] = this->node_id(3)+1;
        conn[5] = this->node_id(4)+1;
        conn[6] = this->node_id(5)+1;
        conn[7] = this->node_id(5)+1;
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}





#ifdef LIBMESH_ENABLE_AMR

const Real InfPrism6::_embedding_matrix[InfPrism6::num_children][InfPrism6::num_nodes][InfPrism6::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //          0           1           2           3           4           5 th parent Node
      {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
    },

    // embedding matrix for child 1
    {
      //          0           1           2           3           4           5 th parent Node
      {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}  // 5
    },

    // embedding matrix for child 2
    {
      //          0           1           2           3           4           5 th parent Node
      {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
      {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
      {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0}  // 5
    },

    // embedding matrix for child 3
    {
      //          0           1           2           3           4           5 th parent Node
      {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
      {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
      {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
      {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
      {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
      {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
    }
  };


#endif


void
InfPrism6::permute(unsigned int perm_num)
{
  libmesh_assert_less (perm_num, 3);

  for (unsigned int i = 0; i != perm_num; ++i)
    {
      swap3nodes(0,1,2);
      swap3nodes(3,4,5);
      swap3neighbors(1,2,3);
    }
}


ElemType
InfPrism6::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 4);
  if (s == 0)
    return TRI3;
  return INFQUAD4;
}

} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
