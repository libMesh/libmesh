// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes

// Local includes cont'd
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/fe_interface.h"
#include "libmesh/side.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{


// ------------------------------------------------------------
// InfPrism6 class static member initializations
const unsigned int InfPrism6::side_nodes_map[4][4] =
  {
    { 0, 1, 2, 99}, // Side 0
    { 0, 1, 3, 4},  // Side 1
    { 1, 2, 4, 5},  // Side 2
    { 2, 0, 5, 3}   // Side 3
  };

const unsigned int InfPrism6::edge_nodes_map[6][2] =
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
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool InfPrism6::is_node_on_edge(const unsigned int n,
                                const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}


std::unique_ptr<Elem> InfPrism6::build_side_ptr (const unsigned int i,
                                                 bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
          // base
        case 0:
          return libmesh_make_unique<Side<Tri3,InfPrism6>>(this,i);

          // ifem sides
        case 1:
        case 2:
        case 3:
          return libmesh_make_unique<Side<InfQuad4,InfPrism6>>(this,i);

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

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(InfPrism6::side_nodes_map[i][n]);

      return face;
    }
}


std::unique_ptr<Elem> InfPrism6::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, n_edges());

  if (i < 3)
    return libmesh_make_unique<SideEdge<Edge2,InfPrism6>>(this,i);

  // infinite edges
  return libmesh_make_unique<SideEdge<InfEdge2,InfPrism6>>(this,i);
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

const float InfPrism6::_embedding_matrix[4][6][6] =
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

} // namespace libMesh

#endif  // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
