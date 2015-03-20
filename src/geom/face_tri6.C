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

// C++ includes

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri6.h"

namespace libMesh
{




// ------------------------------------------------------------
// Tri6 class static member initializations
const unsigned int Tri6::side_nodes_map[3][3] =
  {
    {0, 1, 3}, // Side 0
    {1, 2, 4}, // Side 1
    {2, 0, 5}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

const float Tri6::_embedding_matrix[4][6][6] =
  {
    // embedding matrix for child 0
    {
      //  0      1      2    3    4    5
      { 1.0,   0.0,   0.0, 0.0, 0.0, 0.0}, // 0
      { 0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 1
      { 0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
      {.375, -.125,   0.0, .75, 0.0, 0.0}, // 3
      { 0.0, -.125, -.125, 0.5, .25, 0.5}, // 4
      {.375,   0.0, -.125, 0.0, 0.0, .75}  // 5
    },

    // embedding matrix for child 1
    {
      //  0      1      2    3    4    5
      {  0.0,  0.0,   0.0, 1.0, 0.0, 0.0}, // 0
      {  0.0,  1.0,   0.0, 0.0, 0.0, 0.0}, // 1
      {  0.0,  0.0,   0.0, 0.0, 1.0, 0.0}, // 2
      {-.125, .375,   0.0, .75, 0.0, 0.0}, // 3
      {  0.0, .375, -.125, 0.0, .75, 0.0}, // 4
      {-.125,  0.0, -.125, 0.5, 0.5, .25}  // 5
    },

    // embedding matrix for child 2
    {
      //  0       1     2    3    4    5
      {  0.0,   0.0,  0.0, 0.0, 0.0, 1.0}, // 0
      {  0.0,   0.0,  0.0, 0.0, 1.0, 0.0}, // 1
      {  0.0,   0.0,  1.0, 0.0, 0.0, 0.0}, // 2
      {-.125, -.125,  0.0, .25, 0.5, 0.5}, // 3
      {  0.0, -.125, .375, 0.0, .75, 0.0}, // 4
      {-.125,   0.0, .375, 0.0, 0.0, .75}  // 5
    },

    // embedding matrix for child 3
    {
      //  0       1      2    3    4    5
      {  0.0,   0.0,   0.0, 1.0, 0.0, 0.0}, // 0
      {  0.0,   0.0,   0.0, 0.0, 1.0, 0.0}, // 1
      {  0.0,   0.0,   0.0, 0.0, 0.0, 1.0}, // 2
      {-.125,   0.0, -.125, 0.5, 0.5, .25}, // 3
      {-.125, -.125,   0.0, .25, 0.5, 0.5}, // 4
      {  0.0, -.125, -.125, 0.5, .25, 0.5}  // 5
    }
  };

#endif



// ------------------------------------------------------------
// Tri6 class member functions

bool Tri6::is_vertex(const unsigned int i) const
{
  if (i < 3)
    return true;
  return false;
}

bool Tri6::is_edge(const unsigned int i) const
{
  if (i < 3)
    return false;
  return true;
}

bool Tri6::is_face(const unsigned int) const
{
  return false;
}

bool Tri6::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}



bool Tri6::has_affine_map() const
{
  // Make sure edges are straight
  if (!this->point(3).relative_fuzzy_equals
      ((this->point(0) + this->point(1))/2.))
    return false;
  if (!this->point(4).relative_fuzzy_equals
      ((this->point(1) + this->point(2))/2.))
    return false;
  if (!this->point(5).relative_fuzzy_equals
      ((this->point(2) + this->point(0))/2.))
    return false;

  return true;
}



dof_id_type Tri6::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:

      return
        this->compute_key (this->node(3));

    case 1:

      return
        this->compute_key (this->node(4));

    case 2:

      return
        this->compute_key (this->node(5));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



UniquePtr<Elem> Tri6::build_side (const unsigned int i,
                                  bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return UniquePtr<Elem>(new Side<Edge3,Tri6>(this,i));

  else
    {
      Elem* edge = new Edge3;
      edge->subdomain_id() = this->subdomain_id();

      switch (i)
        {
        case 0:
          {
            edge->set_node(0) = this->get_node(0);
            edge->set_node(1) = this->get_node(1);
            edge->set_node(2) = this->get_node(3);
            break;
          }
        case 1:
          {
            edge->set_node(0) = this->get_node(1);
            edge->set_node(1) = this->get_node(2);
            edge->set_node(2) = this->get_node(4);
            break;
          }
        case 2:
          {
            edge->set_node(0) = this->get_node(2);
            edge->set_node(1) = this->get_node(0);
            edge->set_node(2) = this->get_node(5);
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      return UniquePtr<Elem>(edge);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}


void Tri6::connectivity(const unsigned int sf,
                        const IOPackage iop,
                        std::vector<dof_id_type>& conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        switch(sf)
          {
          case 0:
            // linear sub-triangle 0
            conn[0] = this->node(0)+1;
            conn[1] = this->node(3)+1;
            conn[2] = this->node(5)+1;
            conn[3] = this->node(5)+1;

            return;

          case 1:
            // linear sub-triangle 1
            conn[0] = this->node(3)+1;
            conn[1] = this->node(1)+1;
            conn[2] = this->node(4)+1;
            conn[3] = this->node(4)+1;

            return;

          case 2:
            // linear sub-triangle 2
            conn[0] = this->node(5)+1;
            conn[1] = this->node(4)+1;
            conn[2] = this->node(2)+1;
            conn[3] = this->node(2)+1;

            return;

          case 3:
            // linear sub-triangle 3
            conn[0] = this->node(3)+1;
            conn[1] = this->node(4)+1;
            conn[2] = this->node(5)+1;
            conn[3] = this->node(5)+1;

            return;

          default:
            libmesh_error_msg("Invalid sf = " << sf);
          }
      }

    case VTK:
      {
        // VTK_QUADRATIC_TRIANGLE has same numbering as libmesh TRI6
        conn.resize(6);
        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        conn[3] = this->node(3);
        conn[4] = this->node(4);
        conn[5] = this->node(5);
        return;

        // Used to write out linear sub-triangles for VTK...
        /*
          conn.resize(3);
          switch(sf)
          {
          case 0:
          // linear sub-triangle 0
          conn[0] = this->node(0);
          conn[1] = this->node(3);
          conn[2] = this->node(5);

          return;

          case 1:
          // linear sub-triangle 1
          conn[0] = this->node(3);
          conn[1] = this->node(1);
          conn[2] = this->node(4);

          return;

          case 2:
          // linear sub-triangle 2
          conn[0] = this->node(5);
          conn[1] = this->node(4);
          conn[2] = this->node(2);

          return;

          case 3:
          // linear sub-triangle 3
          conn[0] = this->node(3);
          conn[1] = this->node(4);
          conn[2] = this->node(5);

          return;

          default:
          libmesh_error_msg("Invalid sf = " << sf);
          }
        */
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}





unsigned short int Tri6::second_order_adjacent_vertex (const unsigned int n,
                                                       const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  libmesh_assert_less (v, 2);
  return _second_order_adjacent_vertices[n-this->n_vertices()][v];
}



const unsigned short int Tri6::_second_order_adjacent_vertices[3][2] =
  {
    {0, 1}, // vertices adjacent to node 3
    {1, 2}, // vertices adjacent to node 4
    {0, 2}  // vertices adjacent to node 5
  };



std::pair<unsigned short int, unsigned short int>
Tri6::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}



const unsigned short int Tri6::_second_order_vertex_child_number[6] =
  {
    99,99,99, // Vertices
    0,1,0     // Edges
  };



const unsigned short int Tri6::_second_order_vertex_child_index[6] =
  {
    99,99,99, // Vertices
    1,2,2     // Edges
  };

} // namespace libMesh
