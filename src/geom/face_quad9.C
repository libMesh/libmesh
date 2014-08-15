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
#include "libmesh/face_quad9.h"

namespace libMesh
{




// ------------------------------------------------------------
// Quad9 class static member initializations
const unsigned int Quad9::side_nodes_map[4][3] =
  {
    {0, 1, 4}, // Side 0
    {1, 2, 5}, // Side 1
    {2, 3, 6}, // Side 2
    {3, 0, 7}  // Side 3
  };


#ifdef LIBMESH_ENABLE_AMR

const float Quad9::_embedding_matrix[4][9][9] =
  {
    // embedding matrix for child 0
    {
      //         0           1           2           3           4           5           6           7           8
      {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 3
      {   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 6
      {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
      {   0.140625, -0.0468750,  0.0156250, -0.0468750,   0.281250, -0.0937500, -0.0937500,   0.281250,   0.562500 }  // 8
    },

    // embedding matrix for child 1
    {
      //         0           1           2           3           4           5           6           7           8
      {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
      {  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,    0.00000,   0.750000 }, // 7
      { -0.0468750,   0.140625, -0.0468750,  0.0156250,   0.281250,   0.281250, -0.0937500, -0.0937500,   0.562500 }  // 8
    },

    // embedding matrix for child 2
    {
      //         0           1           2           3           4           5           6           7           8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,   0.750000 }, // 4
      {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 5
      {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
      {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 7
      { -0.0468750,  0.0156250, -0.0468750,   0.140625, -0.0937500, -0.0937500,   0.281250,   0.281250,   0.562500 }  // 8
    },

    // embedding matrix for child 3
    {
      //         0           1           2           3           4           5           6           7           8
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 3
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,   0.375000,    0.00000,  -0.125000,   0.750000 }, // 4
      {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 6
      {    0.00000,    0.00000,    0.00000,    0.00000,  -0.125000,    0.00000,   0.375000,    0.00000,   0.750000 }, // 7
      {  0.0156250, -0.0468750,   0.140625, -0.0468750, -0.0937500,   0.281250,   0.281250, -0.0937500,   0.562500 }  // 8
    }
  };

#endif



// ------------------------------------------------------------
// Quad9 class member functions

bool Quad9::is_vertex(const unsigned int i) const
{
  if (i < 4)
    return true;
  return false;
}

bool Quad9::is_edge(const unsigned int i) const
{
  if (i < 4)
    return false;
  if (i > 7)
    return false;
  return true;
}

bool Quad9::is_face(const unsigned int i) const
{
  if (i > 7)
    return true;
  return false;
}

bool Quad9::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}



bool Quad9::has_affine_map() const
{
  // make sure corners form a parallelogram
  Point v = this->point(1) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(2) - this->point(3)))
    return false;
  // make sure "horizontal" sides are straight
  v /= 2;
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(6) - this->point(3)))
    return false;
  // make sure "vertical" sides are straight
  // and the center node is centered
  v = (this->point(3) - this->point(0))/2;
  if (!v.relative_fuzzy_equals(this->point(7) - this->point(0)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(8) - this->point(4)))
    return false;
  return true;
}



dof_id_type Quad9::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  switch (s)
    {
    case 0:

      return
        this->compute_key (this->node(4));

    case 1:

      return
        this->compute_key (this->node(5));

    case 2:

      return
        this->compute_key (this->node(6));

    case 3:

      return
        this->compute_key (this->node(7));

    default:
      libmesh_error_msg("Invalid side s = " << s);
    }

  libmesh_error_msg("We'll never get here!");
  return 0;
}



AutoPtr<Elem> Quad9::build_side (const unsigned int i,
                                 bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return AutoPtr<Elem>(new Side<Edge3,Quad9>(this,i));

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
            edge->set_node(2) = this->get_node(4);
            break;
          }
        case 1:
          {
            edge->set_node(0) = this->get_node(1);
            edge->set_node(1) = this->get_node(2);
            edge->set_node(2) = this->get_node(5);
            break;
          }
        case 2:
          {
            edge->set_node(0) = this->get_node(2);
            edge->set_node(1) = this->get_node(3);
            edge->set_node(2) = this->get_node(6);
            break;
          }
        case 3:
          {
            edge->set_node(0) = this->get_node(3);
            edge->set_node(1) = this->get_node(0);
            edge->set_node(2) = this->get_node(7);
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      return AutoPtr<Elem>(edge);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<Elem>();
}







void Quad9::connectivity(const unsigned int sf,
                         const IOPackage iop,
                         std::vector<dof_id_type>& conn) const
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
            conn[0] = this->node(0)+1;
            conn[1] = this->node(4)+1;
            conn[2] = this->node(8)+1;
            conn[3] = this->node(7)+1;
            return;

          case 1:
            // linear sub-quad 1
            conn[0] = this->node(4)+1;
            conn[1] = this->node(1)+1;
            conn[2] = this->node(5)+1;
            conn[3] = this->node(8)+1;
            return;

          case 2:
            // linear sub-quad 2
            conn[0] = this->node(7)+1;
            conn[1] = this->node(8)+1;
            conn[2] = this->node(6)+1;
            conn[3] = this->node(3)+1;
            return;

          case 3:
            // linear sub-quad 3
            conn[0] = this->node(8)+1;
            conn[1] = this->node(5)+1;
            conn[2] = this->node(2)+1;
            conn[3] = this->node(6)+1;
            return;

          default:
            libmesh_error_msg("Invalid sf = " << sf);
          }
      }

    case VTK:
      {
        conn.resize(9);
        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        conn[3] = this->node(3);
        conn[4] = this->node(4);
        conn[5] = this->node(5);
        conn[6] = this->node(6);
        conn[7] = this->node(7);
        conn[8] = this->node(8);
        return;

        /*
          switch(sf)
          {
          case 0:
          // linear sub-quad 0
          conn[0] = this->node(0);
          conn[1] = this->node(4);
          conn[2] = this->node(8);
          conn[3] = this->node(7);

          return;

          case 1:
          // linear sub-quad 1
          conn[0] = this->node(4);
          conn[1] = this->node(1);
          conn[2] = this->node(5);
          conn[3] = this->node(8);

          return;

          case 2:
          // linear sub-quad 2
          conn[0] = this->node(7);
          conn[1] = this->node(8);
          conn[2] = this->node(6);
          conn[3] = this->node(3);

          return;

          case 3:
          // linear sub-quad 3
          conn[0] = this->node(8);
          conn[1] = this->node(5);
          conn[2] = this->node(2);
          conn[3] = this->node(6);

          return;

          default:
          libmesh_error_msg("Invalid sf = " << sf);
          }*/
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}




unsigned int Quad9::n_second_order_adjacent_vertices (const unsigned int n) const
{
  switch (n)
    {
    case 4:
    case 5:
    case 6:
    case 7:
      return 2;

    case 8:
      return 4;

    default:
      libmesh_error_msg("Invalid n = " << n);
    }

  libmesh_error_msg("We'll never get here!");
  return libMesh::invalid_uint;
}



unsigned short int Quad9::second_order_adjacent_vertex (const unsigned int n,
                                                        const unsigned int v) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());

  switch (n)
    {
    case 8:
      {
        libmesh_assert_less (v, 4);
        return static_cast<unsigned short int>(v);
      }

    default:
      {
        libmesh_assert_less (v, 2);
        // use the matrix that we inherited from \p Quad
        return _second_order_adjacent_vertices[n-this->n_vertices()][v];
      }
    }
}



std::pair<unsigned short int, unsigned short int>
Quad9::second_order_child_vertex (const unsigned int n) const
{
  libmesh_assert_greater_equal (n, this->n_vertices());
  libmesh_assert_less (n, this->n_nodes());
  /*
   * the _second_order_vertex_child_* vectors are
   * stored in face_quad.C, since they are identical
   * for Quad8 and Quad9 (for the first 4 higher-order nodes)
   */
  return std::pair<unsigned short int, unsigned short int>
    (_second_order_vertex_child_number[n],
     _second_order_vertex_child_index[n]);
}

} // namespace libMesh
