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

// C++ includes

// Local includes
#include "libmesh/face_quad.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"


namespace libMesh
{



// ------------------------------------------------------------
// Quad class static member initializations


// We need to require C++11...
const Real Quad::_master_points[9][3] =
  {
    {-1, -1},
    {1, -1},
    {1, 1},
    {-1, 1},
    {0, -1},
    {1, 0},
    {0, 1},
    {-1, 0},
    {0, 0}
  };



// ------------------------------------------------------------
// Quad class member functions
dof_id_type Quad::key (const unsigned int s) const
{
  libmesh_assert_less (s, this->n_sides());

  return this->compute_key(this->node_id(Quad4::side_nodes_map[s][0]),
                           this->node_id(Quad4::side_nodes_map[s][1]));
}



unsigned int Quad::which_node_am_i(unsigned int side,
                                   unsigned int side_node) const
{
  libmesh_assert_less (side, this->n_sides());
  libmesh_assert_less (side_node, 2);

  return Quad4::side_nodes_map[side][side_node];
}



dof_id_type Quad::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2),
                           this->node_id(3));
}



UniquePtr<Elem> Quad::side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  Elem * edge = new Edge2;

  for (unsigned n=0; n<edge->n_nodes(); ++n)
    edge->set_node(n) = this->node_ptr(Quad4::side_nodes_map[i][n]);

  return UniquePtr<Elem>(edge);
}



bool Quad::is_child_on_side(const unsigned int c,
                            const unsigned int s) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (s, this->n_sides());

  // A quad's children and nodes don't share the same ordering:
  // child 2 and 3 are swapped;
  unsigned int n = (c < 2) ? c : 5-c;
  return (n == s || n == (s+1)%4);
}



unsigned int Quad::opposite_side(const unsigned int side_in) const
{
  libmesh_assert_less (side_in, 4);

  return (side_in + 2) % 4;
}



unsigned int Quad::opposite_node(const unsigned int node_in,
                                 const unsigned int side_in) const
{
  libmesh_assert_less (node_in, 8);
  libmesh_assert_less (node_in, this->n_nodes());
  libmesh_assert_less (side_in, this->n_sides());
  libmesh_assert(this->is_node_on_side(node_in, side_in));

  static const unsigned char side02_nodes_map[] =
    {3, 2, 1, 0, 6, 255, 4, 255};
  static const unsigned char side13_nodes_map[] =
    {1, 0, 3, 2, 255, 7, 255, 5};

  switch (side_in)
    {
    case 0:
    case 2:
      return side02_nodes_map[node_in];
    case 1:
    case 3:
      return side13_nodes_map[node_in];
    default:
      libmesh_error_msg("Unsupported side_in = " << side_in);
    }

  libmesh_error_msg("We'll never get here!");
  return 255;
}


Real Quad::quality (const ElemQuality q) const
{
  switch (q)
    {

      /**
       * Compue the min/max diagonal ratio.
       * This is modeled after the Hex element
       */
    case DISTORTION:
    case DIAGONAL:
    case STRETCH:
      {
        // Diagonal between node 0 and node 2
        const Real d02 = this->length(0,2);

        // Diagonal between node 1 and node 3
        const Real d13 = this->length(1,3);

        // Find the biggest and smallest diagonals
        if ((d02 > 0.) && (d13 >0.))
          if (d02 < d13) return d02 / d13;
          else return d13 / d02;
        else
          return 0.;
        break;
      }

    default:
      return Elem::quality(q);
    }

  /**
   * I don't know what to do for this metric.
   * Maybe the base class knows.  We won't get
   * here because of the default case above.
   */
  return Elem::quality(q);
}






std::pair<Real, Real> Quad::qual_bounds (const ElemQuality q) const
{
  std::pair<Real, Real> bounds;

  switch (q)
    {

    case ASPECT_RATIO:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case SKEW:
      bounds.first  = 0.;
      bounds.second = 0.5;
      break;

    case TAPER:
      bounds.first  = 0.;
      bounds.second = 0.7;
      break;

    case WARP:
      bounds.first  = 0.9;
      bounds.second = 1.;
      break;

    case STRETCH:
      bounds.first  = 0.25;
      bounds.second = 1.;
      break;

    case MIN_ANGLE:
      bounds.first  = 45.;
      bounds.second = 90.;
      break;

    case MAX_ANGLE:
      bounds.first  = 90.;
      bounds.second = 135.;
      break;

    case CONDITION:
      bounds.first  = 1.;
      bounds.second = 4.;
      break;

    case JACOBIAN:
      bounds.first  = 0.5;
      bounds.second = 1.;
      break;

    case SHEAR:
    case SHAPE:
    case SIZE:
      bounds.first  = 0.3;
      bounds.second = 1.;
      break;

    case DISTORTION:
      bounds.first  = 0.6;
      bounds.second = 1.;
      break;

    default:
      libMesh::out << "Warning: Invalid quality measure chosen." << std::endl;
      bounds.first  = -1;
      bounds.second = -1;
    }

  return bounds;
}




const unsigned short int Quad::_second_order_adjacent_vertices[4][2] =
  {
    {0, 1}, // vertices adjacent to node 4
    {1, 2}, // vertices adjacent to node 5
    {2, 3}, // vertices adjacent to node 6
    {0, 3}  // vertices adjacent to node 7
  };



const unsigned short int Quad::_second_order_vertex_child_number[9] =
  {
    99,99,99,99, // Vertices
    0,1,2,0,     // Edges
    0            // Interior
  };



const unsigned short int Quad::_second_order_vertex_child_index[9] =
  {
    99,99,99,99, // Vertices
    1,2,3,3,     // Edges
    2            // Interior
  };



#ifdef LIBMESH_ENABLE_AMR

// We number 25 "possible node locations" for a 2x2 refinement of
// quads with up to 3x3 nodes each
const int Quad::_child_node_lookup[4][9] =
  {
    // node lookup for child 0 (near node 0)
    { 0, 2, 12, 10,  1, 7, 11, 5,  6},

    // node lookup for child 1 (near node 1)
    { 2, 4, 14, 12,  3, 9, 13, 7,  8},

    // node lookup for child 2 (near node 3)
    { 10, 12, 22, 20,  11, 17, 21, 15,  16},

    // node lookup for child 3 (near node 2)
    { 12, 14, 24, 22,  13, 19, 23, 17,  18}
  };


#endif

} // namespace libMesh
