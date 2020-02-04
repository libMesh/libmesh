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



#ifndef LIBMESH_CELL_HEX_H
#define LIBMESH_CELL_HEX_H

// Local includes
#include "libmesh/cell.h"

namespace libMesh
{

template <typename>
class Hex8Templ;

/**
 * The \p Hex is an element in 3D with 6 sides.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all hexahedral element types.
 */
template <typename RealType = Real>
class HexTempl : public CellTempl<RealType>
{
public:
  typedef HexTempl<RealType> Hex;
  typedef Hex8Templ<RealType> Hex8;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;

  /**
   * Default brick element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  HexTempl(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
      Cell(nn, HexTempl<RealType>::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
      this->set_interior_parent(nullptr);
  }

  HexTempl (Hex &&) = delete;
  HexTempl (const Hex &) = delete;
  Hex & operator= (const Hex &) = delete;
  Hex & operator= (Hex &&) = delete;
  virtual ~HexTempl() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override final
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * \returns 6.
   */
  virtual unsigned int n_sides() const override final { return 6; }

  /**
   * \returns 8.  All hexahedra have 8 vertices.
   */
  virtual unsigned int n_vertices() const override final { return 8; }

  /**
   * \returns 12.  All hexahedra have 12 edges.
   */
  virtual unsigned int n_edges() const override final { return 12; }

  /**
   * \returns 6.  All hexahedra have 6 faces.
   */
  virtual unsigned int n_faces() const override final { return 6; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_children() const override final { return 8; }

  /**
   * \returns \p true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override final;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const override final;

  /**
   * \returns The side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const override final;

  /**
   * \returns The local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const override final;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using ElemTempl<RealType>::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Hex8::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns A primitive (4-noded) quad for face i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a primitive (4-noded) quad for face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override final;

  /**
   * \returns A quantitative assessment of element quality based on
   * the quality metric \p q specified by the user.
   */
  virtual Real quality (const ElemQuality q) const override;

  /**
   * \returns The suggested quality bounds for the hex based on
   * quality measure \p q.  These are the values suggested by the
   * CUBIT User's Manual.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const override;


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[7+(LIBMESH_DIM>3)];

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  This matrix
   * is kept here, since the matrix (for the first 12
   * higher-order nodes) is identical for \p Hex20 and
   * \p Hex27.
   */
  static const unsigned short int _second_order_adjacent_vertices[12][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[27];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[27];

  /**
   * Master element node locations
   */
  static const Real _master_points[27][3];

  /**
   * Lookup table from child id, child node id to "possible node
   * location" (a simple dictionary-index in a 5x5x5 grid)
   */
  static const int _child_node_lookup[8][27];
};

// ------------------------------------------------------------
// Hex class static member initializations

template <typename RealType>
const Real HexTempl<RealType>::_master_points[27][3] =
  {
    {-1, -1, -1},
    {1, -1, -1},
    {1, 1, -1},
    {-1, 1, -1},
    {-1, -1, 1},
    {1, -1, 1},
    {1, 1, 1},
    {-1, 1, 1},
    {0, -1, -1},
    {1, 0, -1},
    {0, 1, -1},
    {-1, 0, -1},
    {-1, -1, 0},
    {1, -1, 0},
    {1, 1, 0},
    {-1, 1, 0},
    {0, -1, 1},
    {1, 0, 1},
    {0, 1, 1},
    {-1, 0, 1},
    {0, 0, -1},
    {0, -1, 0},
    {1, 0, 0},
    {0, 1, 0},
    {-1, 0, 0},
    {0, 0, 1}
  };

template <typename RealType>
const unsigned short int HexTempl<RealType>::_second_order_vertex_child_number[27] =
  {
    99,99,99,99,99,99,99,99, // Vertices
    0,1,2,0,0,1,2,3,4,5,6,5, // Edges
    0,0,1,2,0,4,             // Faces
    0                        // Interior
  };



template <typename RealType>
const unsigned short int HexTempl<RealType>::_second_order_vertex_child_index[27] =
  {
    99,99,99,99,99,99,99,99, // Vertices
    1,2,3,3,4,5,6,7,5,6,7,7, // Edges
    2,5,6,7,7,6,             // Faces
    6                        // Interior
  };


template <typename RealType>
const unsigned short int HexTempl<RealType>::_second_order_adjacent_vertices[12][2] =
  {
    { 0,  1}, // vertices adjacent to node 8
    { 1,  2}, // vertices adjacent to node 9
    { 2,  3}, // vertices adjacent to node 10
    { 0,  3}, // vertices adjacent to node 11

    { 0,  4}, // vertices adjacent to node 12
    { 1,  5}, // vertices adjacent to node 13
    { 2,  6}, // vertices adjacent to node 14
    { 3,  7}, // vertices adjacent to node 15

    { 4,  5}, // vertices adjacent to node 16
    { 5,  6}, // vertices adjacent to node 17
    { 6,  7}, // vertices adjacent to node 18
    { 4,  7}  // vertices adjacent to node 19
  };


#ifdef LIBMESH_ENABLE_AMR

// We number 125 "possible node locations" for a 2x2x2 refinement of
// hexes with up to 3x3x3 nodes each
template <typename RealType>
const int HexTempl<RealType>::_child_node_lookup[8][27] =
  {
    // node lookup for child 0 (near node 0)
    { 0, 2, 12, 10,  50, 52, 62, 60,  1, 7, 11, 5,  25, 27, 37, 35,
      51, 57, 61, 55,  6,  26, 32, 36, 30,  56,  31},

    // node lookup for child 1 (near node 1)
    { 2, 4, 14, 12,  52, 54, 64, 62,  3, 9, 13, 7,  27, 29, 39, 37,
      53, 59, 63, 57,  8,  28, 34, 38, 32,  58,  33},

    // node lookup for child 2 (near node 3)
    { 10, 12, 22, 20,  60, 62, 72, 70,  11, 17, 21, 15,  35, 37, 47, 45,
      61, 67, 71, 65,  16,  36, 42, 46, 40,  66,  41},

    // node lookup for child 3 (near node 2)
    { 12, 14, 24, 22,  62, 64, 74, 72,  13, 19, 23, 17,  37, 39, 49, 47,
      63, 69, 73, 67,  18,  38, 44, 48, 42,  68,  43},

    // node lookup for child 4 (near node 4)
    { 50, 52, 62, 60,  100, 102, 112, 110,  51, 57, 61, 55,  75, 77, 87, 85,
      101, 107, 111, 105,  56,  76, 82, 86, 80,  106,  81},

    // node lookup for child 5 (near node 5)
    { 52, 54, 64, 62,  102, 104, 114, 112,  53, 59, 63, 57,  77, 79, 89, 87,
      103, 109, 113, 107,  58,  78, 84, 88, 82,  108,  93},

    // node lookup for child 6 (near node 7)
    { 60, 62, 72, 70,  110, 112, 122, 120,  61, 67, 71, 65,  85, 87, 97, 95,
      111, 117, 121, 115,  66,  86, 92, 96, 90,  116,  91},

    // node lookup for child 7 (near node 6)
    { 62, 64, 74, 72,  112, 114, 124, 122,  63, 69, 73, 67,  87, 89, 99, 97,
      113, 119, 123, 117,  68,  88, 94, 98, 92,  118,  103}
  };

#endif // LIBMESH_ENABLE_AMR

typedef HexTempl<Real> Hex;

} // namespace libMesh

#endif // LIBMESH_CELL_HEX_H
