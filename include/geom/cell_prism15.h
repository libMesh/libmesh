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



#ifndef LIBMESH_CELL_PRISM15_H
#define LIBMESH_CELL_PRISM15_H

// Local includes
#include "libmesh/cell_prism.h"

namespace libMesh
{
template <typename> class Tri6Templ;
template <typename> class Quad8Templ;
template <typename> class Edge3Templ;
template <typename> class ElemTempl;
template <typename> class NodeTempl;
template <typename> class PointTempl;
template <typename> class BoundingBoxTempl;

/**
 * The \p Prism15 is an element in 3D composed of 15 nodes.
 * It is numbered like this:
 * \verbatim
 *   PRISM15:
 *                5
 *                o
 *               /:\
 *              / : \
 *             /  :  \
 *            /   :   \
 *        14 o    :    o 13
 *          /     :     \
 *         /      :      \
 *        /       o 11    \
 *     3 /        :        \4
 *      o---------o---------o
 *      |         :12       |
 *      |         :         |
 *      |         :         |            zeta
 *      |         o         |             ^   eta (into page)
 *      |        .2.        |             | /
 *      |       .   .       |             |/
 *    9 o      .     .      o 10          o---> xi
 *      |     .       .     |
 *      |  8 o         o 7  |
 *      |   .           .   |
 *      |  .             .  |
 *      | .               . |
 *      |.                 .|
 *      o---------o---------o
 *      0         6         1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief A 3D prismatic element with 15 nodes.
 */
template <typename RealType = Real>
class Prism15Templ final : public PrismTempl<RealType>
{
public:
  typedef PrismTempl<RealType> Prism;
  typedef Prism15Templ<RealType> Prism15;
  typedef Quad8Templ<RealType> Quad8;
  typedef Tri6Templ<RealType> Tri6;
  typedef Edge3Templ<RealType> Edge3;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Prism15Templ (Elem * p=nullptr) :
    Prism(Prism15::n_nodes(), p, _nodelinks_data)
  {}

  Prism15Templ (Prism15 &&) = delete;
  Prism15Templ (const Prism15 &) = delete;
  Prism15 & operator= (const Prism15 &) = delete;
  Prism15 & operator= (Prism15 &&) = delete;
  virtual ~Prism15Templ() = default;

  /**
   * \returns \p PRISM15.
   */
  virtual ElemType type () const override { return PRISM15; }

  /**
   * \returns 15.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 1.
   */
  virtual unsigned int n_sub_elem() const override { return 1; }

  /**
   * \returns \p true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * \returns \p Prism15::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * Builds a \p QUAD8 or \p TRI6 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p QUAD8 or \p TRI6 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE3 or \p INFEDGE2 coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;
  /**
   * \returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { return 2; }

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 6 \le n < 15 \f$.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const override;

  /**
   * \returns The child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  See
   * elem.h for further details.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * Geometric constants for Prism15.
   */
  static const int num_nodes = 15;
  static const int num_sides = 5;
  static const int num_edges = 9;
  static const int num_children = 8;
  static const int nodes_per_side = 8;
  static const int nodes_per_edge = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[num_edges][nodes_per_edge];

  /**
   * A specialization for computing the volume of a Prism15.
   */
  virtual RealType volume () const override;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[num_nodes];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const override
  { return _embedding_matrix[i][j][k]; }

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[num_children][num_nodes][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif

};

template <typename RealType>
const int Prism15Templ<RealType>::nodes_per_side;

// ------------------------------------------------------------
// Prism15 class static member initializations
template <typename RealType>
const unsigned int Prism15Templ<RealType>::side_nodes_map[Prism15::num_sides][Prism15::nodes_per_side] =
  {
    {0, 2, 1,  8,  7,  6, 99, 99}, // Side 0
    {0, 1, 4,  3,  6, 10, 12,  9}, // Side 1
    {1, 2, 5,  4,  7, 11, 13, 10}, // Side 2
    {2, 0, 3,  5,  8,  9, 14, 11}, // Side 3
    {3, 4, 5, 12, 13, 14, 99, 99}  // Side 4
  };

template <typename RealType>
const unsigned int Prism15Templ<RealType>::edge_nodes_map[Prism15::num_edges][Prism15::nodes_per_edge] =
  {
    {0, 1,  6}, // Edge 0
    {1, 2,  7}, // Edge 1
    {0, 2,  8}, // Edge 2
    {0, 3,  9}, // Edge 3
    {1, 4, 10}, // Edge 4
    {2, 5, 11}, // Edge 5
    {3, 4, 12}, // Edge 6
    {4, 5, 13}, // Edge 7
    {3, 5, 14}  // Edge 8
  };

#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Prism15Templ<RealType>::_embedding_matrix[Prism15::num_children][Prism15::num_nodes][Prism15::num_nodes] =
  {
    // Embedding matrix for child 0
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0 }, //  3
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  4
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  5
      {   0.375,  -0.125,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,  -0.125,  -0.125,       0,       0,       0,     0.5,    0.25,     0.5,       0,       0,       0,       0,       0,       0 }, //  7
      {   0.375,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0 }, //  8
      {   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0 }, //  9
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, // 11
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0,       0 }, // 12
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, // 13
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.75,       0,    0.25,       0,       0,   0.375 }  // 14
    },

    // Embedding matrix for child 1
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0 }, //  4
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  5
      {  -0.125,   0.375,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,   0.375,  -0.125,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0 }, //  7
      {  -0.125,       0,  -0.125,       0,       0,       0,     0.5,     0.5,    0.25,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, //  9
      {       0,   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0 }, // 10
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 11
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0,       0 }, // 12
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0 }, // 13
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }  // 14
    },

    // Embedding matrix for child 2
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  3
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0 }, //  5
      {  -0.125,  -0.125,       0,       0,       0,       0,    0.25,     0.5,     0.5,       0,       0,       0,       0,       0,       0 }, //  6
      {       0,  -0.125,   0.375,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0,       0 }, //  7
      {  -0.125,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 10
      {       0,       0,   0.375,       0,       0,  -0.125,       0,       0,       0,       0,       0,    0.75,       0,       0,       0 }, // 11
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, // 12
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0 }, // 13
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.25,       0,    0.75,       0,       0,   0.375 }  // 14
    },

    // Embedding matrix for child 3
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0 }, //  2
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  3
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  4
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  5
      {  -0.125,       0,  -0.125,       0,       0,       0,     0.5,     0.5,    0.25,       0,       0,       0,       0,       0,       0 }, //  6
      {  -0.125,  -0.125,       0,       0,       0,       0,    0.25,     0.5,     0.5,       0,       0,       0,       0,       0,       0 }, //  7
      {       0,  -0.125,  -0.125,       0,       0,       0,     0.5,    0.25,     0.5,       0,       0,       0,       0,       0,       0 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0,       0 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.75,       0,       0,   0.375,   0.375,       0,    0.25,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.75,   0.375,       0,   0.375,       0,       0,    0.25 }, // 11
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, // 12
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, // 13
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }  // 14
    },

    // Embedding matrix for child 4
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0 }, //  0
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  1
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  2
      {       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  5
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0,       0 }, //  6
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, //  7
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.75,       0,    0.25,       0,       0,   0.375 }, //  8
      {  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0,       0 }, //  9
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, // 11
      {       0,       0,       0,   0.375,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.75,       0,       0 }, // 12
      {       0,       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,    0.25,     0.5 }, // 13
      {       0,       0,       0,   0.375,       0,  -0.125,       0,       0,       0,       0,       0,       0,       0,       0,    0.75 }  // 14
    },

    // Embedding matrix for child 5
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  0
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0 }, //  1
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  3
      {       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  5
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0,       0 }, //  6
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.75,    0.25,       0,   0.375,       0 }, //  7
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, //  9
      {       0,  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0,       0 }, // 10
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 11
      {       0,       0,       0,  -0.125,   0.375,       0,       0,       0,       0,       0,       0,       0,    0.75,       0,       0 }, // 12
      {       0,       0,       0,       0,   0.375,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.75,       0 }, // 13
      {       0,       0,       0,  -0.125,       0,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,     0.5,    0.25 }  // 14
    },

    // Embedding matrix for child 6
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  0
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  1
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  4
      {       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0 }, //  5
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, //  6
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,   0.375,       0,       0,    0.25,    0.75,       0,   0.375,       0 }, //  7
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,   0.375,    0.25,       0,    0.75,       0,       0,   0.375 }, //  8
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 10
      {       0,       0,  -0.125,       0,       0,   0.375,       0,       0,       0,       0,       0,    0.75,       0,       0,       0 }, // 11
      {       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.25,     0.5,     0.5 }, // 12
      {       0,       0,       0,       0,  -0.125,   0.375,       0,       0,       0,       0,       0,       0,       0,    0.75,       0 }, // 13
      {       0,       0,       0,  -0.125,       0,   0.375,       0,       0,       0,       0,       0,       0,       0,       0,    0.75 }  // 14
    },

    // Embedding matrix for child 7
    {
      //       0        1        2        3        4        5        6        7        8        9       10       11       12       13       14
      {   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0 }, //  0
      {       0,   -0.25,   -0.25,       0,   -0.25,   -0.25,       0,     0.5,       0,       0,     0.5,     0.5,       0,     0.5,       0 }, //  1
      {   -0.25,       0,   -0.25,   -0.25,       0,   -0.25,       0,       0,     0.5,     0.5,       0,     0.5,       0,       0,     0.5 }, //  2
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0 }, //  3
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0 }, //  4
      {       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1 }, //  5
      { -0.1875,   -0.25, -0.1875, -0.1875,   -0.25, -0.1875,    0.25,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125 }, //  6
      { -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,   -0.25,   0.125,    0.25,    0.25,    0.25,    0.25,     0.5,   0.125,    0.25,    0.25 }, //  7
      {   -0.25, -0.1875, -0.1875,   -0.25, -0.1875, -0.1875,    0.25,   0.125,    0.25,     0.5,    0.25,    0.25,    0.25,   0.125,    0.25 }, //  8
      { -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0,       0 }, //  9
      {       0, -0.1875, -0.1875,       0, -0.1875, -0.1875,       0,    0.25,       0,       0,   0.375,   0.375,       0,    0.75,       0 }, // 10
      { -0.1875,       0, -0.1875, -0.1875,       0, -0.1875,       0,       0,    0.25,   0.375,       0,   0.375,       0,       0,    0.75 }, // 11
      {       0,       0,       0,  -0.125,       0,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,     0.5,    0.25 }, // 12
      {       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,       0,    0.25,     0.5,     0.5 }, // 13
      {       0,       0,       0,       0,  -0.125,  -0.125,       0,       0,       0,       0,       0,       0,     0.5,    0.25,     0.5 }  // 14
    }
  };

#endif

typedef Prism15Templ<Real> Prism15;

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM15_H
