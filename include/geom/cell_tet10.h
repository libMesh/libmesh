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



#ifndef LIBMESH_CELL_TET10_H
#define LIBMESH_CELL_TET10_H

// Local includes
#include "libmesh/cell_tet.h"

namespace libMesh
{
template <typename> class Tri6Templ;
template <typename> class Edge3Templ;

/**
 * The \p Tet10 is an element in 3D composed of 10 nodes.
 * It is numbered like this:
 * \verbatim
 *               3
 *   TET10:      o
 *              /|\
 *             / | \
 *         7  /  |  \9
 *           o   |   o              zeta
 *          /    |8   \              ^
 *         /     o     \             |
 *        /    6 |      \            |
 *     0 o.....o.|.......o 2         o---> eta
 *        \      |      /             \
 *         \     |     /               \
 *          \    |    /                 xi (out of page)
 *         4 o   |   o 5
 *            \  |  /
 *             \ | /
 *              \|/
 *               o
 *               1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D tetrahedral element with 10 nodes.
 */
template <typename RealType = Real>
class Tet10Templ final : public TetTempl<RealType>
{
public:
  typedef Tet10Templ<RealType> Tet10;
  typedef TetTempl<RealType> Tet;
  typedef CellTempl<RealType> Cell;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;
  typedef Tri6Templ<RealType> Tri6;
  typedef Edge3Templ<RealType> Edge3;
  using typename Tet::Diagonal;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tet10Templ (Elem * p=nullptr) :
    Tet(Tet10::n_nodes(), p, _nodelinks_data)
  {}

  Tet10Templ (Tet10 &&) = delete;
  Tet10Templ (const Tet10 &) = delete;
  Tet10 & operator= (const Tet10 &) = delete;
  Tet10 & operator= (Tet10 &&) = delete;
  virtual ~Tet10Templ() = default;

  /**
   * \returns \p TET10.
   */
  virtual ElemType type () const override { return TET10; }

  /**
   * \returns 10.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_sub_elem() const override { return 8; }

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
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

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
   * \returns \p Tet10::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * Builds a \p TRI6 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a TRI6 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE3 built coincident with edge i.
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
   * \note \p n is counted as depicted above, \f$ 4 \le n < 10 \f$.
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
   * Geometric constants for Tet10.
   */
  static const int num_nodes = 10;
  static const int num_sides = 4;
  static const int num_edges = 6;
  static const int num_children = 8;
  static const int nodes_per_side = 6;
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
   * A specialization for computing the volume of a Tet10.
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
                                  const unsigned int k) const override;

  /**
   * Matrix that computes new nodal locations/solution values
   * from current nodes/solution.
   */
  static const float _embedding_matrix[num_children][num_nodes][num_nodes];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[6][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[10];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[10];
};


// ------------------------------------------------------------
// Tet10 class static member initializations

template <typename RealType>
const int Tet10Templ<RealType>::nodes_per_side;

template <typename RealType>
const unsigned int Tet10Templ<RealType>::side_nodes_map[Tet10::num_sides][Tet10::nodes_per_side] =
  {
    {0, 2, 1, 6, 5, 4}, // Side 0
    {0, 1, 3, 4, 8, 7}, // Side 1
    {1, 2, 3, 5, 9, 8}, // Side 2
    {2, 0, 3, 6, 7, 9}  // Side 3
  };

template <typename RealType>
const unsigned int Tet10Templ<RealType>::edge_nodes_map[Tet10::num_edges][Tet10::nodes_per_edge] =
  {
    {0, 1, 4}, // Edge 0
    {1, 2, 5}, // Edge 1
    {0, 2, 6}, // Edge 2
    {0, 3, 7}, // Edge 3
    {1, 3, 8}, // Edge 4
    {2, 3, 9}  // Edge 5
  };

#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Tet10Templ<RealType>::_embedding_matrix[Tet10::num_children][Tet10::num_nodes][Tet10::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      { 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 5
      { 0.375,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
      { 0.375,    0.,    0.,-0.125,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
      {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
    },

    // embedding matrix for child 1
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.,    0.}, // 4
      {    0., 0.375,-0.125,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
      {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
      {    0., 0.375,    0.,-0.125,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}  // 9
    },

    // embedding matrix for child 2
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
      {    0.,-0.125, 0.375,    0.,    0.,  0.75,    0.,    0.,    0.,    0.}, // 5
      {-0.125,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 7
      {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 8
      {    0.,    0., 0.375,-0.125,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
    },

    // embedding matrix for child 3
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 4
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
      {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}, // 6
      {-0.125,    0.,    0., 0.375,    0.,    0.,    0.,  0.75,    0.,    0.}, // 7
      {    0.,-0.125,    0., 0.375,    0.,    0.,    0.,    0.,  0.75,    0.}, // 8
      {    0.,    0.,-0.125, 0.375,    0.,    0.,    0.,    0.,    0.,  0.75}  // 9
    },

    // embedding matrix for child 4
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 4
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 5
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
      {    0.,-0.125,    0.,-0.125,   0.5,    0.,    0.,   0.5,  0.25,    0.}, // 7
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}  // 9
    },

    // embedding matrix for child 5
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125,    0.,-0.125,    0.,   0.5,   0.5,  0.25,    0.,    0.,    0.}, // 4
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 5
      {    0.,-0.125,-0.125,    0.,   0.5,  0.25,   0.5,    0.,    0.,    0.}, // 6
      {-0.125,    0.,    0.,-0.125,   0.5,    0.,    0.,  0.25,   0.5,    0.}, // 7
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}  // 9
    },

    // embedding matrix for child 6
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 3
      {-0.125,-0.125,    0.,    0.,  0.25,   0.5,   0.5,    0.,    0.,    0.}, // 4
      {    0.,-0.125,    0.,-0.125,    0.,   0.5,    0.,    0.,  0.25,   0.5}, // 5
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 7
      {    0.,    0.,-0.125,-0.125,    0.,   0.5,    0.,    0.,   0.5,  0.25}, // 8
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}  // 9
    },

    // embedding matrix for child 7
    {
      //    0      1      2      3      4      5      6      7      8      9
      {    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.}, // 0
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.}, // 1
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.}, // 2
      {    0.,    0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.}, // 3
      {-0.125,-0.125,-0.125,-0.125,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25}, // 4
      {    0.,-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5}, // 5
      {-0.125,    0.,    0.,-0.125,    0.,    0.,   0.5,  0.25,    0.,   0.5}, // 6
      {    0.,    0.,-0.125,-0.125,    0.,    0.,   0.5,   0.5,    0.,  0.25}, // 7
      {-0.125,-0.125,    0.,    0.,  0.25,    0.,    0.,   0.5,   0.5,    0.}, // 8
      {-0.125,    0.,-0.125,    0.,    0.,    0.,  0.25,   0.5,    0.,   0.5}  // 9
    }
  };

#endif // LIBMESH_ENABLE_AMR

template <typename RealType>
const unsigned short int Tet10Templ<RealType>::_second_order_vertex_child_number[10] =
  {
    99,99,99,99, // Vertices
    0,1,0,0,1,2  // Edges
  };



template <typename RealType>
const unsigned short int Tet10Templ<RealType>::_second_order_vertex_child_index[10] =
  {
    99,99,99,99, // Vertices
    1,2,2,3,3,3  // Edges
  };

template <typename RealType>
const unsigned short int Tet10Templ<RealType>::_second_order_adjacent_vertices[6][2] =
  {
    {0, 1}, // vertices adjacent to node 4
    {1, 2}, // vertices adjacent to node 5
    {0, 2}, // vertices adjacent to node 6
    {0, 3}, // vertices adjacent to node 7
    {1, 3}, // vertices adjacent to node 8
    {2, 3}  // vertices adjacent to node 9
  };

typedef Tet10Templ<Real> Tet10;

} // namespace libMesh


#endif // LIBMESH_CELL_TET10_H
