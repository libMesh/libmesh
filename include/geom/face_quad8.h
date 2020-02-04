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



#ifndef LIBMESH_FACE_QUAD8_H
#define LIBMESH_FACE_QUAD8_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_quad.h"

namespace libMesh
{
template <typename> class Edge3Templ;

/**
 * The \p QUAD8 is an element in 2D composed of 8 nodes.
 * It is numbered like this:
 * \verbatim
 *          3     6     2
 *   QUAD8: o-----o-----o
 *          |           |           eta
 *          |           |            ^
 *        7 o           o 5          |
 *          |           |            |
 *          |           |            o---> xi
 *          o-----o-----o
 *          0     4     1
 * \endverbatim
 * (xi, eta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D quadrilateral element with 8 nodes.
 */
template <typename RealType = Real>
class Quad8Templ : public QuadTempl<RealType>
{
public:
  typedef QuadTempl<RealType> Quad;
  typedef Quad8Templ<RealType> Quad8;
  typedef Edge3Templ<RealType> Edge3;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor. By default this element has no parent.
   */
  explicit
  Quad8Templ (Elem * p=nullptr) :
    Quad(Quad8::n_nodes(), p, _nodelinks_data) {}

  Quad8Templ (Quad8 &&) = delete;
  Quad8Templ (const Quad8 &) = delete;
  Quad8 & operator= (const Quad8 &) = delete;
  Quad8 & operator= (Quad8 &&) = delete;
  virtual ~Quad8Templ() = default;

  /**
   * \returns \p QUAD8.
   */
  virtual ElemType type () const override { return QUAD8; }

  /**
   * \returns 8.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 5.
   */
  virtual unsigned int n_sub_elem() const override { return 5; }

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
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

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
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   *
   * We reimplement this method here for the \p Quad8 since we can
   * use the center node of each edge to provide a perfect (unique)
   * key.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns \p Quad8::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds an EDGE3 coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  virtual void connectivity(const unsigned int sf,
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
   * \note \p n is counted as depicted above, \f$ 4 \le n < 8 \f$.
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
   * Geometric constants for Quad8.
   */
  static const int num_nodes = 8;
  static const int num_sides = 4;
  static const int num_children = 4;
  static const int nodes_per_side = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * An optimized method for approximating the area of a
   * QUAD8 using quadrature.
   */
  virtual RealType volume () const override;

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const override;

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

#endif // LIBMESH_ENABLE_AMR

};

template <typename RealType>
const int Quad8Templ<RealType>::nodes_per_side;

template <typename RealType>
const unsigned int Quad8Templ<RealType>::side_nodes_map[Quad8::num_sides][Quad8::nodes_per_side] =
  {
    {0, 1, 4}, // Side 0
    {1, 2, 5}, // Side 1
    {2, 3, 6}, // Side 2
    {3, 0, 7}  // Side 3
  };


#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Quad8Templ<RealType>::_embedding_matrix[Quad8::num_children][Quad8::num_nodes][Quad8::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //         0           1           2           3           4           5           6           7
      {    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 3
      {   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 4
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.750000,   0.375000,   0.250000,   0.375000 }, // 5
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.250000,   0.375000,   0.750000 }, // 6
      {   0.375000,    0.00000,    0.00000,  -0.125000,    0.00000,    0.00000,    0.00000,   0.750000 }  // 7
    },

    // embedding matrix for child 1
    {
      //         0           1           2           3           4           5           6           7
      {    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000 }, // 0
      {    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 2
      {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 3
      {  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000,    0.00000 }, // 4
      {    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 5
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.750000,   0.375000,   0.250000 }, // 6
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.750000,   0.375000,   0.250000,   0.375000 }  // 7
    },

    // embedding matrix for child 2
    {
      //         0           1           2           3           4           5           6           7
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000 }, // 0
      {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 1
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 3
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.250000,   0.375000,   0.750000 }, // 4
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.250000,   0.375000,   0.750000,   0.375000 }, // 5
      {    0.00000,    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 6
      {  -0.125000,    0.00000,    0.00000,   0.375000,    0.00000,    0.00000,    0.00000,   0.750000 }  // 7
    },

    // embedding matrix for child 3
    {
      //         0           1           2           3           4           5           6           7
      {  -0.250000,  -0.250000,  -0.250000,  -0.250000,   0.500000,   0.500000,   0.500000,   0.500000 }, // 0
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000,    0.00000 }, // 1
      {    0.00000,    0.00000,    1.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000 }, // 2
      {    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    0.00000,    1.00000,    0.00000 }, // 3
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.375000,   0.750000,   0.375000,   0.250000 }, // 4
      {    0.00000,  -0.125000,   0.375000,    0.00000,    0.00000,   0.750000,    0.00000,    0.00000 }, // 5
      {    0.00000,    0.00000,   0.375000,  -0.125000,    0.00000,    0.00000,   0.750000,    0.00000 }, // 6
      {  -0.187500,  -0.187500,  -0.187500,  -0.187500,   0.250000,   0.375000,   0.750000,   0.375000 }  // 7
    }
  };


#endif

typedef Quad8Templ<Real> Quad8;

} // namespace libMesh

#endif // LIBMESH_FACE_QUAD8_H
