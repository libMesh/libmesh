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



#ifndef LIBMESH_FACE_TRI6_H
#define LIBMESH_FACE_TRI6_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_tri.h"

namespace libMesh
{
template <typename> class Edge3Templ;

/**
 * The \p Tri6 is an element in 2D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *   TRI6:
 *    2
 *    o
 *    |\
 *    | \
 *    |  \            eta
 *    |   \            ^
 *    |    \           |
 *  5 o     o 4        |
 *    |      \         o---> xi
 *    |       \
 *    |        \
 *    o----o----o
 *    0    3    1
 * \endverbatim
 * (xi, eta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D triangular element with 6 nodes.
 */
template <typename RealType = Real>
class Tri6Templ : public TriTempl<RealType>
{
public:
  typedef Tri6Templ<RealType> Tri6;
  typedef TriTempl<RealType> Tri;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;
  typedef Edge3Templ<RealType> Edge3;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tri6Templ (Elem * p=nullptr) :
    Tri(Tri6::n_nodes(), p, _nodelinks_data) {}

  Tri6Templ (Tri6 &&) = delete;
  Tri6Templ (const Tri6 &) = delete;
  Tri6 & operator= (const Tri6 &) = delete;
  Tri6 & operator= (Tri6 &&) = delete;
  virtual ~Tri6Templ() = default;

  /**
   * \returns \p TRI6.
   */
  virtual ElemType type () const override { return TRI6; }

  /**
   * \returns 6.
   */
  virtual unsigned int n_nodes() const override { return num_nodes; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_sub_elem() const override { return 4; }

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
   * \returns \p Tri6::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const override;

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds an EDGE2 coincident with face i.
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
   * \returns The element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.
   *
   * \note \p n is counted as depicted above, \f$ 3 \le n < 6 \f$.
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
   * Geometric constants for Tri6.
   */
  static const int num_nodes = 6;
  static const int num_sides = 3;
  static const int num_children = 4;
  static const int nodes_per_side = 3;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * An optimized method for approximating the area of a
   * TRI6 using quadrature.
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

private:

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes
   */
  static const unsigned short int _second_order_adjacent_vertices[num_sides][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[num_nodes];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[num_nodes];
};

template <typename RealType>
const int Tri6Templ<RealType>::nodes_per_side;

template <typename RealType>
const unsigned int Tri6Templ<RealType>::side_nodes_map[Tri6::num_sides][Tri6::nodes_per_side] =
  {
    {0, 1, 3}, // Side 0
    {1, 2, 4}, // Side 1
    {2, 0, 5}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Tri6Templ<RealType>::_embedding_matrix[Tri6::num_children][Tri6::num_nodes][Tri6::num_nodes] =
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


template <typename RealType>
const unsigned short int Tri6Templ<RealType>::_second_order_vertex_child_number[Tri6::num_nodes] =
  {
    99,99,99, // Vertices
    0,1,0     // Edges
  };



template <typename RealType>
const unsigned short int Tri6Templ<RealType>::_second_order_vertex_child_index[Tri6::num_nodes] =
  {
    99,99,99, // Vertices
    1,2,2     // Edges
  };

template <typename RealType>
const unsigned short int Tri6Templ<RealType>::_second_order_adjacent_vertices[Tri6::num_sides][2] =
  {
    {0, 1}, // vertices adjacent to node 3
    {1, 2}, // vertices adjacent to node 4
    {0, 2}  // vertices adjacent to node 5
  };

typedef Tri6Templ<Real> Tri6;


} // namespace libMesh

#endif // LIBMESH_FACE_TRI6_H
