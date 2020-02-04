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



#ifndef LIBMESH_CELL_HEX8_H
#define LIBMESH_CELL_HEX8_H

// Local includes
#include "libmesh/cell_hex.h"

namespace libMesh
{

template <typename>
class Quad4Templ;
template <typename>
class Edge2Templ;

/**
 * The \p Hex8 is an element in 3D composed of 8 nodes.
 * It is numbered like this:
 * \verbatim
 *   HEX8: 7        6
 *         o--------z
 *        /:       /|         zeta
 *       / :      / |          ^   eta (into page)
 *    4 /  :   5 /  |          | /
 *     o--------o   |          |/
 *     |   o....|...o 2        o---> xi
 *     |  .3    |  /
 *     | .      | /
 *     |.       |/
 *     o--------o
 *     0        1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D hexahedral element with 8 nodes.
 */
template <typename RealType = Real>
class Hex8Templ final : public HexTempl<RealType>
{
public:
  typedef Hex8Templ<RealType> Hex8;
  typedef Quad4Templ<RealType> Quad4;
  typedef Edge2Templ<RealType> Edge2;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Hex8Templ (Elem * p=nullptr) :
    Hex(Hex8::n_nodes(), p, _nodelinks_data)
  {}

  Hex8Templ (Hex8 &&) = delete;
  Hex8Templ (const Hex8 &) = delete;
  Hex8 & operator= (const Hex8 &) = delete;
  Hex8 & operator= (Hex8 &&) = delete;
  virtual ~Hex8Templ() = default;

  /**
   * \returns \p HEX8.
   */
  virtual ElemType type () const override { return HEX8; }

  /**
   * \returns 8.
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
   * \returns FIRST.
   */
  virtual Order default_order() const override;

  /**
   * Builds a QUAD4 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p QUAD4 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a EDGE2 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Geometric constants for Hex8.
   */
  static const int num_nodes = 8;
  static const int num_sides = 6;
  static const int num_edges = 12;
  static const int num_children = 8;
  static const int nodes_per_side = 4;
  static const int nodes_per_edge = 2;

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
   * A specialization for computing the area of a hexahedron
   * with flat sides.
   */
  virtual RealType volume () const override;

  /**
   * Builds a bounding box out of the nodal positions
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

// ------------------------------------------------------------
// Hex8 class static member initializations

template <typename RealType>
const unsigned int Hex8Templ<RealType>::
side_nodes_map[Hex8Templ<RealType>::num_sides][Hex8Templ<RealType>::nodes_per_side] =
  {
    {0, 3, 2, 1}, // Side 0
    {0, 1, 5, 4}, // Side 1
    {1, 2, 6, 5}, // Side 2
    {2, 3, 7, 6}, // Side 3
    {3, 0, 4, 7}, // Side 4
    {4, 5, 6, 7}  // Side 5
  };

template <typename RealType>
const unsigned int Hex8Templ<RealType>::
edge_nodes_map[Hex8Templ<RealType>::num_edges][Hex8Templ<RealType>::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 5}, // Edge 5
    {2, 6}, // Edge 6
    {3, 7}, // Edge 7
    {4, 5}, // Edge 8
    {5, 6}, // Edge 9
    {6, 7}, // Edge 10
    {4, 7}  // Edge 11
  };

#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Hex8Templ<RealType>::
_embedding_matrix[Hex8Templ<RealType>::num_children]
                 [Hex8Templ<RealType>::num_nodes]
                 [Hex8Templ<RealType>::num_nodes] =
  {
    // The 8 children of the Hex-type elements can be thought of as being
    // associated with the 8 vertices of the Hex.  Some of the children are
    // numbered the same as their corresponding vertex, while some are
    // not.  The children which are numbered differently have been marked
    // with ** in the comments below.

    // embedding matrix for child 0 (child 0 is associated with vertex 0)
    {
      //  0     1     2     3     4     5     6     7
      { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
      { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 4
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 5
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 6
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}  // 7
    },

    // embedding matrix for child 1 (child 1 is associated with vertex 1)
    {
      //  0     1     2     3     4     5     6     7
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 3
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 4
      { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 5
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 6
      {.125, .125, .125, .125, .125, .125, .125, .125}  // 7
    },

    // embedding matrix for child 2 (child 2 is associated with vertex 3**)
    {
      //  0      1    2     3     4     5     6     7
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 3
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 4
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 5
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 6
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}  // 7
    },

    // embedding matrix for child 3 (child 3 is associated with vertex 2**)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  .25,  .25,  .25,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 2
      { 0.0,  0.0,  0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 3
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 4
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 5
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 6
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}  // 7
    },

    // embedding matrix for child 4 (child 4 is associated with vertex 4)
    {
      //  0      1    2     3     4     5     6     7
      { 0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 1
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 2
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 3
      { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 7
    },

    // embedding matrix for child 5 (child 5 is associated with vertex 5)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0}, // 1
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 2
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 6
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}  // 7
    },

    // embedding matrix for child 6 (child 6 is associated with vertex 7**)
    {
      //  0      1    2     3     4     5     6     7
      { .25,  0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25}, // 0
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 1
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 7
    },

    // embedding matrix for child 7 (child 7 is associated with vertex 6**)
    {
      //  0      1    2     3     4     5     6     7
      {.125, .125, .125, .125, .125, .125, .125, .125}, // 0
      { 0.0,  .25,  .25,  0.0,  0.0,  .25,  .25,  0.0}, // 1
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.0}, // 2
      { 0.0,  0.0,  .25,  .25,  0.0,  0.0,  .25,  .25}, // 3
      { 0.0,  0.0,  0.0,  0.0,  .25,  .25,  .25,  .25}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 5
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 6
      { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 7
    }
  };

#endif


typedef Hex8Templ<Real> Hex8;

} // namespace libMesh

#endif // LIBMESH_CELL_HEX8_H
