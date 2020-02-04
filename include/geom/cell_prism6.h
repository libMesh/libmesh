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



#ifndef LIBMESH_CELL_PRISM6_H
#define LIBMESH_CELL_PRISM6_H

// Local includes
#include "libmesh/cell_prism.h"

namespace libMesh
{
template <typename> class Quad4Templ;
template <typename> class Tri3Templ;
template <typename> class Edge2Templ;

/**
 * The \p Prism6 is an element in 3D composed of 6 nodes.
 * It is numbered like this:
 * \verbatim
 *   PRISM6:
 *           5
 *           o
 *          /:\
 *         / : \             zeta
 *        /  o  \             ^   eta (into page)
 *     3 o-------o 4          | /
 *       | . 2 . |            |/
 *       |.     .|            o---> xi
 *       o-------o
 *       0       1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D prismatic element with 6 nodes.
 */
template <typename RealType = Real>
class Prism6Templ final : public PrismTempl<RealType>
{
public:
  typedef PrismTempl<RealType> Prism;
  typedef Prism6Templ<RealType> Prism6;
  typedef Quad4Templ<RealType> Quad4;
  typedef Tri3Templ<RealType> Tri3;
  typedef Edge2Templ<RealType> Edge2;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Prism6Templ (Elem * p=nullptr) :
    Prism(Prism6::n_nodes(), p, _nodelinks_data)
  {}

  Prism6Templ (Prism6 &&) = delete;
  Prism6Templ (const Prism6 &) = delete;
  Prism6 & operator= (const Prism6 &) = delete;
  Prism6 & operator= (Prism6 &&) = delete;
  virtual ~Prism6Templ() = default;

  /**
   * \returns \p PRISM6.
   */
  virtual ElemType type () const override { return PRISM6; }

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
   * Builds a \p QUAD4 or \p TRI3 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=true) override;

  /**
   * Rebuilds a \p QUAD4 or \p TRI3 built coincident with face i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  /**
   * Builds a \p EDGE2 or \p INFEDGE2 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Geometric constants for Prism6.
   */
  static const int num_nodes = 6;
  static const int num_sides = 5;
  static const int num_edges = 9;
  static const int num_children = 8;
  static const int nodes_per_side = 4;
  static const int nodes_per_edge = 2;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * This maps the child elements with the associated side of the parent element
   */
  static const unsigned int side_elems_map[num_sides][nodes_per_side];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[num_edges][nodes_per_edge];

  /**
   * Specialized function for computing the element volume.
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
// Prism6 class static member initializations
template <typename RealType>
const unsigned int Prism6Templ<RealType>::side_nodes_map[Prism6::num_sides][Prism6::nodes_per_side] =
  {
    {0, 2, 1, 99}, // Side 0
    {0, 1, 4,  3}, // Side 1
    {1, 2, 5,  4}, // Side 2
    {2, 0, 3,  5}, // Side 3
    {3, 4, 5, 99}  // Side 4
  };

template <typename RealType>
const unsigned int Prism6Templ<RealType>::side_elems_map[Prism6::num_sides][Prism6::nodes_per_side] =
  {
    {0, 1, 2, 3}, // Side 0
    {0, 1, 4, 5}, // Side 1
    {1, 2, 5, 6}, // Side 2
    {0, 2, 4, 6}, // Side 3
    {4, 5, 6, 7}  // Side 4
  };

template <typename RealType>
const unsigned int Prism6Templ<RealType>::edge_nodes_map[Prism6::num_edges][Prism6::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {0, 2}, // Edge 2
    {0, 3}, // Edge 3
    {1, 4}, // Edge 4
    {2, 5}, // Edge 5
    {3, 4}, // Edge 6
    {4, 5}, // Edge 7
    {3, 5}  // Edge 8
  };


#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Prism6Templ<RealType>::_embedding_matrix[Prism6::num_children][Prism6::num_nodes][Prism6::num_nodes] =
  {
    // embedding matrix for child 0
    {
      //  0     1     2     3     4     5
      { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 3
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 1
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 4
      { 0.0,  .25,  .25,  0.0,  .25,  .25}  // 5
    },

    // embedding matrix for child 2
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 2
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 3
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 4
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 0
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 5
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 1
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 5
    },

    // embedding matrix for child 6
    {
      //  0     1     2     3     4     5
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 5
    },

    // embedding matrix for child 7
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    }
  };

#endif

typedef Prism6Templ<Real> Prism6;

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM6_H
