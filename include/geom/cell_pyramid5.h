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



#ifndef LIBMESH_CELL_PYRAMID5_H
#define LIBMESH_CELL_PYRAMID5_H

// Local includes
#include "libmesh/cell_pyramid.h"

namespace libMesh
{
template <typename> class Quad4Templ;
template <typename> class Tri3Templ;
template <typename> class Edge2Templ;

/**
 * The \p Pyramid5 is an element in 3D composed of 5 nodes.
 * It is numbered with a counter-clockwise base like this:
 * \verbatim
 *   PYRAMID5:
 *             o 4
 *           //|\
 *          // | \          zeta
 *         //  |  \          ^   eta (into page)
 *      3 o/...|...o 2       | /
 *       ./    |  /          |/
 *      ./     | /           o---> xi
 *     ./      |/
 *    o--------o
 *    0        1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D pyramid element with 5 nodes.
 */
template <typename RealType = Real>
class Pyramid5Templ final : public PyramidTempl<RealType>
{
public:
  typedef PyramidTempl<RealType> Pyramid;
  typedef Pyramid5Templ<RealType> Pyramid5;
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
  Pyramid5Templ (Elem * p=nullptr) :
    Pyramid(Pyramid5::n_nodes(), p, _nodelinks_data)
  {}

  Pyramid5Templ (Pyramid5 &&) = delete;
  Pyramid5Templ (const Pyramid5 &) = delete;
  Pyramid5 & operator= (const Pyramid5 &) = delete;
  Pyramid5 & operator= (Pyramid5 &&) = delete;
  virtual ~Pyramid5Templ() = default;

  /**
   * \returns \p PRYAMID.
   */
  virtual ElemType type () const override { return PYRAMID5; }

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
   * Builds a \p EDGE2 built coincident with edge i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Geometric constants for Pyramid5.
   */
  static const int num_nodes = 5;
  static const int num_sides = 5;
  static const int num_edges = 8;
  static const int num_children = 0; // not implemented
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
   * Specialization for computing the volume of a pyramid.
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
  virtual float embedding_matrix (const unsigned int,
                                  const unsigned int,
                                  const unsigned int) const override
  { libmesh_not_implemented(); return 0.; }

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif // LIBMESH_ENABLE_AMR

};

// ------------------------------------------------------------
// Pyramid5 class static member initializations

template <typename RealType>
const unsigned int Pyramid5Templ<RealType>::side_nodes_map[Pyramid5::num_sides][Pyramid5::nodes_per_side] =
  {
    {0, 1, 4, 99}, // Side 0
    {1, 2, 4, 99}, // Side 1
    {2, 3, 4, 99}, // Side 2
    {3, 0, 4, 99}, // Side 3
    {0, 3, 2,  1}  // Side 4
  };

template <typename RealType>
const unsigned int Pyramid5Templ<RealType>::edge_nodes_map[Pyramid5::num_edges][Pyramid5::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 4}, // Edge 5
    {2, 4}, // Edge 6
    {3, 4}  // Edge 7
  };


typedef Pyramid5Templ<Real> Pyramid5;

} // namespace libMesh


#endif // LIBMESH_CELL_PYRAMID5_H
