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



#ifndef LIBMESH_FACE_QUAD4_H
#define LIBMESH_FACE_QUAD4_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_quad.h"

namespace libMesh
{

template <typename>
class Edge2Templ;
template <typename>
class BoundingBoxTempl;

/**
 * The \p QUAD4 is an element in 2D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *          3           2
 *   QUAD4: o-----------o
 *          |           |           eta
 *          |           |            ^
 *          |           |            |
 *          |           |            |
 *          |           |            o---> xi
 *          o-----------o
 *          0           1
 * \endverbatim
 * (xi, eta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D quadrilateral element with 4 nodes.
 */
template <typename RealType>
class Quad4Templ : public QuadTempl<RealType>
{
public:
  typedef Quad4Templ<RealType> Quad4;
  typedef Edge2Templ<RealType> Edge2;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Quad4Templ (Elem * p=nullptr) :
    Quad(Quad4::n_nodes(), p, _nodelinks_data) {}

  Quad4Templ (Quad4 &&) = delete;
  Quad4Templ (const Quad4 &) = delete;
  Quad4 & operator= (const Quad4 &) = delete;
  Quad4 & operator= (Quad4 &&) = delete;
  virtual ~Quad4Templ() = default;

  /**
   * \returns \p QUAD4.
   */
  virtual ElemType type () const override { return QUAD4; }

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
   * \returns FIRST.
   */
  virtual Order default_order() const override;

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
   * Geometric constants for Quad4.
   */
  static const int num_nodes = 4;
  static const int num_sides = 4;
  static const int num_children = 4;
  static const int nodes_per_side = 2;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * An optimized method for computing the area of a
   * 4-node quad with straight sides, but not necessarily a
   * parallelogram.
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
// Quad4 class static member initialization

template <typename RealType>
const unsigned int Quad4Templ<RealType>::
side_nodes_map[Quad4Templ<RealType>::num_sides][Quad4Templ<RealType>::nodes_per_side] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 3}, // Side 2
    {3, 0}  // Side 3
  };


#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Quad4Templ<RealType>::
_embedding_matrix[Quad4Templ<RealType>::num_children]
                 [Quad4Templ<RealType>::num_nodes]
                 [Quad4Templ<RealType>::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2    3
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {.25, .25, .25, .25}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 1
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {.25, .25, .25, .25}  // 3
    },

    // embedding matrix for child 2
    {
      // 0    1    2    3
      {0.5, 0.0, 0.0, 0.5}, // 0
      {.25, .25, .25, .25}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    },

    // embedding matrix for child 3
    {
      // 0    1    2    3
      {.25, .25, .25, .25}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    }
  };

#endif

typedef Quad4Templ<Real> Quad4;

} // namespace libMesh

#endif // LIBMESH_FACE_QUAD4_H
