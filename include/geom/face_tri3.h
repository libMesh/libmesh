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



#ifndef LIBMESH_FACE_TRI3_H
#define LIBMESH_FACE_TRI3_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_tri.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

template <typename> class Edge2Templ;

/**
 * The \p Tri3 is an element in 2D composed of 3 nodes.
 * It is numbered like this:
 * \verbatim
 *   TRI3:
 *    2
 *    o
 *    |\            eta
 *    | \            ^
 *    |  \           |
 *    |   \          |
 *    |    \         o---> xi
 *    o-----o
 *    0      1
 * \endverbatim
 * (xi, eta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 2D triangular element with 3 nodes.
 */
template <typename RealType = Real>
class Tri3Templ : public TriTempl<RealType>
{
public:
  typedef Tri3Templ<RealType> Tri3;
  typedef TriTempl<RealType> Tri;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;
  typedef Edge2Templ<RealType> Edge2;
  typedef BoundingBoxTempl<RealType> BoundingBox;

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tri3Templ (Elem * p=nullptr) :
    Tri(Tri3::n_nodes(), p, _nodelinks_data) {}

  Tri3Templ (Tri3 &&) = delete;
  Tri3Templ (const Tri3 &) = delete;
  Tri3 & operator= (const Tri3 &) = delete;
  Tri3 & operator= (Tri3 &&) = delete;
  virtual ~Tri3Templ() = default;

  /**
   * \returns \p TRI3.
   */
  virtual ElemType type () const override { return TRI3; }

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
  virtual bool has_affine_map () const override { return true; }

  /**
   * \returns \p true if the Lagrange shape functions on this element
   * are linear.
   */
  virtual bool is_linear () const override { return true; }

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
   * Geometric constants for Tri3.
   */
  static const int num_nodes = 3;
  static const int num_sides = 3;
  static const int num_children = 4;
  static const int nodes_per_side = 2;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[num_sides][nodes_per_side];

  /**
   * An optimized method for computing the area of a 3-node triangle.
   */
  virtual RealType volume () const override;

  /**
   * \returns The minimum and maximum angles for the triangle
   * (in radians) in a std::pair.  The first entry in the pair
   * is the minimum angle, the second entry is the max angle.
   */
  std::pair<Real, Real> min_and_max_angle() const;

  /**
   * Specialization for tri3 elements. These elements are guaranteed to be planar
   * so a simple linear geometric test can be used.
   */
  virtual bool contains_point (const Point & p, Real tol) const override;

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
// Tri3 class static member initializations

template <typename RealType>
const unsigned int Tri3Templ<RealType>::side_nodes_map[Tri3::num_sides][Tri3::nodes_per_side] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 0}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

template <typename RealType>
const float Tri3Templ<RealType>::_embedding_matrix[Tri3::num_children][Tri3::num_nodes][Tri3::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5}  // 2
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 1.0, 0.0}, // 1
      {0.0, 0.5, 0.5}  // 2
    },

    // embedding matrix for child 2
    {
      // 0    1    2
      {0.5, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.0, 0.0, 1.0}  // 2
    },

    // embedding matrix for child 3
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.5, 0.0, 0.5}  // 2
    }
  };

#endif


typedef Tri3Templ<Real> Tri3;

} // namespace libMesh

#endif // LIBMESH_FACE_TRI3_H
