// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_TET4_H
#define LIBMESH_CELL_TET4_H

// Local includes
#include "libmesh/cell_tet.h"

namespace libMesh
{

/**
 * The \p Tet4 is an element in 3D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 *   TET4:
 *         3
 *         o                 zeta
 *        /|\                 ^
 *       / | \                |
 *      /  |  \               |
 *   0 o...|...o 2            o---> eta
 *      \  |  /                \
 *       \ | /                  \
 *        \|/                    xi (out of page)
 *         o
 *         1
 * \endverbatim
 * (xi, eta, zeta) are the reference element coordinates associated with
 * the given numbering.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief A 3D tetrahedral element with 4 nodes.
 */
class Tet4 libmesh_final : public Tet
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  Tet4 (Elem * p=libmesh_nullptr) :
    Tet(Tet4::n_nodes(), p, _nodelinks_data)
  {}

  Tet4 (Tet4 &&) = delete;
  Tet4 (const Tet4 &) = delete;
  Tet4 & operator= (const Tet4 &) = delete;
  Tet4 & operator= (Tet4 &&) = delete;
  virtual ~Tet4() = default;

  /**
   * \returns \p TET4.
   */
  virtual ElemType type () const override { return TET4; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_nodes() const override { return 4; }

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

  /**
   * Builds a \p TRI3 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy) override;

  /**
   * Builds a \p EDGE2 built coincident with face i.
   * The \p std::unique_ptr<Elem> handles the memory aspect.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
   * element node numbers.
   */
  static const unsigned int side_nodes_map[4][3];

  /**
   * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ edge to
   * element node numbers.
   */
  static const unsigned int edge_nodes_map[6][2];

  /**
   * An optimized method for computing the area of a
   * 4-node tetrahedron.
   */
  virtual Real volume () const override;

  /**
   * \returns The min and max *dihedral* angles for the tetrahedron.
   *
   * \note There are 6 dihedral angles (angles between the planar
   * faces) for the Tet4.  Dihedral angles near 180 deg. are generally
   * bad for interpolation.  Small dihedral angles are not necessarily
   * bad for interpolation, but they can affect the stiffness matrix
   * condition number.
   */
  std::pair<Real, Real> min_and_max_angle() const;

  /**
   * Don't hide Tet::key(side) defined in the base class.
   */
  using Tet::key;

  /**
   * \returns An id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const override;

  /**
   * Uses simple geometric tests to determine if the point p is inside
   * the tetrahedron.
   */
  virtual bool contains_point (const Point & p, Real tol) const override;

protected:

  /**
   * Data for links to nodes.
   */
  Node * _nodelinks_data[4];



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
  static const float _embedding_matrix[8][4][4];

  LIBMESH_ENABLE_TOPOLOGY_CACHES;

#endif

};

} // namespace libMesh


#endif // LIBMESH_CELL_TET4_H
