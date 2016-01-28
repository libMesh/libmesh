// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes
#include <cstddef>

namespace libMesh
{




/**
 * The \p Tet4 is an element in 3D composed of 4 nodes.
 * It is numbered like this:
 * \verbatim
 * TET4:
 *       3
 *       o
 *      /|\
 *     / | \
 *    /  |  \
 * 0 o...|...o 2
 *    \  |  /
 *     \ | /
 *      \|/
 *       o
 *       1
 * \endverbatim
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

  /**
   * @returns \p TET4
   */
  virtual ElemType type () const libmesh_override { return TET4; }

  /**
   * @returns 4
   */
  virtual unsigned int n_nodes() const libmesh_override { return 4; }

  /**
   * @returns 1
   */
  virtual unsigned int n_sub_elem() const libmesh_override { return 1; }

  /**
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const libmesh_override;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const libmesh_override;

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the element map is definitely affine within
   * numerical tolerances
   */
  virtual bool has_affine_map () const libmesh_override { return true; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const libmesh_override { return true; }

  /**
   * @returns FIRST
   */
  virtual Order default_order() const libmesh_override { return FIRST; }

  /**
   * Builds a \p TRI3 built coincident with face i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_side (const unsigned int i,
                                      bool proxy) const libmesh_override;

  /**
   * Builds a \p EDGE2 built coincident with face i.
   * The \p UniquePtr<Elem> handles the memory aspect.
   */
  virtual UniquePtr<Elem> build_edge (const unsigned int i) const libmesh_override;

  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const libmesh_override;

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
  virtual Real volume () const libmesh_override;

  /**
   * Returns the min and max *dihedral* angles for the tetrahedron.
   * Note there are 6 dihedral angles (angles between the planar
   * faces) for the Tet4.  Dihedral angles near 180 deg. are generally
   * bad for interplation.  Small dihedral angles are not necessarily
   * bad for interplation, but they can effect the stiffness matrix
   * condition number.
   */
  std::pair<Real, Real> min_and_max_angle() const;

  /**
   * Don't hide Tet::key(side) defined in the base class.
   */
  using Tet::key;

  /**
   * @returns an id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const libmesh_override;

protected:

  /**
   * Data for links to nodes
   */
  Node * _nodelinks_data[4];



#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual float embedding_matrix (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k) const libmesh_override;

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
