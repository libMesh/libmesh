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



#ifndef LIBMESH_FACE_TRI_H
#define LIBMESH_FACE_TRI_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face.h"

// C++ includes

namespace libMesh
{


// Forward declarations





/**
 * The \p Tri is an element in 2D composed of 3 sides.
 * It looks like this:
 * \verbatim
 *
 *          ^
 *         / \
 *        /   \
 *       /     \
 *      /       \
 *     /         \
 *     -----------
 *
 * \endverbatim
 */
class Tri : public Face
{
public:

  /**
   * Default triangular element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Tri (const unsigned int nn,
       Elem * p,
       Node ** nodelinkdata) :
    Face(nn, Tri::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 2)
      this->set_interior_parent(libmesh_nullptr);
  }

  /**
   * @returns the \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const libmesh_override
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * @returns 3.  All tri-derivatives are guaranteed to have at
   * least 3 nodes.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 3; }

  /**
   * @returns 3
   */
  virtual unsigned int n_sides() const libmesh_override { return 3; }

  /**
   * @returns 3.  All triangles have 3 vertices.
   */
  virtual unsigned int n_vertices() const libmesh_override { return 3; }

  /**
   * @returns 3.  All triangles have 3 edges.
   */
  virtual unsigned int n_edges() const libmesh_override { return 3; }

  /**
   * @returns 4
   */
  virtual unsigned int n_children() const libmesh_override { return 4; }

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const libmesh_override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const libmesh_override;

  /**
   * @returns an id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const libmesh_override;

  /**
   * @returns a primitive (2-noded) edge for
   * edge i.
   */
  virtual UniquePtr<Elem> side (const unsigned int i) const libmesh_override;

  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  virtual Real quality (const ElemQuality q) const libmesh_override;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const libmesh_override;


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[4+(LIBMESH_DIM>2)];

  /**
   * Master element node locations
   */
  static const Real _master_points[6][3];
};

} // namespace libMesh

#endif // LIBMESH_FACE_TRI_H
