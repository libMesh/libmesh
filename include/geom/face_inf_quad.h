// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FACE_INF_QUAD_H
#define LIBMESH_FACE_INF_QUAD_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/elem.h"

namespace libMesh
{

/**
 * The \p InfQuad is an abstract element type that lives in two
 * dimensions.  Here, an infinite face is always a quadrilateral,
 * so this class is directly derived from \p Elem, without an
 * intermediate \p InfFace class.
 *
 * It looks like this:
 * \verbatim
 *                                   closer to infinity
 *          |           |
 *          |           |
 *   side 2 |           | side 1
 *          |           |
 *          |           |
 *           -----------             base side
 *
 *             side 0
 * \endverbatim
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief The base class for all 2D infinite quadrilateral element types.
 */
class InfQuad : public Elem
{
public:

  /**
   * Constructor.  Derived classes implement 'true' elements.
   */
  explicit
  InfQuad (const unsigned int nn,
           Elem * p,
           Node ** nodelinkdata) :
    Elem(nn, InfQuad::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 2)
      this->set_interior_parent(libmesh_nullptr);
  }

  /**
   * \returns The \p Point associated with local \p Node \p i,
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
   * \returns 2, the dimensionality of the object.
   */
  virtual unsigned int dim() const libmesh_override { return 2; }

  /**
   * \returns 3.  Infinite faces have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  virtual unsigned int n_sides() const libmesh_override { return 3; }

  /**
   * \returns 4.  All infinite quads (in our setting) have 4 vertices.
   */
  virtual unsigned int n_vertices() const libmesh_override { return 4; }

  /**
   * \returns 3.  All infinite quads have 1 edge in the
   * base, and 2 perpendicular to the base.
   */
  virtual unsigned int n_edges() const libmesh_override { return 3; }

  /**
   * \returns 0.  All 2D elements have no faces, just edges.
   */
  virtual unsigned int n_faces() const libmesh_override { return 0; }

  /**
   * \returns 2.
   */
  virtual unsigned int n_children() const libmesh_override { return 2; }

  /**
   * \returns \p true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int i) const
    libmesh_override { return (i > 2 && i < 4); }

  /**
   * \returns \p true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const libmesh_override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const libmesh_override;

  /**
   * \returns \p InfQuad4::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int which_node_am_i(unsigned int side,
                                       unsigned int side_node) const libmesh_override;

  /**
   * \returns A primitive (2-noded) edge or infedge for edge \p i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) libmesh_override;

  /**
   * build_edge_ptr() and build_side_ptr() are identical in 2D.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) libmesh_override
  { return build_side_ptr(i); }

  /**
   * is_edge_on_side is trivial in 2D.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const libmesh_override
  { return (e == s); }

  /**
   * \returns A quantitative assessment of element quality based on
   * the quality metric \p q specified by the user.
   */
  virtual Real quality (const ElemQuality q) const libmesh_override;

  /**
   * \returns The suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const libmesh_override;

  /**
   * \returns \p true.  All classes derived from \p InfQuad
   * are infinite elements.
   */
  virtual bool infinite () const libmesh_override { return true; }

  /**
   * \returns The origin of this infinite element.
   */
  virtual Point origin () const libmesh_override
  {
    return ( this->point(0)*2 - this->point(this->n_vertices()/2) );
  }


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


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_FACE_INF_QUAD_H
