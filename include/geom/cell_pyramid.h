// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_PYRAMID_H
#define LIBMESH_CELL_PYRAMID_H

// Local includes
#include "libmesh/cell.h"

namespace libMesh
{




/**
 * The \p Pyramid is an element in 3D with 5 sides.
 */
class Pyramid : public Cell
{
public:

  /**
   * Default pyramid, one quad face, four triangular faces,
   * takes number of nodes and parent.
   * Derived classes implement 'true' elements.
   */
  Pyramid(const unsigned int nn, Elem* p, Node** nodelinkdata) :
    Cell(nn, Pyramid::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
      this->set_interior_parent(NULL);
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
   * @returns 5.  All pyramid-derivatives are guaranteed to have at
   * least 5 nodes.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 5; }

  /**
   * @returns 5
   */
  virtual unsigned int n_sides() const libmesh_override { return 5; }

  /**
   * @returns 5.  All pyramids have 5 vertices.
   */
  virtual unsigned int n_vertices() const libmesh_override { return 5; }

  /**
   * @returns 8.  All pyramids have 8 edges.
   */
  virtual unsigned int n_edges() const libmesh_override { return 8; }

  /**
   * @returns 5.  All pyramids have 5 faces.
   */
  virtual unsigned int n_faces() const libmesh_override { return 5; }

  /**
   * @returns 10
   */
  virtual unsigned int n_children() const libmesh_override { return 10; }

  /*
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const libmesh_override;

  /*
   * @returns true iff the specified edge is on the specified side
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const libmesh_override;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const libmesh_override;

  /**
   * @returns a primitive triangle or quad for
   * face i.
   */
  virtual UniquePtr<Elem> side (const unsigned int i) const libmesh_override;


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem* _elemlinks_data[6+(LIBMESH_DIM>3)];

  /**
   * Master element node locations
   */
  static const Real _master_points[14][3];

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  unsigned int side_children_matrix (const unsigned int,
                                     const unsigned int) const
  { libmesh_not_implemented(); return 0; }

#endif

};

} // namespace libMesh

#endif // LIBMESH_CELL_PYRAMID_H
