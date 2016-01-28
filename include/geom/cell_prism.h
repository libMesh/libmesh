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



#ifndef LIBMESH_CELL_PRISM_H
#define LIBMESH_CELL_PRISM_H

// Local includes
#include "libmesh/cell.h"

namespace libMesh
{




/**
 * The \p Prism is an element in 3D with 5 sides.
 */
class Prism : public Cell
{
public:

  /**
   * Default prismatic element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Prism(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    Cell(nn, Prism::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
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
   * @returns 6.  All prism-derivatives are guaranteed to have at
   * least 6 nodes.
   */
  virtual unsigned int n_nodes() const libmesh_override { return 6; }

  /**
   * @returns 5
   */
  virtual unsigned int n_sides() const libmesh_override { return 5; }

  /**
   * @returns 6.  All prisms have 6 vertices.
   */
  virtual unsigned int n_vertices() const libmesh_override { return 6; }

  /**
   * @returns 9.  All prisms have 9 edges.
   */
  virtual unsigned int n_edges() const libmesh_override { return 9; }

  /**
   * @returns 5.  All prisms have 5 faces.
   */
  virtual unsigned int n_faces() const libmesh_override { return 5; }

  /**
   * @returns 8
   */
  virtual unsigned int n_children() const libmesh_override { return 8; }

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
   * @returns a primitive triangle or quad for
   * face i.
   */
  virtual UniquePtr<Elem> side (const unsigned int i) const libmesh_override;


protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[6+(LIBMESH_DIM>3)];

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  for the first
   * 9 second-order nodes, this matrix is identical for
   * \p Prism15 and \p Prism18, therefore store it here.
   */
  static const unsigned short int _second_order_adjacent_vertices[9][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[18];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[18];

  /**
   * Master element node locations
   */
  static const Real _master_points[18][3];
};

} // namespace libMesh

#endif // LIBMESH_CELL_PRISM_H
