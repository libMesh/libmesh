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



#ifndef LIBMESH_CELL_HEX_H
#define LIBMESH_CELL_HEX_H

// Local includes
#include "libmesh/cell.h"

namespace libMesh
{




/**
 * The \p Hex is an element in 3D with 6 sides.
 */
class Hex : public Cell
{
public:

  /**
   * Default brick element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  Hex(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    Cell(nn, Hex::n_sides(), p, _elemlinks_data, nodelinkdata)
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
   * @returns 6
   */
  virtual unsigned int n_sides() const libmesh_override { return 6; }

  /**
   * @returns 8.  All hexahedrals have 8 vertices.
   */
  virtual unsigned int n_vertices() const libmesh_override { return 8; }

  /**
   * @returns 12.  All hexahedrals have 12 edges.
   */
  virtual unsigned int n_edges() const libmesh_override { return 12; }

  /**
   * @returns 6.  All hexahedrals have 6 faces.
   */
  virtual unsigned int n_faces() const libmesh_override { return 6; }

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
   * @returns the side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const libmesh_override;

  /**
   * @returns the local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
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
   * @returns a primitive (4-noded) quad for
   * face i.
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
  Elem * _elemlinks_data[7+(LIBMESH_DIM>3)];

  /**
   * Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.  This matrix
   * is kept here, since the matrix (for the first 12
   * higher-order nodes) is identical for \p Hex20 and
   * \p Hex27.
   */
  static const unsigned short int _second_order_adjacent_vertices[12][2];

  /**
   * Vector that names a child sharing each second order node.
   */
  static const unsigned short int _second_order_vertex_child_number[27];

  /**
   * Vector that names the child vertex index for each second order node.
   */
  static const unsigned short int _second_order_vertex_child_index[27];

  /**
   * Master element node locations
   */
  static const Real _master_points[27][3];

  /**
   * Lookup table from child id, child node id to "possible node
   * location" (a simple dictionary-index in a 5x5x5 grid)
   */
  static const int _child_node_lookup[8][27];
};

} // namespace libMesh

#endif // LIBMESH_CELL_HEX_H
