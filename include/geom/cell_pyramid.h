// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all pyramid element types.
 */
class Pyramid : public Cell
{
public:

  /**
   * Default pyramid, one quad face, four triangular faces,
   * takes number of nodes and parent.
   * Derived classes implement 'true' elements.
   */
  Pyramid(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    Cell(nn, num_sides, p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
      this->set_interior_parent(nullptr);
  }

  Pyramid (Pyramid &&) = delete;
  Pyramid (const Pyramid &) = delete;
  Pyramid & operator= (const Pyramid &) = delete;
  Pyramid & operator= (Pyramid &&) = delete;
  virtual ~Pyramid() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * Geometric constants for all Pyramids.
   */
  static const int num_sides = 5;
  static const int num_edges = 8;
  static const int num_children = 0; // not implemented

  /**
   * \returns 5.  All pyramid-derivatives are guaranteed to have at
   * least 5 nodes.
   */
  virtual unsigned int n_nodes() const override { return 5; }

  /**
   * \returns 5
   */
  virtual unsigned int n_sides() const override { return 5; }

  /**
   * \returns 5.  All pyramids have 5 vertices.
   */
  virtual unsigned int n_vertices() const override { return 5; }

  /**
   * \returns 8.  All pyramids have 8 edges.
   */
  virtual unsigned int n_edges() const override { return 8; }

  /**
   * \returns 5.  All pyramids have 5 faces.
   */
  virtual unsigned int n_faces() const override { return 5; }

  /**
   * \returns 10.
   */
  virtual unsigned int n_children() const override { return 10; }

  /**
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

   /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns An id associated with the \p s side of this element, as
   * defined solely by element vertices.  The id is not necessarily
   * unique, but should be close.  This is particularly useful in the
   * \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type low_order_key (const unsigned int s) const override;

 /**
   * \returns \p Pyramid5::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p Pyramid5::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A primitive triangle or quad for face i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override;

  /**
   * Rebuilds a primitive triangle or quad for face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override;

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int e) const override final;

  /**
  * \returns Node 4 (the apex) if \p is at the apex within the tolerance
  * \p tol. If \p is not at the apex within the tolerance \p tol, will
  * return invalid_uint.
  *
  * The Jacobian at the apex is singular because directional derivatives
  * are defined in four directions, independently (four edges meet).
  * Therefore, we can check for this in a caught, failed inverse mapping
  * and return the correct point without failure if it is the apex.
  * This behavior _does_ require exceptions to be enabled.
  */
  unsigned int local_singular_node(const Point & p, const Real tol = TOLERANCE*TOLERANCE) const override final;

  /**
   * \returns true iff the node at the given index has a singular
   * mapping; i.e. is the degree-4 node on a Pyramid.
   */
  virtual bool is_singular_node(unsigned int node_idx) const override final { return (node_idx == 4); }

  /**
   * One quad side, four orientations.
   */
  virtual unsigned int n_permutations() const override final { return 4; }

  virtual bool is_flipped() const override final;

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int n) const override;

  /**
   * This maps each edge to the sides that contain said edge.
   */
  static const unsigned int edge_sides_map[8][2];

  virtual bool on_reference_element(const Point & p,
                                    const Real eps = TOLERANCE) const override final;

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[6+(LIBMESH_DIM>3)];

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

  /**
   * This maps the \f$ j^{th} \f$ node to the 3 or 4 edge ids adjacent
   * to the node. The edge numbering matches the ones used in the
   * derived classes' edge_nodes_map. An edge index of 99 is used to
   * indicate that there is no adjacent edge. This data structure is
   * used in the Pyramid::edges_adjacent_to_node() override and is
   * shared by all the derived Hex types.
   */
  static const unsigned int adjacent_edges_map[/*num_vertices*/5][/*max_adjacent_edges*/4];
};

} // namespace libMesh

#endif // LIBMESH_CELL_PYRAMID_H
