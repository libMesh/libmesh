// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_INF_PRISM_H
#define LIBMESH_CELL_INF_PRISM_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf.h"

namespace libMesh
{

/**
 * The \p InfPrism is an element in 3D with 4 sides.
 * The \f$ 5^{th} \f$ side is theoretically located at infinity,
 * and therefore not accounted for.
 * However, one could say that the \f$ 5^{th} \f$ side actually
 * does exist in the mesh, since the outer nodes are located
 * at a specific distance from the mesh origin (and therefore
 * define a side).  Still, this face is not to be used!
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief The base class for all 3D infinite prismatic element types.
 */
class InfPrism : public InfCell
{
public:

  /**
   * Default infinite prism element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  InfPrism(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    InfCell(nn, InfPrism::n_sides(), p, _elemlinks_data, nodelinkdata)
  {}

  InfPrism (InfPrism &&) = delete;
  InfPrism (const InfPrism &) = delete;
  InfPrism & operator= (const InfPrism &) = delete;
  InfPrism & operator= (InfPrism &&) = delete;
  virtual ~InfPrism() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const override final
  {
    libmesh_assert_less(i, this->n_nodes());
    return Point(_master_points[i][0],
                 _master_points[i][1],
                 _master_points[i][2]);
  }

  /**
   * \returns 4.  Infinite elements have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  virtual unsigned int n_sides() const override final { return 4; }

  /**
   * \returns 6.  All infinite prisms (in our
   * setting) have 6 vertices.
   */
  virtual unsigned int n_vertices() const override final { return 6; }

  /**
   * \returns 6.  All infinite prisms have 6 edges,
   * 3 lying in the base, and 3 perpendicular to the base.
   */
  virtual unsigned int n_edges() const override final { return 6; }

  /**
   * \returns 4.  All prisms have 4 faces.
   */
  virtual unsigned int n_faces() const override final { return 4; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_children() const override final { return 4; }

  /**
   * \returns \p true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int i) const
    override final { return (i > 2 && i < 6); }

  /**
   * \returns \p true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override final;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const override final;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const override;

  /**
   * \returns InfPrism6::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns InfPrism6::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A primitive (3-noded) tri or (4-noded) infquad for
   * face i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a primitive (3-noded) tri or (4-noded) infquad for face
   * i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override final;

  /**
   * @returns \p true when this element contains the point
   * \p p.  Customized for infinite elements, since knowledge
   * about the envelope can be helpful.
   */
  virtual bool contains_point (const Point & p, Real tol=TOLERANCE) const override;

  /**
   * One non-infinite side, three orientations.
   */
  virtual unsigned int n_permutations() const override final { return 3; }

  std::vector<unsigned int> sides_on_edge(const unsigned int e) const override final;

  /**
   * This maps each edge to the sides that contain said edge.
   */
  static const unsigned int edge_sides_map[6][2];

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[5+(LIBMESH_DIM>3)];

  /**
   * Master element node locations
   */
  static const Real _master_points[12][3];
};

} // namespace libMesh

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_CELL_INF_PRISM_H
