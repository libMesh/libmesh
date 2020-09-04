// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_CELL_INF_HEX_H
#define LIBMESH_CELL_INF_HEX_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/cell_inf.h"

namespace libMesh
{

/**
 * The \p InfHex is an element in 3D with 5 sides.
 * The \f$ 6^{th} \f$ side is theoretically located at infinity,
 * and therefore not accounted for.
 * However, one could say that the \f$ 6^{th} \f$ side actually
 * does exist in the mesh, since the outer nodes are located
 * at a specific distance from the mesh origin (and therefore
 * define a side).  Still, this face is not to be used!
 *
 * \author Daniel Dreyer
 * \date 2003
 * \brief The base class for all 3D infinite hexahedral element types.
 */
class InfHex : public InfCell
{
public:

  /**
   * Default infinite brick element, takes number of nodes and
   * parent. Derived classes implement 'true' elements.
   */
  InfHex(const unsigned int nn, Elem * p, Node ** nodelinkdata) :
    InfCell(nn, InfHex::n_sides(), p, _elemlinks_data, nodelinkdata)
  {
    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 3)
      this->set_interior_parent(nullptr);
  }

  InfHex (InfHex &&) = delete;
  InfHex (const InfHex &) = delete;
  InfHex & operator= (const InfHex &) = delete;
  InfHex & operator= (InfHex &&) = delete;
  virtual ~InfHex() = default;

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
   * \returns 5.  Infinite elements have one side less
   * than their conventional counterparts, since one
   * side is supposed to be located at infinity.
   */
  virtual unsigned int n_sides() const override final { return 5; }

  /**
   * \returns 8.  All infinite hexahedra (in our
   * setting) have 8 vertices.
   */
  virtual unsigned int n_vertices() const override final { return 8; }

  /**
   * \returns \p true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int i) const
    override final { return (i > 3 && i < 8); }

  /**
   * \returns 8.  All infinite hexahedra have 8 edges,
   * 4 lying in the base, and 4 perpendicular to the base.
   */
  virtual unsigned int n_edges() const override final { return 8; }

  /**
   * \returns 5.  All infinite hexahedra have 5 faces.
   */
  virtual unsigned int n_faces() const override final { return 5; }

  /**
   * \returns 4.
   */
  virtual unsigned int n_children() const override final { return 4; }

  /**
   * \returns \p true if the specified child is on the
   * specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override final;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const override final;

  virtual std::vector<unsigned int> sides_on_edge(const unsigned int e) const override final;

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
   * \returns \p InfHex8::side_nodes_map[side][side_node] after doing some range checking.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * \returns \p InfHex8::edge_nodes_map[edge][edge_node] after doing some range checking.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A primitive (4-noded) quad or infquad for face i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a primitive (4-noded) quad or infquad for face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) override final;

  /**
   * \returns A quantitative assessment of element quality based on
   * the metric \p q specified by the user.
   */
  virtual Real quality (const ElemQuality q) const override;

  /**
   * \returns The suggested quality bounds for the hex based on
   * quality measure \p q.  These are the values suggested by the
   * CUBIT User's Manual.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const override;

  /**
   * \returns \p true when this element contains the point
   * \p p.  Customized for infinite elements, since knowledge
   * about the envelope can be helpful.
   */
  virtual bool contains_point (const Point & p, Real tol=TOLERANCE) const override;

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   */
  Elem * _elemlinks_data[6+(LIBMESH_DIM>3)];



  /**
   * For higher-order elements, namely \p InfHex16 and
   * \p InfHex18, the matrices for adjacent vertices
   * of second order nodes are quite similar (apart from
   * the face nodes, which are directly handled by \p InfHex18).
   * Therefore hold this matrix here, so that both can re-use
   * this.  Matrix that tells which vertices define the location
   * of mid-side (or second-order) nodes.
   */
  static const unsigned short int _second_order_adjacent_vertices[8][2];

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

#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#endif // LIBMESH_CELL_INF_HEX_H
