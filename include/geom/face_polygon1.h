// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FACE_POLYGON1_H
#define LIBMESH_FACE_POLYGON1_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face_polygon.h"

namespace libMesh
{

/**
 * The \p POLYGON1 is an element in 2D with an arbitrary (but fixed)
 * number of first-order (EDGE2) sides.
 *
 * Sides and vertices are numbered clockwise starting from 0.
 *
 * In master space, a POLYGON1 has point 0 at (0,0) and is a regular
 * polygon of side length 1.
 *
 * \author Roy H. Stogner
 * \date 2025
 * \brief A 2D element with an arbitrary number of first-order sides.
 */
class Polygon1 : public Polygon
{
public:

  /**
   * Constructor.  Takes the number of sides as an input, and
   * allocates the same number of nodes, one node for the vertex
   * between each two sides.
   *
   * By default this element has no parent.
   *
   * We'll set up a simple default triangulation here, but if users
   * want to make a non-convex polygon they'll want to retriangulate()
   * after setting up its nodes.
   */
  explicit
  Polygon1 (const unsigned int num_sides, Elem * p=nullptr);

  Polygon1 (Polygon1 &&) = delete;
  Polygon1 (const Polygon1 &) = delete;
  Polygon1 & operator= (const Polygon1 &) = delete;
  Polygon1 & operator= (Polygon1 &&) = delete;
  virtual ~Polygon1() = default;

  /**
   * \returns \p POLYGON1.
   */
  virtual ElemType type () const override final { return POLYGON1; }

  /**
   * \returns the number of triangles to break this into for
   * visualization.
   */
  virtual unsigned int n_sub_elem() const override { return this->n_sides(); }

  /**
   * \returns The local node number for the node opposite to node n on
   * side \p opposite_side(s) (for n_sides() even), or throws an error
   * otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const override final;

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
   * Polygon1 polygons have only vertex nodes, two per side/edge
   */
  virtual unsigned int n_nodes_per_side() const override { return 2;}

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge (== is_node_on_side in 2D).
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override
  { return this->is_node_on_side(n,e); }

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.
   */
  virtual bool has_affine_map () const override;

  /**
   * \returns SECOND.
   */
  virtual Order default_order() const override;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i,
                                                bool proxy=false) override;

  /**
   * Rebuilds an EDGE2 coincident with side i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * \returns 2 for all \p n.
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int) const override
  { return 2; }

  /**
   * Element refinement is not implemented for polygons.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * An optimized method for calculating the area of a POLYGON1.
   */
  virtual Real volume () const override;

  /**
   * An optimized method for calculating the centroid of a POLYGON1.
   */
  virtual Point true_centroid () const override;

  virtual void permute(unsigned int perm_num) override final;

  virtual void flip(BoundaryInfo *) override final;

  ElemType side_type (const unsigned int s) const override final;

  /**
   * Create a triangulation from the current node locations.
   */
  virtual void retriangulate() override final;

protected:

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual Real embedding_matrix (const unsigned int /*i*/,
                                 const unsigned int /*j*/,
                                 const unsigned int /*k*/) const override
  { libmesh_not_implemented(); return 0; }

#endif // LIBMESH_ENABLE_AMR

};


} // namespace libMesh

#endif // LIBMESH_FACE_POLYGON_H
