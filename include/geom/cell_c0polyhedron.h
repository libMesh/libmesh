// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FACE_C0POLYHEDRON_H
#define LIBMESH_FACE_C0POLYHEDRON_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/cell_polyhedron.h"
#include "libmesh/enum_order.h"

namespace libMesh
{

class Node;

/**
 * The \p C0Polyhedron is an element in 3D with an arbitrary (but fixed)
 * number of polygonal first-order (C0Polygon) sides.
 *
 * In master space ... I give up.  We've practically got to have a
 * definition distinct for every physical element anyway, and we've
 * got to give FE and QBase objects access to the physical element, so
 * we'll just do all our FE and Quadrature calculations there.
 *
 * \author Roy H. Stogner
 * \date 2025
 * \brief A 3D element with an arbitrary number of polygonal
 * first-order sides.
 */
class C0Polyhedron : public Polyhedron
{
public:

  /**
   * Constructor.  Takes the sides as an input, and allocates nodes
   * accordingly.
   *
   * By default this element has no parent.
   *
   * We'll set up a simple default tetrahedralization here, but if users
   * want to adjust node positions later they'll want to retriangulate()
   * after.
   * @param mid_elem_node a reference to the unique pointer set to the mid-element node,
   *        only if the default (optimal) tetrahedralization algorithm fails, and the backup
   *        trivial algorithm is used to tetrahedralize the element.
   *        It's the caller's job to add this node to the mesh after!
   * @param p pointer to the parent element of this one (useful for h-refinement)
   */
  explicit
  C0Polyhedron (const std::vector<std::shared_ptr<Polygon>> & sides,
                std::unique_ptr<Node> & mid_elem_node,
                Elem * p=nullptr);

  C0Polyhedron (C0Polyhedron &&) = delete;
  C0Polyhedron (const C0Polyhedron &) = delete;
  C0Polyhedron & operator= (const C0Polyhedron &) = delete;
  C0Polyhedron & operator= (C0Polyhedron &&) = delete;
  virtual ~C0Polyhedron();

  /**
   * \returns An Elem of the specified type, which should be the same
   * type as \p this, wrapped in a smart pointer.
   *
   * This is not a complete clone() method (since e.g. it does not set
   * node pointers; the standard use case reassigns node pointers from
   * a different mesh), but it allows us to more easily copy polyhedra
   * objects that depend on side objects and require existing nodes
   * and so cannot be constructed with Elem::build().
   */
  std::unique_ptr<Elem> disconnected_clone() const override;

  /**
   * \returns \p C0POLYHEDRON.
   */
  virtual ElemType type () const override final { return C0POLYHEDRON; }

  /**
   * \returns the number of vertices.  For the simple C0Polyhedron
   * every node is a vertex.
   */
  virtual unsigned int n_vertices() const override final;

  /**
   * \returns the number of tetrahedra to break this into for
   * visualization.
   */
  virtual unsigned int n_sub_elem() const override
  { return this->n_subelements(); }

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
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_side(const unsigned int s) const override;

  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int e) const override;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const override;

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int n) const override;

  /**
   * \returns \p true if the element map is definitely affine within
   * numerical tolerances.  We don't even share a master element from
   * element to element, so we're going to return false here.
   */
  virtual bool has_affine_map () const override { return false; }

  /**
   * \returns FIRST
   */
  virtual Order default_order() const override { return FIRST; }

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  virtual void connectivity(const unsigned int sf,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const override;

  /**
   * Element refinement is not implemented for polyhedra.
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const override;

  /**
   * An optimized method for calculating the area of a C0Polyhedron.
   */
  virtual Real volume () const override;

  /**
   * An optimized method for calculating the centroid of a C0Polyhedron.
   */
  virtual Point true_centroid () const override;

  ElemType side_type (const unsigned int s) const override final;

  Point side_vertex_average_normal(const unsigned int s) const override final;

  /**
   * \returns the local side-indices of subelement sides of the
   * polyhedron
   * Each subelement here is a tetrahedron
   * @param i the index of the sub-element
   * The return array is indexed by sides of the subelement, and contains the
   * index of the side of the polyhedron if the sub-element side is on that side,
   * or libMesh::invalid_uint if the subelement side is internal to the polyhedron
   * NOTE: this routine involves a loop over the sides of the polyhedron
   */
  virtual std::array<int, 4> subelement_sides_to_poly_sides(unsigned int i) const;

  /**
   * Create a triangulation (tetrahedralization) based on the current
   * sides' triangulations.
   */
  virtual void retriangulate() override final;

protected:

  /**
   * Add to our triangulation (tetrahedralization).
   */
  void add_tet(int n1, int n2, int n3, int n4);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix used to create the elements children.
   */
  virtual Real embedding_matrix (const unsigned int /*i*/,
                                 const unsigned int /*j*/,
                                 const unsigned int /*k*/) const override
  { libmesh_not_implemented(); return 0; }

#endif // LIBMESH_ENABLE_AMR

  /**
   * Whether we have a mid element node
   */
  bool _has_mid_elem_node;
};


} // namespace libMesh

#endif // LIBMESH_FACE_C0POLYHEDRON_H
