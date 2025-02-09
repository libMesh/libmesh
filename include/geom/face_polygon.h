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



#ifndef LIBMESH_FACE_POLYGON_H
#define LIBMESH_FACE_POLYGON_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/face.h"

namespace libMesh
{

/**
 * The \p Polygon is an element in 2D with an arbitrary (but fixed)
 * number of sides.
 *
 * \author Roy H. Stogner
 * \date 2025
 * \brief A 2D element with an arbitrary number of sides.
 */
class Polygon : public Face
{
public:

  /**
   * Arbitrary polygonal element, takes number of nodes and
   * sides, and a (probably unused) parent link. Derived classes
   * implement 'true' elements.
   *
   * We heap-allocate our element and node links, but we can't do that
   * until *after* we've intitialized our parent class, which means
   * we need to do a little initialization of those links manually too.
   */
  Polygon (const unsigned int nn,
           const unsigned int ns,
           Elem * p) :
    Face(nn, ns, p, nullptr, nullptr),
    _elemlinks_data(ns+2), // neighbors + parent + interior_parent
    _nodelinks_data(nn)
  {
    // Do the manual initialization that Elem::Elem couldn't.  No need
    // to manually set nullptr, though, since std::vector does that.
    this->_elemlinks = _elemlinks_data.data();
    this->_nodes = _nodelinks_data.data();
    this->_elemlinks[0] = p;

    // Is this likely to ever be used?  We may do refinement with
    // polygons but it's probably not going to have a hierarchy...
    if (p)
      {
        this->subdomain_id() = p->subdomain_id();
        this->processor_id() = p->processor_id();
        _map_type = p->mapping_type();
        _map_data = p->mapping_data();

#ifdef LIBMESH_ENABLE_AMR
        this->set_p_level(p->p_level());
#endif
      }

    // Make sure the interior parent isn't undefined
    if (LIBMESH_DIM > 2)
      this->set_interior_parent(nullptr);
  }

  Polygon (Polygon &&) = delete;
  Polygon (const Polygon &) = delete;
  Polygon & operator= (const Polygon &) = delete;
  Polygon & operator= (Polygon &&) = delete;
  virtual ~Polygon() = default;

  /**
   * \returns true - polygon subclasses can have numbers of sides and
   * nodes which vary at runtime.
   */
  virtual bool runtime_topology() const { return true; }

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   *
   * This implementation returns the master vertices; subclasses with
   * more points will need to further override.
   */
  virtual Point master_point (const unsigned int i) const override;

  static const int num_children = 0; // Refinement not yet supported

  /**
   * \returns the number of nodes in the polygon.
   */
  virtual unsigned int n_nodes() const override { return this->_nodelinks_data.size(); }

  /**
   * \returns the number of sides, which is the number of links we had
   * to save minus one parent link and minus one interior_parent link.
   */
  virtual unsigned int n_sides() const override final { return _elemlinks_data.size()-2; }

  /**
   * \returns the number of vertices.  All polygons have 1 vertex per side.
   */
  virtual unsigned int n_vertices() const override final { return _elemlinks_data.size()-2; }

  /**
   * \returns the number of edges.  All polygons have 1 edge per side.
   */
  virtual unsigned int n_edges() const override final { return _elemlinks_data.size()-2; }

  /**
   * \returns the number of nodes on each side.  All polygons have the
   * same number of nodes on each side.
   */
  virtual unsigned int n_nodes_per_side() const = 0;

  /**
   * \returns num_children.
   */
  virtual unsigned int n_children() const override final { return num_children; }

  /**
   * \returns \p true if the specified child is on the
   * specified side.
   *
   * Not implemented ... indefinitely?  I don't think we'll be doing
   * hierarchic refinement for general polygons.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

  /**
   * \returns The side number opposite to \p s (for n_sides() even),
   * or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const override final;

  /**
   * Don't hide Elem::key() defined in the base class.
   */
  using Elem::key;

  /**
   * \returns An id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close.
   */
  virtual dof_id_type key () const override;

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
   * \returns The local node id for node \p side_node on side \p side of
   * this Elem.
   *
   * This implementation assumes a particular node numbering "style"
   * and may be overridden in subclasses.
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * Calls local_side_node(edge, edge_node). For 2D elements, there is an implied
   * equivalence between edges and sides, e.g. n_edges() == n_sides(), so we treat
   * these two functions the same.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A primitive (2-noded) edge for edge i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds an EDGE2 coincident with face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & elem,
                         const unsigned int i) override final;

  /**
   * \returns The suggested quality bounds for
   * the hex based on quality measure q.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const override;

  /**
   * n_sides sides, one orientation each.  Not marked final because
   * subclasses may have interior nodes with less symmetry
   */
  virtual unsigned int n_permutations() const override { return this->n_sides(); }

  virtual bool is_flipped() const override final;

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int n) const override;

  /**
   * Create a triangulation from the current node locations.
   *
   * If this is not called by the user, it will be called
   * automatically when needed.
   *
   * If the user moves the polygon nodes, the triangulation may need
   * to be regenerated manually, or the changed mapping may be lower
   * quality or even inverted.
   *
   * Pure virtual because this strongly depends on what mid-edge or
   * mid-face nodes the subclass might have.
   */
  virtual void retriangulate() = 0;

  /**
   * \returns the size of the triangulation of the polygon
   */
  unsigned int n_subtriangles() const { return cast_int<unsigned int>(this->_triangulation.size()); }

  /**
   * \returns the local indices of points on a subtriangle of the polygon
   */
  virtual std::array<int, 3> subtriangle (unsigned int i) const
  {
    libmesh_assert_less(i, this->_triangulation.size());
    return this->_triangulation[i];
  }

  /**
   * \returns the master-space points of a subtriangle of the polygon
   */
  virtual std::array<Point, 3> master_subtriangle (unsigned int i) const;

  /**
   * \returns the index of a subtriangle containing the master-space
   * point \p p, along with barycentric coordinates for \p, or return
   * invalid_uint if no subtriangle contains p to within tolerance \p
   * tol.
   */
  std::tuple<unsigned int, Real, Real>
    subtriangle_coordinates (const Point & p,
                             Real tol = TOLERANCE*TOLERANCE) const;

  virtual bool on_reference_element(const Point & p,
                                    const Real eps = TOLERANCE) const override final;

protected:

  /**
   * Data for links to parent/neighbor/interior_parent elements.
   *
   * There should be num_sides neighbors, preceded by any parent link
   * and followed by any interior_parent link.
   *
   * We're stuck with a heap allocation here since the number of sides
   * in a subclass is chosen at runtime.
   */
  std::vector<Elem *> _elemlinks_data;

  /**
   * Data for links to nodes.  We're stuck with a heap allocation here
   * since the number of nodes in a subclass is chosen at runtime.
   */
  std::vector<Node *> _nodelinks_data;

  /**
   * Data for a triangulation of the polygon.  Our mapping from master
   * space to physical space will be based on these subelements.
   *
   * Positive int values should correspond to local node numbers.
   *
   * Negative int values may correspond to "special" points, e.g. a
    * centroid or skeleton point.
   */
  std::vector<std::array<int, 3>> _triangulation;
};

} // namespace libMesh

#endif // LIBMESH_FACE_POLYGON_H
