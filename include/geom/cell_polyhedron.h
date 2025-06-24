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



#ifndef LIBMESH_CELL_POLYHEDRON_H
#define LIBMESH_CELL_POLYHEDRON_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/cell.h"

namespace libMesh
{

// Forward declarations
class Polygon;

/**
 * The \p Polyhedron is an element in 3D with an arbitrary
 * number of polygonal faces.
 *
 * \author Roy H. Stogner
 * \date 2025
 * \brief A 3D element with an arbitrary number of polygonal faces
 */
class Polyhedron : public Cell
{
public:

  /**
   * Arbitrary polyhedral element, takes a vector of shared pointers
   * to sides (which should already be constructed), and a (probably
   * unused) parent link. Derived classes implement 'true' elements.
   */
  Polyhedron (const std::vector<std::shared_ptr<Polygon>> & sides,
              Elem * p);

  Polyhedron (Polyhedron &&) = delete;
  Polyhedron (const Polyhedron &) = delete;
  Polyhedron & operator= (const Polyhedron &) = delete;
  Polyhedron & operator= (Polyhedron &&) = delete;
  virtual ~Polyhedron() = default;

  /**
   * \returns true - polyhedron subclasses can have numbers of sides
   * and nodes which vary at runtime.
   */
  virtual bool runtime_topology() const override { return true; }

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   *
   * This implementation returns the master vertices; subclasses with
   * more points will need to further override.
   *
   * Trying to come up with a "master" shape for an arbitrary
   * polyhedron is *too* arbitrary, so we're just returning physical
   * points here.
   */
  virtual Point master_point (const unsigned int i) const override;

  static const int num_children = 0; // Refinement not yet supported

  /**
   * \returns the number of nodes in the polyhedron.
   */
  virtual unsigned int n_nodes() const override { return this->_nodelinks_data.size(); }

  /**
   * \returns the number of sides, which is the number of links we had
   * to save minus one parent link and minus one interior_parent link.
   */
  virtual unsigned int n_sides() const override final { return _elemlinks_data.size()-2; }

  /**
   * \returns the number of edges.
   */
  virtual unsigned int n_edges() const override final { return _edge_lookup.size(); }

  /**
   * \returns the number of faces, which is the number of links we had
   * to save minus one parent link and minus one interior_parent link.
   */
  virtual unsigned int n_faces() const override final { return _elemlinks_data.size()-2; }

  /**
   * \returns num_children.
   */
  virtual unsigned int n_children() const override final { return num_children; }

  /**
   * \returns \p true if the specified child is on the
   * specified side.
   *
   * Not implemented ... indefinitely?  I don't think we'll be doing
   * hierarchic refinement for general polyhedra.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const override;

  /**
   * Throws an error.  I'm not sure how to define the "opposite
   * side" of a polyhedron.
   */
  virtual unsigned int opposite_side(const unsigned int s) const override final;

  /**
   * Throws an error - \p opposite_side(s) is too hard to define in
   * general on polyhedra.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const override final;

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
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const override;

  /**
   * Similar to Elem::local_side_node(), but instead of a side id, takes
   * an edge id and a node id on that edge and returns a local node number
   * for the Elem.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const override;

  /**
   * \returns A simple Polygon face for side i.
   */
  virtual std::unique_ptr<Elem> side_ptr (const unsigned int i) override final;

  /**
   * Rebuilds a simple Polygon coincident with face i.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & elem,
                         const unsigned int i) override final;

  /**
   * Copies the Polygon side coincident with side i.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) override;

  /**
   * Copies the Polygon side coincident with side i.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & elem,
                               const unsigned int i) override;

  // Avoid hiding deprecated version with different signature
  using Elem::build_side_ptr;

  /**
   * \returns An element coincident with edge \p i wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) override final;

  /**
   * Resets the loose element \p edge, which may currently point to a
   * different edge than \p i or even a different element than \p
   * this, to point to edge \p i on \p this.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) override final;

  /**
   * \returns The suggested quality bounds for
   * the polyhedron based on quality measure q.
   */
  virtual std::pair<Real, Real> qual_bounds (const ElemQuality q) const override;

  /**
   * \returns the (local) side numbers that touch the specified edge.
   */
  virtual std::vector<unsigned int> sides_on_edge(const unsigned int e) const override final;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const override final;

  /**
   * Maybe we have non-identity permutations, but trying to figure out
   * how many is an exercise in applied group theory, which is a bit
   * much for just expanded unit test coverage.
   */
  virtual unsigned int n_permutations() const override { return 1; }

  virtual void permute(unsigned int libmesh_dbg_var(perm_num)) override final
  { libmesh_assert_equal_to(perm_num, 0); }

  virtual bool is_flipped() const override final;

  /**
   * A flip is one of those general non-identity permutations we can't
   * handle.
   *
   * But we'll call it "not implemented" just in cases someone someday
   * wants to write an implementation to handle polyhedra with
   * topological symmetry planes.
   */
  virtual void flip(BoundaryInfo *) override final { libmesh_not_implemented(); };

  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int n) const override;

  /**
   * Create a triangulation (tetrahedralization) from the current node
   * locations and face triangulations.
   *
   * If this is not called by the user, it will be called
   * automatically when needed.
   *
   * If the user moves the polyhedron nodes, the triangulation may
   * need to be regenerated manually, or the changed mapping may be
   * lower quality or even inverted.
   *
   * Pure virtual because this strongly depends on what mid-edge or
   * mid-face nodes the subclass might have.
   */
  virtual void retriangulate() = 0;

  /**
   * \returns true iff the polyhedron is convex.  Some Polyhedron
   * methods currently assume (or in debugging modes, test for and
   * throw errors if not finding) convexity.
   */
  bool convex();

  /**
   * \returns the size of the triangulation of the polyhedron
   */
  unsigned int n_subelements() const { return cast_int<unsigned int>(this->_triangulation.size()); }

  /**
   * \returns the local indices of points on a subelement of the
   * polyhedron
   *
   * Each subelement here is a tetrahedron
   */
  virtual std::array<int, 4> subelement (unsigned int i) const
  {
    libmesh_assert_less(i, this->_triangulation.size());
    return this->_triangulation[i];
  }

  /**
   * \returns the master-space points of a subelement of the
   * polyhedron
   */
  virtual std::array<Point, 4> master_subelement (unsigned int i) const;

  /**
   * \returns the index of a subelement containing the master-space
   * point \p p, along with (the first three) barycentric coordinates
   * for \p; or return invalid_uint if no subelement contains p to
   * within tolerance \p tol.
   */
  std::tuple<unsigned int, Real, Real, Real>
    subelement_coordinates (const Point & p,
                            Real tol = TOLERANCE*TOLERANCE) const;

  virtual bool on_reference_element(const Point & p,
                                    const Real eps = TOLERANCE) const override final;

protected:

  /**
   * \returns Clones of the sides of \p this, wrapped in smart
   * pointers.
   *
   * This factors out much of the code required by the
   * disconnected_clone() method of Polyhedron subclasses.
   */
  std::vector<std::shared_ptr<Polygon>> side_clones() const;

  /**
   * Helper method for finding the non-cached side that shares an
   * edge, by examining the local node ids there.
   */
  bool side_has_edge_nodes(unsigned int side,
                           unsigned int min_node,
                           unsigned int max_node) const;

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
   * Data for links to sides.  We use shared_ptr so we can keep sides
   * consistent between neighbors.  We also have a bool for each
   * polygon to indicate whether the zeta (xi cross eta) direction for
   * the polygon points into this polyhedron (true) or out of it
   * (false), and a map from side-local node number to element-local
   * node number.
   */
  std::vector<std::tuple<std::shared_ptr<Polygon>,
                         bool,
                         std::vector<unsigned int>>> _sidelinks_data;

  /**
   * One entry for each polyhedron edge, a pair indicating the side
   * number and the edge-of-side number which identifies the edge.
   *
   * Although each edge can be reached from two sides, only the
   * lower-numbered side is in this lookup table.  For edges accessed
   * via the same side, the side edge numbering matches our edge
   * numbering.
   */
  std::vector<std::pair<unsigned int, unsigned int>> _edge_lookup;

  /**
   * Data for a triangulation (tetrahedralization) of the polyhedron.
   *
   * Positive int values should correspond to local node numbers.
   *
   * Negative int values may correspond to "special" points, e.g. a
   * centroid or skeleton point.
   */
  std::vector<std::array<int, 4>> _triangulation;
};


} // namespace libMesh

#endif // LIBMESH_CELL_POLYHEDRON_H
