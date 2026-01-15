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



#ifndef LIBMESH_ELEM_H
#define LIBMESH_ELEM_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/bounding_box.h"
#include "libmesh/dof_object.h"
#include "libmesh/id_types.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/node.h"
#include "libmesh/enum_elem_type.h" // INVALID_ELEM
#include "libmesh/multi_predicates.h"
#include "libmesh/pointer_to_pointer_iter.h"
#include "libmesh/int_range.h"
#include "libmesh/simple_range.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/hashword.h" // Used in compute_key() functions

// C++ includes
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <limits.h> // CHAR_BIT
#include <set>
#include <vector>
#include <memory>
#include <array>

namespace libMesh
{

// Forward declarations
class BoundaryInfo;
class Elem;
class MeshBase;
class MeshRefinement;
#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
#endif
template <class SideType, class ParentType>
class Side;
enum ElemQuality : int;
enum IOPackage : int;
enum Order : int;


/**
 * This is the base class from which all geometric element types are
 * derived.  The \p Elem class provides standard information such as
 * the number of nodes, edges, faces, vertices, children, and
 * neighbors it has, as well as access to (or the ability to
 * construct) these entities.
 *
 * An \p Elem has pointers to its \p Node objects.  Some of these
 * nodes live at the vertices of the element, while others may live on
 * edges (and faces in 3D), or interior to the element.  The number of
 * nodes in a given element, \p n_nodes(), is encoded into the name of
 * the class.  For example, a \p Tri3 has three nodes which correspond
 * to the vertices, while a \p Tri6 has six nodes, three of which are
 * located at vertices, and three which are located at the midpoint of
 * each edge.  Nodes on edges, faces, and element interiors are called
 * second-order nodes.
 *
 * A 1D Elem is an \p Edge, a 2D Elem is a \p Face, and a 3D Elem is a
 * \p Cell.  An \p Elem is composed of a number of sides, which can
 * be accessed as dim-1 dimensional \p Elem types.  For example, a \p
 * Hex8 is a 3D hexahedral element. A \p Hex8 has 6 sides, which are
 * \p Faces of type Quad4.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2007
 * \brief The base class for all geometric element types.
 */
class Elem : public ReferenceCountedObject<Elem>,
             public DofObject
{
protected:

  /**
   * Constructor.  Creates an element with \p n_nodes nodes,
   * \p n_sides sides, \p n_children possible children, and
   * parent \p p.  The constructor allocates the memory necessary
   * to support this data.
   */
  Elem (const unsigned int n_nodes,
        const unsigned int n_sides,
        Elem * parent,
        Elem ** elemlinkdata,
        Node ** nodelinkdata);

public:

  /**
   * Elems are responsible for allocating and deleting space for
   * storing pointers to their children during refinement, so they
   * cannot currently be (default) copy-constructed or copy-
   * assigned. We therefore explicitly delete these operations. In
   * addition, the DofObject base class currently has private copy
   * construction and assignment operators, so that prevents us from
   * copying Elems as well.
   */
  Elem (Elem &&) = delete;
  Elem (const Elem &) = delete;
  Elem & operator= (const Elem &) = delete;
  Elem & operator= (Elem &&) = delete;

  /**
   * Destructor.
   */
  virtual ~Elem() = default;

  /**
   * \returns The \p Point associated with local \p Node \p i.
   */
  const Point & point (const unsigned int i) const;

  /**
   * \returns The \p Point associated with local \p Node \p i
   * as a writable reference.
   */
  Point & point (const unsigned int i);

  /**
   * \returns The \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const = 0;

  /**
   * \returns The global id number of local \p Node \p i.
   */
  dof_id_type node_id (const unsigned int i) const;

  /**
   * \returns The local id number of global \p Node id \p i,
   * or \p invalid_uint if Node id \p i is not local.
   */
  unsigned int local_node (const dof_id_type i) const;

  /**
   * \returns The local index for the \p Node pointer \p node_ptr,
   * or \p invalid_uint if \p node_ptr is not a local node.
   */
  unsigned int get_node_index (const Node * node_ptr) const;

  /**
   * \returns A pointer to an array of local node pointers.
   */
  const Node * const * get_nodes () const;

  /**
   * \returns A const pointer to local \p Node \p i.
   */
  const Node * node_ptr (const unsigned int i) const;

  /**
   * \returns A non-const pointer to local \p Node \p i.
   */
  Node * node_ptr (const unsigned int i);

  /**
   * \returns A const reference to local \p Node \p i.
   */
  const Node & node_ref (const unsigned int i) const;

  /**
   * \returns A writable reference to local \p Node \p i.
   */
  Node & node_ref (const unsigned int i);

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * \returns The pointer to the \p Node with local number \p i as a
   * writable reference.
   *
   * \deprecated This setter cannot update the multiple node pointers
   * used in a general polyhedron; use the \p set_node overload that
   * takes an argument.
   */
  virtual Node * & set_node (const unsigned int i);
#endif // LIBMESH_ENABLE_DEPRECATED

  /**
   * Sets local \p Node \p i to refer to \p node.
   */
  virtual void set_node (const unsigned int i,
                         Node * node);

  /**
   * Nested classes for use iterating over all nodes of an element.
   */
  class NodeRefIter;
  class ConstNodeRefIter;

  /**
   * Returns a range with all nodes of an element, usable in
   * range-based for loops.  The exact type of the return value here
   * may be subject to change in future libMesh releases, but the
   * iterators will always dereference to produce a reference to a
   * Node.
   */
  SimpleRange<NodeRefIter> node_ref_range();

  SimpleRange<ConstNodeRefIter> node_ref_range() const;

  /**
   * \returns The subdomain that this element belongs to.
   */
  subdomain_id_type subdomain_id () const;

  /**
   * \returns The subdomain that this element belongs to as a
   * writable reference.
   */
  subdomain_id_type & subdomain_id ();

  /**
   * A static integral constant representing an invalid subdomain id.
   * See also DofObject::{invalid_id, invalid_unique_id, invalid_processor_id}.
   *
   * \note We don't use the static_cast(-1) trick here since
   * \p subdomain_id_type is sometimes a *signed* integer for
   * compatibility reasons (see libmesh/id_types.h).
   */
  static constexpr subdomain_id_type invalid_subdomain_id
    = std::numeric_limits<subdomain_id_type>::max();

  /**
   * \returns true iff this element type can vary in topology (e.g.
   * have different numbers of sides and/or nodes) at runtime.  For
   * such general polygons or polyhedra, APIs which assume a fixed
   * topology are not safe to use.
   */
  virtual bool runtime_topology() const { return false; }

  /**
   * \returns A pointer to the "reference element" associated
   * with this element.  The reference element is the image of this
   * element in reference parametric space. Importantly, it is *not*
   * an actual element in the mesh, but rather a Singleton-type
   * object, so for example all \p Quad4 elements share the same
   * \p reference_elem().
   *
   * If the element is of a type that can admit multiple topologies,
   * such as a Polygon subtype, then there is no reference element;
   * for such types this method should not be used.
   */
  const Elem * reference_elem () const;

  /**
   * \returns An id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.
   */
  virtual dof_id_type key (const unsigned int s) const = 0;

  /**
   * \returns An id associated with the \p s side of this element, as
   * defined solely by element vertices.  The id is not necessarily
   * unique, but should be close.  This is particularly useful in the
   * \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type low_order_key (const unsigned int s) const = 0;

  /**
   * \returns An id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close. Uses the same hash as the key(s) function, so for example
   * if "tri3" is side 0 of "tet4", then tri3->key()==tet4->key(0).
   */
  virtual dof_id_type key () const;

  /**
   * \returns \p true if two elements are equivalent, \p false
   * otherwise.  This is true if the elements are connected to
   * identical global nodes, regardless of how those nodes might be
   * numbered local to the elements.
   */
  bool operator == (const Elem & rhs) const;

  /**
   * \returns \p false if two elements are equivalent, \p true
   * otherwise.
   */
  bool operator != (const Elem & rhs) const;

  /**
   * \returns \p true if two elements have equal topologies, false
   * otherwise.
   * This is true if the elements connect to nodes of the same id in
   * the same order, and neighbors of the same id on each side, the
   * same id on any parent and/or interior_parent link, etc.
   */
  bool topologically_equal (const Elem & rhs) const;

  /**
   * \returns A const pointer to the \f$ i^{th} \f$ neighbor of this
   * element, or \p nullptr if \p MeshBase::find_neighbors() has not been
   * called.
   *
   * \note If \p MeshBase::find_neighbors() has been called and this
   * function still returns \p nullptr, then the side is on a boundary of
   * the domain.
   */
  const Elem * neighbor_ptr (unsigned int i) const;

  /**
   * \returns A non-const pointer to the \f$ i^{th} \f$ neighbor of this element.
   */
  Elem * neighbor_ptr (unsigned int i);

  /**
   * Nested "classes" for use iterating over all neighbors of an element.
   */
  typedef Elem * const *       NeighborPtrIter;
  typedef const Elem * const * ConstNeighborPtrIter;

  /**
   * Returns a range with all neighbors of an element, usable in
   * range-based for loops.  The exact type of the return value here
   * may be subject to change in future libMesh releases, but the
   * iterators will always dereference to produce a pointer to a
   * neighbor element (or a null pointer, for sides which have no
   * neighbors).
   */
  SimpleRange<NeighborPtrIter> neighbor_ptr_range();

  SimpleRange<ConstNeighborPtrIter> neighbor_ptr_range() const;

#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * \returns A pointer to the \f$ i^{th} \f$ neighbor of this element
   * for interior elements.  If an element is on a periodic
   * boundary, it will return a corresponding element on the opposite
   * side.
   */
  const Elem * topological_neighbor (const unsigned int i,
                                     const MeshBase & mesh,
                                     const PointLocatorBase & point_locator,
                                     const PeriodicBoundaries * pb) const;

  /**
   * \returns A writable pointer to the \f$ i^{th} \f$ neighbor of
   * this element for interior elements.  If an element is on a
   * periodic boundary, it will return a corresponding element on the
   * opposite side.
   */
  Elem * topological_neighbor (const unsigned int i,
                               MeshBase & mesh,
                               const PointLocatorBase & point_locator,
                               const PeriodicBoundaries * pb);

  /**
   * \returns \p true if the element \p elem in question is a neighbor or
   * topological neighbor of this element, \p false otherwise.
   */
  bool has_topological_neighbor (const Elem * elem,
                                 const MeshBase & mesh,
                                 const PointLocatorBase & point_locator,
                                 const PeriodicBoundaries * pb) const;
#endif

  /**
   * Assigns \p n as the \f$ i^{th} \f$ neighbor.
   */
  void set_neighbor (const unsigned int i, Elem * n);

  /**
   * \returns \p true if the element \p elem in question is a neighbor
   * of this element, \p false otherwise.
   */
  bool has_neighbor (const Elem * elem) const;

  /**
   * \returns If \p elem is a neighbor of a child of this element, a
   * pointer to that child, otherwise \p nullptr.
   */
  Elem * child_neighbor (Elem * elem);

  /**
   * \returns If \p elem is a neighbor of a child of this element, a
   * pointer to that child, otherwise \p nullptr.
   */
  const Elem * child_neighbor (const Elem * elem) const;

  /**
   * \returns \p true if this element has a side coincident
   * with a boundary (indicated by a \p nullptr neighbor), \p false
   * otherwise.
   */
  bool on_boundary () const;

  /**
   * \returns \p true if this element is "semilocal" to the calling
   * processor, which must specify its rank.
   *
   * This method is discouraged, as it uses the *old* definition of
   * semilocal (elements which are not local but which are point
   * neighbors of something local) rather than any of the new
   * definitions discussed in ghosting_functor.h
   */
  bool is_semilocal (const processor_id_type my_pid) const;

  /**
   * This function tells you which neighbor \p e is.
   * I.e. if s = a->which_neighbor_am_i(e); then
   * a->neighbor(s) will be an ancestor of e.
   */
  unsigned int which_neighbor_am_i(const Elem * e) const;

  /**
   * This function tells you which side the boundary element \p e is.
   * I.e. if e = a->build_side_ptr(s) or e = a->side_ptr(s); then
   * a->which_side_am_i(e) will be s.
   *
   * \note An \e exact floating point comparison of the nodal
   * positions of \p e is made with the nodal positions of \p this in
   * order to perform this test. The idea is that the test will return
   * a valid side id if \p e either directly shares Node pointers with
   * \p this, or was created by exactly copying some of the nodes of
   * \p this (e.g. through BoundaryMesh::sync()). In these
   * circumstances, non-fuzzy floating point equality is expected.
   *
   * \returns The side of \p this the element which \p e is, otherwise
   * \p invalid_uint.
   */
  unsigned int which_side_am_i(const Elem * e) const;

  /**
   * \returns The local node id for node \p side_node on side \p side of
   * this Elem. Simply relies on the \p side_nodes_map for each of the
   * derived types. For example,
   * Tri3::local_side_node(0, 0) -> 0
   * Tri3::local_side_node(0, 1) -> 1
   * Tri3::local_side_node(1, 0) -> 1
   * Tri3::local_side_node(1, 1) -> 2
   * etc...
   */
  virtual unsigned int local_side_node(unsigned int side,
                                       unsigned int side_node) const = 0;

  /**
   * Similar to Elem::local_side_node(), but instead of a side id, takes
   * an edge id and a node id on that edge and returns a local node number
   * for the Elem. The implementation relies on the "edge_nodes_map" tables
   * for 3D elements. For 2D elements, calls local_side_node(). Throws an
   * error if called on 1D elements.
   */
  virtual unsigned int local_edge_node(unsigned int edge,
                                       unsigned int edge_node) const = 0;

  /**
   * \returns \p true if a vertex of \p e is contained
   * in this element.  If \p mesh_connection is true, looks
   * specifically for containment possibilities of an element \p e
   * that is connected to \p this via membership in the same manifold
   * of the same mesh.
   */
  bool contains_vertex_of(const Elem * e, bool mesh_connection=false) const;

  /**
   * \returns \p true if an edge of \p e is contained in
   * this element.  (Internally, this is done by checking whether at
   * least two vertices of \p e are contained in this element).
   */
  bool contains_edge_of(const Elem * e) const;

  /**
   * This function finds all active elements (including this one)
   * which are in the same manifold as this element and which touch
   * the current active element at the specified point, which should
   * be a point in the current element.
   *
   * Elements which are not "in the same manifold" (e.g. the
   * interior_parent of a boundary element) will not be found with
   * this method.
   *
   * Elements which overlap the specified point but which are only
   * connected to the current element via elements which do not
   * overlap that point (e.g. in a folded or tangled mesh) are not
   * considered to "touch" the current element and will not be found
   * with this method.
   */
  void find_point_neighbors(const Point & p,
                            std::set<const Elem *> & neighbor_set) const;

  /**
   * This function finds all active elements (including this one) in
   * the same manifold as this element which touch this active element
   * at any point.
   */
  void find_point_neighbors(std::set<const Elem *> & neighbor_set) const;

  /**
   * This function finds all active elements (including this one) in
   * the same manifold as start_elem (which must be active and must
   * touch this element) which touch this element at any point.
   */
  void find_point_neighbors(std::set<const Elem *> & neighbor_set,
                            const Elem * start_elem) const;

  /**
   * Non-const version of function above. Fills a set of non-const Elem pointers.
   */
  void find_point_neighbors(std::set<Elem *> & neighbor_set,
                            Elem * start_elem);

  /**
   * This function finds all active elements in the same manifold as
   * this element which touch the current active element along the
   * whole edge defined by the two points \p p1 and \p p2.
   */
  void find_edge_neighbors(const Point & p1,
                           const Point & p2,
                           std::set<const Elem *> & neighbor_set) const;

  /**
   * This function finds all active elements in the same manifold as
   * this element which touch the current active element along any
   * edge (more precisely, at at least two points).
   *
   * In this case, elements are included even if they do not touch a
   * *whole* edge of this element.
   */
  void find_edge_neighbors(std::set<const Elem *> & neighbor_set) const;

  /**
   * This function finds all active elements (*not* including this
   * one) in the parent manifold of this element whose intersection
   * with this element has non-zero measure.
   */
  void find_interior_neighbors(std::set<const Elem *> & neighbor_set) const;

  /**
   * Non-const version of function above that fills up a vector of
   * non-const Elem pointers instead.
   */
  void find_interior_neighbors(std::set<Elem *> & neighbor_set);

  /**
   * Resets this element's neighbors' appropriate neighbor pointers
   * and its parent's and children's appropriate pointers
   * to point to null instead of to this.
   *
   * To be used before an element is deleted from a mesh.
   */
  void remove_links_to_me ();

  /**
   * Resets this element's neighbors' appropriate neighbor pointers
   * and its parent's and children's appropriate pointers
   * to point to the global remote_elem instead of this.
   * Used by the library before an element becomes remote on the
   * local processor.
   */
  void make_links_to_me_remote ();

  /**
   * Resets the \p neighbor_side pointers of our nth neighbor (and
   * its descendants, if appropriate) to point to this Elem instead of
   * to the global remote_elem.  Used by the library when a formerly
   * remote element is being added to the local processor.
   */
  void make_links_to_me_local (unsigned int n, unsigned int neighbor_side);

  /**
   * \returns \p true if this element is remote, false otherwise.
   *
   * A remote element (see \p RemoteElem) is a syntactic convenience --
   * it is a placeholder for an element which exists on some other
   * processor.  Local elements are required to have valid neighbors,
   * and these ghost elements may have remote neighbors for data
   * structure consistency.  The use of remote elements helps ensure
   * that any element we may access has a \p nullptr neighbor only if it
   * lies on the physical boundary of the domain.
   */
  virtual bool is_remote () const
  { return false; }

  /**
   * \returns The connectivity for this element in a specific
   * format, which is specified by the IOPackage tag.
   */
  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const = 0;

  /**
   * Writes the element connectivity for various IO packages
   * to the passed ostream "out".  Not virtual, since it is
   * implemented in the base class.
   */
  void write_connectivity (std::ostream & out,
                           const IOPackage iop) const;

  /**
   * \returns The type of element that has been derived from this
   * base class.
   */
  virtual ElemType type () const = 0;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the geometric dimension of the element.
   *
   * This is currently usable even for complicated subclasses with
   * runtime-varying topology.
   */
  static const unsigned int type_to_dim_map[INVALID_ELEM];

  /**
   * \returns The dimensionality of the object.
   */
  virtual unsigned short dim () const = 0;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of nodes in the element.
   *
   * This is only usable for simple types for which the node number
   * is fixed; for more general types like Polygon subclasses an actual
   * instantiated Elem must be queried.
   */
  static const unsigned int type_to_n_nodes_map[INVALID_ELEM];

  /**
   * \returns The number of nodes this element contains.
   */
  virtual unsigned int n_nodes () const = 0;

  /**
   * The maximum number of nodes *any* element can contain.
   * This is useful for replacing heap vectors with stack arrays.
   */
  static const unsigned int max_n_nodes = 27;

  /**
   * \returns An integer range from 0 up to (but not including)
   * the number of nodes this element contains.
   */
  IntRange<unsigned short> node_index_range () const;

  /**
   * \returns The number of nodes the given child of this element
   * contains.  Except in odd cases like pyramid refinement this will
   * be the same as the number of nodes in the parent element.
   */
  virtual unsigned int n_nodes_in_child (unsigned int /*c*/) const
  { return this->n_nodes(); }

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of sides on the element.
   *
   * This is only usable for simple types for which the node number
   * is fixed; for more general types like Polygon subclasses an actual
   * instantiated Elem must be queried.
   */
  static const unsigned int type_to_n_sides_map[INVALID_ELEM];

  /**
   * \returns The number of sides the element that has been derived
   * from this class has. In 2D the number of sides is the number
   * of edges, in 3D the number of sides is the number of faces.
   */
  virtual unsigned int n_sides () const = 0;

  /**
   * \returns The type of element for side \p s.
   */
  virtual ElemType side_type (const unsigned int s) const = 0;

  /**
   * \returns the normal (outwards-facing) of the side of the element at the vertex-average of the side
   * @param s the side of interest
   */
  virtual Point side_vertex_average_normal(const unsigned int s) const;

  /**
   * \returns An integer range from 0 up to (but not including)
   * the number of sides this element has.
   */
  IntRange<unsigned short> side_index_range () const;

  /**
   * \returns The number of neighbors the element that has been derived
   * from this class has.
   *
   * Only face (or edge in 2D) neighbors are stored, so this method
   * returns n_sides().  At one point we intended to allow derived
   * classes to override this, but too much current libMesh code
   * assumes n_neighbors==n_sides.
   */
  unsigned int n_neighbors () const
  { return this->n_sides(); }

  /**
   * \returns The number of vertices the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_vertices () const = 0;

  /**
   * \returns The number of edges the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_edges () const = 0;

  /**
   * \returns An integer range from 0 up to (but not including)
   * the number of edges this element has.
   */
  IntRange<unsigned short> edge_index_range () const;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of edges on the element.
   *
   * This is only usable for simple types for which the node number
   * is fixed; for more general types like Polygon subclasses an actual
   * instantiated Elem must be queried.
   */
  static const unsigned int type_to_n_edges_map[INVALID_ELEM];

  /**
   * \returns The number of faces the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_faces () const = 0;

  /**
   * \returns An integer range from 0 up to (but not including)
   * the number of faces this element has.
   */
  IntRange<unsigned short> face_index_range () const;

  /**
   * \returns The number of children the element that has been derived
   * from this class may have.
   */
  virtual unsigned int n_children () const = 0;

  /**
   * \returns \p true if the specified (local) node number is a vertex node.
   */
  virtual bool is_vertex(const unsigned int i) const = 0;

  /**
   * \returns \p true if the specified child has a vertex at the
   * specified (child-local) node number.
   * Except in odd cases like pyramid refinement the child will have
   * the same local structure as the parent element.
   */
  virtual bool is_vertex_on_child (unsigned int /*c*/,
                                   unsigned int n) const
  { return this->is_vertex(n); }

  /**
   * \returns \p true if this element has a vertex at the specified
   * (child-local) node number \p n of the specified child \p c.
   */
  virtual bool is_vertex_on_parent(unsigned int c,
                                   unsigned int n) const;

  /**
   * \returns \p true if the specified (local) node number is an edge node.
   * For 1D elements, is_edge() is equivalent to is_internal().
   */
  virtual bool is_edge(const unsigned int i) const = 0;

  /**
   * \returns \p true if the specified (local) node number is a face node.
   * For 2D elements, is_face() is equivalent to is_internal().
   * For 1D elements, is_face() == false.
   */
  virtual bool is_face(const unsigned int i) const = 0;

  /**
   * \returns \p true if the specified (local) node number is an internal node.
   */
  bool is_internal(const unsigned int i) const;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const = 0;

  /**
   * \returns the (local) node numbers on the specified side
   */
  virtual std::vector<unsigned int> nodes_on_side(const unsigned int /*s*/) const = 0;

  /**
   * \returns the (local) node numbers on the specified edge
   */
  virtual std::vector<unsigned int> nodes_on_edge(const unsigned int /*e*/) const = 0;

  /**
   * \returns the (local) side numbers that touch the specified edge
   */
  virtual std::vector<unsigned int> sides_on_edge(const unsigned int /*e*/) const = 0;

  /**
   * \returns the (local) edge numbers that touch the specified node
   */
  virtual std::vector<unsigned int> edges_adjacent_to_node(const unsigned int /*n*/) const = 0;

  /**
   * \returns \p true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const = 0;

  /**
   * \returns \p true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const = 0;

  /**
   * \returns The side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const;

  /**
   * \returns The local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const;

  /**
   * \returns The number of sub-elements this element may be broken
   * down into for visualization purposes.  For example, 1 for a
   * linear triangle, 4 for a quadratic (6-noded) triangle, etc...
   */
  virtual unsigned int n_sub_elem () const = 0;

  /**
   * \returns A temporary element coincident with side \p i.
   *
   * This method returns the _minimum_ element necessary to uniquely
   * identify the side.  For example, the side of a hexahedron is
   * always returned as a 4-noded quadrilateral, regardless of what
   * type of hex you are dealing with.  Important data like subdomain
   * id, p level, or mapping type may be omitted from the temporary
   * element.  If you want a first-class full-ordered face (i.e. a
   * 9-noded quad face for a 27-noded hexahedron), use the
   * build_side_ptr method.
   *
   * \note The const version of this function is non-virtual; it
   * simply calls the virtual non-const version and const_casts the
   * return type.
   */
  virtual std::unique_ptr<Elem> side_ptr (unsigned int i) = 0;
  std::unique_ptr<const Elem> side_ptr (unsigned int i) const;

  /**
   * Resets the loose element \p side, which may currently point to a
   * different side than \p i or even a different element than \p
   * this, to point to side \p i on \p this.  If \p side is currently
   * an element of the wrong type, it will be freed and a new element
   * allocated; otherwise no memory allocation will occur.
   *
   * This will cause \p side to be a minimum-ordered element, even if
   * it is handed a higher-ordered element that must be replaced.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual void side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) = 0;
  void side_ptr (std::unique_ptr<const Elem> & side, const unsigned int i) const;

  /**
   * \returns An temporary element coincident with side \p i wrapped
   * in a smart pointer.
   *
   * The element returned is full-ordered and full-featured, in
   * contrast to the side method.  For example, calling
   * build_side_ptr(0) on a 20-noded hex in subdomain 5 will build a
   * 8-noded quadrilateral coincident with face 0, assign it subdomain
   * id 5, and pass back the pointer.
   *
   * The side element's id() is undefined; it is a temporary element
   * not added to any mesh.
   *
   * A \p std::unique_ptr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i) = 0;
  std::unique_ptr<const Elem> build_side_ptr (const unsigned int i) const;

#ifdef LIBMESH_ENABLE_DEPRECATED
  /*
   * Older versions of libMesh supported a "proxy" option here.
   */
  virtual std::unique_ptr<Elem> build_side_ptr (const unsigned int i, bool proxy)
  { if (proxy) libmesh_error(); libmesh_deprecated(); return this->build_side_ptr(i); }

  std::unique_ptr<const Elem> build_side_ptr (const unsigned int i, bool proxy) const
  { if (proxy) libmesh_error(); libmesh_deprecated(); return this->build_side_ptr(i); }
#endif

  /**
   * Resets the loose element \p side, which may currently point to a
   * different side than \p i or even a different element than \p
   * this, to point to side \p i on \p this.  If \p side is currently
   * an element of the wrong type, it will be freed and a new element
   * allocated; otherwise no memory allocation will occur.
   *
   * This will cause \p side to be a full-ordered element, even if it
   * is handed a lower-ordered element that must be replaced.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual void build_side_ptr (std::unique_ptr<Elem> & side, const unsigned int i) = 0;
  void build_side_ptr (std::unique_ptr<const Elem> & side, const unsigned int i) const;

  /**
   * \returns An element coincident with edge \p i wrapped in a smart pointer.
   *
   * The element returned is full-ordered.  For example, calling
   * build_edge_ptr(0) on a 20-noded hex will build a 3-noded edge
   * coincident with edge 0 and pass back the pointer.  A \p
   * std::unique_ptr<Elem> is returned to prevent a memory leak.  This way
   * the user need not remember to delete the object.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual std::unique_ptr<Elem> build_edge_ptr (const unsigned int i) = 0;
  std::unique_ptr<const Elem> build_edge_ptr (const unsigned int i) const;

  /**
   * Resets the loose element \p edge, which may currently point to a
   * different edge than \p i or even a different element than \p
   * this, to point to edge \p i on \p this.  If \p edge is currently
   * an element of the wrong type, it will be freed and a new element
   * allocated; otherwise no memory allocation will occur.
   *
   * This will cause \p edge to be a full-ordered element, even if it
   * is handed a lower-ordered element that must be replaced.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual void build_edge_ptr (std::unique_ptr<Elem> & edge, const unsigned int i) = 0;
  void build_edge_ptr (std::unique_ptr<const Elem> & edge, const unsigned int i) const;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the default approximation order of elements of that type.
   *
   * This is currently usable even for complicated subclasses with
   * runtime-varying topology.
   */
  static const Order type_to_default_order_map[INVALID_ELEM];

  /**
   * \returns The default approximation order for this element.  This
   * is the order that will be used to compute the map to the
   * reference element.
   */
  virtual Order default_order () const = 0;

  /**
   * \returns The maximum supported approximation order for nodal
   * (Lagrange or Rational Bezier-Bernstein) variables on this element
   * type.  This is usually the same as the default order.
   */
  virtual Order supported_nodal_order() const { return default_order(); }

  /**
   * \returns The default approximation order for side elements of
   * this element type.  This may be lower for elements with 'bubble
   * functions' in the Lagrange basis.
   */
  virtual Order default_side_order () const { return default_order(); }

  /**
   * \returns The "true" geometric centroid of the element, c=(cx, cy,
   * cz), where:
   *
   * [cx]            [\int x dV]
   * [cy] := (1/V) * [\int y dV]
   * [cz]            [\int z dV]
   *
   * This method is virtual since some derived elements might want to
   * use shortcuts to compute their centroid. For most element types,
   * this method is more expensive than calling vertex_average(), so
   * if you only need a point which is located "somewhere" in the
   * interior of the element, consider calling vertex_average() instead.
   */
  virtual Point true_centroid () const;

  /**
   * \returns A Point at the average of the elment's vertices.
   *
   * \note This used to be the base class centroid() implementation, but
   * the centroid is only equal to the vertex average in some special cases.
   * The centroid() implementation now returns the "true" centroid of the
   * element (up to quadrature error).
   */
  Point vertex_average () const;

  /**
   * \returns The "circumcenter of mass" (area-weighted average of
   * triangulation circumcenters) of the element.
   *
   * Not implemented for infinite elements, not currently implemented
   * for 3D elements, currently ignores curvature of element edges.
   */
  virtual Point quasicircumcenter () const
  { libmesh_not_implemented(); }

  /**
   * \returns The minimum vertex separation for the element.
   */
  virtual Real hmin () const;

  /**
   * \returns The maximum vertex separation for the element.
   */
  virtual Real hmax () const;

  /**
   * \returns The (length/area/volume) of the geometric element.
   *
   * If the element is twisted or inverted such that the mapping
   * Jacobian is singular at any point, implementations of this method
   * may return a "net" volume or may simply return NaN.
   */
  virtual Real volume () const;

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   *
   * The base class implementation determines a bounding box for the
   * element *nodes*, which should be sufficient for first order
   * finite elements.  Higher order geometric elements will need to
   * override with an implementation which takes curved elements into
   * account.
   */
  virtual BoundingBox loose_bounding_box () const;

  /**
   * \returns A quantitative assessment of element quality based on
   * the quality metric \p q specified by the user. Not all ElemQuality
   * metrics are supported for all Elem types; consult the Elem::quality()
   * overrides for specific Elem types to determine which quality metrics
   * are supported. The ElemQuality metrics with generic support for all
   * Elems with dimension > 1 are:
   * .) EDGE_LENGTH_RATIO - ratio of maximum to minimum edge (in 2D,
   *    side) length, where the min/max is taken over all Elem edges.
   * .) MIN,MAX_ANGLE - The minimum (respectively maximum) angle
   *    between all pairs of adjacent Elem edges, in degrees. Note
   *    that, in 3D, these are *not* the dihedral angles (angle
   *    between adjacent planar faces of the element), which we plan to
   *    add support for in the future. In 2D, we compute the angle
   *    between adjacent sides for this metric.
   */
  virtual Real quality (const ElemQuality q) const;

  /**
   * \returns The suggested quality bounds for the Elem based on
   * quality measure \p q.
   *
   * These are the values suggested by the CUBIT User's Manual.  Since
   * this function can have no possible meaning for an abstract Elem,
   * it is an error in the base class.
   */
  virtual std::pair<Real,Real> qual_bounds (const ElemQuality) const
  { libmesh_not_implemented(); return std::make_pair(0.,0.); }

  /**
   * \returns \p true if the physical point p is contained in this
   * element, false otherwise.
   *
   * For linear elements, performs an initial tight bounding box check
   * (as an optimization step) and (if that passes) then uses the
   * user-defined tolerance "tol" in a call to inverse_map() to actually
   * test if the point is in the element.  For quadratic elements, the
   * bounding box optimization is skipped, and only the inverse_map()
   * steps are performed.
   *
   * \note This routine should not be used to determine if a point
   * is merely "nearby" an element to within some tolerance. For that,
   * use Elem::close_to_point() instead.
   */
  virtual bool contains_point (const Point & p, Real tol=TOLERANCE) const;

  /**
   * \returns \p true if the master-space point p is contained in the
   * reference element corresponding to this element, false otherwise.
   *
   * Since we are doing floating point comparisons here the parameter
   * \p eps can be specified to indicate a tolerance.  For example,
   * \f$ x \le 1 \f$  becomes \f$ x \le 1 + \epsilon \f$.
   */
  virtual bool on_reference_element(const Point & p,
                                    const Real eps = TOLERANCE) const = 0;

  /**
   * \returns \p true if this element is "close" to the point p, where
   * "close" is determined by the tolerance tol.
   */
  virtual bool close_to_point(const Point & p, Real tol) const;

  /**
   * \returns \p true if edge \p i is positively oriented. An edge is
   * positively oriented iff its first vertex (i.e. zeroth node) is
   * lexicographically greater than its second vertex (i.e. first node).
   */
  bool positive_edge_orientation(const unsigned int i) const;

  /**
   * \returns \p true if face \p i is positively oriented. A face is
   * positively oriented iff the triangle defined by the lexicographically
   * least vertex and its two adjacent vertices on the same face is
   * positively oriented. Said triangle is positively oriented iff its
   * vertices are an odd permutation of their lexicographic ordering.
   */
  bool positive_face_orientation(const unsigned int i) const;

  /**
   * \returns \p true iff, for an edge \p e on side \p s, the node map for
   * side \p s is such that the first vertex (i.e. zeroth node) of \p e is
   * lower positioned than the second vertex (i.e. first node) of \p e.
   */
  bool relative_edge_face_order(const unsigned int e, const unsigned int s) const;

  /**
   * A helper function for copying generic element data (mapping,
   * subdomain, processor) from an element to a derived (child, side,
   * edge) element.  Useful for forwards compatibility when new data
   * is added.
   */
  void inherit_data_from(const Elem & src);

private:
  /**
   * Shared private implementation used by the contains_point()
   * and close_to_point() routines.  The box_tol tolerance is
   * used in the bounding box optimization, the map_tol tolerance is used
   * in the calls to inverse_map() and on_reference_element().
   */
  bool point_test(const Point & p, Real box_tol, Real map_tol) const;

public:
  /**
   * \returns \p true if the element map is definitely affine (i.e. the same at
   * every quadrature point) within numerical tolerances.
   */
  virtual bool has_affine_map () const { return false; }

  /**
   * \returns \p true if the element map is invertible everywhere on
   * the element, to within a user-specified tolerance. The tolerance
   * is generally used in comparisons against zero, so it should be an
   * absolute rather than a relative tolerance. Throws a
   * libmesh_not_implemented() error unless specialized by derived
   * classes.
   */
  virtual bool has_invertible_map(Real tol = TOLERANCE*TOLERANCE) const;

  /**
   * \returns \p true if the Lagrange shape functions on this element
   * are linear.
   */
  virtual bool is_linear () const { return false; }

  /**
   * Prints relevant information about the element.
   */
  void print_info (std::ostream & os=libMesh::out) const;

  /**
   * Prints relevant information about the element to a string.
   */
  std::string get_info () const;

  /**
   * \returns \p true if the element is active (i.e. has no active
   * descendants) or AMR is disabled, \p false otherwise.
   *
   * \note It suffices to check the first child only.
   */
  bool active () const;

  /**
   * \returns \p true if the element is an ancestor (i.e. has an
   * active child or ancestor child), \p false otherwise or when AMR
   * is disabled.
   */
  bool ancestor () const;

  /**
   * \returns \p true if the element is subactive (i.e. has no active
   * descendants), \p false otherwise or if AMR is disabled.
   */
  bool subactive () const;

  /**
   * \returns \p true if the element has any children (active or not),
   * \p false otherwise, or if AMR is disabled.
   */
  bool has_children () const;

  /**
   * \returns \p true if the element has any descendants other than
   * its immediate children, \p false otherwise, or if AMR is disabled.
   */
  bool has_ancestor_children () const;

  /**
   * \returns \p true if \p descendant is a child of \p this, or a
   * child of a child of \p this, etc., \p false otherwise or if AMR
   * is disabled.
   */
  bool is_ancestor_of(const Elem * descendant) const;

  /**
   * \returns A const pointer to the element's parent, or \p nullptr if
   * the element was not created via refinement.
   */
  const Elem * parent () const;

  /**
   * \returns A pointer to the element's parent, or \p nullptr if
   * the element was not created via refinement.
   */
  Elem * parent ();

  /**
   * Sets the pointer to the element's parent.
   * Dangerous! Only use this if you know what you are doing!
   */
  void set_parent (Elem * p);

  /**
   * \returns A pointer to the element's top-most (i.e. level-0) parent.
   *
   * That is, \p this if this is a level-0 element, this element's parent
   * if this is a level-1 element, this element's grandparent if this is
   * a level-2 element, etc...
   */
  const Elem * top_parent () const;

  /**
   * \returns The higher-dimensional Elem for which this Elem is a face.
   *
   * In some cases it is desirable to extract the boundary (or a subset thereof)
   * of a D-dimensional mesh as a (D-1)-dimensional manifold.  In this case
   * we may want to know the 'parent' element from which the manifold elements
   * were extracted.  We can easily do that for the level-0 manifold elements
   * by storing the D-dimensional parent.  This method provides access to that
   * element.
   *
   * This method returns nullptr if this->dim() == LIBMESH_DIM; in
   * such cases no data storage for an interior parent pointer has
   * been allocated.
   */
  const Elem * interior_parent () const;

  Elem * interior_parent ();

  /**
   * Sets the pointer to the element's interior_parent.
   * Dangerous! Only use this if you know what you are doing!
   */
  void set_interior_parent (Elem * p);

  /**
   * \returns The distance between nodes n1 and n2.
   *
   * Useful for computing the lengths of the sides of elements.
   */
  Real length (const unsigned int n1,
               const unsigned int n2) const;

  /**
   * \returns The number of adjacent vertices that uniquely define the
   * location of the \f$ n^{th} \f$ second-order node, or 0 for linear
   * elements.
   *
   * This method is useful when converting linear elements to quadratic
   * elements.
   *
   * \note \p n has to be greater than or equal to \p this->n_vertices().
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const;

  /**
   * \returns The element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node, or 0 for
   * linear elements.
   *
   * \note The value is always less than \p this->n_vertices(), while
   * \p n has to be greater than or equal to \p this->n_vertices().
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const;

  /**
   * \returns A pair (c,v), where
   * c == child index, and
   * v == element-local index of the \p \f$ n^{th} \f$
   *      second-order node on the parent element.
   * For linear elements, (0,0) is returned.
   *
   * \note The return values are always less than \p this->n_children()
   * and \p this->child_ptr(c)->n_vertices().
   *
   * \note \p n has to be greater than or equal to \p this->n_vertices().
   *
   * \note On refined second-order elements, the return value will
   * satisfy \p this->node_ptr(n) == this->child_ptr(c)->node_ptr(v).
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const;

  /**
   * \returns The ElemType of the associated second-order element
   * (which will be the same as the input if the input is already a
   * second-order ElemType) or INVALID_ELEM for elements that cannot be
   * converted into higher order equivalents.
   *
   * For example, when \p this is a \p TET4, then \p TET10 is returned.
   *
   * For some elements, there exist two second-order equivalents, e.g.
   * for \p Quad4 there is \p Quad8 and \p Quad9.  When the optional
   * \p full_ordered is \p true, then \p QUAD9 is returned.  When
   * \p full_ordered is \p false, then \p QUAD8 is returned.
   */
  static ElemType second_order_equivalent_type (const ElemType et,
                                                const bool full_ordered=true);

  /**
   * \returns The element type of the associated first-order element,
   * or \p INVALID_ELEM for first-order or other elements that cannot be
   * converted into lower order equivalents.
   *
   * For example, when \p this is a \p TET10, then \p TET4 is returned.
   */
  static ElemType first_order_equivalent_type (const ElemType et);

  /**
   * \returns The ElemType of the associated "complete" order element
   * (which will be the same as the input if the input is already a
   * complete-order ElemType), or INVALID_ELEM for elements that cannot be
   * converted into complete-order equivalents.
   *
   * The "complete" version of an element is an element which can
   * represent the same geometry but which has nodes available to
   * restore degrees of freedom on any vertex, edge, or face.
   *
   * For example, when \p this is a \p TET4, then \p TET14 is returned.
   */
  static ElemType complete_order_equivalent_type (const ElemType et);

  /**
   * \returns The refinement level of the current element.
   *
   * If the element's parent is \p nullptr then by convention it is at
   * level 0, otherwise it is simply at one level greater than its
   * parent.
   */
  unsigned int level () const;

  /**
   * \returns The value of the p refinement level of an active
   * element, or the minimum value of the p refinement levels
   * of an ancestor element's descendants.
   */
  unsigned int p_level () const;

  /**
   * \returns \p true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const = 0;

  /**
   * \returns The value of the mapping type for the element.
   */
  ElemMappingType mapping_type () const;

  /**
   * Sets the value of the mapping type for the element.
   */
  void set_mapping_type (const ElemMappingType type);

  /**
   * \returns The value of the mapping data for the element.
   */
  unsigned char mapping_data () const;

  /**
   * Sets the value of the mapping data for the element.
   */
  void set_mapping_data (const unsigned char data);


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Enumeration of possible element refinement states.
   */
  enum RefinementState { COARSEN = 0,
                         DO_NOTHING,
                         REFINE,
                         JUST_REFINED,
                         JUST_COARSENED,
                         INACTIVE,
                         COARSEN_INACTIVE,
                         INVALID_REFINEMENTSTATE };

  /**
   * \returns A constant pointer to the \f$ i^{th} \f$ child for this element.
   * For internal use only - skips assertions about null pointers.
   */
  const Elem * raw_child_ptr (unsigned int i) const;

  /**
   * \returns A constant pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  const Elem * child_ptr (unsigned int i) const;

  /**
   * \returns A non-constant pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  Elem * child_ptr (unsigned int i);

  /**
   * Nested classes for use iterating over all children of a parent
   * element.
   */
  class ChildRefIter;
  class ConstChildRefIter;

  /**
   * Returns a range with all children of a parent element, usable in
   * range-based for loops.  The exact type of the return value here
   * may be subject to change in future libMesh releases, but the
   * iterators will always dereference to produce a reference to a
   * child element.
   */
  SimpleRange<ChildRefIter> child_ref_range();

  SimpleRange<ConstChildRefIter> child_ref_range() const;

private:
  /**
   * Sets the pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  void set_child (unsigned int c, Elem * elem);

public:
  /**
   * \returns The child index which \p e corresponds to.
   *
   * I.e. if c = a->which_child_am_i(e); then a->child_ptr(c) will be
   * e.
   */
  unsigned int which_child_am_i(const Elem * e) const;

  /**
   * \returns \p true if the specified child is on the specified edge.
   */
  virtual bool is_child_on_edge(const unsigned int c,
                                const unsigned int e) const;

  /**
   * Adds a child pointer to the array of children of this element.
   * If this is the first child to be added, this method allocates
   * memory in the parent's _children array, otherwise, it just sets
   * the pointer.
   */
  void add_child (Elem * elem);

  /**
   * Adds a new child pointer to the specified index in the array of
   * children of this element.  If this is the first child to be added,
   * this method allocates memory in the parent's _children array,
   * otherwise, it just sets the pointer.
   */
  void add_child (Elem * elem, unsigned int c);

  /**
   * Replaces the child pointer at the specified index in the child array.
   */
  void replace_child (Elem * elem, unsigned int c);

  /**
   * Fills the vector \p family with the children of this element,
   * recursively.  Calling this method on a twice-refined element
   * will give you the element itself, its direct children, and their
   * children, etc...  When the optional parameter \p reset is
   * true, the vector will be cleared before the element and its
   * descendants are added.
   *
   * The family tree only includes ancestor and active elements. To
   * include subactive elements as well, use total_family_tree().
   */
  void family_tree (std::vector<const Elem *> & family,
                    bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void family_tree (std::vector<Elem *> & family,
                    bool reset = true);

  /**
   * Same as the \p family_tree() member, but also adds any subactive
   * descendants.
   */
  void total_family_tree (std::vector<const Elem *> & family,
                          bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void total_family_tree (std::vector<Elem *> & family,
                          bool reset = true);

  /**
   * Same as the \p family_tree() member, but only adds the active
   * children.  Can be thought of as removing all the inactive
   * elements from the vector created by \p family_tree, but is
   * implemented more efficiently.
   */
  void active_family_tree (std::vector<const Elem *> & active_family,
                           bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void active_family_tree (std::vector<Elem *> & active_family,
                           bool reset = true);

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void family_tree_by_side (std::vector<const Elem *> & family,
                            unsigned int side,
                            bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void family_tree_by_side (std::vector<Elem *> & family,
                            unsigned int side,
                            bool reset = true);

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void active_family_tree_by_side (std::vector<const Elem *> & family,
                                   unsigned int side,
                                   bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void active_family_tree_by_side (std::vector<Elem *> & family,
                                   unsigned int side,
                                   bool reset = true);

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void family_tree_by_neighbor (std::vector<const Elem *> & family,
                                const Elem * neighbor,
                                bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void family_tree_by_neighbor (std::vector<Elem *> & family,
                                Elem * neighbor,
                                bool reset = true);

  /**
   * Same as the \p family_tree_by_neighbor() member, but also adds
   * any subactive descendants.
   */
  void total_family_tree_by_neighbor (std::vector<const Elem *> & family,
                                      const Elem * neighbor,
                                      bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void total_family_tree_by_neighbor (std::vector<Elem *> & family,
                                      Elem * neighbor,
                                      bool reset = true);

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p subneighbor.  Only applicable when
   * \p this->has_neighbor(neighbor) and
   * \p neighbor->is_ancestor(subneighbor)
   */
  void family_tree_by_subneighbor (std::vector<const Elem *> & family,
                                   const Elem * neighbor,
                                   const Elem * subneighbor,
                                   bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void family_tree_by_subneighbor (std::vector<Elem *> & family,
                                   Elem * neighbor,
                                   Elem * subneighbor,
                                   bool reset = true);

  /**
   * Same as the \p family_tree_by_subneighbor() member, but also adds
   * any subactive descendants.
   */
  void total_family_tree_by_subneighbor (std::vector<const Elem *> & family,
                                         const Elem * neighbor,
                                         const Elem * subneighbor,
                                         bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void total_family_tree_by_subneighbor (std::vector<Elem *> & family,
                                         Elem * neighbor,
                                         Elem * subneighbor,
                                         bool reset = true);

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void active_family_tree_by_neighbor (std::vector<const Elem *> & family,
                                       const Elem * neighbor,
                                       bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void active_family_tree_by_neighbor (std::vector<Elem *> & family,
                                       Elem * neighbor,
                                       bool reset = true);

  /**
   * Same as the \p active_family_tree_by_neighbor() member, but the
   * \p neighbor here may be a topological (e.g. periodic boundary
   * condition) neighbor, not just a local neighbor.
   */
  void active_family_tree_by_topological_neighbor (std::vector<const Elem *> & family,
                                                   const Elem * neighbor,
                                                   const MeshBase & mesh,
                                                   const PointLocatorBase & point_locator,
                                                   const PeriodicBoundaries * pb,
                                                   bool reset = true) const;

  /**
   * Non-const version of function above; fills a vector of non-const pointers.
   */
  void active_family_tree_by_topological_neighbor (std::vector<Elem *> & family,
                                                   Elem * neighbor,
                                                   const MeshBase & mesh,
                                                   const PointLocatorBase & point_locator,
                                                   const PeriodicBoundaries * pb,
                                                   bool reset = true);

  /**
   * \returns The value of the refinement flag for the element.
   */
  RefinementState refinement_flag () const;

  /**
   * Sets the value of the refinement flag for the element.
   */
  void set_refinement_flag (const RefinementState rflag);

  /**
   * \returns The value of the p-refinement flag for the element.
   */
  RefinementState p_refinement_flag () const;

  /**
   * Sets the value of the p-refinement flag for the element.
   */
  void set_p_refinement_flag (const RefinementState pflag);

  /**
   * \returns The maximum value of the p-refinement levels of
   * an ancestor element's descendants.
   */
  unsigned int max_descendant_p_level () const;

  /**
   * \returns The minimum p-refinement level of elements which are
   * descended from this element, and which share a side with the
   * active \p neighbor.
   */
  unsigned int min_p_level_by_neighbor (const Elem * neighbor,
                                        unsigned int current_min) const;

  /**
   * \returns The minimum new p-refinement level (i.e. after refinement
   * and coarsening is done) of elements which are descended from this
   * element and which share a side with the active \p neighbor.
   */
  unsigned int min_new_p_level_by_neighbor (const Elem * neighbor,
                                            unsigned int current_min) const;

  /**
   * Sets the value of the p-refinement level for the element.
   *
   * \note The maximum p-refinement level is currently 255.
   */
  void set_p_level (const unsigned int p);

  /**
   * Sets the value of the p-refinement level for the element
   * without altering the p-level of its ancestors
   */
  void hack_p_level (const unsigned int p);

  /**
   * Sets the value of the p-refinement level for the element
   * without altering the p-level of its ancestors; also sets the
   * p_refinement_flag, simultaneously so that they can be safely
   * checked for mutual consistency
   */
  void hack_p_level_and_refinement_flag (const unsigned int p,
                                         RefinementState pflag);

  /**
   * Refine the element.
   */
  virtual void refine (MeshRefinement & mesh_refinement);

  /**
   * Coarsen the element.  This function is non-virtual since it is the same
   * for all element types.
   */
  void coarsen ();

  /**
   * Contract an active element, i.e. remove pointers to any
   * subactive children.  This should only be called via
   * MeshRefinement::contract, which will also remove subactive
   * children from the mesh.
   */
  void contract ();

#endif

#ifndef NDEBUG
  /**
   * Checks for consistent neighbor links on this element.
   */
  void libmesh_assert_valid_neighbors() const;

  /**
   * Checks for a valid id and pointers to nodes with valid ids on
   * this element.
   */
  void libmesh_assert_valid_node_pointers() const;
#endif // !NDEBUG

  /**
   * \returns The local node index of the given point IF said node
   * has a singular Jacobian for this element. If the given point
   * is not a node or is a node and does not have a singular Jacobian,
   * this will return invalid_uint.
   *
   * The intention is for this to be overridden in derived element
   * classes that do have nodes that have singular Jacobians. When
   * mapping failures are caught, we can check this to see if the
   * failed physical point is actually a singular point and
   * return the correct master point.
   */
  virtual unsigned int local_singular_node(const Point & /* p */, const Real /* tol */ = TOLERANCE*TOLERANCE) const
  { return invalid_uint; }

  /**
   * \returns true iff the node at the given index has a singular
   * mapping; i.e. is the degree-4 node on a Pyramid.
   */
  virtual bool is_singular_node(unsigned int /* node_i */) const { return false; }

  /**
   * \returns The local index of the center node on the side \p side.
   *
   * A center node is a node that is located at the centroid of the given side.
   * If the given side does not have a center node, this will return invalid_uint.
   */
  virtual unsigned int center_node_on_side(const unsigned short side) const;

protected:

  /**
   * The protected nested SideIter class is used to iterate over the
   * sides of this Elem.  It is a specially-designed class since
   * no sides are actually stored by the element.  This iterator-like
   * class has to provide the following three operations
   * 1) operator*
   * 2) operator++
   * 3) operator==
   * The definition can be found at the end of this header file.
   */
  class SideIter;

public:
  /**
   * Useful iterator typedefs
   */
  typedef Predicates::multi_predicate Predicate;

  /**
   * Data structure for iterating over sides.  Defined at the end of
   * this header file.
   */
  struct side_iterator;

  /**
   * Iterator accessor functions
   */
  side_iterator boundary_sides_begin();
  side_iterator boundary_sides_end();

private:
  /**
   * Side iterator helper functions.  Used to replace the begin()
   * and end() functions of the STL containers.
   */
  SideIter _first_side();
  SideIter _last_side();

public:

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p true if the element is an infinite element,
   * \p false otherwise.
   */
  virtual bool infinite () const = 0;

  /**
   * \returns \p true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.
   *
   * This is false for all nodes on non-infinite elements, so we won't
   * make it pure virtual, to simplify their code.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int /* n */) const
  { libmesh_assert (!this->infinite()); return false; }

  /**
   * \returns The origin for an infinite element.
   *
   * Currently, all infinite elements used in a mesh share the same
   * origin.  Override this in infinite element classes.
   */
  virtual Point origin () const { libmesh_not_implemented(); return Point(); }

#else

  static constexpr bool infinite () { return false; }

#endif

  /**
   * \returns An Elem of type \p type wrapped in a smart pointer.
   */
  static std::unique_ptr<Elem> build (const ElemType type,
                                      Elem * p=nullptr);

  /**
   * Calls the build() method above with a nullptr parent, and
   * additionally sets the newly-created Elem's id. This can be useful
   * when adding pre-numbered Elems to a Mesh via add_elem() calls.
   */
  static std::unique_ptr<Elem> build_with_id (const ElemType type,
                                              dof_id_type id);

  /**
   * \returns An Elem of the same type as \p this, wrapped in a smart
   * pointer.
   *
   * This is not a complete clone() method (since e.g. it does not set
   * node pointers; the standard use case reassigns node pointers from
   * a different mesh), but it is necessary to use this instead of
   * build() for runtime-polymorphic elements like Polygon subtypes
   * whose "type" depends on more than their type(), and it is useful
   * to use this for elements whose id, unique_id, extra integers,
   * etc. should be preserved in the near-clone.
   */
  virtual std::unique_ptr<Elem> disconnected_clone () const;

  /**
   * Returns the number of independent permutations of element nodes -
   * e.g. a cube can be reoriented to put side 0 where side N is (for
   * 0 <= N < 6) and then rotated in one of four ways, giving 24
   * possible permutations.
   *
   * Permutations which change the mapping Jacobian of an element
   * (i.e. flipping the element) are not allowed in this definition.
   */
  virtual unsigned int n_permutations() const = 0;

  /**
   * Permutes the element (by swapping node and neighbor pointers)
   * according to the specified index.
   *
   * This is useful for regression testing, by making it easy to make
   * a structured mesh behave more like an arbitrarily unstructured
   * mesh.
   *
   * This is so far *only* used for regression testing, so we do
   * not currently provide a way to permute any boundary side/edge ids
   * along with the element permutation.
   */
  virtual void permute(unsigned int perm_num) = 0;

  /**
   * Flips the element (by swapping node and neighbor pointers) to
   * have a mapping Jacobian of opposite sign.
   *
   * This is useful for automatically fixing up elements that have
   * been newly created (e.g. from extrusions) with a negative
   * Jacobian.
   *
   * If \p boundary_info is not null, swap boundary side/edge ids
   * consistently.
   */
  virtual void flip(BoundaryInfo * boundary_info) = 0;

  /**
   * \returns Whether the element is flipped compared to standard
   * libMesh (e.g. clockwise for 2D elements) node orientations.
   *
   * Always returns \p false if a 2D element is not in the XY plane or
   * a 1D element is not on the X axis; user code designed to work for
   * embedded manifolds should handle any consistent orientation, and
   * determining whether an orientation is consistent is not a local
   * operation.
   */
   virtual bool is_flipped() const = 0;

  /**
   * Flips the element (by swapping node and neighbor pointers) to
   * have a mapping Jacobian of opposite sign, iff we find a negative
   * orientation.  This only fixes flipped elements; for tangled
   * elements the only fixes possible are non-local.
   */
  void orient(BoundaryInfo * boundary_info);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * \returns The local node id on the parent which corresponds to node
   * \p n of child \p c, or \p invalid_uint if no such parent
   * node exists.
   */
  virtual unsigned int as_parent_node (unsigned int c,
                                       unsigned int n) const;

  /**
   * \returns All the pairs of nodes (indexed by local node id) which
   * should bracket node \p n of child \p c.
   */
  virtual
  const std::vector<std::pair<unsigned char, unsigned char>> &
  parent_bracketing_nodes(unsigned int c,
                          unsigned int n) const;

  /**
   * \returns All the pairs of nodes (indexed by global node id) which
   * should bracket node \p n of child \p c.
   */
  virtual
  const std::vector<std::pair<dof_id_type, dof_id_type>>
  bracketing_nodes(unsigned int c,
                   unsigned int n) const;


  /**
   * \returns The embedding matrix entry for the requested child.
   */
  virtual Real embedding_matrix (const unsigned int child_num,
                                 const unsigned int child_node_num,
                                 const unsigned int parent_node_num) const = 0;

  /**
   * \returns A "version number" that identifies which embedding
   * matrix is in use.
   *
   * Some element types may use a different embedding matrix depending
   * on their geometric characteristics.
   */
  virtual unsigned int embedding_matrix_version () const { return 0; }

#endif // LIBMESH_ENABLE_AMR


protected:

  /**
   * Default tolerance to use in has_affine_map().
   */
  static constexpr Real affine_tol = TOLERANCE*TOLERANCE;

  /**
   * \returns A hash key computed from a single node id.
   */
  static dof_id_type compute_key (dof_id_type n0);

  /**
   * \returns A hash key computed from two node ids.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1);

  /**
   * \returns A hash key computed from three node ids.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1,
                                  dof_id_type n2);

  /**
   * \returns A hash key computed from four node ids.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1,
                                  dof_id_type n2,
                                  dof_id_type n3);

  /**
   * Swaps two node_ptrs
   */
  void swap2nodes(unsigned int n1, unsigned int n2)
  {
    Node * temp = this->node_ptr(n1);
    this->set_node(n1, this->node_ptr(n2));
    this->set_node(n2, temp);
  }

  /**
   * Swaps two neighbor_ptrs
   */
  void swap2neighbors(unsigned int n1, unsigned int n2)
  {
    Elem * temp = this->neighbor_ptr(n1);
    this->set_neighbor(n1, this->neighbor_ptr(n2));
    this->set_neighbor(n2, temp);
  }

  /**
   * Swaps two sides in \p boundary_info, if it is non-null.
   */
  void swap2boundarysides(unsigned short s1, unsigned short s2,
                          BoundaryInfo * boundary_info) const;

  /**
   * Swaps two edges in \p boundary_info, if it is non-null.
   */
  void swap2boundaryedges(unsigned short e1, unsigned short e2,
                          BoundaryInfo * boundary_info) const;

  /**
   * Swaps three node_ptrs, "rotating" them.
   */
  void swap3nodes(unsigned int n1, unsigned int n2, unsigned int n3)
  {
    swap2nodes(n1, n2);
    swap2nodes(n2, n3);
  }

  /**
   * Swaps three neighbor_ptrs, "rotating" them.
   */
  void swap3neighbors(unsigned int n1, unsigned int n2,
                      unsigned int n3)
  {
    swap2neighbors(n1, n2);
    swap2neighbors(n2, n3);
  }

  /**
   * Swaps four node_ptrs, "rotating" them.
   */
  void swap4nodes(unsigned int n1, unsigned int n2, unsigned int n3,
                  unsigned int n4)
  {
    swap3nodes(n1, n2, n3);
    swap2nodes(n3, n4);
  }

  /**
   * Swaps four neighbor_ptrs, "rotating" them.
   */
  void swap4neighbors(unsigned int n1, unsigned int n2,
                      unsigned int n3, unsigned int n4)
  {
    swap3neighbors(n1, n2, n3);
    swap2neighbors(n3, n4);
  }


  /**
   * An implementation for simple (all sides equal) elements
   */
  template <typename Sideclass, typename Subclass>
  std::unique_ptr<Elem>
  simple_build_side_ptr(const unsigned int i);

  /**
   * An implementation for simple (all sides equal) elements
   */
  template <typename Subclass>
  void simple_build_side_ptr(std::unique_ptr<Elem> & side,
                             const unsigned int i,
                             ElemType sidetype);

  /**
   * An implementation for simple (all sides equal) elements
   */
  template <typename Subclass, typename Mapclass>
  void simple_side_ptr(std::unique_ptr<Elem> & side,
                       const unsigned int i,
                       ElemType sidetype);

  /**
   * An implementation for simple (all edges equal) elements
   */
  template <typename Edgeclass, typename Subclass>
  std::unique_ptr<Elem>
  simple_build_edge_ptr(const unsigned int i);

  /**
   * An implementation for simple (all edges equal) elements
   */
  template <typename Subclass>
  void simple_build_edge_ptr(std::unique_ptr<Elem> & edge,
                             const unsigned int i,
                             ElemType edgetype);


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Elem subclasses which don't do their own bracketing node
   * calculations will need to supply a static cache, since the
   * default calculation is slow.
   */
  virtual
  std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char>>>>> &
  _get_bracketing_node_cache() const
  {
    static std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char>>>>> c;
    libmesh_error();
    return c;
  }

  /**
   * Elem subclasses which don't do their own child-to-parent node
   * calculations will need to supply a static cache, since the
   * default calculation is slow.
   */
  virtual
  std::vector<std::vector<std::vector<signed char>>> &
  _get_parent_indices_cache() const
  {
    static std::vector<std::vector<std::vector<signed char>>> c;
    libmesh_error();
    return c;
  }

#endif // LIBMESH_ENABLE_AMR

public:

  /**
   * Replaces this element with \p nullptr for all of its neighbors.
   * This is useful when deleting an element.
   */
  void nullify_neighbors ();

protected:

  /**
   * Pointers to the nodes we are connected to.
   */
  Node ** _nodes;

  /**
   * Pointers to this element's parent and neighbors, and for
   * lower-dimensional elements' interior_parent.
   */
  Elem ** _elemlinks;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * unique_ptr to array of this element's children.
   *
   * A Mesh ultimately owns the child Elems so we are not responsible
   * for deleting them, but we are responsible for cleaning up the
   * array allocated to hold those Elems, hence the unique_ptr.
   */
  std::unique_ptr<Elem *[]> _children;
#endif

  /**
   * The subdomain to which this element belongs.
   */
  subdomain_id_type _sbd_id;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * h refinement flag. This is stored as an unsigned char
   * to save space.
   */
  unsigned char _rflag;

  /**
   * p refinement flag. This is stored as an unsigned char
   * to save space.
   */
  unsigned char _pflag;

  /**
   * p refinement level - the difference between the
   * polynomial degree on this element and the minimum
   * polynomial degree on the mesh.
   * This is stored as an unsigned char to save space.
   * In theory, these last four bytes might have
   * been padding anyway.
   */
  unsigned char _p_level;
#endif

  /**
   * Mapping function type; currently either 0 (LAGRANGE) or 1
   * (RATIONAL_BERNSTEIN).
   */
  unsigned char _map_type;

  /**
   * Mapping function data; currently used when needed to store the
   * RATIONAL_BERNSTEIN nodal weight data index.
   */
  unsigned char _map_data;
};



// ------------------------------------------------------------
// Elem helper classes
//
class
Elem::NodeRefIter : public PointerToPointerIter<Node>
{
public:
  NodeRefIter (Node * const * nodepp) : PointerToPointerIter<Node>(nodepp) {}
};


class
Elem::ConstNodeRefIter : public PointerToPointerIter<const Node>
{
public:
  ConstNodeRefIter (const Node * const * nodepp) : PointerToPointerIter<const Node>(nodepp) {}
};


#ifdef LIBMESH_ENABLE_AMR
class
Elem::ChildRefIter : public PointerToPointerIter<Elem>
{
public:
  ChildRefIter (Elem * const * childpp) : PointerToPointerIter<Elem>(childpp) {}
};


class
Elem::ConstChildRefIter : public PointerToPointerIter<const Elem>
{
public:
  ConstChildRefIter (const Elem * const * childpp) : PointerToPointerIter<const Elem>(childpp) {}
};



inline
SimpleRange<Elem::ChildRefIter> Elem::child_ref_range()
{
  libmesh_assert(_children);
  return {_children.get(), _children.get() + this->n_children()};
}


inline
SimpleRange<Elem::ConstChildRefIter> Elem::child_ref_range() const
{
  libmesh_assert(_children);
  return {_children.get(), _children.get() + this->n_children()};
}
#endif // LIBMESH_ENABLE_AMR




// ------------------------------------------------------------
// global Elem functions

inline
std::ostream & operator << (std::ostream & os, const Elem & e)
{
  e.print_info(os);
  return os;
}


// ------------------------------------------------------------
// Elem class member functions
inline
Elem::Elem(const unsigned int nn,
           const unsigned int ns,
           Elem * p,
           Elem ** elemlinkdata,
           Node ** nodelinkdata) :
  _nodes(nodelinkdata),
  _elemlinks(elemlinkdata),
  _sbd_id(0),
#ifdef LIBMESH_ENABLE_AMR
  _rflag(Elem::DO_NOTHING),
  _pflag(Elem::DO_NOTHING),
  _p_level(0),
#endif
  _map_type(p ? p->mapping_type() : 0),
  _map_data(p ? p->mapping_data() : 0)
{
  this->processor_id() = DofObject::invalid_processor_id;

  // If this ever legitimately fails we need to increase max_n_nodes
  libmesh_assert_less_equal(nn, max_n_nodes);

  // We currently only support refinement of elements into child
  // elements of the same type.  We can't test elem->type() here,
  // because that's virtual and we're still in the base class
  // constructor, but we can at least usually verify constency with
  // the arguments we were handed.
#ifndef NDEBUG
  if (p && !p->runtime_topology())
    {
      libmesh_assert_equal_to(nn, p->n_nodes());
      libmesh_assert_equal_to(ns, p->n_sides());
    }
#endif

  // Initialize the nodes data structure if we're given a pointer to
  // memory for it.
  if (_nodes)
    {
      for (unsigned int n=0; n<nn; n++)
        _nodes[n] = nullptr;
    }

  // Initialize the neighbors/parent data structure
  // _elemlinks = new Elem *[ns+1];

  // Initialize the elements data structure if we're given a pointer
  // to memory for it.  If we *weren't* given memory for it, e.g.
  // because a subclass like an arbitrary Polygon needs to
  // heap-allocate this memory, then that subclass will have to handle
  // this initialization too.
  if (_elemlinks)
    {
      _elemlinks[0] = p;

      for (unsigned int n=1; n<ns+1; n++)
        _elemlinks[n] = nullptr;

      // Optionally initialize data from the parent
      if (this->parent())
        {
          this->subdomain_id() = this->parent()->subdomain_id();
          this->processor_id() = this->parent()->processor_id();
          _map_type = this->parent()->_map_type;
          _map_data = this->parent()->_map_data;

#ifdef LIBMESH_ENABLE_AMR
          this->set_p_level(this->parent()->p_level());
#endif
        }
    }
}



inline
const Point & Elem::point (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_not_equal_to (_nodes[i]->id(), Node::invalid_id);

  return *_nodes[i];
}



inline
Point & Elem::point (const unsigned int i)
{
  libmesh_assert_less (i, this->n_nodes());

  return *_nodes[i];
}



inline
dof_id_type Elem::node_id (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_not_equal_to (_nodes[i]->id(), Node::invalid_id);

  return _nodes[i]->id();
}



inline
unsigned int Elem::local_node (const dof_id_type i) const
{
  for (auto n : make_range(this->n_nodes()))
    if (this->node_id(n) == i)
      return n;

  return libMesh::invalid_uint;
}



inline
const Node * const * Elem::get_nodes () const
{
  return _nodes;
}



inline
const Node * Elem::node_ptr (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);

  return _nodes[i];
}



inline
Node * Elem::node_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);

  return _nodes[i];
}



inline
const Node & Elem::node_ref (const unsigned int i) const
{
  return *this->node_ptr(i);
}



inline
Node & Elem::node_ref (const unsigned int i)
{
  return *this->node_ptr(i);
}



inline
unsigned int Elem::get_node_index (const Node * node_ptr) const
{
  for (auto n : make_range(this->n_nodes()))
    if (this->_nodes[n] == node_ptr)
      return n;

  return libMesh::invalid_uint;
}



#ifdef LIBMESH_ENABLE_DEPRECATED
inline
Node * & Elem::set_node (const unsigned int i)
{
  libmesh_assert_less (i, this->n_nodes());

  libmesh_deprecated();

  return _nodes[i];
}
#endif // LIBMESH_ENABLE_DEPRECATED



inline
void Elem::set_node (const unsigned int i,
                     Node * node)
{
  libmesh_assert_less (i, this->n_nodes());

  _nodes[i] = node;
}



inline
subdomain_id_type Elem::subdomain_id () const
{
  return _sbd_id;
}



inline
subdomain_id_type & Elem::subdomain_id ()
{
  return _sbd_id;
}



inline
bool Elem::operator != (const Elem & rhs) const
{
  return !(*this == rhs);
}



inline
const Elem * Elem::neighbor_ptr (unsigned int i) const
{
  libmesh_assert_less (i, this->n_neighbors());

  return _elemlinks[i+1];
}



inline
Elem * Elem::neighbor_ptr (unsigned int i)
{
  libmesh_assert_less (i, this->n_neighbors());

  return _elemlinks[i+1];
}



inline
void Elem::set_neighbor (const unsigned int i, Elem * n)
{
  libmesh_assert_less (i, this->n_neighbors());

  _elemlinks[i+1] = n;
}



inline
bool Elem::has_neighbor (const Elem * elem) const
{
  for (auto n : this->neighbor_ptr_range())
    if (n == elem)
      return true;

  return false;
}



inline
Elem * Elem::child_neighbor (Elem * elem)
{
  for (auto n : elem->neighbor_ptr_range())
    if (n && n->parent() == this)
      return n;

  return nullptr;
}



inline
const Elem * Elem::child_neighbor (const Elem * elem) const
{
  for (auto n : elem->neighbor_ptr_range())
    if (n && n->parent() == this)
      return n;

  return nullptr;
}



inline
SimpleRange<Elem::NodeRefIter>
Elem::node_ref_range()
{
  return {_nodes, _nodes+this->n_nodes()};
}



inline
SimpleRange<Elem::ConstNodeRefIter>
Elem::node_ref_range() const
{
  return {_nodes, _nodes+this->n_nodes()};
}



inline
IntRange<unsigned short>
Elem::node_index_range() const
{
  return {0, cast_int<unsigned short>(this->n_nodes())};
}



inline
IntRange<unsigned short>
Elem::edge_index_range() const
{
  return {0, cast_int<unsigned short>(this->n_edges())};
}



inline
IntRange<unsigned short>
Elem::face_index_range() const
{
  return {0, cast_int<unsigned short>(this->n_faces())};
}



inline
IntRange<unsigned short>
Elem::side_index_range() const
{
  return {0, cast_int<unsigned short>(this->n_sides())};
}




inline
std::unique_ptr<const Elem> Elem::side_ptr (unsigned int i) const
{
  // Call the non-const version of this function, return the result as
  // a std::unique_ptr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  return me->side_ptr(i);
}



inline
void
Elem::side_ptr (std::unique_ptr<const Elem> & elem,
                const unsigned int i) const
{
  // Hand off to the non-const version of this function
  Elem * me = const_cast<Elem *>(this);
  std::unique_ptr<Elem> e {const_cast<Elem *>(elem.release())};
  me->side_ptr(e, i);
  elem = std::move(e);
}



inline
std::unique_ptr<const Elem>
Elem::build_side_ptr (const unsigned int i) const
{
  // Call the non-const version of this function, return the result as
  // a std::unique_ptr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  return me->build_side_ptr(i);
}



inline
void
Elem::build_side_ptr (std::unique_ptr<const Elem> & elem,
                      const unsigned int i) const
{
  // Hand off to the non-const version of this function
  Elem * me = const_cast<Elem *>(this);
  std::unique_ptr<Elem> e {const_cast<Elem *>(elem.release())};
  me->build_side_ptr(e, i);
  elem = std::move(e);
}



template <typename Sideclass, typename Subclass>
inline
std::unique_ptr<Elem>
Elem::simple_build_side_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> face = std::make_unique<Sideclass>();
  for (auto n : face->node_index_range())
    face->set_node(n, this->node_ptr(Subclass::side_nodes_map[i][n]));

  face->set_interior_parent(this);
  face->inherit_data_from(*this);

  return face;
}



template <typename Subclass>
inline
void
Elem::simple_build_side_ptr (std::unique_ptr<Elem> & side,
                             const unsigned int i,
                             ElemType sidetype)
{
  libmesh_assert_less (i, this->n_sides());

  if (!side.get() || side->type() != sidetype)
    {
      Subclass & real_me = cast_ref<Subclass&>(*this);
      side = real_me.Subclass::build_side_ptr(i);
    }
  else
    {
      side->set_interior_parent(this);
      side->inherit_data_from(*this);
      for (auto n : side->node_index_range())
        side->set_node(n, this->node_ptr(Subclass::side_nodes_map[i][n]));
    }
}



template <typename Subclass, typename Mapclass>
inline
void
Elem::simple_side_ptr (std::unique_ptr<Elem> & side,
                       const unsigned int i,
                       ElemType sidetype)
{
  libmesh_assert_less (i, this->n_sides());

  if (!side.get() || side->type() != sidetype)
    {
      Subclass & real_me = cast_ref<Subclass&>(*this);
      side = real_me.Subclass::side_ptr(i);
    }
  else
    {
      side->subdomain_id() = this->subdomain_id();

      for (auto n : side->node_index_range())
        side->set_node(n, this->node_ptr(Mapclass::side_nodes_map[i][n]));
    }
}



inline
std::unique_ptr<const Elem>
Elem::build_edge_ptr (const unsigned int i) const
{
  // Call the non-const version of this function, return the result as
  // a std::unique_ptr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  return me->build_edge_ptr(i);
}



inline
void
Elem::build_edge_ptr (std::unique_ptr<const Elem> & elem,
                      const unsigned int i) const
{
  // Hand off to the non-const version of this function
  Elem * me = const_cast<Elem *>(this);
  std::unique_ptr<Elem> e {const_cast<Elem *>(elem.release())};
  me->build_edge_ptr(e, i);
  elem = std::move(e);
}


template <typename Edgeclass, typename Subclass>
inline
std::unique_ptr<Elem>
Elem::simple_build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  std::unique_ptr<Elem> edge = std::make_unique<Edgeclass>();

  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(Subclass::edge_nodes_map[i][n]));

  edge->set_interior_parent(this);
  edge->inherit_data_from(*this);

  return edge;
}




template <typename Subclass>
inline
void
Elem::simple_build_edge_ptr (std::unique_ptr<Elem> & edge,
                             const unsigned int i,
                             ElemType edgetype)
{
  libmesh_assert_less (i, this->n_edges());

  if (!edge.get() || edge->type() != edgetype)
    {
      Subclass & real_me = cast_ref<Subclass&>(*this);
      edge = real_me.Subclass::build_edge_ptr(i);
    }
  else
    {
      edge->inherit_data_from(*this);
      for (auto n : edge->node_index_range())
        edge->set_node(n, this->node_ptr(Subclass::edge_nodes_map[i][n]));
    }
}



inline
bool Elem::on_boundary () const
{
  // By convention, the element is on the boundary
  // if it has a nullptr neighbor.
  return this->has_neighbor(nullptr);
}



inline
unsigned int Elem::which_neighbor_am_i (const Elem * e) const
{
  libmesh_assert(e);

  const Elem * eparent = e;

  while (eparent->level() > this->level())
    {
      eparent = eparent->parent();
      libmesh_assert(eparent);
    }

  for (auto s : make_range(this->n_sides()))
    if (this->neighbor_ptr(s) == eparent)
      return s;

  return libMesh::invalid_uint;
}



inline
bool Elem::active() const
{
#ifdef LIBMESH_ENABLE_AMR
  if ((this->refinement_flag() == INACTIVE) ||
      (this->refinement_flag() == COARSEN_INACTIVE))
    return false;
  else
    return true;
#else
  return true;
#endif
}





inline
bool Elem::subactive() const
{
#ifdef LIBMESH_ENABLE_AMR
  if (this->active())
    return false;
  if (!this->has_children())
    return true;
  for (const Elem * my_ancestor = this->parent();
       my_ancestor != nullptr;
       my_ancestor = my_ancestor->parent())
    if (my_ancestor->active())
      return true;
#endif

  return false;
}



inline
bool Elem::has_children() const
{
#ifdef LIBMESH_ENABLE_AMR
  if (!_children)
    return false;
  else
    return true;
#else
  return false;
#endif
}


inline
bool Elem::has_ancestor_children() const
{
#ifdef LIBMESH_ENABLE_AMR
  if (!_children)
    return false;
  else
    for (auto & c : child_ref_range())
      if (c.has_children())
        return true;
#endif
  return false;
}



inline
bool Elem::is_ancestor_of(const Elem *
#ifdef LIBMESH_ENABLE_AMR
                          descendant
#endif
                          ) const
{
#ifdef LIBMESH_ENABLE_AMR
  const Elem * e = descendant;
  while (e)
    {
      if (this == e)
        return true;
      e = e->parent();
    }
#endif
  return false;
}



inline
const Elem * Elem::parent () const
{
  return _elemlinks[0];
}



inline
Elem * Elem::parent ()
{
  return _elemlinks[0];
}



inline
void Elem::set_parent (Elem * p)
{
  // We no longer support using parent() as interior_parent()
  libmesh_assert_equal_to(this->dim(), p ? p->dim() : this->dim());
  _elemlinks[0] = p;
}



inline
const Elem * Elem::top_parent () const
{
  const Elem * tp = this;

  // Keep getting the element's parent
  // until that parent is at level-0
  while (tp->parent() != nullptr)
    tp = tp->parent();

  libmesh_assert(tp);
  libmesh_assert_equal_to (tp->level(), 0);

  return tp;
}



inline
unsigned int Elem::level() const
{
#ifdef LIBMESH_ENABLE_AMR

  // if I don't have a parent I was
  // created directly from file
  // or by the user, so I am a
  // level-0 element
  if (this->parent() == nullptr)
    return 0;

  // if the parent and this element are of different
  // dimensionality we are at the same level as
  // the parent (e.g. we are the 2D side of a
  // 3D element)
  if (this->dim() != this->parent()->dim())
    return this->parent()->level();

  // otherwise we are at a level one
  // higher than our parent
  return (this->parent()->level() + 1);

#else

  // Without AMR all elements are
  // at level 0.
  return 0;

#endif
}



inline
unsigned int Elem::p_level() const
{
#ifdef LIBMESH_ENABLE_AMR
  return _p_level;
#else
  return 0;
#endif
}



inline
ElemMappingType Elem::mapping_type () const
{
  return static_cast<ElemMappingType>(_map_type);
}



inline
void Elem::set_mapping_type(const ElemMappingType type)
{
  _map_type = cast_int<unsigned char>(type);
}



inline
unsigned char Elem::mapping_data () const
{
  return _map_data;
}



inline
void Elem::set_mapping_data(const unsigned char data)
{
  _map_data = data;
}



#ifdef LIBMESH_ENABLE_AMR

inline
const Elem * Elem::raw_child_ptr (unsigned int i) const
{
  if (!_children)
    return nullptr;

  return _children[i];
}

inline
const Elem * Elem::child_ptr (unsigned int i) const
{
  libmesh_assert(_children);
  libmesh_assert(_children[i]);

  return _children[i];
}

inline
Elem * Elem::child_ptr (unsigned int i)
{
  libmesh_assert(_children);
  libmesh_assert(_children[i]);

  return _children[i];
}


inline
void Elem::set_child (unsigned int c, Elem * elem)
{
  libmesh_assert (this->has_children());

  _children[c] = elem;
}



inline
unsigned int Elem::which_child_am_i (const Elem * e) const
{
  libmesh_assert(e);
  libmesh_assert (this->has_children());

  unsigned int nc = this->n_children();
  for (unsigned int c=0; c != nc; c++)
    if (this->child_ptr(c) == e)
      return c;

  libmesh_error_msg("ERROR:  which_child_am_i() was called with a non-child!");

  return libMesh::invalid_uint;
}



inline
Elem::RefinementState Elem::refinement_flag () const
{
  return static_cast<RefinementState>(_rflag);
}



inline
void Elem::set_refinement_flag(RefinementState rflag)
{
  _rflag = cast_int<unsigned char>(rflag);
}



inline
Elem::RefinementState Elem::p_refinement_flag () const
{
  return static_cast<RefinementState>(_pflag);
}



inline
void Elem::set_p_refinement_flag(RefinementState pflag)
{
  if (this->p_level() == 0)
    libmesh_assert_not_equal_to
      (pflag, Elem::JUST_REFINED);

  _pflag = cast_int<unsigned char>(pflag);
}



inline
unsigned int Elem::max_descendant_p_level () const
{
  // This is undefined for subactive elements,
  // which have no active descendants
  libmesh_assert (!this->subactive());
  if (this->active())
    return this->p_level();

  unsigned int max_p_level = _p_level;
  for (auto & c : child_ref_range())
    max_p_level = std::max(max_p_level,
                           c.max_descendant_p_level());
  return max_p_level;
}



inline
void Elem::hack_p_level(unsigned int p)
{
  if (p == 0)
    libmesh_assert_not_equal_to
      (this->p_refinement_flag(), Elem::JUST_REFINED);

  _p_level = cast_int<unsigned char>(p);
}


inline
void Elem::hack_p_level_and_refinement_flag (unsigned int p,
                                             RefinementState pflag)
{
  _pflag = cast_int<unsigned char>(pflag);
  this->hack_p_level(p);
}

#endif // ifdef LIBMESH_ENABLE_AMR


inline
void Elem::orient(BoundaryInfo * boundary_info)
{
  if (this->is_flipped())
    this->flip(boundary_info);
}


inline
dof_id_type Elem::compute_key (dof_id_type n0)
{
  return n0;
}



inline
dof_id_type Elem::compute_key (dof_id_type n0,
                               dof_id_type n1)
{
  // Order the two so that n0 < n1
  if (n0 > n1) std::swap (n0, n1);

  return Utility::hashword2(n0, n1);
}



inline
dof_id_type Elem::compute_key (dof_id_type n0,
                               dof_id_type n1,
                               dof_id_type n2)
{
  std::array<dof_id_type, 3> array = {{n0, n1, n2}};
  std::sort(array.begin(), array.end());
  return Utility::hashword(array);
}



inline
dof_id_type Elem::compute_key (dof_id_type n0,
                               dof_id_type n1,
                               dof_id_type n2,
                               dof_id_type n3)
{
  std::array<dof_id_type, 4> array = {{n0, n1, n2, n3}};
  std::sort(array.begin(), array.end());
  return Utility::hashword(array);
}



inline
void Elem::inherit_data_from (const Elem & src)
{
  this->set_mapping_type(src.mapping_type());
  this->set_mapping_data(src.mapping_data());
  this->subdomain_id() = src.subdomain_id();
  this->processor_id(src.processor_id());
#ifdef LIBMESH_ENABLE_AMR
  this->set_p_level(src.p_level());
#endif
}



/**
 * The definition of the protected nested SideIter class.
 */
class Elem::SideIter
{
public:
  // Constructor with arguments.
  SideIter(const unsigned int side_number,
           Elem * parent)
    : _side(),
      _side_ptr(nullptr),
      _parent(parent),
      _side_number(side_number)
  {}


  // Empty constructor.
  SideIter()
    : _side(),
      _side_ptr(nullptr),
      _parent(nullptr),
      _side_number(libMesh::invalid_uint)
  {}


  // Copy constructor
  SideIter(const SideIter & other)
    : _side(),
      _side_ptr(nullptr),
      _parent(other._parent),
      _side_number(other._side_number)
  {}


  // op=
  SideIter & operator=(const SideIter & other)
  {
    this->_parent      = other._parent;
    this->_side_number = other._side_number;
    return *this;
  }

  // unary op*
  Elem *& operator*() const
  {
    // Set the std::unique_ptr
    this->_update_side_ptr();

    // Return a reference to _side_ptr
    return this->_side_ptr;
  }

  // op++
  SideIter & operator++()
  {
    ++_side_number;
    return *this;
  }

  // op==  Two side iterators are equal if they have
  // the same side number and the same parent element.
  bool operator == (const SideIter & other) const
  {
    return (this->_side_number == other._side_number &&
            this->_parent      == other._parent);
  }


  // Consults the parent Elem to determine if the side
  // is a boundary side.  Note: currently side N is a
  // boundary side if neighbor N is nullptr.  Be careful,
  // this could possibly change in the future?
  bool side_on_boundary() const
  {
    return this->_parent->neighbor_ptr(_side_number) == nullptr;
  }

private:
  // Update the _side pointer by building the correct side.
  // This has to be called before dereferencing.
  void _update_side_ptr() const
  {
    // Construct new side, store in std::unique_ptr
    this->_side = this->_parent->build_side_ptr(this->_side_number);

    // Also set our internal naked pointer.  Memory is still owned
    // by the std::unique_ptr.
    this->_side_ptr = _side.get();
  }

  // std::unique_ptr to the actual side, handles memory management for
  // the sides which are created during the course of iteration.
  mutable std::unique_ptr<Elem> _side;

  // Raw pointer needed to facilitate passing back to the user a
  // reference to a non-temporary raw pointer in order to conform to
  // the variant_filter_iterator interface.  It points to the same
  // thing the std::unique_ptr "_side" above holds.  What happens if the user
  // calls delete on the pointer passed back?  Well, this is an issue
  // which is not addressed by the iterators in libMesh.  Basically it
  // is a bad idea to ever call delete on an iterator from the library.
  mutable Elem * _side_ptr;

  // Pointer to the parent Elem class which generated this iterator
  Elem * _parent;

  // A counter variable which keeps track of the side number
  unsigned int _side_number;
};






// Private implementation functions in the Elem class for the side iterators.
// They have to come after the definition of the SideIter class.
inline
Elem::SideIter Elem::_first_side()
{
  return SideIter(0, this);
}



inline
Elem::SideIter Elem::_last_side()
{
  return SideIter(this->n_neighbors(), this);
}




/**
 * The definition of the struct used for iterating over sides.
 */
struct
Elem::side_iterator : variant_filter_iterator<Elem::Predicate, Elem *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  side_iterator (const IterType & d,
                 const IterType & e,
                 const PredType & p ) :
    variant_filter_iterator<Elem::Predicate, Elem *>(d,e,p) {}
};



inline
SimpleRange<Elem::NeighborPtrIter> Elem::neighbor_ptr_range()
{
  return {_elemlinks+1, _elemlinks + 1 + this->n_neighbors()};
}


inline
SimpleRange<Elem::ConstNeighborPtrIter> Elem::neighbor_ptr_range() const
{
  return {_elemlinks+1, _elemlinks + 1 + this->n_neighbors()};
}

} // namespace libMesh


// Helper function for default caches in Elem subclasses

#define LIBMESH_ENABLE_TOPOLOGY_CACHES                                  \
  virtual                                                               \
  std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char>>>>> & \
  _get_bracketing_node_cache() const override                   \
  {                                                                     \
    static std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char>>>>> c; \
    return c;                                                           \
  }                                                                     \
                                                                        \
  virtual                                                               \
  std::vector<std::vector<std::vector<signed char>>> &                  \
  _get_parent_indices_cache() const override                    \
  {                                                                     \
    static std::vector<std::vector<std::vector<signed char>>> c;        \
    return c;                                                           \
  }






#endif // LIBMESH_ELEM_H
