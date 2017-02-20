// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_elem_quality.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/multi_predicates.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/hashword.h" // Used in compute_key() functions

// C++ includes
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <limits.h> // CHAR_BIT
#include <set>
#include <vector>

namespace libMesh
{

// Forward declarations
class MeshBase;
class MeshRefinement;
class Elem;
#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
#endif

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
 * \e second-order nodes.
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
   * Destructor.  Frees all the memory associated with the element.
   */
  virtual ~Elem();

  /**
   * @returns the \p Point associated with local \p Node \p i.
   */
  const Point & point (const unsigned int i) const;

  /**
   * @returns the \p Point associated with local \p Node \p i
   * as a writeable reference.
   */
  Point & point (const unsigned int i);

  /**
   * @returns the \p Point associated with local \p Node \p i,
   * in master element rather than physical coordinates.
   */
  virtual Point master_point (const unsigned int i) const = 0;

  /**
   * @returns the global id number of local \p Node \p i.
   */
  dof_id_type node_id (const unsigned int i) const;

  /**
   * @returns the global id number of local \p Node \p i.
   *
   * This method is deprecated; use the less ambiguously named
   * node_id() instead.
   */
  dof_id_type node (const unsigned int i) const;

  /**
   * @returns the local id number of global \p Node id \p i,
   * or \p invalid_uint if Node id \p i is not local.
   */
  unsigned int local_node (const dof_id_type i) const;

  /**
   * @returns the local index for the \p Node pointer \p node_ptr,
   * or \p invalid_id if \p node_ptr is not a local node.
   */
  unsigned int get_node_index (const Node * node_ptr) const;

  /**
   * @returns a pointer to an array of local node pointers.
   */
  const Node * const * get_nodes () const;

  /**
   * @returns a const pointer to local \p Node \p i.
   */
  const Node * node_ptr (const unsigned int i) const;

  /**
   * @returns a non-const pointer to local \p Node \p i.
   */
  Node * node_ptr (const unsigned int i);

  /**
   * @returns a const reference to local \p Node \p i.
   */
  const Node & node_ref (const unsigned int i) const;

  /**
   * @returns a writable reference to local \p Node \p i.
   */
  Node & node_ref (const unsigned int i);

  /**
   * @returns the pointer to local \p Node \p i.
   *
   * This method is deprecated.  Use the more consistently and less
   * confusingly named node_ptr() instead.
   */
  Node * get_node (const unsigned int i) const;

  /**
   * @returns the pointer to local \p Node \p i as a writeable reference.
   */
  virtual Node * & set_node (const unsigned int i);

  /**
   * @returns the subdomain that this element belongs to.
   */
  subdomain_id_type subdomain_id () const;

  /**
   * @returns the subdomain that this element belongs to as a
   * writeable reference.
   */
  subdomain_id_type & subdomain_id ();

  /**
   * A static integral constant representing an invalid subdomain id.
   * See also DofObject::{invalid_id, invalid_unique_id, invalid_processor_id}.
   *
   * Note 1: we don't use the static_cast(-1) trick here since
   * subdomain_id_type is sometimes a *signed* integer for
   * compatibility reasons (see libmesh/id_types.h).
   *
   * Note 2: Normally you can declare static const integral types
   * directly in the header file (C++ standard, 9.4.2/4) but
   * std::numeric_limits<T>::max() is not considered a "constant
   * expression".  This one is therefore defined in elem.C.
   * http://stackoverflow.com/questions/2738435/using-numeric-limitsmax-in-constant-expressions
   */
  static const subdomain_id_type invalid_subdomain_id;

  /**
   * @returns a pointer to the "reference element" associated
   * with this element.  The reference element is the image of this
   * element in reference parametric space. Importantly, it is *not*
   * an actual element in the mesh, but rather a Singleton-type
   * object, so for example all \p Quad4 elements share the same
   * \p reference_elem().
   */
  const Elem * reference_elem () const;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessarily unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const = 0;

  /**
   * @returns an id associated with the global node ids of this
   * element.  The id is not necessarily unique, but should be
   * close. Uses the same hash as the key(s) function, so for example
   * if "tri3" is side 0 of "tet4", then tri3->key()==tet4->key(0).
   */
  virtual dof_id_type key () const;

  /**
   * @returns true if two elements are identical, false otherwise.
   * This is true if the elements are connected to identical global
   * nodes, regardless of how those nodes might be numbered local
   * to the elements.
   */
  bool operator == (const Elem & rhs) const;

  /**
   * @returns a const pointer to the \f$ i^{th} \f$ neighbor of this element.
   * If \p MeshBase::find_neighbors() has not been called this
   * simply returns \p NULL.  If \p MeshBase::find_neighbors()
   * has been called and this returns \p NULL then the side is on
   * a boundary of the domain.
   */
  const Elem * neighbor_ptr (unsigned int i) const;

  /**
   * @returns a non-const pointer to the \f$ i^{th} \f$ neighbor of this element.
   */
  Elem * neighbor_ptr (unsigned int i);

  /**
   * This function is deprecated.  Use the more specifically named and
   * const-correct neighbor_ptr() functions instead.
   */
  Elem * neighbor (const unsigned int i) const;


#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * @returns a pointer to the \f$ i^{th} \f$ neighbor of this element
   * for interior elements.  If an element is on a periodic
   * boundary, it will return a corresponding element on the opposite
   * side.
   */
  const Elem * topological_neighbor (const unsigned int i,
                                     const MeshBase & mesh,
                                     const PointLocatorBase & point_locator,
                                     const PeriodicBoundaries * pb) const;

  /**
   * @returns a writeable pointer to the \f$ i^{th} \f$ neighbor of
   * this element for interior elements.  If an element is on a
   * periodic boundary, it will return a corresponding element on the
   * opposite side.
   */
  Elem * topological_neighbor (const unsigned int i,
                               MeshBase & mesh,
                               const PointLocatorBase & point_locator,
                               const PeriodicBoundaries * pb);

  /**
   * @return \p true if the element \p elem in question is a neighbor or
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
   * @returns \p true if the element \p elem in question is a neighbor
   * of this element, \p false otherwise.
   */
  bool has_neighbor (const Elem * elem) const;

  /**
   * If the element \p elem in question is a neighbor
   * of a child of this element, this returns a pointer
   * to that child.  Otherwise it returns NULL.
   */
  Elem * child_neighbor (Elem * elem);

  /**
   * If the element \p elem in question is a neighbor
   * of a child of this element, this returns a pointer
   * to that child.  Otherwise it returns NULL.
   */
  const Elem * child_neighbor (const Elem * elem) const;

  /**
   * @returns \p true if this element has a side coincident
   * with a boundary (indicated by a \p NULL neighbor), \p false
   * otherwise.
   */
  bool on_boundary () const;

  /**
   * @returns \p true if this element is semilocal to the calling
   * processor, which must specify its rank.
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
   * Returns \p invalid_uint if \p e is not a side of \p this.
   */
  unsigned int which_side_am_i(const Elem * e) const;

  /**
   * This function returns true if a vertex of \p e is contained
   * in this element.
   */
  bool contains_vertex_of(const Elem * e) const;

  /**
   * This function returns true if an edge of \p e is contained in
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
   * Resets this element's neighbors' appropriate neighbor pointers
   * and its parent's and children's appropriate pointers
   * to point to the global remote_elem instead of this.
   * Used by the library before a remote element is deleted on the
   * local processor.
   */
  void make_links_to_me_remote ();

  /**
   * Resets the appropriate neighbor pointers of our nth neighbor (and
   * its descendants, if appropriate) to point to this Elem instead of
   * to the global remote_elem.  Used by the library when a formerly
   * remote element is being added to the local processor.
   */
  void make_links_to_me_local (unsigned int n);

  /**
   * Returns true if this element is remote, false otherwise.  A
   * remote element (see \p RemoteElem) is a syntactic convenience --
   * it is a placeholder for an element which exists on some other
   * processor.  Local elements are required to have valid neighbors,
   * and these ghost elements may have remote neighbors for data
   * structure consistency.  The use of remote elements helps ensure
   * that any element we may access has a NULL neighbor only if it
   * lies on the physical boundary of the domain.
   */
  virtual bool is_remote () const
  { return false; }

  /**
   * Returns the connectivity for this element in a specific
   * format, which is specified by the IOPackage tag.  This
   * method supercedes the tecplot_connectivity(...) and vtk_connectivity(...)
   * routines.
   */
  virtual void connectivity(const unsigned int sc,
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const = 0;

  /**
   * Writes the element connectivity for various IO packages
   * to the passed ostream "out".  Not virtual, since it is
   * implemented in the base class.  This function supercedes the
   * write_tecplot_connectivity(...) and write_ucd_connectivity(...)
   * routines.
   */
  void write_connectivity (std::ostream & out,
                           const IOPackage iop) const;

  /**
   * @returns the type of element that has been derived from this
   * base class.
   */
  virtual ElemType type () const = 0;

  /**
   * @returns the dimensionality of the object.
   */
  virtual unsigned int dim () const = 0;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of nodes in the element.
   */
  static const unsigned int type_to_n_nodes_map[INVALID_ELEM];

  /**
   * @returns the number of nodes this element contains.
   */
  virtual unsigned int n_nodes () const = 0;

  /**
   * @returns the number of nodes the given child of this element
   * contains.  Except in odd cases like pyramid refinement this will
   * be the same as the number of nodes in the parent element.
   */
  virtual unsigned int n_nodes_in_child (unsigned int /*c*/) const
  { return this->n_nodes(); }

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of sides on the element.
   */
  static const unsigned int type_to_n_sides_map[INVALID_ELEM];

  /**
   * @returns the number of sides the element that has been derived
   * from this class has. In 2D the number of sides is the number
   * of edges, in 3D the number of sides is the number of faces.
   */
  virtual unsigned int n_sides () const = 0;

  /**
   * @returns the number of neighbors the element that has been derived
   * from this class has.  By default only face (or edge in 2D)
   * neighbors are stored, so this method returns n_sides(),
   * however it may be overloaded in a derived class.
   */
  virtual unsigned int n_neighbors () const
  { return this->n_sides(); }

  /**
   * @returns the number of vertices the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_vertices () const = 0;

  /**
   * @returns the number of edges the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_edges () const = 0;

  /**
   * This array maps the integer representation of the \p ElemType enum
   * to the number of edges on the element.
   */
  static const unsigned int type_to_n_edges_map[INVALID_ELEM];

  /**
   * @returns the number of faces the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_faces () const = 0;

  /**
   * @returns the number of children the element that has been derived
   * from this class may have.
   */
  virtual unsigned int n_children () const = 0;

  /**
   * @returns true if the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const = 0;

  /**
   * @returns true if the specified child has a vertex at the
   * specified (child-local) node number.
   * Except in odd cases like pyramid refinement the child will have
   * the same local structure as the parent element.
   */
  virtual unsigned int is_vertex_on_child (unsigned int /*c*/,
                                           unsigned int n) const
  { return this->is_vertex(n); }

  /**
   * @returns true if this element has a vertex at the specified
   * (child-local) node number \p n of the specified child \p c.
   */
  virtual bool is_vertex_on_parent(unsigned int c,
                                   unsigned int n) const;

  /**
   * @returns true if the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const = 0;

  /**
   * @returns true if the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const = 0;

  /**
   * @returns true if the specified (local) node number is on the
   * specified side.
   */
  virtual bool is_node_on_side(const unsigned int n,
                               const unsigned int s) const = 0;

  /**
   * @returns true if the specified (local) node number is on the
   * specified edge.
   */
  virtual bool is_node_on_edge(const unsigned int n,
                               const unsigned int e) const = 0;

  /**
   * @returns true if the specified edge is on the specified side.
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const = 0;

  /**
   * @returns the side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const;

  /**
   * @returns the local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const;

  /**
   * @returns the number of sub-elements this element may be broken
   * down into for visualization purposes.  For example, this returns
   * 1 for a linear triangle, 4 for a quadratic (6-noded) triangle, etc...
   */
  virtual unsigned int n_sub_elem () const = 0;

  /**
   * @returns a proxy element coincident with side \p i.  This method returns
   * the _minimum_ element necessary to uniquely identify the side.  So,
   * for example, the side of a hexahedron is always returned as a 4-noded
   * quadrilateral, regardless of what type of hex you are dealing with.  If
   * you want the full-ordered face (i.e. a 9-noded quad face for a 27-noded
   * hexahedron) use the build_side method.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual UniquePtr<Elem> side_ptr (unsigned int i) = 0;
  UniquePtr<const Elem> side_ptr (unsigned int i) const;

  /**
   * @returns a proxy element coincident with side \p i.
   *
   * This method is deprecated and will be removed, since it
   * returns a non-const pointer to a side that could be used to
   * indirectly modify this.  Please use the the const-correct
   * side_ptr() function instead.
   */
  UniquePtr<Elem> side (const unsigned int i) const;

  /**
   * Creates an element coincident with side \p i. The element returned is
   * full-ordered, in contrast to the side method.  For example, calling
   * build_side_ptr(0) on a 20-noded hex will build a 8-noded quadrilateral
   * coincident with face 0 and pass back the pointer.
   *
   * A \p UniquePtr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   *
   * The second argument, which is true by default, specifies that a
   * "proxy" element (of type Side) will be returned.  This type of
   * return value is useful because it does not allocate additional
   * memory, and is usually sufficient for FE calculation purposes.
   * If you really need a full-ordered, non-proxy side object, call
   * this function with proxy=false.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual UniquePtr<Elem> build_side_ptr (const unsigned int i, bool proxy=true) = 0;
  UniquePtr<const Elem> build_side_ptr (const unsigned int i, bool proxy=true) const;

  /**
   * @returns a proxy element coincident with side \p i.
   *
   * This method is deprecated and will be removed, since it
   * returns a non-const pointer to a side that could be used to
   * indirectly modify this.  Please use the the const-correct
   * build_side_ptr() function instead.
   */
  UniquePtr<Elem> build_side (const unsigned int i, bool proxy=true) const;

  /**
   * Creates an element coincident with edge \p i. The element returned is
   * full-ordered.  For example, calling build_edge_ptr(0) on a 20-noded hex will
   * build a 3-noded edge coincident with edge 0 and pass back the pointer.
   *
   * A \p UniquePtr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   *
   * The const version of this function is non-virtual; it simply
   * calls the virtual non-const version and const_casts the return
   * type.
   */
  virtual UniquePtr<Elem> build_edge_ptr (const unsigned int i) = 0;
  UniquePtr<const Elem> build_edge_ptr (const unsigned int i) const;

  /**
   * Creates an element coincident with edge \p i.
   *
   * This is deprecated and will be removed, since it
   * returns a non-const pointer to an edge that could be used to
   * indirectly modify this.  Please use the the const-correct
   * build_edge_ptr() function instead.
   */
  UniquePtr<Elem> build_edge (const unsigned int i) const;

  /**
   * @returns the default approximation order for this element type.
   * This is the order that will be used to compute the map to the
   * reference element.
   */
  virtual Order default_order () const = 0;

  /**
   * @returns the centriod of the element. The centroid is
   * computed as the average of all the element vertices.
   * This method is overloadable since some derived elements
   * might want to use shortcuts to compute their centroid.
   */
  virtual Point centroid () const;

  /**
   * @returns the minimum vertex separation for the element.
   */
  virtual Real hmin () const;

  /**
   * @returns the maximum vertex separation for the element.
   */
  virtual Real hmax () const;

  /**
   * @return the (length/area/volume) of the geometric element.
   */
  virtual Real volume () const;

  /**
   * @return a bounding box (not necessarily the minimal bounding box)
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
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  virtual Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for the Elem based on
   * quality measure q.  These are the values suggested by the CUBIT
   * User's Manual.  Since this function can have no possible meaning
   * for an abstract Elem, it is an error in the base class.
   */
  virtual std::pair<Real,Real> qual_bounds (const ElemQuality) const
  { libmesh_not_implemented(); return std::make_pair(0.,0.); }

  /**
   * @returns true if the point p is contained in this element,
   * false otherwise.
   *
   * For linear elements, performs an initial tight bounding box check
   * (as an optimization step) and (if that passes) then uses the
   * user-defined tolerance "tol" in a call to inverse_map() to actually
   * test if the point is in the element.  For quadratic elements, the
   * bounding box optimization is skipped, and only the inverse_map()
   * steps are performed.
   *
   * Note that this routine should not be used to determine if a point
   * is merely "nearby" an element to within some tolerance. For that,
   * use Elem::close_to_point() instead.
   */
  virtual bool contains_point (const Point & p, Real tol=TOLERANCE) const;

  /**
   * @returns true if this element is "close" to the point p, where
   * "close" is determined by the tolerance tol.
   */
  virtual bool close_to_point(const Point & p, Real tol) const;

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
   * @returns true if the element map is definitely affine (i.e. the same at
   * every quadrature point) within numerical tolerances.
   */
  virtual bool has_affine_map () const { return false; }

  /**
   * @returns true if the Lagrange shape functions on this element
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
   * @returns \p true if the element is active (i.e. has no active
   * descendants), \p false otherwise. Note that it suffices to check the
   * first child only. Always returns \p true if AMR is disabled.
   */
  bool active () const;

  /**
   * @returns \p true if the element is an ancestor (i.e. has an
   * active child or ancestor child), \p false otherwise. Always
   * returns \p false if AMR is disabled.
   */
  bool ancestor () const;

  /**
   * @returns \p true if the element is subactive (i.e. has no active
   * descendants), \p false otherwise. Always returns \p false if AMR
   * is disabled.
   */
  bool subactive () const;

  /**
   * @returns \p true if the element has any children (active or not),
   * \p false otherwise. Always returns \p false if AMR is disabled.
   */
  bool has_children () const;

  /**
   * @returns \p true if the element has any descendants other than
   * its immediate children, \p false otherwise. Always returns \p
   * false if AMR is disabled.
   */
  bool has_ancestor_children () const;

  /**
   * @returns \p true if \p descendant is a child of \p this, or a
   * child of a child of \p this, etc. Always returns \p false if AMR
   * is disabled.
   */
  bool is_ancestor_of(const Elem * descendant) const;

  /**
   * @returns a const pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement.
   */
  const Elem * parent () const;

  /**
   * @returns a pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement.
   */
  Elem * parent ();

  /**
   * Sets the pointer to the element's parent.
   * Dangerous! Only use this if you know what you are doing!
   */
  void set_parent (Elem * p);

  /**
   * @returns a pointer to the element's top-most (i.e. level-0) parent.
   * Returns \p this if this is a level-0 element, this element's parent
   * if this is a level-1 element, this element's grandparent if this is
   * a level-2 element, etc...
   */
  const Elem * top_parent () const;

  /**
   * In some cases it is desireable to extract the boundary (or a subset thereof)
   * of a D-dimensional mesh as a (D-1)-dimensional manifold.  In this case
   * we may want to know the 'parent' element from which the manifold elements
   * were extracted.  We can easily do that for the level-0 manifold elements
   * by storing the D-dimensional parent.  This method provides access to that
   * element.
   *
   * This method is not safe to call if this->dim() == LIBMESH_DIM; in
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
   * @returns the distance between nodes n1 and n2.
   * Useful for computing the lengths of the sides of elements.
   */
  Real length (const unsigned int n1,
               const unsigned int n2) const;

  /**
   * @returns the number of adjacent vertices that uniquely define
   * the location of the \f$ n^{th} \f$ @e second-order node.  For linear
   * elements (\p default_order()==FIRST), this returns 0.
   * This method is useful when converting linear elements to quadratic
   * elements.  Note that \p n has to be greater than or equal to
   * \p this->n_vertices().
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const;

  /**
   * @returns the element-local number of the \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.  Note that the
   * return value is always less than \p this->n_vertices(), while \p
   * n has to be greater than or equal to \p this->n_vertices().  For
   * linear elements, this returns 0.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
                                                           const unsigned int v) const;

  /**
   * @returns a pair (c,v), where
   * c == child index, and
   * v == element-local index of the \p \f$ n^{th} \f$
   *      second-order node on the parent element.
   * For linear elements, this returns (0,0).
   *
   * Notes:
   * .) The return values are always less than \p this->n_children()
   *    and \p this->child_ptr(c)->n_vertices().
   * .) \p n has to be greater than or equal to \p this->n_vertices().
   * .) On refined second-order elements, the return value will
   *    satisfy \p this->node_ptr(n) == this->child_ptr(c)->node_ptr(v).
   */
  virtual std::pair<unsigned short int, unsigned short int>
  second_order_child_vertex (const unsigned int n) const;

  /**
   * @returns the element type of the associated second-order element,
   * e.g. when \p this is a \p TET4, then \p TET10 is returned.
   * Returns \p INVALID_ELEM for second-order or other elements that
   * cannot be converted into higher order equivalents.
   *
   * For some elements, there exist two second-order equivalents, e.g.
   * for \p Quad4 there is \p Quad8 and \p Quad9.  When the optional
   * \p full_ordered is \p true, then \p QUAD9 is returned.  When
   * \p full_ordered is \p false, then \p QUAD8 is returned.
   */
  static ElemType second_order_equivalent_type (const ElemType et,
                                                const bool full_ordered=true);

  /**
   * @returns the element type of the associated first-order element,
   * e.g. when \p this is a \p TET10, then \p TET4 is returned.  Returns
   * \p INVALID_ELEM for first-order or other elements that cannot be
   * converted into lower order equivalents.
   */
  static ElemType first_order_equivalent_type (const ElemType et);


  /**
   * @returns the refinement level of the current element.  If the
   * element's parent is \p NULL then by convention it is at
   * level 0, otherwise it is simply at one level greater than
   * its parent.
   */
  unsigned int level () const;

  /**
   * Returns the value of the p refinement level of an active
   * element, or the minimum value of the p refinement levels
   * of an ancestor element's descendants.
   */
  unsigned int p_level () const;

  /**
   * @returns true if the specified child is on the specified side.
   */
  virtual bool is_child_on_side(const unsigned int c,
                                const unsigned int s) const = 0;

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
   * @returns a constant pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  const Elem * child_ptr (unsigned int i) const;

  /**
   * @returns a non-constant pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  Elem * child_ptr (unsigned int i);

  /**
   * This function is now deprecated, use the more accurately-named and
   * const correct child_ptr() function instead.
   */
  Elem * child (const unsigned int i) const;


private:
  /**
   * Sets the pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  void set_child (unsigned int c, Elem * elem);

public:
  /**
   * This function tells you which child \p e is.
   * I.e. if c = a->which_child_am_i(e); then
   * a->child_ptr(c) will be e.
   */
  unsigned int which_child_am_i(const Elem * e) const;

  /**
   * @returns true if the specified child is on the specified edge.
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
                    const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but also adds any subactive
   * descendants.
   */
  void total_family_tree (std::vector<const Elem *> & active_family,
                          const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds the active
   * children.  Can be thought of as removing all the inactive
   * elements from the vector created by \p family_tree, but is
   * implemented more efficiently.
   */
  void active_family_tree (std::vector<const Elem *> & active_family,
                           const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void family_tree_by_side (std::vector<const Elem *> & family,
                            const unsigned int side,
                            const bool reset=true) const;

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void active_family_tree_by_side (std::vector<const Elem *> & family,
                                   const unsigned int side,
                                   const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void family_tree_by_neighbor (std::vector<const Elem *> & family,
                                const Elem * neighbor,
                                const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p subneighbor.  Only applicable when
   * \p this->has_neighbor(neighbor) and
   * \p neighbor->is_ancestor(subneighbor)
   */
  void family_tree_by_subneighbor (std::vector<const Elem *> & family,
                                   const Elem * neighbor,
                                   const Elem * subneighbor,
                                   const bool reset=true) const;

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void active_family_tree_by_neighbor (std::vector<const Elem *> & family,
                                       const Elem * neighbor,
                                       const bool reset=true) const;

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
                                                   const bool reset=true) const;

  /**
   * Returns the value of the refinement flag for the element.
   */
  RefinementState refinement_flag () const;

  /**
   * Sets the value of the refinement flag for the element.
   */
  void set_refinement_flag (const RefinementState rflag);

  /**
   * Returns the value of the p-refinement flag for the element.
   */
  RefinementState p_refinement_flag () const;

  /**
   * Sets the value of the p-refinement flag for the element.
   */
  void set_p_refinement_flag (const RefinementState pflag);

  /**
   * Returns the maximum value of the p-refinement levels of
   * an ancestor element's descendants.
   */
  unsigned int max_descendant_p_level () const;

  /**
   * Returns the minimum p-refinement level of elements which are
   * descended from this element, and which share a side with the
   * active \p neighbor.
   */
  unsigned int min_p_level_by_neighbor (const Elem * neighbor,
                                        unsigned int current_min) const;

  /**
   * Returns the minimum new p-refinement level (i.e. after refinement
   * and coarsening is done) of elements which are descended from this
   * element and which share a side with the active \p neighbor.
   */
  unsigned int min_new_p_level_by_neighbor (const Elem * neighbor,
                                            unsigned int current_min) const;

  /**
   * Sets the value of the p-refinement level for the element.
   * Note that the maximum p-refinement level is currently 255.
   */
  void set_p_level (const unsigned int p);

  /**
   * Sets the value of the p-refinement level for the element
   * without altering the p-level of its ancestors
   */
  void hack_p_level (const unsigned int p);

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

#ifdef DEBUG
  /**
   * Checks for consistent neighbor links on this element.
   */
  void libmesh_assert_valid_neighbors() const;

  /**
   * Checks for a valid id and pointers to nodes with valid ids on
   * this element.
   */
  void libmesh_assert_valid_node_pointers() const;
#endif // DEBUG

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
   * @returns \p true if the element is an infinite element,
   * \p false otherwise.
   */
  virtual bool infinite () const = 0;

  /**
   * @returns true if the specified (local) node number is a
   * "mid-edge" node on an infinite element edge.  This is false for
   * all nodes on non-infinite elements, so we won't make it pure
   * virtual, to simplify their code.
   */
  virtual bool is_mid_infinite_edge_node(const unsigned int /* n */) const
  { libmesh_assert (!this->infinite()); return false; }

  /**
   * @returns the origin for an infinite element.  Currently, \e all
   * infinite elements used in a mesh share the same origin.  Overload
   * this in infinite element classes.
   */
  virtual Point origin () const { libmesh_not_implemented(); return Point(); }

#endif




  /**
   * Build an element of type \p type.  Since this method allocates
   * memory, the new \p Elem is returned in a \p UniquePtr<>.
   */
  static UniquePtr<Elem> build (const ElemType type,
                                Elem * p=libmesh_nullptr);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Returns the local node id on the parent which corresponds to node
   * \p n of child \p c, or returns invalid_uint if no such parent
   * node exists.
   */
  virtual unsigned int as_parent_node (unsigned int c,
                                       unsigned int n) const;

  /**
   * Returns all the pairs of nodes (indexed by local node id) which
   * should bracket node \p n of child \p c.
   */
  virtual
  const std::vector<std::pair<unsigned char, unsigned char> > &
  parent_bracketing_nodes(unsigned int c,
                          unsigned int n) const;

  /**
   * Returns all the pairs of nodes (indexed by global node id) which
   * should bracket node \p n of child \p c.
   */
  virtual
  const std::vector<std::pair<dof_id_type, dof_id_type> >
  bracketing_nodes(unsigned int c,
                   unsigned int n) const;


  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes.
   */
  virtual float embedding_matrix (const unsigned int child_num,
                                  const unsigned int child_node_num,
                                  const unsigned int parent_node_num) const = 0;

  /**
   * Some element types may use a different embedding matrix for
   * different elements.  But we may want to cache data based on that
   * matrix.  So we return a "version number" that can be used to
   * identify which matrix is in use.
   */
  virtual unsigned int embedding_matrix_version () const { return 0; }

#endif // LIBMESH_ENABLE_AMR


protected:

  /**
   * Compute a hash key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0);

  /**
   * Compute a hash key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1);

  /**
   * Compute a hash key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1,
                                  dof_id_type n2);

  /**
   * Compute a hash key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
                                  dof_id_type n1,
                                  dof_id_type n2,
                                  dof_id_type n3);


#ifdef LIBMESH_ENABLE_AMR

  /**
   * Elem subclasses which don't do their own bracketing node
   * calculations will need to supply a static cache, since the
   * default calculation is slow.
   */
  virtual
  std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char> > > > > &
  _get_bracketing_node_cache() const
  {
    static std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char> > > > > c;
    libmesh_error();
    return c;
  }

  /**
   * Elem subclasses which don't do their own child-to-parent node
   * calculations will need to supply a static cache, since the
   * default calculation is slow.
   */
  virtual
  std::vector<std::vector<std::vector<signed char> > > &
  _get_parent_indices_cache() const
  {
    static std::vector<std::vector<std::vector<signed char> > > c;
    libmesh_error();
    return c;
  }

#endif // LIBMESH_ENABLE_AMR

public:

  /**
   * Replaces this element with \p NULL for all of its neighbors.
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
   * Pointers to this element's children.
   */
  Elem ** _children;
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
};



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
#ifdef LIBMESH_ENABLE_AMR
  _children(libmesh_nullptr),
#endif
  _sbd_id(0)
#ifdef LIBMESH_ENABLE_AMR
  ,
  _rflag(Elem::DO_NOTHING),
  _pflag(Elem::DO_NOTHING),
  _p_level(0)
#endif
{
  this->processor_id() = DofObject::invalid_processor_id;

  // Initialize the nodes data structure
  if (_nodes)
    {
      for (unsigned int n=0; n<nn; n++)
        _nodes[n] = libmesh_nullptr;
    }

  // Initialize the neighbors/parent data structure
  // _elemlinks = new Elem *[ns+1];

  if (_elemlinks)
    {
      _elemlinks[0] = p;

      for (unsigned int n=1; n<ns+1; n++)
        _elemlinks[n] = libmesh_nullptr;
    }

  // Optionally initialize data from the parent
  if (this->parent() != libmesh_nullptr)
    {
      this->subdomain_id() = this->parent()->subdomain_id();
      this->processor_id() = this->parent()->processor_id();
    }

#ifdef LIBMESH_ENABLE_AMR
  if (this->parent())
    this->set_p_level(this->parent()->p_level());
#endif
}



inline
Elem::~Elem()
{
  // Deleting my parent/neighbor/nodes storage isn't necessary since it's
  // handled by the subclass

  // if (_nodes != libmesh_nullptr)
  //   delete [] _nodes;
  // _nodes = libmesh_nullptr;

  // delete [] _elemlinks;

#ifdef LIBMESH_ENABLE_AMR

  // Delete my children's storage
  if (_children != libmesh_nullptr)
    delete [] _children;
  _children = libmesh_nullptr;

#endif
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
dof_id_type Elem::node (const unsigned int i) const
{
  libmesh_deprecated();
  return this->node_id(i);
}



inline
unsigned int Elem::local_node (const dof_id_type i) const
{
  for (unsigned int n=0; n != this->n_nodes(); ++n)
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
Node * Elem::get_node (const unsigned int i) const
{
  // This const function has incorrectly returned a non-const pointer
  // for years.  Now that it is reimplemented in terms of the new
  // interface which does return a const pointer, we need to use a
  // const_cast to mimic the old (incorrect) behavior.  This function
  // is now deprecated and eventually will be removed entirely,
  // obviating the need for this ugly cast.
  libmesh_deprecated();
  return const_cast<Node *>(this->node_ptr(i));
}



inline
unsigned int Elem::get_node_index (const Node * node_ptr) const
{
  for (unsigned int n=0; n != this->n_nodes(); ++n)
    if (this->_nodes[n] == node_ptr)
      return n;

  return libMesh::invalid_uint;
}



inline
Node * & Elem::set_node (const unsigned int i)
{
  libmesh_assert_less (i, this->n_nodes());

  return _nodes[i];
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
Elem * Elem::neighbor (const unsigned int i) const
{
  // Support the deprecated interface by calling the new,
  // const-correct interface and casting the result to an Elem *.
  libmesh_deprecated();
  return const_cast<Elem *>(this->neighbor_ptr(i));
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
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor_ptr(n) == elem)
      return true;

  return false;
}



inline
Elem * Elem::child_neighbor (Elem * elem)
{
  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    if (elem->neighbor_ptr(n) &&
        elem->neighbor_ptr(n)->parent() == this)
      return elem->neighbor_ptr(n);

  return libmesh_nullptr;
}



inline
const Elem * Elem::child_neighbor (const Elem * elem) const
{
  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    if (elem->neighbor_ptr(n) &&
        elem->neighbor_ptr(n)->parent() == this)
      return elem->neighbor_ptr(n);

  return libmesh_nullptr;
}



inline
UniquePtr<const Elem> Elem::side_ptr (unsigned int i) const
{
  // Call the non-const version of this function, return the result as
  // a UniquePtr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  const Elem * s = const_cast<const Elem *>(me->side_ptr(i).release());
  return UniquePtr<const Elem>(s);
}



inline
UniquePtr<Elem> Elem::side (const unsigned int i) const
{
  // Call the const version of side_ptr(), and const_cast the result.
  libmesh_deprecated();
  Elem * s = const_cast<Elem *>(this->side_ptr(i).release());
  return UniquePtr<Elem>(s);
}



inline
UniquePtr<const Elem>
Elem::build_side_ptr (const unsigned int i, bool proxy) const
{
  // Call the non-const version of this function, return the result as
  // a UniquePtr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  const Elem * s = const_cast<const Elem *>(me->build_side_ptr(i, proxy).release());
  return UniquePtr<const Elem>(s);
}



inline
UniquePtr<Elem>
Elem::build_side (const unsigned int i, bool proxy) const
{
  // Call the const version of build_side_ptr(), and const_cast the result.
  libmesh_deprecated();
  Elem * s = const_cast<Elem *>(this->build_side_ptr(i, proxy).release());
  return UniquePtr<Elem>(s);
}



inline
UniquePtr<const Elem>
Elem::build_edge_ptr (const unsigned int i) const
{
  // Call the non-const version of this function, return the result as
  // a UniquePtr<const Elem>.
  Elem * me = const_cast<Elem *>(this);
  const Elem * e = const_cast<const Elem *>(me->build_edge_ptr(i).release());
  return UniquePtr<const Elem>(e);
}



inline
UniquePtr<Elem>
Elem::build_edge (const unsigned int i) const
{
  // Call the const version of build_edge_ptr(), and const_cast the result.
  libmesh_deprecated();
  Elem * e = const_cast<Elem *>(this->build_edge_ptr(i).release());
  return UniquePtr<Elem>(e);
}



inline
bool Elem::on_boundary () const
{
  // By convention, the element is on the boundary
  // if it has a NULL neighbor.
  return this->has_neighbor(libmesh_nullptr);
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

  for (unsigned int s=0; s<this->n_neighbors(); s++)
    if (this->neighbor_ptr(s) == eparent)
      return s;

  return libMesh::invalid_uint;
}



inline
unsigned int Elem::which_side_am_i (const Elem * e) const
{
  libmesh_assert(e);

  const unsigned int ns = this->n_sides();
  const unsigned int nn = this->n_nodes();

  const unsigned int en = e->n_nodes();

  // e might be on any side until proven otherwise
  std::vector<bool> might_be_side(ns, true);

  for (unsigned int i=0; i != en; ++i)
    {
      Point side_point = e->point(i);
      unsigned int local_node_id = libMesh::invalid_uint;

      // Look for a node of this that's contiguous with node i of e
      for (unsigned int j=0; j != nn; ++j)
        if (this->point(j) == side_point)
          local_node_id = j;

      // If a node of e isn't contiguous with some node of this, then
      // e isn't a side of this.
      if (local_node_id == libMesh::invalid_uint)
        return libMesh::invalid_uint;

      // If a node of e isn't contiguous with some node on side s of
      // this, then e isn't on side s.
      for (unsigned int s=0; s != ns; ++s)
        if (!this->is_node_on_side(local_node_id, s))
          might_be_side[s] = false;
    }

  for (unsigned int s=0; s != ns; ++s)
    if (might_be_side[s])
      {
#ifdef DEBUG
        for (unsigned int s2=s+1; s2 < ns; ++s2)
          libmesh_assert (!might_be_side[s2]);
#endif
        return s;
      }

  // Didn't find any matching side
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
       my_ancestor != libmesh_nullptr;
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
  if (_children == libmesh_nullptr)
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
  if (_children == libmesh_nullptr)
    return false;
  else
    for (unsigned int c=0; c != this->n_children(); c++)
      if (this->child_ptr(c)->has_children())
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
  _elemlinks[0] = p;
}



inline
const Elem * Elem::top_parent () const
{
  const Elem * tp = this;

  // Keep getting the element's parent
  // until that parent is at level-0
  while (tp->parent() != libmesh_nullptr)
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
  if (this->parent() == libmesh_nullptr)
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



#ifdef LIBMESH_ENABLE_AMR

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
Elem * Elem::child (const unsigned int i) const
{
  // Support the deprecated interface by calling the new,
  // const-correct interface and casting the result to an Elem *.
  libmesh_deprecated();
  return const_cast<Elem *>(this->child_ptr(i));
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

  for (unsigned int c=0; c<this->n_children(); c++)
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
  _rflag = cast_int<RefinementState>(rflag);
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
  for (unsigned int c=0; c != this->n_children(); c++)
    max_p_level = std::max(max_p_level,
                           this->child_ptr(c)->max_descendant_p_level());
  return max_p_level;
}



inline
void Elem::set_p_level(unsigned int p)
{
  // Maintain the parent's p level as the minimum of it's children
  if (this->parent() != libmesh_nullptr)
    {
      unsigned int parent_p_level = this->parent()->p_level();

      // If our new p level is less than our parents, our parents drops
      if (parent_p_level > p)
        {
          this->parent()->set_p_level(p);

          // And we should keep track of the drop, in case we need to
          // do a projection later.
          this->parent()->set_p_refinement_flag(Elem::JUST_COARSENED);
        }
      // If we are the lowest p level and it increases, so might
      // our parent's, but we have to check every other child to see
      else if (parent_p_level == _p_level && _p_level < p)
        {
          _p_level = cast_int<unsigned char>(p);
          parent_p_level = cast_int<unsigned char>(p);
          for (unsigned int c=0; c != this->parent()->n_children(); c++)
            parent_p_level = std::min(parent_p_level,
                                      this->parent()->child_ptr(c)->p_level());

          // When its children all have a higher p level, the parent's
          // should rise
          if (parent_p_level > this->parent()->p_level())
            {
              this->parent()->set_p_level(parent_p_level);

              // And we should keep track of the rise, in case we need to
              // do a projection later.
              this->parent()->set_p_refinement_flag(Elem::JUST_REFINED);
            }

          return;
        }
    }

  _p_level = cast_int<unsigned char>(p);
}



inline
void Elem::hack_p_level(unsigned int p)
{
  if (p == 0)
    libmesh_assert_not_equal_to
      (this->p_refinement_flag(), Elem::JUST_REFINED);

  _p_level = cast_int<unsigned char>(p);
}



#endif // ifdef LIBMESH_ENABLE_AMR


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
  // Order the numbers such that n0 < n1 < n2.
  // We'll do it in 3 steps like this:
  //
  //     n0         n1                n2
  //     min(n0,n1) max(n0,n1)        n2
  //     min(n0,n1) min(n2,max(n0,n1) max(n2,max(n0,n1)
  //           |\   /|                  |
  //           | \ / |                  |
  //           |  /  |                  |
  //           | /  \|                  |
  //  gb min= min   max              gb max

  // Step 1
  if (n0 > n1) std::swap (n0, n1);

  // Step 2
  if (n1 > n2) std::swap (n1, n2);

  // Step 3
  if (n0 > n1) std::swap (n0, n1);

  libmesh_assert ((n0 < n1) && (n1 < n2));

  dof_id_type array[3] = {n0, n1, n2};
  return Utility::hashword(array, 3);
}



inline
dof_id_type Elem::compute_key (dof_id_type n0,
                               dof_id_type n1,
                               dof_id_type n2,
                               dof_id_type n3)
{
  // Sort first
  // Step 1
  if (n0 > n1) std::swap (n0, n1);

  // Step 2
  if (n2 > n3) std::swap (n2, n3);

  // Step 3
  if (n0 > n2) std::swap (n0, n2);

  // Step 4
  if (n1 > n3) std::swap (n1, n3);

  // Finally sort step 5
  if (n1 > n2) std::swap (n1, n2);

  libmesh_assert ((n0 < n1) && (n1 < n2) && (n2 < n3));

  dof_id_type array[4] = {n0, n1, n2, n3};
  return Utility::hashword(array, 4);
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
      _side_ptr(libmesh_nullptr),
      _parent(parent),
      _side_number(side_number)
  {}


  // Empty constructor.
  SideIter()
    : _side(),
      _side_ptr(libmesh_nullptr),
      _parent(libmesh_nullptr),
      _side_number(libMesh::invalid_uint)
  {}


  // Copy constructor
  SideIter(const SideIter & other)
    : _side(),
      _side_ptr(libmesh_nullptr),
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
    // Set the UniquePtr
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
  // boundary side if nieghbor N is NULL.  Be careful,
  // this could possibly change in the future?
  bool side_on_boundary() const
  {
    return this->_parent->neighbor_ptr(_side_number) == libmesh_nullptr;
  }

private:
  // Update the _side pointer by building the correct side.
  // This has to be called before dereferencing.
  void _update_side_ptr() const
  {
    // Construct new side, store in UniquePtr
    this->_side = this->_parent->build_side_ptr(this->_side_number);

    // Also set our internal naked pointer.  Memory is still owned
    // by the UniquePtr.
    this->_side_ptr = _side.get();
  }

  // UniquePtr to the actual side, handles memory management for
  // the sides which are created during the course of iteration.
  mutable UniquePtr<Elem> _side;

  // Raw pointer needed to facilitate passing back to the user a
  // reference to a non-temporary raw pointer in order to conform to
  // the variant_filter_iterator interface.  It points to the same
  // thing the UniquePtr "_side" above holds.  What happens if the user
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


} // namespace libMesh


// Helper function for default caches in Elem subclases

#define LIBMESH_ENABLE_TOPOLOGY_CACHES                                  \
  virtual                                                               \
  std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char> > > > > & \
  _get_bracketing_node_cache() const libmesh_override                   \
  {                                                                     \
    static std::vector<std::vector<std::vector<std::vector<std::pair<unsigned char, unsigned char> > > > > c; \
    return c;                                                           \
  }                                                                     \
                                                                        \
  virtual                                                               \
  std::vector<std::vector<std::vector<signed char> > > &                \
  _get_parent_indices_cache() const libmesh_override                    \
  {                                                                     \
    static std::vector<std::vector<std::vector<signed char> > > c;      \
    return c;                                                           \
  }


#endif // LIBMESH_ELEM_H
