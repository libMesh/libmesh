// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
 * This is the base class from which all geometric entities
 * (elements) are derived.  The \p Elem class contains information
 * that every entity might need, such as its number of nodes and
 * pointers to the nodes to which it is connected.  This class
 * also provides virtual functions that will be overloaded by
 * derived classes.  These functions provide information such as
 * the number of sides the element has, who its neighbors are,
 * how many children it might have, and who they are.
 *
 * In an \p Elem becomes an \p Edge in 1D, a \p Face in 2D, and a \p
 * Cell in 3D.  An \p Elem is composed of a number of sides, which you
 * may access as \p Elem types in dimension \p D-1.  For example, a
 * concrete element type in 3D is a \p Hex8, which is a hexahedral. A
 * \p Hex8 has 6 sides, which are \p Faces.  You may access these
 * sides.
 *
 * An \p Elem is composed of a number of \p Node objects.  Some of
 * these nodes live at the vertices of the element, and others may
 * live on edges (and faces in 3D) or interior to the element.  The
 * number of vertices an element contains \p n_vertices() is
 * determined strictly by the type of geometric object it corresponds
 * to.  For example, a \p Tri is a type of \p Face that always
 * contains 3 vertices.  A \p Tri3 is a specific triangular element
 * type with three 3 nodes, all located at the vertices.  A \p Tri6 is
 * another triangular element with 6 nodes, 3 of which are located at
 * vertices and another 3 that live on the edges.
 * In all that follows, nodes that live either on edges, faces or the
 * interior are named @e second-order nodes.
 *
 * \author Benjamin S. Kirk, 2002-2007
 */

// ------------------------------------------------------------
// Elem class definition
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
	Elem* parent,
	Elem** elemlinkdata,
	Node** nodelinkdata);

 public:

  /**
   * Destructor.  Frees all the memory associated with the element.
   */
  virtual ~Elem();

  /**
   * @returns the \p Point associated with local \p Node \p i.
   */
  virtual const Point & point (const unsigned int i) const;

  /**
   * @returns the \p Point associated with local \p Node \p i
   * as a writeable reference.
   */
  virtual Point & point (const unsigned int i);

  /**
   * @returns the global id number of local \p Node \p i.
   */
  virtual dof_id_type node (const unsigned int i) const;

  /**
   * @returns the local id number of global \p Node id \p i,
   * or \p invalid_uint if Node id \p i is not local.
   */
  virtual unsigned int local_node (const dof_id_type i) const;

  /**
   * @returns the local index for the \p Node pointer \p node_ptr,
   * or \p invalid_id if \p node_ptr is not a local node.
   */
  unsigned int get_node_index (const Node* node_ptr) const;

  /**
   * @returns the pointer to local \p Node \p i.
   */
  virtual Node* get_node (const unsigned int i) const;

  /**
   * @returns the pointer to local \p Node \p i as a writeable reference.
   */
  virtual Node* & set_node (const unsigned int i);

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
   * @returns a pointer to the "reference element" associated
   * with this element.  The reference element is the image of this
   * element in reference parametric space. Importantly, it is *not*
   * an actual element in the mesh, but rather a Singleton-type
   * object, so for example all \p Quad4 elements share the same
   * \p reference_elem().
   */
  const Elem* reference_elem () const;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual dof_id_type key (const unsigned int s) const = 0;

  /**
   * @returns true if two elements are identical, false otherwise.
   * This is true if the elements are connected to identical global
   * nodes, regardless of how those nodes might be numbered local
   * to the elements.
   */
  bool operator == (const Elem& rhs) const;

  /**
   * @returns a pointer to the \f$ i^{th} \f$ neighbor of this element.
   * If \p MeshBase::find_neighbors() has not been called this
   * simply returns \p NULL.  If \p MeshBase::find_neighbors()
   * has been called and this returns \p NULL then the side is on
   * a boundary of the domain.
   */
  Elem* neighbor (const unsigned int i) const;

#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * @returns a pointer to the \f$ i^{th} \f$ neighbor of this element
   * for interior elements.  If an element is on a periodic
   * boundary, it will return a corresponding element on the opposite
   * side.
   */
  const Elem* topological_neighbor (const unsigned int i,
                                    const MeshBase& mesh,
                                    const PointLocatorBase& point_locator,
                                    const PeriodicBoundaries* pb) const;

  /**
   * @returns a writeable pointer to the \f$ i^{th} \f$ neighbor of
   * this element for interior elements.  If an element is on a
   * periodic boundary, it will return a corresponding element on the
   * opposite side.
   */
  Elem* topological_neighbor (const unsigned int i,
                              MeshBase& mesh,
                              const PointLocatorBase& point_locator,
                              const PeriodicBoundaries* pb);

  /**
   * @return \p true if the element \p elem in question is a neighbor or
   * topological neighbor of this element, \p false otherwise.
   */
  bool has_topological_neighbor (const Elem* elem,
                                 const MeshBase& mesh,
                                 const PointLocatorBase& point_locator,
                                 PeriodicBoundaries* pb) const;
#endif

  /**
   * Assigns \p n as the \f$ i^{th} \f$ neighbor.
   */
  void set_neighbor (const unsigned int i, Elem* n);

  /**
   * @returns \p true if the element \p elem in question is a neighbor
   * of this element, \p false otherwise.
   */
  bool has_neighbor (const Elem* elem) const;

  /**
   * If the element \p elem in question is a neighbor
   * of a child of this element, this returns a pointer
   * to that child.  Otherwise it returns NULL.
   */
  Elem* child_neighbor (Elem* elem) const;

  /**
   * If the element \p elem in question is a neighbor
   * of a child of this element, this returns a pointer
   * to that child.  Otherwise it returns NULL.
   */
  const Elem* child_neighbor (const Elem* elem) const;

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
   * This function tells you which neighbor you \p (e) are.
   * I.e. if s = a->which_neighbor_am_i(e); then
   * a->neighbor(s) will be an ancestor of e;
   */
  unsigned int which_neighbor_am_i(const Elem *e) const;

  /**
   * This function tells you which side the boundary element \p e is.
   * I.e. if e = a->build_side(s) or e = a->side(s); then
   * a->which_side_am_i(e) will be s.
   *
   * Returns \p invalid_uint if \p e is not a side of \p this
   */
  unsigned int which_side_am_i(const Elem *e) const;

  /**
   * This function returns true iff a vertex of e is contained
   * in this element
   */
  bool contains_vertex_of(const Elem *e) const;

  /**
   * This function returns true iff an edge of \p e is contained in
   * this element.  (Internally, this is done by checking whether at
   * least two vertices of \p e are contained in this element).
   */
  bool contains_edge_of(const Elem *e) const;

  /**
   * This function finds all elements (including this one) which
   * touch the current active element at the specified point, which
   * should be a point in the current element.
   */
  void find_point_neighbors(const Point &p,
                            std::set<const Elem *> &neighbor_set) const;

  /**
   * This function finds all elements (including this one) which
   * touch the current element at any point
   */
  void find_point_neighbors(std::set<const Elem *> &neighbor_set) const;

  /**
   * This function finds all active elements which touch the current
   * active element along the specified edge defined by the two points
   * \p p1 and \p p2
   */
  void find_edge_neighbors(const Point& p1,
                           const Point& p2,
			   std::set<const Elem *> &neighbor_set) const;

  /**
   * This function finds all active elements which touch the current
   * active element along any edge (more precisely, at at least two
   * points).
   */
  void find_edge_neighbors(std::set<const Elem *> &neighbor_set) const;

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
   * Returns true if this element is remote, false otherwise.
   * A remote element (see \p RemoteElem) is a syntactic convenience --
   * it is a placeholder for an element which exists on some other
   * processor.  Local elements are required to have valid neighbors,
   * and these ghost elements may have remote neighbors for data
   * structure consistency.  The use of remote elements helps assure
   * that any element we may access has a NULL neighbor if and only if
   * it lies on the physical boundary of the domain.
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
			    std::vector<dof_id_type>& conn) const = 0;

  /**
   * Writes the element connectivity for various IO packages
   * to the passed ostream "out".  Not virtual, since it is
   * implemented in the base class.  This function supercedes the
   * write_tecplot_connectivity(...) and write_ucd_connectivity(...)
   * routines.
   */
  void write_connectivity (std::ostream& out,
			   const IOPackage iop) const;

//   /**
//    * @returns the VTK element type of the sc-th sub-element.
//    */
//   virtual unsigned int vtk_element_type (const unsigned int sc) const = 0;

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
   * however it may be overloaded in a derived class
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
   * @returns true iff the specified (local) node number is a vertex.
   */
  virtual bool is_vertex(const unsigned int i) const = 0;

  /**
   * @returns true iff the specified (local) node number is an edge.
   */
  virtual bool is_edge(const unsigned int i) const = 0;

  /**
   * @returns true iff the specified (local) node number is a face.
   */
  virtual bool is_face(const unsigned int i) const = 0;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified side
   */
  virtual bool is_node_on_side(const unsigned int n,
			       const unsigned int s) const = 0;

  /*
   * @returns true iff the specified (local) node number is on the
   * specified edge
   */
  virtual bool is_node_on_edge(const unsigned int n,
			       const unsigned int e) const = 0;

  /*
   * @returns true iff the specified edge is on the specified side
   */
  virtual bool is_edge_on_side(const unsigned int e,
			       const unsigned int s) const = 0;

  /*
   * @returns the side number opposite to \p s (for a tensor product
   * element), or throws an error otherwise.
   */
  virtual unsigned int opposite_side(const unsigned int s) const;

  /*
   * @returns the local node number for the node opposite to node n
   * on side \p opposite_side(s) (for a tensor product element), or
   * throws an error otherwise.
   */
  virtual unsigned int opposite_node(const unsigned int n,
                                     const unsigned int s) const;

//   /**
//    * @returns the number of children this element has that
//    * share side \p s
//    */
//   virtual unsigned int n_children_per_side (const unsigned int) const = 0;

  /**
   * @returns the number of sub-elements this element may be broken
   * down into for visualization purposes.  For example, this returns
   * 1 for a linear triangle, 4 for a quadratic (6-noded) triangle, etc...
   */
  virtual unsigned int n_sub_elem () const = 0;

  /**
   * @returns a proxy element coincident with side \p i.  This method returns
   * the _minimum_ element necessary to uniquely identify the side.  So,
   * for example, the side of a hexahedral is always returned as a 4-noded
   * quadrilateral, regardless of what type of hex you are dealing with.  If
   * you want the full-ordered face (i.e. a 9-noded quad face for a 27-noded
   * hexahedral) use the build_side method.
   */
  virtual AutoPtr<Elem> side (const unsigned int i) const = 0;

  /**
   * Creates an element coincident with side \p i. The element returned is
   * full-ordered, in contrast to the side method.  For example, calling
   * build_side(0) on a 20-noded hex will build a 8-noded quadrilateral
   * coincident with face 0 and pass back the pointer.
   *
   * A \p AutoPtr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   *
   * The second argument, which is true by default, specifies that a
   * "proxy" element (of type Side) will be returned.  This type of
   * return value is useful because it does not allocate additional
   * memory, and is usually sufficient for FE calculation purposes.
   * If you really need a full-ordered, non-proxy side object, call
   * this function with proxy=false.
   */
  virtual AutoPtr<Elem> build_side (const unsigned int i,
				    bool proxy=true) const = 0;

  /**
   * Creates an element coincident with edge \p i. The element returned is
   * full-ordered.  For example, calling build_edge(0) on a 20-noded hex will
   * build a 3-noded edge coincident with edge 0 and pass back the pointer.
   *
   * A \p AutoPtr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   */
  virtual AutoPtr<Elem> build_edge (const unsigned int i) const = 0;

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
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  virtual Real quality (const ElemQuality q) const;

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   * Since this function can have no possible meaning
   * for an abstract Elem, it is an error.
   */
  virtual std::pair<Real,Real> qual_bounds (const ElemQuality) const
  { libmesh_error(); return std::make_pair(0.,0.); }

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
  virtual bool contains_point (const Point& p, Real tol=TOLERANCE) const;

  /**
   * @returns true if this element is "close" to the point p, where
   * "close" is determined by the tolerance tol.
   */
  virtual bool close_to_point(const Point& p, Real tol) const;

private:
  /**
   * Shared private implementation used by the contains_point()
   * and close_to_point() routines.  The box_tol tolerance is
   * used in the bounding box optimization, the map_tol tolerance is used
   * in the calls to inverse_map() and on_reference_element().
   */
  bool point_test(const Point& p, Real box_tol, Real map_tol) const;

public:
  /**
   * @returns true iff the element map is definitely affine (i.e. the same at
   * every quadrature point) within numerical tolerances
   */
  virtual bool has_affine_map () const { return false; }

  /**
   * @returns true iff the Lagrange shape functions on this element
   * are linear
   */
  virtual bool is_linear () const { return false; }

  /**
   * Prints relevant information about the element.
   */
  void print_info (std::ostream& os=libMesh::out) const;

  /**
   * Prints relevant information about the element to a string.
   */
  std::string get_info () const;

  /**
   * @returns \p true if the element is active (i.e. has no active
   * descendants), \p false  otherwise. Note that it suffices to check the
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
   * \p false  otherwise. Always returns \p false if AMR is disabled.
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
   * child of a child of \p this, etc.
   * Always returns \p false if AMR is disabled.
   */
  bool is_ancestor_of(const Elem *descendant) const;

  /**
   * @returns a const pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement, i.e. was read from file.
   */
  const Elem* parent () const;

  /**
   * @returns a pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement, i.e. was read from file.
   */
  Elem* parent ();

  /**
   * Sets the pointer to the element's parent.
   * Dangerous to use in high-level code.
   */
  void set_parent (Elem *p);

  /**
   * @returns a pointer to the element's top-most (i.e. level-0) parent.
   * Returns \p this if this is a level-0 element, this element's parent
   * if this is a level-1 element, this element's grandparent if this is
   * a level-2 element, etc...
   */
  const Elem* top_parent () const;

  /**
   * In some cases it is desireable to extract the boundary (or a subset thereof)
   * of a D-dimensional mesh as a (D-1)-dimensional manifold.  In this case
   * we may want to know the 'parent' element from which the manifold elements
   * were extracted.  We can easily do that for the level-0 manifold elements
   * by storing the D-dimensional parent.  This method provides access to that
   * element.
   */
  const Elem* interior_parent () const;

  /**
   * Sets the pointer to the element's interior_parent.
   * Dangerous to use in high-level code.
   */
  void set_interior_parent (Elem *p);

  /**
   * @returns the magnitude of the distance between nodes n1 and n2.
   * Useful for computing the lengths of the sides of elements.
   */
  Real length (const unsigned int n1,
	       const unsigned int n2) const;

  /**
   * @returns the number of adjacent vertices, that uniquely define
   * the location of the \f$ n^{th} \f$ @e second-order node.  For linear
   * elements ( \p default_order()==FIRST ), this returns 0.
   * This method is useful when converting linear elements to quadratic
   * elements.  Note that \p n has to be greater or equal
   * \p this->n_vertices().
   */
  virtual unsigned int n_second_order_adjacent_vertices (const unsigned int n) const;

  /**
   * @returns the element-local number of the  \f$ v^{th} \f$ vertex
   * that defines the \f$ n^{th} \f$ second-order node.  Note that
   * the return value is always less \p this->n_vertices(), while
   * \p n has to be greater or equal \p this->n_vertices().  For
   * linear elements this returns 0.
   */
  virtual unsigned short int second_order_adjacent_vertex (const unsigned int n,
							   const unsigned int v) const;

  /**
   * @returns the child number \p c and element-local index \p v of the
   * \f$ n^{th} \f$ second-order node on the parent element.  Note that
   * the return values are always less \p this->n_children() and
   * \p this->child(c)->n_vertices(), while \p n has to be greater or equal
   * to \p * this->n_vertices().  For linear elements this returns 0,0.
   * On refined second order elements, the return value will satisfy
   * \p this->get_node(n)==this->child(c)->get_node(v)
   */
  virtual std::pair<unsigned short int, unsigned short int>
	  second_order_child_vertex (const unsigned int n) const;

  /**
   * @returns the element type of the associated second-order element,
   * e.g. when \p this is a \p TET4, then \p TET10 is returned.  Returns
   * \p INVALID_ELEM for second order or other elements that should not
   * or cannot be converted into higher order equivalents.
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
   * \p INVALID_ELEM for first order or other elements that should not
   * or cannot be converted into lower order equivalents.
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
   * of an ancestor element's descendants
   */
  unsigned int p_level () const;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Useful ENUM describing the refinement state of
   * an element.
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
   * @returns a pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  Elem* child (const unsigned int i) const;

private:
  /**
   * @sets the pointer to the \f$ i^{th} \f$ child for this element.
   * Do not call if this element has no children, i.e. is active.
   */
  void set_child (unsigned int c, Elem *elem);

public:
  /**
   * This function tells you which child you \p (e) are.
   * I.e. if c = a->which_child_am_i(e); then
   * a->child(c) will be e;
   */
  unsigned int which_child_am_i(const Elem *e) const;

  /**
   * @returns true iff the specified child is on the
   * specified side
   */
  virtual bool is_child_on_side(const unsigned int c,
			        const unsigned int s) const = 0;

  /**
   * @returns true iff the specified child is on the
   * specified edge
   */
  virtual bool is_child_on_edge(const unsigned int c,
			        const unsigned int e) const;

  /**
   * Adds a child pointer to the array of children of this element.
   * If this is the first child to be added, this method allocates
   * memory in the parent's _children array, otherwise, it just sets
   * the pointer.
   */
  void add_child (Elem* elem);

  /**
   * Adds a new child pointer to the specified index in the array of
   * children of this element.  If this is the first child to be added,
   * this method allocates memory in the parent's _children array,
   * otherwise, it just sets the pointer.
   */
  void add_child (Elem* elem, unsigned int c);

  /**
   * Replaces the child pointer at the specified index in the array of
   * children of this element.
   */
  void replace_child (Elem* elem, unsigned int c);

  /**
   * Fills the vector \p family with the children of this element,
   * recursively.  So, calling this method on a twice-refined element
   * will give you the element itself, its direct children, and their
   * children, etc...  When the optional parameter \p reset is
   * true then the vector will be cleared before the element and its
   * descendants are added.
   *
   * The family tree only includes ancestor and active elements; for
   * subactive elements as well, use total_family_tree.
   */
  void family_tree (std::vector<const Elem*>& family,
		    const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but also adds any subactive
   * descendants.
   */
  void total_family_tree (std::vector<const Elem*>& active_family,
			  const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds the active
   * children.  Can be thought of as removing all the inactive
   * elements from the vector created by \p family_tree, but is
   * implemented more efficiently.
   */
  void active_family_tree (std::vector<const Elem*>& active_family,
			   const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void family_tree_by_side (std::vector<const Elem*>& family,
                            const unsigned int side,
                            const bool reset=true) const;

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p side.
   */
  void active_family_tree_by_side (std::vector<const Elem*>& family,
                                   const unsigned int side,
                                   const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void family_tree_by_neighbor (std::vector<const Elem*>& family,
		                const Elem *neighbor,
		                const bool reset=true) const;

  /**
   * Same as the \p family_tree() member, but only adds elements
   * which are next to \p subneighbor.  Only applicable when
   * \p this->has_neighbor(neighbor) and
   * \p neighbor->is_ancestor(subneighbor)
   */
  void family_tree_by_subneighbor (std::vector<const Elem*>& family,
		                   const Elem *neighbor,
		                   const Elem *subneighbor,
		                   const bool reset=true) const;

  /**
   * Same as the \p active_family_tree() member, but only adds elements
   * which are next to \p neighbor.
   */
  void active_family_tree_by_neighbor (std::vector<const Elem*>& family,
		                       const Elem *neighbor,
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
   * Returns the value of the p refinement flag for the element.
   */
  RefinementState p_refinement_flag () const;

  /**
   * Sets the value of the p refinement flag for the element.
   */
  void set_p_refinement_flag (const RefinementState pflag);

  /**
   * Returns the maximum value of the p refinement levels of
   * an ancestor element's descendants
   */
  unsigned int max_descendant_p_level () const;

  /**
   * Returns the minimum p refinement level of elements which
   * are descended from this and which share a side with the
   * active \p neighbor
   */
  unsigned int min_p_level_by_neighbor (const Elem* neighbor,
					unsigned int current_min) const;

  /**
   * Returns the minimum new p refinement level (i.e. after
   * refinement and coarsening is done) of elements which are
   * descended from this and which share a side with the
   * active \p neighbor
   */
  unsigned int min_new_p_level_by_neighbor (const Elem* neighbor,
					    unsigned int current_min) const;

  /**
   * Sets the value of the p refinement level for the element
   * Note that the maximum p refinement level is currently 255
   */
  void set_p_level (const unsigned int p);

  /**
   * Sets the value of the p refinement level for the element
   * without altering the p level of its ancestors
   */
  void hack_p_level (const unsigned int p);

  /**
   * Refine the element.
   */
  virtual void refine (MeshRefinement& mesh_refinement);

  /**
   * Coarsen the element.  This is not
   * virtual since it is the same for all
   * element types.
   */
  void coarsen ();

  /**
   * Contract an active element, i.e. remove pointers to any
   * subactive children.  This should only be called via
   * MeshRefinement::contract, which will also remove subactive
   * children from the mesh
   */
  void contract ();

#endif

#ifdef DEBUG
  /**
   * This function checks for consistent neighbor links at this
   * element.
   */
  void libmesh_assert_valid_neighbors() const;

  /**
   * This function checks for a valid id and for pointers to nodes
   * with valid ids at this element.
   */
  void libmesh_assert_valid_node_pointers() const;
#endif // DEBUG

protected:

  /**
   * The protected nested SideIter class is used to iterate over the
   * sides of this Elem.  It is a specially designed class since
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
  //typedef variant_filter_iterator<Elem*, Predicate> side_iterator;

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
   * \p false  otherwise.
   */
  virtual bool infinite () const = 0;

  /**
   * @returns the origin for an infinite element.  Currently,
   * @e all infinite elements used in a mesh share the same
   * origin.  Overload this in infinite element classes.
   * By default, issues an error, because returning the
   * all zero point would very likely lead to unexpected
   * behavior.
   */
  virtual Point origin () const { libmesh_error(); return Point(); }

#endif




  /**
   * Build an element of type \p type.  Since this method
   * allocates memory the new \p Elem is returned in a
   * \p AutoPtr<>
   */
  static AutoPtr<Elem> build (const ElemType type,
			      Elem* p=NULL);

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes
   */
  virtual float embedding_matrix (const unsigned int i,
				  const unsigned int j,
				  const unsigned int k) const = 0;

#endif


  //--------------------------------------------------------
  /**
   * Convenient way to communicate elements.  This struct
   * packes up an element so that it can easily be communicated through
   * an MPI array.
   *
   * \author Benjamin S. Kirk
   * \date 2008
   */
  class PackedElem;

  unsigned int packed_size() const;

 protected:

  //-------------------------------------------------------
  // These methods compute has keys from the specified
  // global node numbers
  //
  /**
   * Compute a key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0);

  /**
   * Compute a key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
				  dof_id_type n1);

  /**
   * Compute a key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
				  dof_id_type n1,
				  dof_id_type n2);

  /**
   * Compute a key from the specified nodes.
   */
  static dof_id_type compute_key (dof_id_type n0,
				  dof_id_type n1,
				  dof_id_type n2,
				  dof_id_type n3);
  //-------------------------------------------------------



  /**
   * Replaces this element with \p NULL for all of
   * its neighbors.  This is useful when deleting an
   * element.
   */
  void nullify_neighbors ();

  /**
   * Pointers to the nodes we are connected to.
   */
  Node** _nodes;

  /**
   * Pointers to this element's parent and neighbors, and for
   * lower-dimensional elements interior_parent.
   */
  Elem** _elemlinks;

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Pointers to this element's children.
   */
  Elem** _children;

  /**
   * h refinement flag. This is stored as an unsigned char
   * to save space.
   */
  unsigned char _rflag;
  //RefinementState _rflag;

  /**
   * p refinement flag. This is stored as an unsigned char
   * to save space.
   */
  unsigned char _pflag;
  //RefinementState _pflag;

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
   * The subdomain to which this element belongs.
   */
  subdomain_id_type _sbd_id;

  /**
   * Make the classes that need to access our build
   * member friends.  These classes do not really fit
   * the profile of what a "friend" should be, but
   * if we are going to protect the constructor and
   * the build method, there's no way around it.
   *
   * Do we *really* need to protect the build member?
   * It would seem that we are just getting around it
   * by using friends!
   */
  friend class MeshRefinement;    // (Elem::nullify_neighbors)
};

// ------------------------------------------------------------
// global Elem functions

inline
std::ostream& operator << (std::ostream& os, const Elem& e)
{
  e.print_info(os);
  return os;
}


// ------------------------------------------------------------
// Elem class member functions
inline
Elem::Elem(const unsigned int nn,
	   const unsigned int ns,
	   Elem* p,
           Elem** elemlinkdata,
           Node** nodelinkdata) :
  _nodes(nodelinkdata),
  _elemlinks(elemlinkdata),
#ifdef LIBMESH_ENABLE_AMR
  _children(NULL),
  _rflag(Elem::DO_NOTHING),
  _pflag(Elem::DO_NOTHING),
  _p_level(0),
#endif
  _sbd_id(0)
{
  this->processor_id() = DofObject::invalid_processor_id;

  // Initialize the nodes data structure
  if (_nodes)
    {
      for (unsigned int n=0; n<nn; n++)
	_nodes[n] = NULL;
    }

  // Initialize the neighbors/parent data structure
  // _elemlinks = new Elem*[ns+1];

  if (_elemlinks)
    {
      _elemlinks[0] = p;

      for (unsigned int n=1; n<ns+1; n++)
        _elemlinks[n] = NULL;
    }

  // Optionally initialize data from the parent
  if (this->parent() != NULL)
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

  // if (_nodes != NULL)
  //   delete [] _nodes;
  // _nodes = NULL;

  // delete [] _elemlinks;

#ifdef LIBMESH_ENABLE_AMR

  // Delete my children's storage
  if (_children != NULL)
    delete [] _children;
  _children = NULL;

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
dof_id_type Elem::node (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);
  libmesh_assert_not_equal_to (_nodes[i]->id(), Node::invalid_id);

  return _nodes[i]->id();
}



inline
unsigned int Elem::local_node (const dof_id_type i) const
{
  for (unsigned int n=0; n != this->n_nodes(); ++n)
    if (this->node(n) == i)
      return n;

  return libMesh::invalid_uint;
}



inline
Node* Elem::get_node (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_nodes());
  libmesh_assert(_nodes[i]);

  return _nodes[i];
}



inline
unsigned int Elem::get_node_index (const Node* node_ptr) const
{
  for (unsigned int n=0; n != this->n_nodes(); ++n)
    if (this->_nodes[n] == node_ptr)
      return n;

  return libMesh::invalid_uint;
}



inline
Node* & Elem::set_node (const unsigned int i)
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
Elem* Elem::neighbor (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_neighbors());

  return _elemlinks[i+1];
}



inline
void Elem::set_neighbor (const unsigned int i, Elem* n)
{
  libmesh_assert_less (i, this->n_neighbors());

  _elemlinks[i+1] = n;
}



inline
bool Elem::has_neighbor (const Elem* elem) const
{
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor(n) == elem)
      return true;

  return false;
}



inline
Elem* Elem::child_neighbor (Elem* elem) const
{
  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    if (elem->neighbor(n) &&
	elem->neighbor(n)->parent() == this)
      return elem->neighbor(n);

  return NULL;
}



inline
const Elem* Elem::child_neighbor (const Elem* elem) const
{
  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    if (elem->neighbor(n) &&
	elem->neighbor(n)->parent() == this)
      return elem->neighbor(n);

  return NULL;
}



inline
bool Elem::on_boundary () const
{
  // By convention, the element is on the boundary
  // if it has a NULL neighbor.
  return this->has_neighbor(NULL);
}



inline
unsigned int Elem::which_neighbor_am_i (const Elem* e) const
{
  libmesh_assert(e);

  const Elem* eparent = e;

  while (eparent->level() > this->level())
    {
      eparent = eparent->parent();
      libmesh_assert(eparent);
    }

  for (unsigned int s=0; s<this->n_neighbors(); s++)
    if (this->neighbor(s) == eparent)
      return s;

  return libMesh::invalid_uint;
}



inline
unsigned int Elem::which_side_am_i (const Elem* e) const
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
  for (const Elem* my_ancestor = this->parent();
       my_ancestor != NULL;
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
  if (_children == NULL)
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
  if (_children == NULL)
    return false;
  else
    for (unsigned int c=0; c != this->n_children(); c++)
      if (this->child(c)->has_children())
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
  const Elem *e = descendant;
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
const Elem* Elem::parent () const
{
  return _elemlinks[0];
}



inline
Elem* Elem::parent ()
{
  return _elemlinks[0];
}



inline
void Elem::set_parent (Elem *p)
{
  _elemlinks[0] = p;
}



inline
const Elem* Elem::top_parent () const
{
  const Elem* tp = this;

  // Keep getting the element's parent
  // until that parent is at level-0
  while (tp->parent() != NULL)
    tp = tp->parent();

  libmesh_assert(tp);
  libmesh_assert_equal_to (tp->level(), 0);

  return tp;
}



inline
const Elem* Elem::interior_parent () const
{
  // interior parents make no sense for full-dimensional elements.
  libmesh_assert_less (this->dim(), LIBMESH_DIM);

  // // and they [USED TO BE] only good for level-0 elements
  // if (this->level() != 0)
    // return this->parent()->interior_parent();

  // We store the interior_parent pointer after both the parent
  // neighbor and neighbor pointers
  Elem *interior_p = _elemlinks[1+this->n_sides()];

  // If we have an interior_parent, it had better be the
  // one-higher-dimensional interior element we are looking for.
  libmesh_assert (!interior_p ||
                  interior_p->dim() == (this->dim()+1));

  return interior_p;
}



inline
void Elem::set_interior_parent (Elem *p)
{
  // interior parents make no sense for full-dimensional elements.
  libmesh_assert_less (this->dim(), LIBMESH_DIM);

  // this had better be a one-higher-dimensional interior element
  libmesh_assert (!p ||
                  p->dim() == (this->dim()+1));

  _elemlinks[1+this->n_sides()] = p;
}



inline
unsigned int Elem::level() const
{
#ifdef LIBMESH_ENABLE_AMR

  // if I don't have a parent I was
  // created directly from file
  // or by the user, so I am a
  // level-0 element
  if (this->parent() == NULL)
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
Elem* Elem::child (const unsigned int i) const
{
  libmesh_assert(_children);
  libmesh_assert(_children[i]);

  return _children[i];
}



inline
void Elem::set_child (unsigned int c, Elem *elem)
{
  libmesh_assert (this->has_children());

  _children[c] = elem;
}



inline
unsigned int Elem::which_child_am_i (const Elem* e) const
{
  libmesh_assert(e);
  libmesh_assert (this->has_children());

  for (unsigned int c=0; c<this->n_children(); c++)
    if (this->child(c) == e)
      return c;

  libMesh::err << "ERROR:  which_child_am_i() was called with a non-child!"
	        << std::endl;

  libmesh_error();

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
  _rflag = libmesh_cast_int<RefinementState>(rflag);
}



inline
Elem::RefinementState Elem::p_refinement_flag () const
{
  return static_cast<RefinementState>(_pflag);
}



inline
void Elem::set_p_refinement_flag(RefinementState pflag)
{
  _pflag = libmesh_cast_int<unsigned char>(pflag);
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
			   this->child(c)->max_descendant_p_level());
  return max_p_level;
}



inline
void Elem::set_p_level(unsigned int p)
{
  // Maintain the parent's p level as the minimum of it's children
  if (this->parent() != NULL)
    {
      unsigned int parent_p_level = this->parent()->p_level();

      // If our new p level is less than our parents, our parents drops
      if (parent_p_level > p)
	{
          this->parent()->set_p_level(p);
	}
      // If we are the lowest p level and it increases, so might
      // our parent's, but we have to check every other child to see
      else if (parent_p_level == _p_level && _p_level < p)
	{
	  _p_level = libmesh_cast_int<unsigned char>(p);
	  parent_p_level = libmesh_cast_int<unsigned char>(p);
	  for (unsigned int c=0; c != this->parent()->n_children(); c++)
	    parent_p_level = std::min(parent_p_level,
				      this->parent()->child(c)->p_level());

	  if (parent_p_level != this->parent()->p_level())
	    this->parent()->set_p_level(parent_p_level);

	  return;
	}
    }

  _p_level = libmesh_cast_int<unsigned char>(p);
}



inline
void Elem::hack_p_level(unsigned int p)
{
  _p_level = libmesh_cast_int<unsigned char>(p);
}



#endif /* ifdef LIBMESH_ENABLE_AMR */


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



//-----------------------------------------------------------------
/**
 * Convenient way to communicate elements.  This class
 * packes up an element so that it can easily be communicated through
 * an MPI array.
 *
 * \author Benjamin S. Kirk
 * \date 2008
 */
class Elem::PackedElem
{
private:

  /**
   * Iterator pointing to the beginning of this packed element's index buffer.
   */
  const std::vector<largest_id_type>::const_iterator _buf_begin;

public:

  /**
   * Constructor.  Takes an input iterator pointing to the
   * beginning of the connectivity block for this element.
   */
  PackedElem (const std::vector<largest_id_type>::const_iterator _buf_in) :
    _buf_begin(_buf_in)
  {}

  /**
   * An \p Elem can be packed into an integer array which is
   * \p header_size + elem->n_nodes() in length.
   */
  static const unsigned int header_size; /* = 10 with AMR, 4 without */

  /**
   * For each element the serialization is of the form
   * [ level p_level r_flag p_r_flag etype processor_id subdomain_id
   *  self_ID parent_ID which_child node_0 node_1 ... node_n
   *  dof_object_buffer_1 ...]
   * We cannot use unsigned int because parent_ID can be negative
   */
  static void pack (std::vector<largest_id_type> &conn, const Elem* elem);

  /**
   * Unpacks this packed element.  Returns a pointer to the new element.
   * Takes a pointer to the parent, which is required unless this packed
   * element is at level 0.
   */
  Elem * unpack (MeshBase &mesh, Elem *parent = NULL) const;

  /**
   * \p return the level of this packed element.
   */
  unsigned int level () const
  {
    return static_cast<unsigned int>(*_buf_begin);
  }

  /**
   * \p return the p-level of this packed element.
   */
  unsigned int p_level () const
  {
    return static_cast<unsigned int>(*(_buf_begin+1));
  }

#ifdef LIBMESH_ENABLE_AMR
  /**
   * \p return the refinement state of this packed element.
   */
  Elem::RefinementState refinement_flag () const
  {
    // libmesh_assert_greater_equal (*(_buf_begin+2), 0);
    libmesh_assert_less (*(_buf_begin+2), INVALID_REFINEMENTSTATE);
    return static_cast<Elem::RefinementState>(*(_buf_begin+2));
  }

  /**
   * \p return the p-refinement state of this packed element.
   */
  Elem::RefinementState p_refinement_flag () const
  {
    // libmesh_assert_greater_equal (*(_buf_begin+3), 0);
    libmesh_assert_less (*(_buf_begin+3), INVALID_REFINEMENTSTATE);
    return static_cast<Elem::RefinementState>(*(_buf_begin+3));
  }
#endif // LIBMESH_ENABLE_AMR

  /**
   * \p return the element type of this packed element.
   */
  ElemType type () const
  {
    // libmesh_assert_greater_equal (*(_buf_begin+4), 0);
    libmesh_assert_less (*(_buf_begin+4), INVALID_ELEM);
    return static_cast<ElemType>(*(_buf_begin+4));
  }

  /**
   * \p return the processor id of this packed element.
   */
  processor_id_type processor_id () const
  {
    // libmesh_assert_greater_equal (*(_buf_begin+5), 0);
    libmesh_assert_less (static_cast<unsigned int>(*(_buf_begin+5)),
                   libMesh::n_processors() ||
                   static_cast<processor_id_type>(*(_buf_begin+5)) == DofObject::invalid_processor_id);
    return static_cast<processor_id_type>(*(_buf_begin+5));
  }

  /**
   * \p return the subdomain id of this packed element.
   */
  subdomain_id_type subdomain_id () const
  {
    return static_cast<subdomain_id_type>(*(_buf_begin+6));
  }

  /**
   * \p return the id of this packed element.
   */
  dof_id_type id () const
  {
    return static_cast<dof_id_type>(*(_buf_begin+7));
  }

  /**
   * \p return the parent id of this packed element.
   */
  int parent_id () const
  {
    return *(_buf_begin+8);
  }

  /**
   * \p return which child this packed element is.
   */
  unsigned int which_child_am_i () const
  {
    return static_cast<unsigned int>(*(_buf_begin+9));
  }

  /**
   * \p return the number of nodes in this packed element
   */
  unsigned int n_nodes () const
  {
    return Elem::type_to_n_nodes_map[this->type()];
  }

  /**
   * \p return the global index of the packed element's nth node.
   */
  unsigned int node (const unsigned int n) const
  {
    return static_cast<unsigned int>(*(_buf_begin+header_size+n));
  }

  /**
   * \p return the number of neighbors of the packed element
   */
  unsigned int n_neighbors() const
  {
    return Elem::type_to_n_sides_map[this->type()];
  }

  /**
   * \p return the global index of the packed element's nth neighbor
   */
  unsigned int neighbor (const unsigned int n) const
  {
    return static_cast<unsigned int>
      (*(_buf_begin + header_size + this->n_nodes() + n));
  }


  std::vector<largest_id_type>::const_iterator indices() const
  {
    return _buf_begin + header_size + this->n_nodes() +
           this->n_neighbors();
  }

  unsigned int packed_size() const
  {
    return this->header_size + this->n_nodes() +
           this->n_neighbors() +
           DofObject::unpackable_indexing_size(this->indices());
  }
}; // end class PackedElem

/**
 * The definition of the protected nested SideIter class.
 */
class Elem::SideIter
{
public:
  // Constructor with arguments.
  SideIter(const unsigned int side_number,
	   Elem* parent)
    : _side(),
      _side_ptr(NULL),
      _parent(parent),
      _side_number(side_number)
  {}


  // Empty constructor.
  SideIter()
    : _side(),
      _side_ptr(NULL),
      _parent(NULL),
      _side_number(libMesh::invalid_uint)
  {}


  // Copy constructor
  SideIter(const SideIter& other)
    : _side(),
      _side_ptr(NULL),
      _parent(other._parent),
      _side_number(other._side_number)
  {}


  // op=
  SideIter& operator=(const SideIter& other)
  {
    this->_parent      = other._parent;
    this->_side_number = other._side_number;
    return *this;
  }

  // unary op*
  Elem*& operator*() const
  {
    // Set the AutoPtr
    this->_update_side_ptr();

    // Return a reference to _side_ptr
    return this->_side_ptr;
  }

  // op++
  SideIter& operator++()
  {
    ++_side_number;
    return *this;
  }

  // op==  Two side iterators are equal if they have
  // the same side number and the same parent element.
  bool operator == (const SideIter& other) const
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
    return this->_parent->neighbor(_side_number) == NULL;
  }

private:
  // Update the _side pointer by building the correct side.
  // This has to be called before dereferencing.
  void _update_side_ptr() const
  {
    // Construct new side, store in AutoPtr
    this->_side = this->_parent->build_side(this->_side_number);

    // Also set our internal naked pointer.  Memory is still owned
    // by the AutoPtr.
    this->_side_ptr = _side.get();
  }

  // AutoPtr to the actual side, handles memory management for
  // the sides which are created during the course of iteration.
  mutable AutoPtr<Elem> _side;

  // Raw pointer needed to facilitate passing back to the user a
  // reference to a non-temporary raw pointer in order to conform to
  // the variant_filter_iterator interface.  It points to the same
  // thing the AutoPtr "_side" above holds.  What happens if the user
  // calls delete on the pointer passed back?  Well, this is an issue
  // which is not addressed by the iterators in libMesh.  Basically it
  // is a bad idea to ever call delete on an iterator from the library.
  mutable Elem* _side_ptr;

  // Pointer to the parent Elem class which generated this iterator
  Elem* _parent;

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
Elem::side_iterator :
variant_filter_iterator<Elem::Predicate,
			Elem*>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  side_iterator (const IterType& d,
		 const IterType& e,
		 const PredType& p ) :
    variant_filter_iterator<Elem::Predicate,
			    Elem*>(d,e,p) {}
};


inline
unsigned int
Elem::packed_size() const
{
  return PackedElem::header_size + this->n_nodes() +
         this->n_neighbors() + this->packed_indexing_size();
}


} // namespace libMesh


#endif // LIBMESH_ELEM_H
