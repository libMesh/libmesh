// $Id: elem.h,v 1.6 2004-07-22 19:54:16 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __elem_h__
#define __elem_h__

// C++ includes
#include <algorithm>
#include <vector>

// Local includes
#include "libmesh_common.h"
#include "dof_object.h"
#include "reference_counted_object.h"
#include "node.h"
#include "enum_elem_type.h"
#include "enum_elem_quality.h"
#include "enum_order.h"
#include "enum_io_package.h"
#include "auto_ptr.h"


// Forward declarations
class MeshBase;
class MeshRefinement;
class MeshCommunication;
class MetisPartitioner;
class Elem;



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
 * \author Benjamin S. Kirk, 2002-2003
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
  Elem (const unsigned int n_nodes=0,
	const unsigned int n_sides=0,
	const Elem* parent=NULL);

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
   * @returns the global id number of local \p Node \p i.
   */
  unsigned int node (const unsigned int i) const;

  /**
   * @returns the pointer to local \p Node \p i.
   */
  Node* get_node (const unsigned int i) const;

  /**
   * @returns the pointer to local \p Node \p i as a writeable reference.
   */
  Node* & set_node (const unsigned int i);
  
  /**
   * @returns the subdomain that this element belongs to.
   * To conserve space this is stored as an unsigned char.
   */
  unsigned char subdomain_id () const;
  
  /**
   * @returns the subdomain that this element belongs to as a
   * writeable reference.
   */
  unsigned char & set_subdomain_id ();

  /**
   * @returns an id assocated with this element.  The id is not
   * guaranteed to be unique, but it should be close.  The id
   * is thus useful, for example, as a key in a hash table
   * data structure.
   */
  unsigned int key () const;

  /**
   * @returns an id associated with the \p s side of this element.
   * The id is not necessariy unique, but should be close.  This is
   * particularly useful in the \p MeshBase::find_neighbors() routine.
   */
  virtual unsigned int key (const unsigned int s) const = 0;
  
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

  /**
   * Assigns \p n as the \f$ i^{th} \f$ neighbor.
   */
  void set_neighbor (const unsigned int i, Elem* n);

  /**
   * @returns \p true if the element \p elem in question is a neighbor
   * of this element, \p false otherwise.
   */
  bool is_neighbor (const Elem* elem) const;
  
  /**
   * @returns \p true if this element has a side coincident
   * with a boundary (indicated by a \p NULL neighbor), \p false
   * otherwise.
   */
  bool on_boundary () const;
  
  /**
   * This function tells you which neighbor you \p (e) are.
   * What is returned is the index of the side _in_the_neighbor_
   * which you share. 
   */
  unsigned int which_neighbor_am_i(const Elem *e) const; 

  /**
   * Returns the connectivity for this element in a specific
   * format, which is specified by the IOPackage tag.  This
   * method supercedes the tecplot_connectivity(...) and vtk_connectivity(...)
   * routines.
   */
  virtual void connectivity(const unsigned int sc,
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const = 0;

  /**
   * Writes the element connectivity for various IO packages
   * to the passed ostream "out".  Not virtual, since it is
   * implemented in the base class.  This function supercedes the
   * write_tecplot_connectivity(...) and write_ucd_connectivity(...)
   * routines.
   */
  void write_connectivity (std::ostream& out,
			   const IOPackage iop) const;

  /**
   * @returns the VTK element type of the sc-th sub-element.
   */
  virtual unsigned int vtk_element_type (const unsigned int sc) const = 0;

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
   * @returns the number of nodes this element contains. 
   */
  virtual unsigned int n_nodes () const = 0;

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
   * @returns the number of children this element has that
   * share side \p s
   */
  virtual unsigned int n_children_per_side (const unsigned int) const = 0;
  
  /**
   * @returns the number of sub-elements this element may be broken
   * down into for visualization purposes.  For example, this returns
   * 1 for a linear triangle, 4 for a quadratic (6-noded) triangle, etc...
   */
  virtual unsigned int n_sub_elem () const = 0;

  /**
   * @returns an element coincident with side \p i.  This method returns
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
   */
  virtual AutoPtr<Elem> build_side (const unsigned int i) const = 0;

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
   * This method is overloadable since some derived elements 
   * might want to use shortcuts to compute their centroid.
   */
  virtual Real hmin () const;
  
  /**
   * @returns the maximum vertex separation for the element.
   * This method is overloadable since some derived elements 
   * might want to use shortcuts to compute their centroid.
   */
  virtual Real hmax () const;
  
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
  { error(); return std::make_pair(0.,0.); }
  
  /**
   * @returns true if the point p is contained in this element, 
   * false otherwise.
   */
  virtual bool contains_point (const Point& p) const;

  /**
   * @returns \p true if the element is active (i.e. has no children),
   * \p false  otherwise. Note that it suffices to check the first
   * child only. Always returns \p true if AMR is disabled. 
   */
  bool active () const;
  
  /**
   * @returns a pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement, i.e. was read from file.
   */
  const Elem* parent () const;

  /**
   * @returns a pointer to the element's top-most (i.e. level-0) parent.
   * Returns \p this if this is a level-0 element, this element's parent
   * if this is a level-1 element, this element's grandparent if this is
   * a level-2 element, etc...
   */
  const Elem* top_parent () const;
  
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
   * @returns the refinement level of the current element.  If the
   * element's parent is \p NULL then by convention it is at
   * level 0, otherwise it is simply at one level greater than
   * its parent.
   */
  unsigned int level () const;
  
#ifdef ENABLE_AMR

  /**
   * Useful ENUM describing the refinement state of
   * an element.
   */
  enum RefinementState { COARSEN = 0,
			 DO_NOTHING,
			 REFINE,
			 JUST_REFINED };
  
  /**
   * @returns a pointer to the \f$ i^{th} \f$ child for this element.
   * Returns \p NULL  if this element has no children, i.e. is active.
   */
  Elem* child (const unsigned int i) const;

  /**
   * Fills the vector \p family with the children of this element,
   * recursively.  So, calling this method on a twice-refined element
   * will give you the element itself, its direct children, and their
   * children, etc...  When the optional parameter \p reset is
   * true then the vector will be cleared before the element and its
   * descendants are added.
   */
  void family_tree (std::vector<const Elem*>& family,
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
   * Returns the value of the refinement flag for the element.
   */
  RefinementState refinement_flag () const;

  /**
   * Sets the value of the refinement flag for the element.
   */     
  void set_refinement_flag (const RefinementState rflag);

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
   * The non-const begin and end accessor functions.
   */
  std::pair<Elem**, Elem**> neighbors_begin ();
  std::pair<Elem**, Elem**> neighbors_end ();

  /**
   * The const begin and end accessor functions. 
   */
  std::pair<const Elem**, const Elem**> neighbors_begin () const;
  std::pair<const Elem**, const Elem**> neighbors_end () const;
  
#endif



#ifdef ENABLE_INFINITE_ELEMENTS

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
  virtual Point origin () const { error(); return Point(); }

#endif


 protected:

  
  /**
   * Build an element of type \p type.  Be careful when using this
   * function: it allocates memory, so you must remember to delete
   * the pointer when you are done with it.  The user should not call this,
   * so it is protected.
   */
  static Elem* build (const ElemType type,
		      const Elem* p=NULL);


  //-------------------------------------------------------
  // These methods compute has keys from the specified
  // global node numbers
  //
  /**
   * Compute a key from the specified nodes.
   */
  static unsigned int compute_key (unsigned int n0);

  /**
   * Compute a key from the specified nodes.
   */
  static unsigned int compute_key (unsigned int n0,
				   unsigned int n1);

  /**
   * Compute a key from the specified nodes.
   */
  static unsigned int compute_key (unsigned int n0,
				   unsigned int n1,
				   unsigned int n2);

  /**
   * Compute a key from the specified nodes.
   */
  static unsigned int compute_key (unsigned int n0,
				   unsigned int n1,
				   unsigned int n2,
				   unsigned int n3);
  //-------------------------------------------------------


  
  /**
   * Replaces this element with \p NULL for all of
   * its neighbors.  This is useful when deleting an
   * element.
   */
  void nullify_neighbors ();
  
#ifdef ENABLE_AMR
  
  /**
   * Matrix that transforms the parents nodes into the children's
   * nodes
   */
  virtual float embedding_matrix (const unsigned int i,
				  const unsigned int j,
				  const unsigned int k) const = 0;
  
  /**
   * Matrix that allows children to inherit boundary conditions.
   */
  virtual unsigned int side_children_matrix (const unsigned int i,
					     const unsigned int j) const = 0;

#endif
  
  /**
   * Pointers to the nodes we are conneted to.
   */
  Node** _nodes;

  /**
   * The subdomain to which this element belongs.
   */
  unsigned char _sbd_id;

  /**
   * Pointers to this element's neighbors.
   */
  Elem** _neighbors;
 
  /**
   * A pointer to this element's parent.  This is constant
   * since there is no way it should be able to change!.
   */
  const Elem* const _parent;

#ifdef ENABLE_AMR
  
  /**
   * Pointers to this element's children.
   */
  Elem** _children;

  /**
   * Refinement flag. This is stored as an unsigned char
   * to save space.
   */
  unsigned char _rflag;
  //RefinementState _rflag;
  
#endif

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
  friend class XdrInterface;      // (Elem::build)
  friend class MeshRefinement;    // (Elem::nullify_neighbors)
  friend class Mesh;              // (Elem::build in mesh_exodus_support)
  friend class MeshBase;          // (Elem::build in mesh_base_modification)
  friend class MeshCommunication; // (Elem::build in mesh_communication)
};




// ------------------------------------------------------------
// Elem class member functions
inline
Elem::Elem(const unsigned int nn,
	   const unsigned int ns,
	   const Elem* p) :
  _parent(p)
{
  assert (nn);

  this->set_subdomain_id() = 0;
  this->set_processor_id() = 0;

  // Initialize the nodes data structure
  {
    _nodes = NULL;

    if (nn != 0)
      {
	_nodes = new Node*[nn]; 
	
	for (unsigned int n=0; n<nn; n++)
	  _nodes[n] = NULL;
      }
  }
  
  // Initialize the neighbors data structure
  {
    _neighbors = NULL;

    if (ns != 0)
      {
	_neighbors = new Elem*[ns]; 
	
	for (unsigned int n=0; n<ns; n++)
	  _neighbors[n] = NULL;
      }
  }

  // Optionally initialize data from the parent
  if (this->parent() != NULL)
    {
      this->set_subdomain_id() = this->parent()->subdomain_id();
      this->set_processor_id() = this->parent()->processor_id();
    }  

#ifdef ENABLE_AMR
  
  _children = NULL;

  this->set_refinement_flag(Elem::DO_NOTHING);

#endif  
}



inline
Elem::~Elem() 
{
  // Delete my node storage
  if (_nodes != NULL)
    delete [] _nodes;
  _nodes = NULL;

  // Delete my neighbor storage
  if (_neighbors != NULL)
    delete [] _neighbors;
  _neighbors = NULL;

#ifdef ENABLE_AMR

  // Delete my children's storage
  if (_children != NULL)
    delete [] _children;
  _children = NULL;
  
#endif
}



inline
const Point & Elem::point (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);

  return *_nodes[i];
}



inline
Point & Elem::point (const unsigned int i)
{
  assert (i < this->n_nodes());

  return *_nodes[i];
}



inline
unsigned int Elem::node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);

  return _nodes[i]->id();
}



inline
Node* Elem::get_node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);

  return _nodes[i];
}



inline
Node* & Elem::set_node (const unsigned int i)
{
  assert (i < this->n_nodes());

  return _nodes[i];
}



inline
unsigned char Elem::subdomain_id () const
{
  return _sbd_id;
}



inline
unsigned char & Elem::set_subdomain_id ()
{
  return _sbd_id;
}



inline
Elem* Elem::neighbor (const unsigned int i) const
{
  assert (i < this->n_neighbors());

  return _neighbors[i];
}



inline
void Elem::set_neighbor (const unsigned int i, Elem* n)
{
  assert (i < this->n_neighbors());
  
  _neighbors[i] = n;
}



inline
bool Elem::is_neighbor (const Elem* elem) const
{
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor(n) == elem)
      return true;

  return false;
}



inline
bool Elem::on_boundary () const
{
  // By convention, the element is on the boundary
  // if it has a NULL neighbor.
  return this->is_neighbor(NULL);
}



inline
unsigned int Elem::which_neighbor_am_i (const Elem* e) const
{
  assert (e != NULL);
  
  for (unsigned int s=0; s<this->n_neighbors(); s++)
    if (this->neighbor(s) == e)
      return s;
    

  std::cerr << "ERROR:  Elements are not neighbors!" 
	    << std::endl;

  error();

  return libMesh::invalid_uint;
}



inline
bool Elem::active() const
{
#ifdef ENABLE_AMR
  
  if (_children == NULL)
    return true;
  else
    return false;
  
#else
  
  return true;
  
#endif
}




inline
const Elem* Elem::parent () const
{
  return _parent;
}



inline
const Elem* Elem::top_parent () const
{
  const Elem* tp = this;

  // Keep getting the element's parent
  // until that parent is at level-0
  while (tp->parent() != NULL)
    tp = tp->parent();
  
  assert (tp != NULL);
  assert (tp->level() == 0);

  return tp;  
}



inline
unsigned int Elem::level() const
{
#ifdef ENABLE_AMR

  // if I don't have a parent I was
  // created directly from file
  // or by the user, so I am a
  // level-0 element
  if (this->parent() == NULL)
    return 0;

  // otherwise we are at a level one
  // higher than our parent
  return (this->parent()->level() + 1);

#else

  // Without AMR all elements are
  // at level 0.
  return 0;
  
#endif
}


#ifdef ENABLE_AMR

inline
Elem* Elem::child (const unsigned int i) const
{
  assert (_children    != NULL);
  assert (_children[i] != NULL);
  
  return _children[i];
}



inline
Elem::RefinementState Elem::refinement_flag () const
{
  return static_cast<RefinementState>(_rflag);
}



inline
void Elem::set_refinement_flag(RefinementState rflag)
{
  if (rflag != static_cast<RefinementState>(static_cast<unsigned char>(rflag)))
    {
      std::cerr << "ERROR: unsigned int too small to hold Elem::_rflag!"
		<< std::endl
		<< "Recompile with Elem:_flag set to something bigger!"
		<< std::endl;
      error();
    }

  _rflag = rflag;
}



inline
std::pair<Elem**, Elem**> Elem::neighbors_begin ()
{
  return std::make_pair (&_neighbors[0],
			 &_neighbors[this->n_neighbors()]);
}



inline
std::pair<Elem**, Elem**> Elem::neighbors_end ()
{
  return std::make_pair (&_neighbors[this->n_neighbors()],
			 &_neighbors[this->n_neighbors()]);
}



inline
std::pair<const Elem**, const Elem**> Elem::neighbors_begin () const
{
  return std::make_pair (const_cast<const Elem**>(&_neighbors[0]),
			 const_cast<const Elem**>(&_neighbors[this->n_neighbors()]));
}



inline
std::pair<const Elem**, const Elem**> Elem::neighbors_end () const
{
  return std::make_pair (const_cast<const Elem**>(&_neighbors[this->n_neighbors()]),
			 const_cast<const Elem**>(&_neighbors[this->n_neighbors()]));
}


#endif /* ifdef ENABLE_AMR */


inline
unsigned int Elem::compute_key (unsigned int n0)
{
  return n0;
}



inline
unsigned int Elem::compute_key (unsigned int n0,
				unsigned int n1)
{
  // big prime number
  const unsigned int bp = 65449;
  
  // Order the two so that n0 < n1
  if (n0 > n1) std::swap (n0, n1);

  return (n0%bp + (n1<<5)%bp);  
}



inline
unsigned int Elem::compute_key (unsigned int n0,
				unsigned int n1,
				unsigned int n2)
{
  // big prime number
  const unsigned int bp = 65449;

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

  assert ((n0 < n1) && (n1 < n2));

  
  return (n0%bp + (n1<<5)%bp + (n2<<10)%bp);
}



inline
unsigned int Elem::compute_key (unsigned int n0,
				unsigned int n1,
				unsigned int n2,
				unsigned int n3)
{
  // big prime number
  const unsigned int bp = 65449;

  // Step 1
  if (n0 > n1) std::swap (n0, n1);

  // Step 2
  if (n2 > n3) std::swap (n2, n3);

  // Step 3
  if (n0 > n2) std::swap (n0, n2);

  // Step 4
  if (n1 > n3) std::swap (n1, n3);

  // Finally step 5
  if (n1 > n2) std::swap (n1, n2);

  assert ((n0 < n1) && (n1 < n2) && (n2 < n3));
  
  return (n0%bp + (n1<<5)%bp + (n2<<10)%bp + (n3<<15)%bp);
}
				


#endif // end #ifndef __elem_h__
