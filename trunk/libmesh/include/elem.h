// $Id: elem.h,v 1.12 2003-02-13 00:16:47 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include <vector>

// Local includes
#include "mesh_common.h"
#include "reference_counted_object.h"
#include "node.h"
#include "enum_elem_type.h"
#include "enum_elem_quality.h"
#include "enum_order.h"
#include "auto_ptr.h"


// Forward declarations
class Mesh;
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
 *
 * \author Benjamin S. Kirk, 2002-2003
 */

// ------------------------------------------------------------
// Elem class definition
class Elem : public ReferenceCountedObject<Elem>
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
	Elem* _parent=NULL);

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
  unsigned char subdomain_id () const
  { return _sbd_id; };
  
  /**
   * @returns the subdomain that this element belongs to as a
   * writeable reference.
   */
  unsigned char & subdomain_id ()
  { return _sbd_id; };

  /**
   * @returns the processor that this element belongs to.
   * To conserve space this is stored as a short integer.
   */
  unsigned short int processor_id () const
  { return _proc_id; };
  
  /**
   * @returns the processor that this element belongs to as a
   * writeable reference.
   */
  unsigned short int & processor_id ()
  { return _proc_id; };

  /**
   * @returns an id assocated with this element.  The id is not
   * guaranteed to be unique, but it should be close.  The id
   * is thus useful, for example, as a key in a hash table
   * data structure.
   */
  unsigned int key() const;
  
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
  Elem* neighbor(const unsigned int i) const
  { assert (i < n_sides()); return _neighbors[i]; };

  /**
   * Assigns \p n as the \f$ i^{th} \f$ neighbor.
   */
  void set_neighbor(const unsigned int i, Elem* n)
  { assert (i < n_sides()); _neighbors[i]=n; return; };

  /**
   * This function tells you which neighbor you \p (e) are.
   * What is returned is the index of the side _in_the_neighbor_
   * which you share. 
   */
  unsigned int which_neighbor_am_i(const Elem *e) const; 

  /**
   * @returns the connectivity in the \p Tecplot format, which is
   * 1-based.  Also, maps all elements to quadrilaterals in 2D and
   * hexahedrals in 3D, which can be useful for writing hybrid meshes
   * for visualization. 
   */
  virtual const std::vector<unsigned int> tecplot_connectivity(const unsigned int sc=0) const = 0;

  /**
   * @returns the connectivity of the \f$ sc^{th} \f$
   * sub-element in the VTK format.
   */
  virtual void vtk_connectivity(const unsigned int sc,
 				std::vector<unsigned int>* conn) const = 0;

  /**
   * @returns the VTK element type of the sc-th sub-element.
   */
  virtual unsigned int vtk_element_type (const unsigned int sc) const = 0;
  
  /**
   * Writes the \p Tecplot connectivity to the \p out stream.
   * This function is actually defined in the Elem base class
   * because it works the same for all element types.
   */
  virtual void write_tecplot_connectivity(std::ostream &out) const;
  
  /**
   * Writes the \p UCD connectivity to the \p out stream.
   * This function is actually defined in the Elem base class
   * because it works the same for all element types.
   */
  virtual void write_ucd_connectivity(std::ostream &out) const;

  /**
   * @returns the type of element that has been derived from this
   * base class.
   */
  virtual ElemType type() const = 0;

  /**
   * @returns the dimensionality of the object.
   */
  virtual unsigned int dim () const = 0;
  
  /**
   * @returns the number of nodes this element contains. 
   */
  virtual unsigned int n_nodes() const = 0;

  /**
   * @returns the number of sides the element that has been derived
   * from this class has. In 2D the number of sides is the number
   * of edges, in 3D the number of sides is the number of faces.
   */
  virtual unsigned int n_sides() const = 0;

  /**
   * @returns the number of neighbors the element that has been derived
   * from this class has.  By default only face (or edge in 2D)
   * neighbors are stored, so this method returns n_sides(),
   * however it may be overloaded in a derived class
   */
  virtual unsigned int n_neighbors() const
  { return n_sides(); };

  /**
   * @returns the number of vertices the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_vertices() const = 0;

  /**
   * @returns the number of edges the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_edges() const = 0;

  /**
   * @returns the number of faces the element that has been derived
   * from this class has.
   */
  virtual unsigned int n_faces() const = 0;
  
  /**
   * @returns the number of children the element that has been derived
   * from this class may have.
   */
  virtual unsigned int n_children() const = 0;

  /**
   * @returns the number of sub-elements this element may be broken
   * down into for visualization purposes.  For example, this returns
   * 1 for a linear triangle, 4 for a quadratic (6-noded) triangle, etc...
   */
  virtual unsigned int n_sub_elem() const = 0;

  /**
   * @returns an element coincident with side \p i.  This method returns
   * the _minimum_ element necessary to uniquely identify the side.  So, 
   * for example, the side of a hexahedral is always returned as a 4-noded
   * quadrilateral, regardless of what type of hex you are dealing with.  If
   * you want the full-ordered face (i.e. a 9-noded quad face for a 27-noded
   * hexahedral) use the build_side method.
   */
  virtual AutoPtr<Elem> side(const unsigned int i) const = 0;
  
  /**
   * Creates an element coincident with side \p i. The element returned is
   * full-ordered, in contrast to the side method.  For example, calling 
   * build_side(0) on a 20-noded hex will build a 8-noded quadrilateral
   * coincident with face 0 and pass back the pointer.
   *
   * A \p AutoPtr<Elem> is returned to prevent a memory leak.
   * This way the user need not remember to delete the object.
   */
  virtual AutoPtr<Elem> build_side(const unsigned int i) const = 0;

  /**
   * @returns the default approximation order for this element type.
   * This is the order that will be used to compute the map to the
   * reference element.
   */
  virtual Order default_order() const = 0;
  
  /**
   * @returns the centriod of the element. The centroid is 
   * computed as the average of all the element vertices. 
   * This method is overloadable since some derived elements 
   * might want to use shortcuts to compute their centroid.
   */
  virtual Point centroid() const;
  
  /**
   * @returns the minimum vertex separation for the element.  
   * This method is overloadable since some derived elements 
   * might want to use shortcuts to compute their centroid.
   */
  virtual Real hmin() const;
  
  /**
   * @returns the maximum vertex separation for the element.
   * This method is overloadable since some derived elements 
   * might want to use shortcuts to compute their centroid.
   */
  virtual Real hmax() const;
  
  /**
   * Based on the quality metric q specified by the user,
   * returns a quantitative assessment of element quality.
   */
  virtual Real quality(const ElemQuality q) const;  

  /**
   * Returns the suggested quality bounds for
   * the hex based on quality measure q.  These are
   * the values suggested by the CUBIT User's Manual.
   * Since this function can have no possible meaning
   * for an abstract Elem, it is an error.
   */
  virtual std::pair<Real,Real> qual_bounds (const ElemQuality) const
  { error(); std::pair<Real,Real> p(0,0); return p; };
  
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
  bool active() const;
  
  /**
   * @returns a pointer to the element's parent.  Returns \p NULL if
   * the element was not created via refinement, i.e. was read from file.
   */
  const Elem* parent() const
  { return _parent; };
  
  /**
   * @returns the magnitude of the distance between nodes n1 and n2.
   * Useful for computing the lengths of the sides of elements.
   */
  Real length (const unsigned int n1, 
	       const unsigned int n2) const;
  
#ifdef ENABLE_AMR

  /**
   * Useful ENUM describing the refinement state of
   * an element.
   */
  enum RefinementState { COARSEN = -1,
			 DO_NOTHING,
			 REFINE,
			 JUST_REFINED };
  
  /**
   * @returns the refinement level of the current element.  If the
   * element's parent is \p NULL then by convention it is at
   * level 0, otherwise it is simply at one level greater than
   * its parent.
   */
  unsigned int level() const;
  
  /**
   * @returns a pointer to the \f$ i^{th} \f$ child for this element.
   * Returns \p NULL  if this element has no children, i.e. is active.
   */
  Elem* child(const unsigned int i) const
  {
    assert (_children    != NULL);
    assert (_children[i] != NULL);
    return _children[i];
  };

  /**
   * Returns the value of the refinement flag for the element.
   */
  RefinementState refinement_flag () const
  { return _rflag; };

  /**
   * Returns the value of the refinement flag for the element
   * as a writeable reference.
   */     
  RefinementState & set_refinement_flag ()
  { return _rflag; };

  /**
   * Refine the element.
   */
  virtual void refine (Mesh&) = 0;
  
  /**
   * Coarsen the element.  This is not
   * virtual since it is the same for all
   * element types.
   */
  void coarsen();

  /**
   * Convenience typedefs for providing iterator access
   * to the _neighbors array.  These make the accessor
   * functions easier to read.
   */
  typedef std::pair<Elem**, Elem**> ElemPair;
  typedef std::pair<const Elem**, const Elem**> ConstElemPair;

  /**
   * The non-const begin and end accessor functions.
   */
  ElemPair neighbors_begin() { return ElemPair (&_neighbors[0],         &_neighbors[n_sides()]); }
  ElemPair neighbors_end()   { return ElemPair (&_neighbors[n_sides()], &_neighbors[n_sides()]); }

  /**
   * The const begin and end accessor functions. 
   */
  ConstElemPair neighbors_begin() const { return ConstElemPair (const_cast<const Elem**>(&_neighbors[0]),
								const_cast<const Elem**>(&_neighbors[n_sides()])); }
  ConstElemPair neighbors_end()   const { return ConstElemPair (const_cast<const Elem**>(&_neighbors[n_sides()]),
								const_cast<const Elem**>(&_neighbors[n_sides()])); }
#endif

  
 protected:

  
  /**
   * Build an element of type \p type.  Be careful when using this
   * function: it allocates memory, so you must remember to delete
   * the pointer when you are done with it.  The user should not call this,
   * so it is protected.
   */
  static Elem* build (const ElemType type);
  
  /**
   * Pointers to the nodes we are conneted to.
   */
  Node** _nodes;

  /**
   * The subdomain to which this element belongs.
   */
  unsigned char _sbd_id;

  /**
   * The processor to which this element belongs.
   */
  unsigned short int _proc_id;

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
   * Refinement flag.
   */
  RefinementState _rflag;
  
#endif

  /**
   * Make the classes that need to access our build
   * member friends.  These classes do not really fit
   * the profile of what a "friend" should be, but
   * if we are going to protect the constructor and
   * the build method, there's no way around it.
   */
  friend class XdrInterface;
  friend class Mesh;
};




// ------------------------------------------------------------
// Elem class member functions
inline
Elem::Elem(const unsigned int nn,
	   const unsigned int ns,
	   Elem* p) :
  _parent(p)
{
  assert (nn);

  subdomain_id() = 0;
  processor_id() = 0;

  // Initialize the nodes data structure
  {
    _nodes = NULL;

    if (nn != 0)
      {
	_nodes = new Node*[nn]; 
	
	for (unsigned int n=0; n<nn; n++)
	  _nodes[n] = NULL;
      };
  };
  
  // Initialize the neighbors data structure
  {
    _neighbors = NULL;

    if (ns != 0)
      {
	_neighbors = new Elem*[ns]; 
	
	for (unsigned int n=0; n<ns; n++)
	  _neighbors[n] = NULL;
      };
  };

  // Optionally initialize data from the parent
  if (parent() != NULL)
    {
      subdomain_id() = parent()->subdomain_id();
      processor_id() = parent()->processor_id();
    };  

#ifdef ENABLE_AMR
  
  _children = NULL;

  set_refinement_flag() = Elem::DO_NOTHING;

#endif  
};



inline
Elem::~Elem() 
{
  if (_nodes != NULL)
    delete [] _nodes;
  _nodes = NULL;
  
  if (_neighbors != NULL)
    delete [] _neighbors;
  _neighbors = NULL;

#ifdef ENABLE_AMR
  
  if (_children != NULL)
    delete [] _children;
  _children = NULL;
  
#endif
};



inline
const Point & Elem::point (const unsigned int i) const
{
  assert (i < n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);

  return *_nodes[i];
};



inline
Point & Elem::point (const unsigned int i)
{
  assert (i < n_nodes());

  return *_nodes[i];
};



inline
unsigned int Elem::node (const unsigned int i) const
{
  assert (i < n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);

  return _nodes[i]->id();
};



inline
Node* Elem::get_node (const unsigned int i) const
{
  assert (i < n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);

  return _nodes[i];
};



inline
Node* & Elem::set_node (const unsigned int i)
{
  assert (i < n_nodes());

  return _nodes[i];
};



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
};


#endif // end #ifndef __elem_h__
