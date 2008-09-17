// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __node_h__
#define __node_h__

// C++ includes

// Local includes
#include "point.h"
#include "dof_object.h"
#include "reference_counted_object.h"
#include "auto_ptr.h"


// forward declarations
class Node;
class MeshRefinement;


/**
 * A \p Node is like a \p Point, but with more information.  A \p Node
 * is located in space and is associated with some \p (x,y,z)
 * coordinates.  Additionally, a \p Node may be enumerated with a
 * global \p id.  Finally, a \p Node may have an arbitrary number of
 * degrees of freedom associated with it.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \version $Revision$
 */

class Node : public Point,
	     public DofObject,      
	     public ReferenceCountedObject<Node>
{
  
public:
  
  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p DIM dimensions with an \p id of \p Node::invalid_id
   */
  Node  (const Real x,
	 const Real y,
	 const Real z,
	 const unsigned int id = invalid_id);

  /**
   * Copy-constructor.
   */
  Node (const Node& n);

  /**
   * Copy-constructor from a \p Point.  Optionally assigned the \p id.
   */
  explicit Node (const Point& p,
	         const unsigned int id = invalid_id);

  /**
   * Destructor.
   */ 
  virtual ~Node ();

  /**
   * Assign to a node from a point
   */
  Node& operator= (const Point& p);
  
  /**
   * Builds a \p Node and returns an \p AutoPtr<Node> to the
   * newly-created object.  The \p id is copied from \p n.id()
   */
  static AutoPtr<Node> build (const Node& n);

  /**
   * Builds a \p Node from \p Point p and returns an \p AutoPtr<Node>
   * to the newly-created object.  Optionally assignes the \p id.
   */
  static AutoPtr<Node> build (const Point& p,
			      const unsigned int id);
  
  /**
   * Builds a \p Node from specified points and returns an \p AutoPtr<Node>
   * to the newly-created object.  Optionally assigned the \p id.
   */
  static AutoPtr<Node> build (const Real x,
			      const Real y,
			      const Real z,
			      const unsigned int id);

  /**
   * @returns \p true if the node is active.  An active node is
   * defined as one for which \p id() is not \p Node::invalid_id.
   * Inactive nodes are nodes that are in the mesh but are not
   * connected to any elements.
   */
  bool active () const;


  /**
   * @returns \p true if this node equals rhs, false otherwise.
   * Note that rhs is a DofObject, so it is possible that rhs is
   * not a Node. If rhs is not a Node, then we return false of course.
   */
  virtual bool operator==(const DofObject& rhs) const;


#ifdef HAVE_MPI
  /**
   * Convenient way to communicate nodes.  This struct defines a
   * packed up node which can be easily communicated through a
   * derived MPI datatype.
   *
   * \author Benjamin S. Kirk
   * \date 2008
   */
  struct PackedNode
  {
    unsigned int id;
    unsigned int pid;
    Real x;
    Real y;
    Real z;

    PackedNode () :
      id(0),
      x(0.),
      y(0.),
      z(0.)
    {}

    PackedNode (const Node &node) :
      id(node.id()),
      pid(node.processor_id()),
      x(node(0)),
      y(node(1)),
      z(node(2))
    {}

    AutoPtr<Node> build_node () const
    {
      AutoPtr<Node> node(new Node(x,y,z,id));
      node->processor_id() = pid;
      return node;
    }

    Point build_point () const
    {
      return Point(x,y,z);
    }
    
    static MPI_Datatype create_mpi_datatype ();
    
  };
#endif // #ifdef HAVE_MPI
  
private:

  /**
   * This class need access to the node key information,
   * but no one else should be able to mess with it.
   */
  friend class MeshRefinement;
  friend class Elem;
};



//------------------------------------------------------
// Inline functions
inline
Node::Node (const Real x,
	    const Real y,
	    const Real z,
	    const unsigned int id) :
  Point(x,y,z)
{
  this->set_id() = id;
}



inline
Node::Node (const Node& n) :
  Point(n),
  DofObject(n),
  ReferenceCountedObject<Node>()
{
}



inline
Node::Node (const Point& p,
	    const unsigned int id) :
  Point(p)
{
  // optionally assign the id.  We have
  // to do it like this otherwise
  // Node n = Point p would erase
  // the id!
  if (id != invalid_id)
    this->set_id() = id;
}



inline
Node::~Node ()
{
}



inline
Node & Node::operator= (const Point& p)
{
  (*this)(0) = p(0);
  (*this)(1) = p(1);
  (*this)(2) = p(2);

  return *this;
}



inline
AutoPtr<Node> Node::build(const Node& n)
{
  AutoPtr<Node> ap(new Node(n));
  return ap;
}



inline
AutoPtr<Node> Node::build(const Point& p,
			  const unsigned int id)
{
  
  AutoPtr<Node> ap(new Node(p,id));
  return ap;
}



inline
AutoPtr<Node> Node::build(const Real x,
			  const Real y,
			  const Real z,
			  const unsigned int id)
{
  AutoPtr<Node> ap(new Node(x,y,z,id));
  return ap;
}



inline
bool Node::active () const
{
  return (this->id() != Node::invalid_id);
}



#endif
