// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_NODE_H
#define LIBMESH_NODE_H

// Local includes
#include "libmesh/point.h"
#include "libmesh/dof_object.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <iostream>
#include <vector>

namespace libMesh
{


// forward declarations
class Node;
class MeshBase;
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
 */

class Node : public Point,
             public DofObject,
             public ReferenceCountedObject<Node>
{

public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p LIBMESH_DIM dimensions with an \p id of \p Node::invalid_id
   */
  explicit
  Node  (const Real x=0,
         const Real y=0,
         const Real z=0,
         const dof_id_type id = invalid_id);

  /**
   * Copy-constructor.
   */
  Node (const Node& n);

  /**
   * Copy-constructor from a \p Point.  Optionally assigned the \p id.
   */
  explicit Node (const Point& p,
                 const dof_id_type id = invalid_id);

  /**
   * Destructor.
   */
  ~Node ();

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
                              const dof_id_type id);

  /**
   * Builds a \p Node from specified points and returns an \p AutoPtr<Node>
   * to the newly-created object.  Optionally assigned the \p id.
   */
  static AutoPtr<Node> build (const Real x,
                              const Real y,
                              const Real z,
                              const dof_id_type id);

  /**
   * @returns \p true if the node is active.  An active node is
   * defined as one for which \p id() is not \p Node::invalid_id.
   * Inactive nodes are nodes that are in the mesh but are not
   * connected to any elements.
   */
  bool active () const;


  /**
   * @returns \p true if this node equals rhs, false otherwise.
   */
  bool operator ==(const Node& rhs) const;

  /**
   * Prints relevant information about the node
   */
  void print_info (std::ostream& os=libMesh::out) const;

  /**
   * Prints relevant information about the node to a string.
   */
  std::string get_info () const;

#ifdef LIBMESH_HAVE_MPI
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
    static const unsigned int header_size = 2;

    dof_id_type id;
    processor_id_type pid;
    Real x;
    // FIXME: We should drop z (and y) if libMesh is built 2D (or 1D) only
    Real y;
    Real z;

    PackedNode () :
      id(0),
      pid(DofObject::invalid_processor_id),
      x(0.),
      y(0.),
      z(0.)
    {}

    explicit
    PackedNode (const Node &node) :
      id(node.id()),
      pid(node.processor_id()),
      x(node(0)),
#if LIBMESH_DIM > 1
      y(node(1)),
#else
      y(0.),
#endif
#if LIBMESH_DIM > 2
      z(node(2))
#else
      z(0.)
#endif
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

    /**
     * For each node the serialization is of the form
     * [ processor_id self_ID x1 x2 y1 y2 z1 z2
     *  dof_object_buffer_1 ...]
     * There may be 1 or 3 or 4 ints per coordinate depending on
     * machine architecture.
     */
    static void pack (std::vector<largest_id_type> &conn, const Node* node);

    static void unpack (std::vector<largest_id_type>::const_iterator start, Node& node);
  };

  unsigned int packed_size() const
  {
    // use "(a+b-1)/b" trick to get a/b to round up
    static const unsigned int idtypes_per_Real =
      (sizeof(Real) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

    return PackedNode::header_size + LIBMESH_DIM*idtypes_per_Real +
      this->packed_indexing_size();
  }

#endif // #ifdef LIBMESH_HAVE_MPI

  /**
   * @returns the number of nodes connected with this node.
   * Currently, this value is invalid (zero) except for
   * subdivision meshes.
   */
  unsigned int valence() const { return _valence; }

  /**
   * Sets the number of nodes connected with this node.
   */
  void set_valence(unsigned int val) { _valence = val; }

private:

  /**
   * The number of nodes connected with this node.
   * Currently, this value is invalid (zero) except for
   * subdivision meshes.
   */
  unsigned int _valence;

  /**
   * This class need access to the node key information,
   * but no one else should be able to mess with it.
   */
  friend class MeshRefinement;
  friend class Elem;
};



// ------------------------------------------------------------
// global Node functions

inline
std::ostream& operator << (std::ostream& os, const Node& n)
{
  n.print_info(os);
  return os;
}



//------------------------------------------------------
// Inline functions
inline
Node::Node (const Real x,
            const Real y,
            const Real z,
            const dof_id_type dofid) :
  Point(x,y,z),
  _valence(0)
{
  this->set_id() = dofid;
}



inline
Node::Node (const Node& n) :
  Point(n),
  DofObject(n),
  ReferenceCountedObject<Node>(),
  _valence(n.valence())
{
}



inline
Node::Node (const Point& p,
            const dof_id_type dofid) :
  Point(p)
{
  // optionally assign the id.  We have
  // to do it like this otherwise
  // Node n = Point p would erase
  // the id!
  if (dofid != invalid_id)
    this->set_id() = dofid;
}



inline
Node::~Node ()
{
}



inline
Node & Node::operator= (const Point& p)
{
  (*this)(0) = p(0);
#if LIBMESH_DIM > 1
  (*this)(1) = p(1);
#endif
#if LIBMESH_DIM > 2
  (*this)(2) = p(2);
#endif

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
                          const dof_id_type id)
{

  AutoPtr<Node> ap(new Node(p,id));
  return ap;
}



inline
AutoPtr<Node> Node::build(const Real x,
                          const Real y,
                          const Real z,
                          const dof_id_type id)
{
  AutoPtr<Node> ap(new Node(x,y,z,id));
  return ap;
}



inline
bool Node::active () const
{
  return (this->id() != Node::invalid_id);
}


} // namespace libMesh


#endif // LIBMESH_NODE_H
