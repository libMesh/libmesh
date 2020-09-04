// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/auto_ptr.h" // libmesh_make_unique

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
 * \brief A geometric point in (x,y,z) space associated with a DOF.
 */
class Node : public Point,
             public DofObject,
             public ReferenceCountedObject<Node>
{

public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p LIBMESH_DIM dimensions with an \p id of \p Node::invalid_id.
   */
  explicit
  Node  (const Real x=0,
         const Real y=0,
         const Real z=0,
         const dof_id_type id = invalid_id);

  /**
   * "Copy"-constructor: deliberately slices the source Node and only
   * copies its' Point information, since copying anything else would
   * likely be a bug.
   *
   * \deprecated - the constructor from Point is what we really want.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  Node (const Node & n);
#else
  Node (const Node & n) = delete;
#endif

  /**
   * Copy-constructor from a \p Point.  Optionally assigned the \p id.
   */
  explicit Node (const Point & p,
                 const dof_id_type id = invalid_id);

  /**
   * Disambiguate constructing from non-Real scalars
   */
  template <typename T,
            typename = typename
              boostcopy::enable_if_c<ScalarTraits<T>::value,void>::type>
  explicit Node (const T x) :
    Point (x,0,0)
  { this->set_id() = invalid_id; }

  /**
   * Destructor.
   */
  ~Node ();

  /**
   * Assign to a node from a point.
   */
  Node & operator= (const Point & p);

  /**
   * \returns A \p Node copied from \p n and wrapped in a smart pointer.
   *
   * \deprecated - anyone copying a Node would almost certainly be
   * better off copying the much cheaper Point or taking a reference
   * to the Node.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  static std::unique_ptr<Node> build (const Node & n);
#endif

  /**
   * \returns A \p Node copied from \p p with id == \p id and wrapped in a smart pointer.
   */
  static std::unique_ptr<Node> build (const Point & p,
                                      const dof_id_type id);

  /**
   * \returns A \p Node created from the specified (x,y,z) positions
   * with id == \p id and wrapped in a smart pointer.
   */
  static std::unique_ptr<Node> build (const Real x,
                                      const Real y,
                                      const Real z,
                                      const dof_id_type id);

  /**
   * \returns \p true if the node is active.  An active node is
   * defined as one for which \p id() is not \p Node::invalid_id.
   * Inactive nodes are nodes that are in the mesh but are not
   * connected to any elements.
   */
  bool active () const;


  /**
   * \returns \p true if this node equals rhs, false otherwise.
   */
  bool operator ==(const Node & rhs) const;

  /**
   * Prints relevant information about the node.
   */
  void print_info (std::ostream & os=libMesh::out) const;

  /**
   * Prints relevant information about the node to a string.
   */
  std::string get_info () const;

#ifdef LIBMESH_HAVE_MPI
  unsigned int packed_size() const
  {
    const unsigned int header_size = 2;

    // use "(a+b-1)/b" trick to get a/b to round up
    static const unsigned int idtypes_per_Real =
      (sizeof(Real) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

    return header_size + LIBMESH_DIM*idtypes_per_Real +
      this->packed_indexing_size();
  }

#endif // #ifdef LIBMESH_HAVE_MPI

  /**
   * \returns The number of nodes connected with this node.
   * Currently, this value is invalid (zero) except for
   * subdivision meshes.
   */
  unsigned int valence() const
  {
#ifdef LIBMESH_ENABLE_NODE_VALENCE
    return _valence;
#else
    libmesh_not_implemented();
    return libMesh::invalid_uint;
#endif
  }

  /**
   * Sets the number of nodes connected with this node.
   */
  void set_valence(unsigned int val);

  /**
   * Return which of pid1 and pid2 would be preferred by the current
   * load-balancing heuristic applied to this node.
   */
  processor_id_type choose_processor_id(processor_id_type pid1, processor_id_type pid2) const;

private:

  /**
   * This class need access to the node key information,
   * but no one else should be able to mess with it.
   */
  friend class MeshRefinement;
  friend class Elem;

#ifdef LIBMESH_ENABLE_NODE_VALENCE
  /**
   * Type used to store node valence.
   */
  typedef unsigned char valence_idx_t;

  /**
   * The number of nodes connected with this node.
   * Currently, this value is invalid (zero) except for
   * subdivision meshes.
   */
  valence_idx_t _valence;
#endif
};



// ------------------------------------------------------------
// Global Node functions

inline
std::ostream & operator << (std::ostream & os, const Node & n)
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
  Point(x,y,z)
#ifdef LIBMESH_ENABLE_NODE_VALENCE
  ,
  _valence(0)
#endif
{
  this->set_id() = dofid;
}


#ifdef LIBMESH_ENABLE_DEPRECATED
inline
Node::Node (const Node & n) :
  Point(n),
  DofObject(), // Deliberately slicing!
  ReferenceCountedObject<Node>()
#ifdef LIBMESH_ENABLE_NODE_VALENCE
  ,
  _valence(0)
#endif
{
  libmesh_deprecated(); // Cast to Point first!
}
#endif // LIBMESH_ENABLE_DEPRECATED



inline
Node::Node (const Point & p,
            const dof_id_type dofid) :
  Point(p)
#ifdef LIBMESH_ENABLE_NODE_VALENCE
  ,
  _valence(0)
#endif
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
Node & Node::operator= (const Point & p)
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



#ifdef LIBMESH_ENABLE_DEPRECATED
inline
std::unique_ptr<Node> Node::build(const Node & n)
{
  libmesh_deprecated();
  return libmesh_make_unique<Node>(n);
}
#endif



inline
std::unique_ptr<Node> Node::build(const Point & p,
                                  const dof_id_type id)
{
  return libmesh_make_unique<Node>(p,id);
}



inline
std::unique_ptr<Node> Node::build(const Real x,
                                  const Real y,
                                  const Real z,
                                  const dof_id_type id)
{
  return libmesh_make_unique<Node>(x,y,z,id);
}



inline
bool Node::active () const
{
  return (this->id() != Node::invalid_id);
}



#ifdef LIBMESH_ENABLE_NODE_VALENCE

inline
void Node::set_valence (unsigned int val)
{
  _valence = cast_int<valence_idx_t>(val);
}

#else

inline
void Node::set_valence (unsigned int)
{
  libmesh_not_implemented();
}

#endif // #ifdef LIBMESH_ENABLE_NODE_VALENCE

} // namespace libMesh


#endif // LIBMESH_NODE_H
