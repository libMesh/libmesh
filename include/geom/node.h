// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <sstream>

namespace libMesh
{

// forward declarations
template <typename> class MeshBaseTempl;
template <typename> class MeshRefinementTempl;
template <typename> class ElemTempl;

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
template <typename RealType>
class NodeTempl : public PointTempl<RealType>,
                  public DofObject,
                  public ReferenceCountedObject<NodeTempl<RealType>>
{

public:
  typedef NodeTempl<RealType> Node;
  typedef PointTempl<RealType> Point;
  typedef ElemTempl<RealType> Elem;
  typedef MeshRefinementTempl<RealType> MeshRefinement;

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p LIBMESH_DIM dimensions with an \p id of \p Node::invalid_id.
   */
  explicit
  NodeTempl  (const RealType x=0,
              const RealType y=0,
              const RealType z=0,
              const dof_id_type id = invalid_id);

  /**
   * Copy-constructor.
   *
   * \deprecated - anyone copying a Node would almost certainly be
   * better off copying the much cheaper Point or taking a reference
   * to the Node.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  NodeTempl (const NodeTempl & n);
#endif

  /**
   * Copy-constructor from a \p Point.  Optionally assigned the \p id.
   */
  explicit NodeTempl (const PointTempl<RealType> & p,
                      const dof_id_type id = invalid_id);

  /**
   * Disambiguate constructing from non-Real scalars
   */
  template <typename T,
            typename = typename
              boostcopy::enable_if_c<ScalarTraits<T>::value,void>::type>
  NodeTempl (const T x) :
    PointTempl<RealType> (x,0,0)
  { this->set_id() = invalid_id; }

  /**
   * Destructor.
   */
  ~NodeTempl ();

  /**
   * Assign to a node from a point.
   */
  NodeTempl<RealType> & operator= (const PointTempl<RealType> & p);

  /**
   * \returns A \p Node copied from \p n and wrapped in a smart pointer.
   *
   * \deprecated - anyone copying a Node would almost certainly be
   * better off copying the much cheaper Point or taking a reference
   * to the Node.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  static std::unique_ptr<NodeTempl<RealType>> build (const NodeTempl<RealType> & n);
#endif

  /**
   * \returns A \p Node copied from \p p with id == \p id and wrapped in a smart pointer.
   */
  static std::unique_ptr<NodeTempl<RealType>> build (const PointTempl<RealType> & p,
                                                     const dof_id_type id);

  /**
   * \returns A \p Node created from the specified (x,y,z) positions
   * with id == \p id and wrapped in a smart pointer.
   */
  static std::unique_ptr<NodeTempl<RealType>> build (const RealType x,
                                                     const RealType y,
                                                     const RealType z,
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
  bool operator ==(const NodeTempl<RealType> & rhs) const;

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
      (sizeof(RealType) + sizeof(largest_id_type) - 1) / sizeof(largest_id_type);

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
  friend class MeshRefinementTempl<RealType>;
  friend class ElemTempl<RealType>;

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
template <typename RealType>
inline
std::ostream & operator << (std::ostream & os, const NodeTempl<RealType> & n)
{
  n.print_info(os);
  return os;
}



//------------------------------------------------------
// Inline functions
template <typename RealType>
inline
NodeTempl<RealType>::NodeTempl (const RealType x,
                                const RealType y,
                                const RealType z,
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
template <typename RealType>
inline
NodeTempl<RealType>::NodeTempl (const NodeTempl<RealType> & n) :
  PointTempl<RealType>(n),
  DofObject(n),
  ReferenceCountedObject<NodeTempl<RealType>>()
#ifdef LIBMESH_ENABLE_NODE_VALENCE
  ,
  _valence(n._valence)
#endif
{
  libmesh_deprecated();
}
#endif


template <typename RealType>
inline
NodeTempl<RealType>::NodeTempl (const PointTempl<RealType> & p,
                                const dof_id_type dofid) :
  PointTempl<RealType>(p)
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


template <typename RealType>
inline
NodeTempl<RealType>::~NodeTempl ()
{
}


template <typename RealType>
inline
NodeTempl<RealType> &
NodeTempl<RealType>::operator= (const PointTempl<RealType> & p)
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
template <typename RealType>
inline
std::unique_ptr<NodeTempl<RealType>>
NodeTempl<RealType>::build(const NodeTempl<RealType> & n)
{
  libmesh_deprecated();
  return libmesh_make_unique<NodeTempl<RealType>>(n);
}
#endif


template <typename RealType>
inline
std::unique_ptr<NodeTempl<RealType>>
NodeTempl<RealType>::build(const Point & p,
                           const dof_id_type id)
{
  return libmesh_make_unique<NodeTempl<RealType>>(p,id);
}



template <typename RealType>
inline
std::unique_ptr<NodeTempl<RealType>>
NodeTempl<RealType>::build(const RealType x,
                           const RealType y,
                           const RealType z,
                           const dof_id_type id)
{
  return libmesh_make_unique<NodeTempl<RealType>>(x,y,z,id);
}



template <typename RealType>
inline
bool NodeTempl<RealType>::active () const
{
  return (this->id() != NodeTempl<RealType>::invalid_id);
}



#ifdef LIBMESH_ENABLE_NODE_VALENCE

template <typename RealType>
inline
void NodeTempl<RealType>::set_valence (unsigned int val)
{
  _valence = cast_int<valence_idx_t>(val);
}

#else

template <typename RealType>
inline
void NodeTempl<RealType>::set_valence (unsigned int)
{
  libmesh_not_implemented();
}

#endif // #ifdef LIBMESH_ENABLE_NODE_VALENCE

template <typename RealType>
bool NodeTempl<RealType>::operator==(const NodeTempl<RealType> & rhs) const
{
  // Explicitly calling the operator== defined in Point
  return this->Point::operator==(rhs);
}



template <typename RealType>
void NodeTempl<RealType>::print_info (std::ostream & os) const
{
  os << this->get_info()
     << std::endl;
}



template <typename RealType>
std::string NodeTempl<RealType>::get_info () const
{
  std::ostringstream oss;

  oss << "  NodeTempl<RealType> id()=";

  if (this->valid_id())
    oss << this->id();
  else
    oss << "invalid";

  oss << ", processor_id()=" << this->processor_id() <<
    ", Point=" << *static_cast<const Point *>(this) << '\n';

  oss << "    DoFs=";
  for (auto s : IntRange<unsigned int>(0, this->n_systems()))
    for (auto v : IntRange<unsigned int>(0, this->n_vars(s)))
      for (auto c : IntRange<unsigned int>(0, this->n_comp(s,v)))
        oss << '(' << s << '/' << v << '/' << this->dof_number(s,v,c) << ") ";

  return oss.str();
}


template <typename RealType>
processor_id_type
NodeTempl<RealType>::choose_processor_id(processor_id_type pid1, processor_id_type pid2) const
{
  if (pid1 == DofObject::invalid_processor_id)
    return pid2;

  // Do we want the new load-balanced node partitioning heuristic
  // instead of the default partitioner-friendlier heuristic?
  static bool load_balanced_nodes =
    libMesh::on_command_line ("--load-balanced-nodes");

  // For better load balancing, we can use the min
  // even-numberered nodes and the max for odd-numbered.
  if (load_balanced_nodes)
    {
      if (this->id() % 2 &&
          pid2 != DofObject::invalid_processor_id)
        return std::max(pid1, pid2);
      else
        return std::min(pid1, pid2);
    }

  // Our default behavior, which puts too many nodes on lower MPI
  // ranks but which keeps elements' nodes on the same partition more
  // often, is simply:
  return std::min(pid1, pid2);
}

typedef NodeTempl<Real> Node;

} // namespace libMesh


#endif // LIBMESH_NODE_H
