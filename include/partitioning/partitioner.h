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



#ifndef LIBMESH_PARTITIONER_H
#define LIBMESH_PARTITIONER_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh_base.h" // for MeshBase::element_iterator

// C++ Includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
class ErrorVector;

/**
 * The \p Partitioner class provides a uniform interface for
 * partitioning algorithms.  It takes a reference to a \p MeshBase
 * object as input, which it will partition into a number of
 * subdomains.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Base class for all concrete Partitioner instantiations.
 */
class Partitioner
{
public:

  /**
   * Constructor.
   */
  Partitioner () : _weights(libmesh_nullptr) {}

  /**
   * Destructor. Virtual so that we can derive from this class.
   */
  virtual ~Partitioner() {}

  /**
   * Creates a new partitioner of this type and returns it in a \p
   * UniquePtr.  This is used when copying meshes, and must be
   * overloaded in the derived classes.
   */
  virtual UniquePtr<Partitioner> clone () const = 0;

  /**
   * Partitions the \p MeshBase into \p n parts by setting
   * processor_id() on Nodes and Elems.
   *
   * NOTE: If you are implementing a new type of Partitioner, you most
   * likely do *not* want to override the partition() function, see
   * instead the protected virtual _do_partition() method below.  The
   * partition() function is responsible for doing a lot of
   * libmesh-internals-specific setup and finalization before and
   * after the _do_partition() function is called.  The only
   * responsibility of the _do_partition() function, on the other
   * hand, is to set the processor IDs of the elements according to a
   * specific partitioning algorithm.  See, e.g. MetisPartitioner for
   * an example.
   */
  virtual void partition (MeshBase & mesh,
                          const unsigned int n);

  /**
   * Partitions the \p MeshBase into \p mesh.n_processors() by setting
   * processor_id() on Nodes and Elems.
   *
   * NOTE: If you are implementing a new type of Partitioner, you most
   * likely do *not* want to override the partition() function, see
   * instead the protected virtual _do_partition() method below.  The
   * partition() function is responsible for doing a lot of
   * libmesh-internals-specific setup and finalization before and
   * after the _do_partition() function is called.  The only
   * responsibility of the _do_partition() function, on the other
   * hand, is to set the processor IDs of the elements according to a
   * specific partitioning algorithm.  See, e.g. MetisPartitioner for
   * an example.
   */
  virtual void partition (MeshBase & mesh);

  /**
   * Partitions elements in the range (it, end) into n parts.
   * The mesh from which the iterators are created must also be passed
   * in, since it is a parallel object and has other useful
   * information in it.
   *
   * Although partition_range() is part of the public Partitioner
   * interface, it should not generally be called by applications.
   * Its main purpose is to support the SubdomainPartitioner, which
   * uses it internally to individually partition ranges of elements
   * before combining them into the final partitioning. Most of the
   * time, the protected _do_partition() function is implemented in
   * terms of partition_range() by passing a range which includes all
   * the elements of the Mesh.
   */
  virtual void partition_range (MeshBase & /*mesh*/,
                                MeshBase::element_iterator /*beg*/,
                                MeshBase::element_iterator /*end*/,
                                const unsigned int /*n_parts*/)
  { libmesh_not_implemented(); }

  /**
   * Repartitions the \p MeshBase into \p n parts. (Some partitoning
   * algorithms can repartition more efficiently than computing a new
   * partitioning from scratch.)  The default behavior is to simply
   * call this->partition(mesh,n).
   */
  void repartition (MeshBase & mesh,
                    const unsigned int n);

  /**
   * Repartitions the \p MeshBase into \p mesh.n_processors() parts.  This
   * is required since some partitoning algorithms can repartition
   * more efficiently than computing a new partitioning from scratch.
   */
  void repartition (MeshBase & mesh);

  /**
   * These functions assign processor IDs to newly-created elements
   * (in parallel) which are currently assigned to processor 0.
   */
  static void partition_unpartitioned_elements (MeshBase & mesh);

  static void partition_unpartitioned_elements (MeshBase & mesh,
                                                const unsigned int n);

  /**
   * This function is called after partitioning to set the processor IDs
   * for the inactive parent elements.  A parent's processor ID is the same
   * as its first child.
   */
  static void set_parent_processor_ids(MeshBase & mesh);

  /**
   * This function is called after partitioning to set the processor IDs
   * for the nodes.  By definition, a Node's processor ID is the minimum
   * processor ID for all of the elements which share the node.
   */
  static void set_node_processor_ids(MeshBase & mesh);

  /**
   * Attach weights that can be used for partitioning.  This ErrorVector should be
   * _exactly_ the same on every processor and should have mesh->max_elem_id()
   * entries.
   */
  virtual void attach_weights(ErrorVector * /*weights*/) { libmesh_not_implemented(); }

protected:

  /**
   * Trivially "partitions" the mesh for one processor.
   * Simply loops through the elements and assigns all of them
   * to processor 0.  Is is provided as a separate function
   * so that derived classes may use it without reimplementing it.
   */
  void single_partition (MeshBase & mesh);

  /**
   * Slightly generalized version of single_partition which acts on a
   * range of elements defined by the pair of iterators (it, end).
   */
  void single_partition_range(MeshBase::element_iterator it,
                              MeshBase::element_iterator end);

  /**
   * This is the actual partitioning method which must be overloaded
   * in derived classes.  It is called via the public partition()
   * method above by the user.
   */
  virtual void _do_partition(MeshBase & mesh,
                             const unsigned int n) = 0;

  /**
   * This is the actual re-partitioning method which can be overloaded
   * in derived classes.  Note that the default behavior is to simply
   * call the partition function.
   */
  virtual void _do_repartition (MeshBase & mesh,
                                const unsigned int n) { this->_do_partition (mesh, n); }

  /**
   * The blocksize to use when doing blocked parallel communication.  This limits the
   * maximum vector size which can be used in a single communication step.
   */
  static const dof_id_type communication_blocksize;

  /**
   * The weights that might be used for partitioning.
   */
  ErrorVector * _weights;
};

} // namespace libMesh

#endif // LIBMESH_PARTITIONER_H
