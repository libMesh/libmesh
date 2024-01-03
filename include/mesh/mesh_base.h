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



#ifndef LIBMESH_MESH_BASE_H
#define LIBMESH_MESH_BASE_H

// Local Includes
#include "libmesh/dof_object.h" // for invalid_processor_id
#include "libmesh/enum_order.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/multi_predicates.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/parallel_object.h"
#include "libmesh/simple_range.h"

// C++ Includes
#include <cstddef>
#include <string>
#include <memory>

#include "libmesh/vector_value.h"

// periodic boundary condition support
// Use forward declarations inside the libMesh namespace
namespace libMesh
{
  class PeriodicBoundary;
  class PeriodicBoundaries;
}

namespace libMesh
{

// forward declarations
class Elem;
class GhostingFunctor;
class Node;
class Point;
class Partitioner;
class BoundaryInfo;

template <typename T>
class SparseMatrix;

enum ElemType : int;
enum ElemMappingType : unsigned char;

template <class MT>
class MeshInput;


/**
 * This is the \p MeshBase class. This class provides all the data necessary
 * to describe a geometric entity.  It allows for the description of a
 * \p dim dimensional object that lives in \p LIBMESH_DIM-dimensional space.
 * \par
 * A mesh is made of nodes and elements, and this class provides data
 * structures to store and access both.  A mesh may be partitioned into a
 * number of subdomains, and this class provides that functionality.
 * Furthermore, this class provides functions for reading and writing a
 * mesh to disk in various formats.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Base class for Mesh.
 */
class MeshBase : public ParallelObject
{
public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  MeshBase (const Parallel::Communicator & comm_in,
            unsigned char dim=1);

  /**
   * Copy-constructor.
   */
  MeshBase (const MeshBase & other_mesh);

  /**
   * Move-constructor - deleted because after a theoretical move-construction
   * and then destruction of the moved-from object, the moved \p BoundaryInfo
   * would hold an invalid reference to the moved-from mesh
   */
  MeshBase(MeshBase &&) = delete;

  /**
   * Copy and move assignment are not allowed because MeshBase
   * subclasses manually manage memory (Elems and Nodes) and therefore
   * the default versions of these operators would leak memory.  Since
   * we don't want to maintain non-default copy and move assignment
   * operators at this time, the safest and most self-documenting
   * approach is to delete them.
   *
   * If you need to copy a Mesh, use the clone() method.
   */
  MeshBase & operator= (const MeshBase &) = delete;
  MeshBase & operator= (MeshBase && other_mesh);

  /**
   * Shim to allow operator = (&&) to behave like a virtual function
   * without having to be one.
   */
  virtual MeshBase & assign(MeshBase && other_mesh) = 0;

  /**
   * This tests for exactly-equal data in all the senses that a
   * mathematician would care about (element connectivity, nodal
   * coordinates), but in the senses a programmer would care about it
   * allows for non-equal equivalence in some ways (we accept
   * different Elem/Node addresses in memory) but not others (we do
   * not accept different subclass types, nor even different Elem/Node
   * ids).
   *
   * Though this method is non-virtual, its implementation calls the
   * virtual function \p subclass_locally_equals() to test for
   * equality of subclass-specific data as well.
   */
  bool operator== (const MeshBase & other_mesh) const;

  bool operator!= (const MeshBase & other_mesh) const
  {
    return !(*this == other_mesh);
  }

  /**
   * This behaves the same as operator==, but only for the local and
   * ghosted aspects of the mesh; i.e. operator== is true iff local
   * equality is true on every rank.
   */
  bool locally_equals (const MeshBase & other_mesh) const;

  /**
   * Virtual "copy constructor".  The copy will be of the same
   * subclass as \p this, and will satisfy "copy == this" when it is
   * created.
   */
  virtual std::unique_ptr<MeshBase> clone() const = 0;

  /**
   * Destructor.
   */
  virtual ~MeshBase ();

  /**
   * A partitioner to use at each partitioning
   */
  virtual std::unique_ptr<Partitioner> & partitioner() { return _partitioner; }

  /**
   * The information about boundary ids on the mesh
   */
  const BoundaryInfo & get_boundary_info() const { return *boundary_info; }

  /**
   * Writable information about boundary ids on the mesh
   */
  BoundaryInfo & get_boundary_info() { return *boundary_info; }

  /**
   * Deletes all the element and node data that is currently stored.
   *
   * elem and node extra_integer data is nevertheless *retained* here,
   * for better compatibility between that feature and older code's
   * use of MeshBase::clear()
   */
  virtual void clear ();

  /**
   * Deletes all the element data that is currently stored.
   *
   * No Node is removed from the mesh, however even NodeElem elements
   * are deleted, so the remaining Nodes will be considered "unused"
   * and cleared unless they are reconnected to new elements before
   * the next preparation step.
   *
   * This does not affect BoundaryInfo data; any boundary information
   * associated elements should already be cleared.
   */
  virtual void clear_elems () = 0;

  /**
   * \returns \p true if the mesh has undergone all the preparation
   * done in a call to \p prepare_for_use, \p false otherwise.
   */
  bool is_prepared () const
  { return _preparation; }

  /**
   * \returns the \p Preparation structure with details about in what
   * ways \p this mesh is currently prepared or unprepared.  This
   * structure may change in the future when cache designs change.
   */
  struct Preparation;

  Preparation preparation () const
  { return _preparation; }

  /**
   * Tells this we have done some operation where we should no longer
   * consider ourself prepared.  This is a very coarse setting; it may
   * be preferable to mark finer-grained settings instead.
   */
  void set_isnt_prepared()
  { _preparation = false; }

  /**
   * Tells this we have done some operation creating unpartitioned
   * elements.
   */
  void set_isnt_partitioned()
  { _preparation.is_partitioned = false; }

  /**
   * Tells this we have done some operation (e.g. adding elements to a
   * distributed mesh on one processor only) which can lose
   * synchronization of id counts.
   */
  void set_hasnt_synched_id_counts()
  { _preparation.has_synched_id_counts = false; }

  /**
   * Tells this we have done some operation (e.g. adding elements
   * without setting their neighbor pointers) which requires neighbor
   * pointers to be found later
   */
  void set_hasnt_neighbor_ptrs()
  { _preparation.has_neighbor_ptrs = false; }

  /**
   * Tells this we have done some operation (e.g. adding elements with
   * a new dimension or subdomain value) which may invalidate cached
   * summaries of element data
   */
  void set_hasnt_cached_elem_data()
  { _preparation.has_cached_elem_data = false; }

  /**
   * Tells this we have done some operation (e.g. refining elements
   * with interior parents) which requires interior parent pointers to
   * be found later
   */
  void set_hasnt_interior_parent_ptrs()
  { _preparation.has_interior_parent_ptrs = false; }

  /**
   * Tells this we have done some operation (e.g. repartitioning)
   * which may have left elements as ghosted which on a distributed
   * mesh should be remote.
   *
   * User code should probably never need to use this; we can set it
   * in Partitioner.
   */
  void set_hasnt_removed_remote_elements()
  { _preparation.has_removed_remote_elements = false; }

  /**
   * Tells this we have done some operation (e.g. coarsening)
   * which may have left orphaned nodes in need of removal.
   *
   * Most user code should probably never need to use this; we can set
   * it in MeshRefinement.
   */
  void set_hasnt_removed_orphaned_nodes()
  { _preparation.has_removed_orphaned_nodes = false; }

  /**
   * Tells this we have done some operation (e.g. adding or removing
   * elements) which may require a reinit() of custom ghosting
   * functors.
   */
  void hasnt_reinit_ghosting_functors()
  { _preparation.has_reinit_ghosting_functors = false; }

  /**
   * Tells this we have done some operation (e.g. removing elements or
   * adding new boundary entries) which may have invalidated our
   * cached boundary id sets.
   */
  void set_hasnt_boundary_id_sets()
  { _preparation.has_boundary_id_sets = false; }

  /**
   * \returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const
  { return true; }

  /**
   * \returns \p true if all elements and nodes of the mesh
   * exist on the processor 0, \p false otherwise
   */
  virtual bool is_serial_on_zero () const
  { return true; }

  /**
   * Asserts that not all elements and nodes of the mesh necessarily
   * exist on the current processor.  Only valid to call on classes
   * which can be created in a distributed form.
   */
  virtual void set_distributed ()
  { libmesh_error(); }

  /**
   * \returns \p true if new elements and nodes can and should be
   * created in synchronization on all processors, \p false otherwise
   */
  virtual bool is_replicated () const
  { return true; }

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor
   */
  virtual void allgather () {}

  /**
   * Gathers all elements and nodes of the mesh onto
   * processor zero
   */
  virtual void gather_to_zero() {}

  /**
   * When supported, deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   */
  virtual void delete_remote_elements () {
    _preparation.has_removed_remote_elements = true;
  }

  /**
   * Loops over ghosting functors and calls mesh_reinit()
   */
  void reinit_ghosting_functors();

  /**
   * \returns The logical dimension of the mesh; i.e. the manifold
   * dimension of the elements in the mesh.  When we have
   * multi-dimensional meshes (e.g. hexes and quads in the same mesh)
   * then this will return the largest such dimension.
   */
  unsigned int mesh_dimension () const;

  /**
   * Resets the logical dimension of the mesh. If the mesh has
   * elements of multiple dimensions, this should be set to the largest
   * dimension. E.g. if the mesh has 1D and 2D elements, this should
   * be set to 2. If the mesh has 2D and 3D elements, this should be
   * set to 3.
   */
  void set_mesh_dimension (unsigned char d)
  { _elem_dims.clear(); _elem_dims.insert(d); }

  /**
   * \returns A const reference to a std::set of element dimensions
   * present in the mesh.
   */
  const std::set<unsigned char> & elem_dimensions() const
  { return _elem_dims; }

  /**
   * \returns A const reference to a std::set of element default
   * orders present in the mesh.
   */
  const std::set<Order> & elem_default_orders() const
  { return _elem_default_orders; }

  /**
   * \returns The smallest supported_nodal_order() of any element
   * present in the mesh, which is thus the maximum supported nodal
   * order on the mesh as a whole.
   */
  Order supported_nodal_order() const
  { return _supported_nodal_order; }

  /**
   * Most of the time you should not need to call this, as the element
   * dimensions will be set automatically by a call to cache_elem_data(),
   * therefore only call this if you know what you're doing.
   *
   * In some specialized situations, for example when adding a single
   * Elem on all procs, it can be faster to skip calling cache_elem_data()
   * and simply specify the element dimensions manually, which is why this
   * setter exists.
   */
  void set_elem_dimensions(std::set<unsigned char> elem_dims);

  /**
   * Typedef for the "set" container used to store elemset ids. The
   * main requirements are that the entries be sorted and unique, so
   * std::set works for this, but there may be more efficient
   * alternatives.
   */
  typedef std::set<elemset_id_type> elemset_type;

  /**
   * Tabulate a user-defined "code" for elements which belong to the element sets
   * specified in \p id_set. For example, suppose that we have two elemsets A and
   * B with the following Elem ids:
   * Elemset A = {1, 3}
   * Elemset B = {2, 3}
   *
   * This implies the following mapping from elem id to elemset id:
   * Elem 1 -> {A}
   * Elem 2 -> {B}
   * Elem 3 -> {A,B}
   *
   * In this case, we would need to tabulate three different elemset codes, e.g.:
   * 0 -> {A}
   * 1 -> {B}
   * 2 -> {A,B}
   *
   * Also sets up the inverse mapping, so that if one knows all the
   * element sets an Elem belongs to, one can look up the
   * corresponding elemset code.
   */
  void add_elemset_code(dof_id_type code, MeshBase::elemset_type id_set);

  /**
   * Returns the number of unique elemset ids which have been added
   * via add_elemset_code(), which is the size of the _all_elemset_ids
   * set.
   */
  unsigned int n_elemsets() const;

  /**
   * Look up the element sets for a given elemset code and
   * vice-versa. The elemset must have been previously stored by
   * calling add_elemset_code(). If no such code/set is found, returns
   * the empty set or DofObject::invalid_id, respectively.
   */
  void get_elemsets(dof_id_type elemset_code, MeshBase::elemset_type & id_set_to_fill) const;
  dof_id_type get_elemset_code(const MeshBase::elemset_type & id_set) const;

  /**
   * Return a vector of all elemset codes defined on the mesh. We get
   * this by looping over the _elemset_codes map.
   */
  std::vector<dof_id_type> get_elemset_codes() const;

  /**
   * Replace elemset code "old_code" with "new_code". This function loops over
   * all elements and changes the extra integer corresponding to the "elemset_code"
   * label, and updates the _elemset_codes and _elemset_codes_inverse_map members.
   * Does not change the elemset ids of any of the sets.
   */
  void change_elemset_code(dof_id_type old_code, dof_id_type new_code);

  /**
   * Replace elemset id "old_id" with "new_id". Does not change any of the
   * elemset codes, so does not need to loop over the elements themselves.
   */
  void change_elemset_id(elemset_id_type old_id, elemset_id_type new_id);

  /**
   * \returns The "spatial dimension" of the mesh.
   *
   * The spatial dimension is defined as:
   *
   *   1 - for an exactly x-aligned mesh of 1D elements
   *   2 - for an exactly x-y planar mesh of 2D elements
   *   3 - otherwise
   *
   * No tolerance checks are performed to determine whether the Mesh
   * is x-aligned or x-y planar, only strict equality with zero in the
   * higher dimensions is checked.  Also, x-z and y-z planar meshes are
   * considered to have spatial dimension == 3.
   *
   * The spatial dimension is updated during prepare_for_use() based
   * on the dimensions of the various elements present in the Mesh,
   * but is *never automatically decreased* by this function.
   *
   * For example, if the user calls set_spatial_dimension(2) and then
   * later inserts 3D elements into the mesh,
   * Mesh::spatial_dimension() will return 3 after the next call to
   * prepare_for_use().  On the other hand, if the user calls
   * set_spatial_dimension(3) and then inserts only x-aligned 1D
   * elements into the Mesh, mesh.spatial_dimension() will remain 3.
   */
  unsigned int spatial_dimension () const;

  /**
   * Sets the "spatial dimension" of the Mesh.  See the documentation
   * for Mesh::spatial_dimension() for more information.
   */
  void set_spatial_dimension(unsigned char d);

  /**
   * \returns The number of nodes in the mesh.
   *
   * This function and others must be defined in derived classes since
   * the MeshBase class has no specific storage for nodes or elements.
   * The standard \p n_nodes() function may return a cached value on
   * distributed meshes, and so can be called by any processor at any
   * time.
   */
  virtual dof_id_type n_nodes () const = 0;

  /**
   * \returns The number of nodes in the mesh.
   *
   * This function and others must be overridden in derived classes since
   * the MeshBase class has no specific storage for nodes or elements.
   * The \p parallel_n_nodes() function computes a parallel-synchronized
   * value on distributed meshes, and so must be called in parallel
   * only.
   */
  virtual dof_id_type parallel_n_nodes () const = 0;

  /**
   * \returns The number of nodes on processor \p proc.
   */
  dof_id_type n_nodes_on_proc (const processor_id_type proc) const;

  /**
   * \returns The number of nodes on the local processor.
   */
  dof_id_type n_local_nodes () const
  { return this->n_nodes_on_proc (this->processor_id()); }

  /**
   * \returns The number of nodes owned by no processor.
   */
  dof_id_type n_unpartitioned_nodes () const
  { return this->n_nodes_on_proc (DofObject::invalid_processor_id); }

  /**
   * \returns A number one greater than the maximum node id in the
   * mesh. A more apt name for this method would be end_node_id
   */
  virtual dof_id_type max_node_id () const = 0;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * \returns The next unique id to be used.
   */
  unique_id_type next_unique_id() const { return _next_unique_id; }

  /**
   * Sets the next available unique id to be used.  On a
   * ReplicatedMesh, or when adding unpartitioned objects to a
   * DistributedMesh, this must be kept in sync on all processors.
   *
   * On a DistributedMesh, other unique_id values (larger than this
   * one) may be chosen next, to allow unique_id assignment without
   * communication.
   */
  virtual void set_next_unique_id(unique_id_type id) = 0;
#endif

  /**
   * Reserves space for a known number of nodes.
   *
   * \note This method may or may not do anything, depending on the
   * actual \p Mesh implementation.  If you know the number of nodes
   * you will add and call this method before repeatedly calling \p
   * add_point() the implementation will be more efficient.
   */
  virtual void reserve_nodes (const dof_id_type nn) = 0;

  /**
   * \returns The number of elements in the mesh.
   *
   * The standard n_elem() function may return a cached value on
   * distributed meshes, and so can be called by any processor at any
   * time.
   */
  virtual dof_id_type n_elem () const = 0;

  /**
   * \returns The number of elements in the mesh.
   *
   * The parallel_n_elem() function computes a parallel-synchronized
   * value on distributed meshes, and so must be called in parallel
   * only.
   */
  virtual dof_id_type parallel_n_elem () const = 0;

  /**
   * \returns A number one greater than the maximum element id in the
   * mesh. A more apt name for this method would be end_elem_id
   */
  virtual dof_id_type max_elem_id () const = 0;

  /**
   * \returns A number greater than or equal to the maximum unique_id in the
   * mesh.
   */
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const = 0;
#endif

  /**
   * Reserves space for a known number of elements.
   *
   * \note This method may or may not do anything, depending on the
   * actual \p Mesh implementation.  If you know the number of
   * elements you will add and call this method before repeatedly
   * calling \p add_point() the implementation will be more efficient.
   */
  virtual void reserve_elem (const dof_id_type ne) = 0;

  /**
   * Updates parallel caches so that methods like n_elem()
   * accurately reflect changes on other processors
   */
  virtual void update_parallel_id_counts () = 0;

  /**
   * \returns The number of active elements in the mesh.
   *
   * Implemented in terms of active_element_iterators.
   */
  virtual dof_id_type n_active_elem () const = 0;

  /**
   * \returns The number of elements on processor \p proc.
   */
  dof_id_type n_elem_on_proc (const processor_id_type proc) const;

  /**
   * \returns The number of elements on the local processor.
   */
  dof_id_type n_local_elem () const
  { return this->n_elem_on_proc (this->processor_id()); }

  /**
   * \returns The number of elements owned by no processor.
   */
  dof_id_type n_unpartitioned_elem () const
  { return this->n_elem_on_proc (DofObject::invalid_processor_id); }

  /**
   * \returns The number of active elements on processor \p proc.
   */
  dof_id_type n_active_elem_on_proc (const processor_id_type proc) const;

  /**
   * \returns The number of active elements on the local processor.
   */
  dof_id_type n_active_local_elem () const
  { return this->n_active_elem_on_proc (this->processor_id()); }

  /**
   * \returns The number of elements that will be written
   * out in certain I/O formats.
   *
   * For example, a 9-noded quadrilateral will be broken into 4 linear
   * sub-elements for plotting purposes.  Thus, for a mesh of 2 \p
   * QUAD9 elements \p n_tecplot_elem() will return 8.  Implemented in
   * terms of element_iterators.
   */
  dof_id_type n_sub_elem () const;

  /**
   * Same as \p n_sub_elem(), but only counts active elements.
   */
  dof_id_type n_active_sub_elem () const;

  /**
   * \returns A constant reference (for reading only) to the
   * \f$ i^{th} \f$ point, which should be present in this processor's
   * subset of the mesh data structure.
   */
  virtual const Point & point (const dof_id_type i) const = 0;

  /**
   * \returns A constant reference (for reading only) to the
   * \f$ i^{th} \f$ node, which should be present in this processor's
   * subset of the mesh data structure.
   */
  virtual const Node & node_ref (const dof_id_type i) const
  {
    return *this->node_ptr(i);
  }

  /**
   * \returns A reference to the \f$ i^{th} \f$ node, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual Node & node_ref (const dof_id_type i)
  {
    return *this->node_ptr(i);
  }

  /**
   * \returns A pointer to the \f$ i^{th} \f$ node, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Node * node_ptr (const dof_id_type i) const = 0;

  /**
   * \returns A writable pointer to the \f$ i^{th} \f$ node, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Node * node_ptr (const dof_id_type i) = 0;

  /**
   * \returns A pointer to the \f$ i^{th} \f$ node, or \p nullptr if no such
   * node exists in this processor's mesh data structure.
   */
  virtual const Node * query_node_ptr (const dof_id_type i) const = 0;

  /**
   * \returns A writable pointer to the \f$ i^{th} \f$ node, or \p nullptr if
   * no such node exists in this processor's mesh data structure.
   */
  virtual Node * query_node_ptr (const dof_id_type i) = 0;

  /**
   * \returns A reference to the \f$ i^{th} \f$ element, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Elem & elem_ref (const dof_id_type i) const
  {
    return *this->elem_ptr(i);
  }

  /**
   * \returns A writable reference to the \f$ i^{th} \f$ element, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Elem & elem_ref (const dof_id_type i)
  {
    return *this->elem_ptr(i);
  }

  /**
   * \returns A pointer to the \f$ i^{th} \f$ element, which should be
   * present in this processor's subset of the mesh data structure.
   */
  virtual const Elem * elem_ptr (const dof_id_type i) const = 0;

  /**
   * \returns A writable pointer to the \f$ i^{th} \f$ element, which
   * should be present in this processor's subset of the mesh data
   * structure.
   */
  virtual Elem * elem_ptr (const dof_id_type i) = 0;

  /**
   * \returns A pointer to the \f$ i^{th} \f$ element, or nullptr if no
   * such element exists in this processor's mesh data structure.
   */
  virtual const Elem * query_elem_ptr (const dof_id_type i) const = 0;

  /**
   * \returns A writable pointer to the \f$ i^{th} \f$ element, or nullptr
   * if no such element exists in this processor's mesh data structure.
   */
  virtual Elem * query_elem_ptr (const dof_id_type i) = 0;

  /**
   * Add a new \p Node at \p Point \p p to the end of the vertex array,
   * with processor_id \p procid.
   * Use DofObject::invalid_processor_id (default) to add a node to all
   * processors, or this->processor_id() to add a node to the local
   * processor only.
   * If adding a node locally, passing an \p id other than
   * DofObject::invalid_id will set that specific node id.  Only
   * do this in parallel if you are manually keeping ids consistent.
   */
  virtual Node * add_point (const Point & p,
                            const dof_id_type id = DofObject::invalid_id,
                            const processor_id_type proc_id =
                            DofObject::invalid_processor_id) = 0;

  /**
   * Add \p Node \p n to the end of the vertex array.
   */
  virtual Node * add_node (Node * n) = 0;

  /**
   * Version of add_node() taking a std::unique_ptr by value. The version
   * taking a dumb pointer will eventually be deprecated in favor of this
   * version. This API is intended to indicate that ownership of the Node
   * is transferred to the Mesh when this function is called, and it should
   * play more nicely with the Node::build() API which has always returned
   * a std::unique_ptr.
   */
  virtual Node * add_node (std::unique_ptr<Node> n) = 0;

  /**
   * Removes the Node n from the mesh.
   */
  virtual void delete_node (Node * n) = 0;

  /**
   * Takes ownership of node \p n on this partition of a distributed
   * mesh, by setting n.processor_id() to this->processor_id(), as
   * well as changing n.id() and moving it in the mesh's internal
   * container to give it a new authoritative id.
   */
  virtual void own_node (Node &) {}

  /**
   * Changes the id of node \p old_id, both by changing node(old_id)->id() and
   * by moving node(old_id) in the mesh's internal container.  No element with
   * the id \p new_id should already exist.
   */
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id) = 0;

  /**
   * Add elem \p e to the end of the element array.
   * To add an element locally, set e->processor_id() before adding it.
   * To ensure a specific element id, call e->set_id() before adding it;
   * only do this in parallel if you are manually keeping ids consistent.
   *
   * Users should call MeshBase::complete_preparation() after elements are
   * added to and/or deleted from the mesh.
   */
  virtual Elem * add_elem (Elem * e) = 0;

  /**
   * Version of add_elem() taking a std::unique_ptr by value. The version
   * taking a dumb pointer will eventually be deprecated in favor of this
   * version. This API is intended to indicate that ownership of the Elem
   * is transferred to the Mesh when this function is called, and it should
   * play more nicely with the Elem::build() API which has always returned
   * a std::unique_ptr.
   */
  virtual Elem * add_elem (std::unique_ptr<Elem> e) = 0;

  /**
   * Insert elem \p e to the element array, preserving its id
   * and replacing/deleting any existing element with the same id.
   *
   * Users should call MeshBase::complete_preparation() after elements are
   * added to and/or deleted from the mesh.
   */
  virtual Elem * insert_elem (Elem * e) = 0;

  /**
   * Version of insert_elem() taking a std::unique_ptr by value. The version
   * taking a dumb pointer will eventually be deprecated in favor of this
   * version. This API is intended to indicate that ownership of the Elem
   * is transferred to the Mesh when this function is called, and it should
   * play more nicely with the Elem::build() API which has always returned
   * a std::unique_ptr.
   */
  virtual Elem * insert_elem (std::unique_ptr<Elem> e) = 0;

  /**
   * Removes element \p e from the mesh. This method must be
   * implemented in derived classes in such a way that it does not
   * invalidate element iterators.  Users should call
   * MeshBase::complete_preparation() after elements are added to
   * and/or deleted from the mesh.
   *
   * \note Calling this method may produce isolated nodes, i.e. nodes
   * not connected to any element.
   */
  virtual void delete_elem (Elem * e) = 0;

  /**
   * Changes the id of element \p old_id, both by changing elem(old_id)->id()
   * and by moving elem(old_id) in the mesh's internal container.  No element
   * with the id \p new_id should already exist.
   */
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id) = 0;

  /**
   * Returns the default master space to physical space mapping basis
   * functions to be used on newly added elements.
   */
  ElemMappingType default_mapping_type () const
  {
    return _default_mapping_type;
  }

  /**
   * Set the default master space to physical space mapping basis
   * functions to be used on newly added elements.
   */
  void set_default_mapping_type (const ElemMappingType type)
  {
    _default_mapping_type = type;
  }

  /**
   * Returns any default data value used by the master space to
   * physical space mapping.
   */
  unsigned char default_mapping_data () const
  {
    return _default_mapping_data;
  }

  /**
   * Set the default master space to physical space mapping basis
   * functions to be used on newly added elements.
   */
  void set_default_mapping_data (const unsigned char data)
  {
    _default_mapping_data = data;
  }

  /**
   * Locate element face (edge in 2D) neighbors.  This is done with the help
   * of a \p std::map that functions like a hash table.
   * After this routine is called all the elements with a \p nullptr neighbor
   * pointer are guaranteed to be on the boundary.  Thus this routine is
   * useful for automatically determining the boundaries of the domain.
   * If reset_remote_elements is left to false, remote neighbor links are not
   * reset and searched for in the local mesh.  If reset_current_list is
   * left as true, then any existing links will be reset before initiating
   * the algorithm, while honoring the value of the reset_remote_elements
   * flag.
   */
  virtual void find_neighbors (const bool reset_remote_elements = false,
                               const bool reset_current_list    = true) = 0;

  /**
   * Removes any orphaned nodes, nodes not connected to any elements.
   * Typically done automatically in a preparation step
   */
  void remove_orphaned_nodes ();

  /**
   * After partitioning a mesh it is useful to renumber the nodes and elements
   * so that they lie in contiguous blocks on the processors.  This method
   * does just that.
   */
  virtual void renumber_nodes_and_elements () = 0;

  /**
   * There is no reason for a user to ever call this function.
   *
   * This function restores a previously broken element/node numbering such that
   * \p mesh.node_ref(n).id() == n.
   */
  virtual void fix_broken_node_and_element_numbering () = 0;


#ifdef LIBMESH_ENABLE_AMR
  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract () = 0;
#endif

  /**
   * Register an integer datum (of type dof_id_type) to be added to
   * each element in the mesh.
   *
   * If the mesh already has elements, data by default is allocated in
   * each of them.  This may be expensive to do repeatedly; use
   * add_elem_integers instead.  Alternatively, the \p allocate_data
   * option can be manually set to false, but if this is done then a
   * manual call to \p size_elem_extra_integers() will need to be done
   * before the new space is usable.
   *
   * Newly allocated values for the new datum will be initialized to
   * \p default_value
   *
   * \returns The index number for the new datum, or for the existing
   * datum if one by the same name has already been added.
   */
  unsigned int add_elem_integer(std::string name,
                                bool allocate_data = true,
                                dof_id_type default_value = DofObject::invalid_id);

  /**
   * Register integer data (of type dof_id_type) to be added to
   * each element in the mesh, one string name for each new integer.
   *
   * If the mesh already has elements, data by default is allocated in
   * each of them.
   *
   * Newly allocated values for the new datum with name \p names[i]
   * will be initialized to \p default_values[i], or to
   * DofObject::invalid_id if \p default_values is null.
   *
   * \returns The index numbers for the new data, and/or for existing
   * data if data by some of the same names has already been added.
   */
  std::vector<unsigned int> add_elem_integers(const std::vector<std::string> & names,
                                              bool allocate_data = true,
                                              const std::vector<dof_id_type> * default_values = nullptr);

  /*
   * \returns The index number for the named extra element integer
   * datum, which must have already been added.
   */
  unsigned int get_elem_integer_index(std::string_view name) const;

  /*
   * \returns Whether or not the mesh has an element integer with its name.
   */
  bool has_elem_integer(std::string_view name) const;

  /*
   * \returns The name for the indexed extra element integer
   * datum, which must have already been added.
   */
  const std::string & get_elem_integer_name(unsigned int i) const
  { return _elem_integer_names[i]; }

  /*
   * \returns The number of extra element integers for which space is
   * being reserved on this mesh.
   *
   * If non-integer data has been associated, each datum of type T
   * counts for sizeof(T)/sizeof(dof_id_type) times in the return
   * value.
   */
  unsigned int n_elem_integers() const { return _elem_integer_names.size(); }

  /**
   * Register a datum (of type T) to be added to each element in the
   * mesh.
   *
   * If the mesh already has elements, data by default is allocated in
   * each of them.  This may be expensive to do repeatedly; use
   * add_elem_data instead.  Alternatively, the \p allocate_data
   * option can be manually set to false, but if this is done then a
   * manual call to \p size_elem_extra_integers() will need to be done
   * before the new space is usable.
   *
   * Newly allocated values for the new datum will be initialized to
   * \p *default_value if \p default_value is not null, or to
   * meaningless memcpy output otherwise.
   *
   * \returns The index numbers for the new data, and/or for existing
   * data if data by some of the same names has already been added.
   *
   * If type T is larger than dof_id_type, its data will end up
   * spanning multiple index values, but will be queried with the
   * starting index number.
   *
   * No type checking is done with this function!  If you add data of
   * type T, don't try to access it with a call specifying type U.
   */
  template <typename T>
  unsigned int add_elem_datum(const std::string & name,
                              bool allocate_data = true,
                              const T * default_value = nullptr);

  /**
   * Register data (of type T) to be added to each element in the
   * mesh.
   *
   * If the mesh already has elements, data is allocated in each.
   *
   * Newly allocated values for the new datum with name \p names[i]
   * will be initialized to \p default_values[i], or to
   * meaningless memcpy output if \p default_values is null.
   *
   * \returns The starting index number for the new data, or for the
   * existing data if one by the same name has already been added.
   *
   * If type T is larger than dof_id_type, each datum will end up
   * spanning multiple index values, but will be queried with the
   * starting index number.
   *
   * No type checking is done with this function!  If you add data of
   * type T, don't try to access it with a call specifying type U.
   */
  template <typename T>
  std::vector<unsigned int> add_elem_data(const std::vector<std::string> & names,
                                          bool allocate_data = true,
                                          const std::vector<T> * default_values = nullptr);

  /**
   * Register an integer datum (of type dof_id_type) to be added to
   * each node in the mesh.
   *
   * If the mesh already has nodes, data by default is allocated in
   * each of them.  This may be expensive to do repeatedly; use
   * add_node_integers instead.  Alternatively, the \p allocate_data
   * option can be manually set to false, but if this is done then a
   * manual call to \p size_node_extra_integers() will need to be done
   * before the new space is usable.
   *
   * Newly allocated values for the new datum will be initialized to
   * \p default_value
   *
   * \returns The index number for the new datum, or for the existing
   * datum if one by the same name has already been added.
   */
  unsigned int add_node_integer(std::string name,
                                bool allocate_data = true,
                                dof_id_type default_value = DofObject::invalid_id);

  /**
   * Register integer data (of type dof_id_type) to be added to
   * each node in the mesh.
   *
   * If the mesh already has nodes, data by default is allocated in
   * each.
   *
   * Newly allocated values for the new datum with name \p names[i]
   * will be initialized to \p default_values[i], or to
   * DofObject::invalid_id if \p default_values is null.
   *
   * \returns The index numbers for the new data, and/or for existing
   * data if data by some of the same names has already been added.
   */
  std::vector<unsigned int> add_node_integers(const std::vector<std::string> & names,
                                              bool allocate_data = true,
                                              const std::vector<dof_id_type> * default_values = nullptr);

  /*
   * \returns The index number for the named extra node integer
   * datum, which must have already been added.
   */
  unsigned int get_node_integer_index(std::string_view name) const;

  /*
   * \returns Whether or not the mesh has a node integer with its name.
   */
  bool has_node_integer(std::string_view name) const;

  /*
   * \returns The name for the indexed extra node integer
   * datum, which must have already been added.
   */
  const std::string & get_node_integer_name(unsigned int i) const
  { return _node_integer_names[i]; }

  /*
   * \returns The number of extra node integers for which space is
   * being reserved on this mesh.
   *
   * If non-integer data has been associated, each datum of type T
   * counts for sizeof(T)/sizeof(dof_id_type) times in the return
   * value.
   */
  unsigned int n_node_integers() const { return _node_integer_names.size(); }

  /**
   * Register a datum (of type T) to be added to each node in the
   * mesh.
   *
   * If the mesh already has nodes, data by default is allocated in
   * each of them.  This may be expensive to do repeatedly; use
   * add_node_data instead.  Alternatively, the \p allocate_data
   * option can be manually set to false, but if this is done then a
   * manual call to \p size_node_extra_integers() will need to be done
   * before the new space is usable.
   *
   * Newly allocated values for the new datum will be initialized to
   * \p *default_value if \p default_value is not null, or to
   * meaningless memcpy output otherwise.
   *
   * \returns The starting index number for the new datum, or for the
   * existing datum if one by the same name has already been added.
   *
   * If type T is larger than dof_id_type, its data will end up
   * spanning multiple index values, but will be queried with the
   * starting index number.
   *
   * No type checking is done with this function!  If you add data of
   * type T, don't try to access it with a call specifying type U.
   */
  template <typename T>
  unsigned int add_node_datum(const std::string & name,
                              bool allocate_data = true,
                              const T * default_value = nullptr);

  /**
   * Register data (of type T) to be added to each node in the
   * mesh.
   *
   * If the mesh already has nodes, data by default is allocated in each.
   *
   * Newly allocated values for the new datum with name \p names[i]
   * will be initialized to \p default_values[i], or to
   * meaningless memcpy output if \p default_values is null.
   *
   * \returns The starting index number for the new data, or for the
   * existing data if one by the same name has already been added.
   *
   * If type T is larger than dof_id_type, its data will end up
   * spanning multiple index values, but will be queried with the
   * starting index number.
   *
   * No type checking is done with this function!  If you add data of
   * type T, don't try to access it with a call specifying type U.
   */
  template <typename T>
  std::vector<unsigned int> add_node_data(const std::vector<std::string> & name,
                                          bool allocate_data = true,
                                          const std::vector<T> * default_values = nullptr);

  /**
   * Prepare a newly created (or read) mesh for use.
   * This involves several steps:
   *  1.) renumbering (if enabled)
   *  2.) removing any orphaned nodes
   *  3.) updating parallel id counts
   *  4.) finding neighbor links
   *  5.) caching summarized element data
   *  6.) finding interior parent links
   *  7.) clearing any old point locator
   *  8.) calling reinit() on ghosting functors
   *  9.) repartitioning (if enabled)
   *  10.) removing any remote elements (if enabled)
   *  11.) regenerating summarized boundary id sets
   *
   * The argument to skip renumbering is now deprecated - to prevent a
   * mesh from being renumbered, set allow_renumbering(false). The argument to skip
   * finding neighbors is also deprecated. To prevent find_neighbors, set
   * allow_find_neighbors(false)
   *
   * If this is a distributed mesh, local copies of remote elements
   * will be deleted here - to keep those elements replicated during
   * preparation, set allow_remote_element_removal(false).
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void prepare_for_use (const bool skip_renumber_nodes_and_elements, const bool skip_find_neighbors);
  void prepare_for_use (const bool skip_renumber_nodes_and_elements);
#endif // LIBMESH_ENABLE_DEPRECATED
  void prepare_for_use ();

  /*
   * Prepare a newly created or modified mesh for use.
   *
   * Unlike \p prepare_for_use(), \p complete_preparation() performs
   * *only* those preparatory steps that have been marked as
   * necessary in the MeshBase::Preparation state.
   */
  void complete_preparation();

  /**
   * Call the default partitioner (currently \p metis_partition()).
   */
  virtual void partition (const unsigned int n_parts);

  void partition ()
  { this->partition(this->n_processors()); }

  /**
   * Redistribute elements between processors.  This gets called
   * automatically by the Partitioner, and merely notifies any
   * GhostingFunctors of redistribution in the case of a
   * ReplicatedMesh or serialized DistributedMesh
   */
  virtual void redistribute ();

  /**
   * Recalculate any cached data after elements and nodes have been
   * repartitioned.
   */
  virtual void update_post_partitioning () {}

  /**
   * If false is passed in then this mesh will no longer be renumbered
   * when being prepared for use.  This may slightly adversely affect
   * performance during subsequent element access, particularly when
   * using a distributed mesh.
   *
   * Important! When allow_renumbering(false) is set,
   * ReplicatedMesh::n_elem() and ReplicatedMesh::n_nodes() will
   * return *wrong* values whenever adaptive refinement is followed by
   * adaptive coarsening. (Uniform refinement followed by uniform
   * coarsening is OK.) This is due to the fact that n_elem() and
   * n_nodes() are currently O(1) functions that just return the size
   * of the respective underlying vectors, and this size is wrong when
   * the numbering includes "gaps" from nodes and elements that have
   * been deleted. We plan to implement a caching mechanism in the
   * near future that will fix this incorrect behavior.
   */
  void allow_renumbering(bool allow) { _skip_renumber_nodes_and_elements = !allow; }
  bool allow_renumbering() const { return !_skip_renumber_nodes_and_elements; }

  /**
   * If \p false is passed then this mesh will no longer work to find element
   * neighbors when being prepared for use
   */
  void allow_find_neighbors(bool allow) { _skip_find_neighbors = !allow; }
  bool allow_find_neighbors() const { return !_skip_find_neighbors; }

  /**
   * If false is passed in then this mesh will no longer have remote
   * elements deleted when being prepared for use; i.e. even a
   * DistributedMesh will remain (if it is already) serialized.
   * This may adversely affect performance and memory use.
   */
  void allow_remote_element_removal(bool allow) { _allow_remote_element_removal = allow; }
  bool allow_remote_element_removal() const { return _allow_remote_element_removal; }

  /**
   * If \p true is passed, then this mesh will no longer require
   * unique_ids to be unique across the set of all DofObjects. That
   * is, although no two Elems (resp. Nodes) will share the same
   * unique_id, a given Elem and Node might share the same unique_id.
   */
  void allow_node_and_elem_unique_id_overlap(bool allow) { _allow_node_and_elem_unique_id_overlap = allow; }
  bool allow_node_and_elem_unique_id_overlap() const { return _allow_node_and_elem_unique_id_overlap; }

  /**
   * If true is passed in then the elements on this mesh will no
   * longer be (re)partitioned, and the nodes on this mesh will only
   * be repartitioned if they are found "orphaned" via coarsening or
   * other removal of the last element responsible for their
   * node/element processor id consistency.
   *
   * \note It would probably be a bad idea to call this on a
   * DistributedMesh _before_ the first partitioning has happened...
   * because no elements would get assigned to your processor pool.
   *
   * \note Skipping partitioning can have adverse effects on your
   * performance when using AMR... i.e. you could get large load
   * imbalances.  However you might still want to use this if the
   * communication and computation of the rebalance and repartition is
   * too high for your application.
   *
   * It is also possible, for backwards-compatibility purposes, to
   * skip noncritical partitioning by resetting the partitioner()
   * pointer for this mesh.
   */
  void skip_noncritical_partitioning(bool skip)
  { _skip_noncritical_partitioning = skip; }

  bool skip_noncritical_partitioning() const
  { return _skip_noncritical_partitioning || _skip_all_partitioning || !_partitioner.get(); }


  /**
   * If true is passed in then nothing on this mesh will be
   * (re)partitioned.
   *
   * \note The caveats for skip_noncritical_partitioning() still
   * apply, and removing elements from a mesh with this setting
   * enabled can leave node processor ids in an inconsistent state
   * (not matching any attached element), causing failures in other
   * library code.  Do not use this setting along with element
   * deletion or coarsening.
   */
  void skip_partitioning(bool skip) { _skip_all_partitioning = skip; }

  bool skip_partitioning() const { return _skip_all_partitioning; }

  /**
   * Adds a functor which can specify ghosting requirements for use on
   * distributed meshes.  Multiple ghosting functors can be added; any
   * element which is required by any functor will be ghosted.
   *
   * GhostingFunctor memory must be managed by the code which calls
   * this function; the GhostingFunctor lifetime is expected to extend
   * until either the functor is removed or the Mesh is destructed.
   */
  void add_ghosting_functor(GhostingFunctor & ghosting_functor);

  /**
   * Adds a functor which can specify ghosting requirements for use on
   * distributed meshes.  Multiple ghosting functors can be added; any
   * element which is required by any functor will be ghosted.
   *
   * GhostingFunctor memory when using this method is managed by the
   * shared_ptr mechanism.
   */
  void add_ghosting_functor(std::shared_ptr<GhostingFunctor> ghosting_functor)
  { _shared_functors[ghosting_functor.get()] = ghosting_functor;
    this->add_ghosting_functor(*ghosting_functor); }

  /**
   * Removes a functor which was previously added to the set of
   * ghosting functors.
   */
  void remove_ghosting_functor(GhostingFunctor & ghosting_functor);

  /**
   * Iterator type for ghosting functor ranges.  This has changed in
   * the past and may change again; code should use auto or the type
   * here.
   */
  typedef std::vector<GhostingFunctor *>::const_iterator GhostingFunctorIterator;

  /**
   * Beginning of range of ghosting functors
   */
  GhostingFunctorIterator ghosting_functors_begin() const
  { return _ghosting_functors.begin(); }

  /**
   * End of range of ghosting functors
   */
  GhostingFunctorIterator ghosting_functors_end() const
  { return _ghosting_functors.end(); }

  /**
   * Default ghosting functor
   */
  GhostingFunctor & default_ghosting() { return *_default_ghosting; }

  /**
   * Constructs a list of all subdomain identifiers in the local mesh if
   * \p global == false, and in the global mesh if \p global == true (default).
   * Subdomains correspond to separate subsets of the mesh which could correspond
   * e.g. to different materials in a solid mechanics application,
   * or regions where different physical processes are important.  The subdomain
   * mapping is independent from the parallel decomposition.
   *
   * Unpartitioned elements are included in the set in the case that \p
   * global == true. If \p global == false, the unpartitioned elements are not
   * included because unpartitioned elements do not have a sense of locality.
   */
  void subdomain_ids (std::set<subdomain_id_type> & ids, const bool global = true) const;

  /**
   * \returns The number of subdomains in the global mesh. Subdomains correspond
   * to separate subsets of the mesh which could correspond e.g. to different
   * materials in a solid mechanics application, or regions where different
   * physical processes are important.  The subdomain mapping is independent
   * from the parallel decomposition.
   */
  subdomain_id_type n_subdomains () const;

  /**
   * \returns The number of subdomains in the local mesh. Subdomains correspond
   * to separate subsets of the mesh which could correspond e.g. to different
   * materials in a solid mechanics application, or regions where different
   * physical processes are important.  The subdomain mapping is independent
   * from the parallel decomposition.
   */
  subdomain_id_type n_local_subdomains () const;

  /**
   * \returns The number of partitions which have been defined via
   * a call to either mesh.partition() or by building a Partitioner
   * object and calling partition.
   *
   * \note The partitioner object is responsible for setting this
   * value.
   */
  unsigned int n_partitions () const
  { return _n_parts; }

  /**
   * \returns A string containing relevant information
   * about the mesh.
   *
   * \p verbosity sets the verbosity, with 0 being the least and 2 being the greatest.
   * 0 - Dimensions, number of nodes, number of elems, number of subdomains, number of
   *     partitions, prepared status.
   * 1 - Adds the mesh bounding box, mesh element types, specific nodesets/edgesets/sidesets
   *     with element types, number of nodes/edges/sides.
   * 2 - Adds volume information and bounding boxes to boundary information.
   *
   * The \p global parameter pertains primarily to verbosity levels 1 and above.
   * When \p global == true, information is only output on rank 0 and the information
   * is reduced. When \p global == false, information is output on all ranks that pertains
   * only to that local partition.
   */
  std::string get_info (const unsigned int verbosity = 0, const bool global = true) const;

  /**
   * Prints relevant information about the mesh.
   *
   * Take note of the docstring for get_info() for more information pretaining to
   * the \p verbosity and \p global parameters.
   */
  void print_info (std::ostream & os=libMesh::out, const unsigned int verbosity = 0, const bool global = true) const;

  /**
   * Equivalent to calling print_info() above, but now you can write:
   * Mesh mesh;
   * libMesh::out << mesh << std::endl;
   */
  friend std::ostream & operator << (std::ostream & os, const MeshBase & m);

  /**
   * Interfaces for reading/writing a mesh to/from a file.  Must be
   * implemented in derived classes.
   */
  virtual void read  (const std::string & name,
                      void * mesh_data=nullptr,
                      bool skip_renumber_nodes_and_elements=false,
                      bool skip_find_neighbors=false) = 0;
  virtual void write (const std::string & name) const = 0;

  /**
   * Converts a mesh with higher-order
   * elements into a mesh with linear elements.  For
   * example, a mesh consisting of \p Tet10 will be converted
   * to a mesh with \p Tet4 etc.
   */
  virtual void all_first_order () = 0;

  /**
   * We need an empty, generic class to act as a predicate for this
   * and derived mesh classes.
   */
  typedef Predicates::multi_predicate Predicate;

  /**
   * structs for the element_iterator's.
   *
   * \note These iterators were designed so that derived mesh classes
   * could use the _same_ base class iterators interchangeably.  Their
   * definition comes later in the header file.
   */
  struct element_iterator;
  struct const_element_iterator;

  /**
   * structs for the node_iterator's.
   *
   * \note These iterators were designed so that derived mesh classes
   * could use the _same_ base class iterators interchangeably.  Their
   * definition comes later in the header file.
   */
  struct node_iterator;
  struct const_node_iterator;

  /**
   * Converts a set of this Mesh's elements defined by \p range from
   * FIRST order to SECOND order. Must be called on conforming,
   * non-refined meshes. For example, a mesh consisting of \p Tet4
   * will be converted to a mesh with \p Tet10 etc.
   *
   * \note For some elements like \p Hex8 there exist two higher order
   * equivalents, \p Hex20 and \p Hex27.  When \p full_ordered is \p
   * true (default), then \p Hex27 is built.  Otherwise, \p Hex20 is
   * built.  The same holds obviously for \p Quad4, \p Prism6, etc.
   */
  virtual void all_second_order_range(const SimpleRange<element_iterator> & range,
                                      const bool full_ordered = true) = 0;

  /**
   * Calls the range-based version of this function with a range
   * consisting of all elements in the mesh.
   */
  void all_second_order (const bool full_ordered = true);

  /**
   * Converts a set of elements in this (conforming, non-refined) mesh
   * into "complete" order elements, i.e. elements which
   * can store degrees of freedom on any vertex, edge, or face.  For
   * example, a mesh consisting of \p Tet4 or \p Tet10 will be
   * converted to a mesh with \p Tet14 etc.
   */
  virtual void all_complete_order_range(const SimpleRange<element_iterator> & range) = 0;

  /**
   * Calls the range-based version of this function with a range
   * consisting of all elements in the mesh.
   */
  virtual void all_complete_order ();

  /**
   * In a few (very rare) cases, the user may have manually tagged the
   * elements with specific processor IDs by hand, without using a
   * partitioner.  In this case, the Mesh will not know that the total
   * number of partitions, _n_parts, has changed, unless you call this
   * function.  This is an O(N active elements) calculation.  The return
   * value is the number of partitions, and _n_parts is also set by
   * this function.
   */
  unsigned int recalculate_n_partitions();

  /**
   * \returns A pointer to a subordinate \p PointLocatorBase object
   * for this mesh, constructing a master PointLocator first if
   * necessary.  This should not be used in threaded or
   * non-parallel_only code unless the master has already been
   * constructed.
   */
  std::unique_ptr<PointLocatorBase> sub_point_locator () const;

  /**
   * Set value used by PointLocatorBase::close_to_point_tol().
   *
   * Defaults to 0.0. If nonzero, calls close_to_point_tol() whenever
   * a new PointLocator is built for use by this Mesh.  Since the Mesh
   * controls the creation and destruction of the PointLocator, if
   * there are any parameters we need to customize on it, the Mesh
   * will need to know about them.
   */
  void set_point_locator_close_to_point_tol(Real val);
  Real get_point_locator_close_to_point_tol() const;

  /**
   * Releases the current \p PointLocator object.
   */
  void clear_point_locator ();

  /**
   * In the point locator, do we count lower dimensional elements
   * when we refine point locator regions? This is relevant in
   * tree-based point locators, for example.
   */
  void set_count_lower_dim_elems_in_point_locator(bool count_lower_dim_elems);

  /**
   * Get the current value of _count_lower_dim_elems_in_point_locator.
   */
  bool get_count_lower_dim_elems_in_point_locator() const;

  /**
   * Verify id and processor_id consistency of our elements and
   * nodes containers.
   * Calls libmesh_assert() on each possible failure.
   * Currently only implemented on DistributedMesh; a serial data
   * structure is much harder to get out of sync.
   */
  virtual void libmesh_assert_valid_parallel_ids() const {}

  /**
   * \returns A writable reference for getting/setting an optional
   * name for a subdomain.
   */
  std::string & subdomain_name(subdomain_id_type id);
  const std::string & subdomain_name(subdomain_id_type id) const;

  /**
   * \returns The id of the named subdomain if it exists,
   * \p Elem::invalid_subdomain_id otherwise.
   */
  subdomain_id_type get_id_by_name(std::string_view name) const;

  /*
   * We have many combinations of iterators that filter on various
   * characteristics; we use macros to make their abstract base class
   * and their subclass declarations more terse.
   */
#define ABSTRACT_ELEM_ITERATORS(TYPE, ARGDECL) \
  virtual element_iterator TYPE##elements_begin(ARGDECL) = 0; \
  virtual element_iterator TYPE##elements_end(ARGDECL) = 0; \
  virtual const_element_iterator TYPE##elements_begin(ARGDECL) const = 0; \
  virtual const_element_iterator TYPE##elements_end(ARGDECL) const = 0; \
  virtual SimpleRange<element_iterator> TYPE##element_ptr_range(ARGDECL) = 0; \
  virtual SimpleRange<const_element_iterator> TYPE##element_ptr_range(ARGDECL) const = 0;

#define DECLARE_ELEM_ITERATORS(TYPE, ARGDECL, ARGS) \
  virtual element_iterator TYPE##elements_begin(ARGDECL) override final; \
  virtual element_iterator TYPE##elements_end(ARGDECL) override final; \
  virtual const_element_iterator TYPE##elements_begin(ARGDECL) const override final; \
  virtual const_element_iterator TYPE##elements_end(ARGDECL) const override final; \
  virtual SimpleRange<element_iterator> TYPE##element_ptr_range(ARGDECL) override final { return {TYPE##elements_begin(ARGS), TYPE##elements_end(ARGS)}; } \
  virtual SimpleRange<const_element_iterator> TYPE##element_ptr_range(ARGDECL) const override final  { return {TYPE##elements_begin(ARGS), TYPE##elements_end(ARGS)}; }

#define ABSTRACT_NODE_ITERATORS(TYPE, ARGDECL) \
  virtual node_iterator TYPE##nodes_begin(ARGDECL) = 0; \
  virtual node_iterator TYPE##nodes_end(ARGDECL) = 0; \
  virtual const_node_iterator TYPE##nodes_begin(ARGDECL) const = 0; \
  virtual const_node_iterator TYPE##nodes_end(ARGDECL) const = 0; \
  virtual SimpleRange<node_iterator> TYPE##node_ptr_range(ARGDECL) = 0; \
  virtual SimpleRange<const_node_iterator> TYPE##node_ptr_range(ARGDECL) const = 0;

#define DECLARE_NODE_ITERATORS(TYPE, ARGDECL, ARGS) \
  virtual node_iterator TYPE##nodes_begin(ARGDECL) override final; \
  virtual node_iterator TYPE##nodes_end(ARGDECL) override final; \
  virtual const_node_iterator TYPE##nodes_begin(ARGDECL) const override final; \
  virtual const_node_iterator TYPE##nodes_end(ARGDECL) const override final; \
  virtual SimpleRange<node_iterator> TYPE##node_ptr_range(ARGDECL) override final { return {TYPE##nodes_begin(ARGS), TYPE##nodes_end(ARGS)}; } \
  virtual SimpleRange<const_node_iterator> TYPE##node_ptr_range(ARGDECL) const override final  { return {TYPE##nodes_begin(ARGS), TYPE##nodes_end(ARGS)}; }

#define LIBMESH_COMMA ,

  /*
   * element_iterator accessors
   *
   * The basic elements_begin() and elements_end() iterators iterate
   * over all elements in a mesh, returning element pointers or const
   * element pointers when dereferenced (depending on whether the mesh
   * reference was const). range-for loops can be written using
   * element_ptr_range()
   *
   * Filtered versions of these iterators, which skip over all
   * elements not matching some predicate, are also available, by
   * adding a prefix to the methods above.  E.g. local_ (in a form
   * like local_elements_begin() or local_element_ptr_range()) will
   * iterate only over elements whose processor_id() is the current
   * processor, or active_ will iterate only over active elements even
   * if the mesh is refined, or active_local_ will iterate over
   * elements that are both active and local.  Negation forms such as
   * not_local_ also exist.
   *
   * For some iterator prefixes, such as type_, an argument is needed
   * for the filter; e.g. the ElemType to select for in that case.
   *
   * All valid prefixes and their corresponding arguments can be found
   * in the macro invocations below.
   */
  ABSTRACT_ELEM_ITERATORS(,)                // elements_begin(), element_ptr_range(): all elements
  ABSTRACT_ELEM_ITERATORS(active_,)         // Elem::active() == true
  ABSTRACT_ELEM_ITERATORS(ancestor_,)       // Elem::ancestor() == true
  ABSTRACT_ELEM_ITERATORS(subactive_,)      // Elem::subactive() == true
  ABSTRACT_ELEM_ITERATORS(local_,)          // Elem::processor_id() == this processor
  ABSTRACT_ELEM_ITERATORS(unpartitioned_,)  // Elem::processor_id() == invalid_processor_id
  ABSTRACT_ELEM_ITERATORS(facelocal_,)      // is on or has a neighbor on this processor
  ABSTRACT_ELEM_ITERATORS(level_,unsigned int level)   // Elem::level() == level
  ABSTRACT_ELEM_ITERATORS(pid_,processor_id_type pid)  // Elem::processor_id() == pid
  ABSTRACT_ELEM_ITERATORS(type_,ElemType type)         // Elem::type() == type

  ABSTRACT_ELEM_ITERATORS(active_subdomain_,subdomain_id_type sid) // active && Elem::subdomain_id() == sid
  ABSTRACT_ELEM_ITERATORS(active_subdomain_set_,std::set<subdomain_id_type> ss) // active && ss.contains(Elem::subdomain_id())

  // Iterators which use negations of filters described above
  ABSTRACT_ELEM_ITERATORS(not_active_,)
  ABSTRACT_ELEM_ITERATORS(not_ancestor_,)
  ABSTRACT_ELEM_ITERATORS(not_subactive_,)
  ABSTRACT_ELEM_ITERATORS(not_local_,)
  ABSTRACT_ELEM_ITERATORS(not_level_,unsigned int level)

  // Iterators which combine multiple of the filters described above
  ABSTRACT_ELEM_ITERATORS(active_local_,)
  ABSTRACT_ELEM_ITERATORS(active_not_local_,)
  ABSTRACT_ELEM_ITERATORS(active_unpartitioned_,)
  ABSTRACT_ELEM_ITERATORS(active_type_,ElemType type)
  ABSTRACT_ELEM_ITERATORS(active_pid_,processor_id_type pid)
  ABSTRACT_ELEM_ITERATORS(local_level_,unsigned int level)
  ABSTRACT_ELEM_ITERATORS(local_not_level_,unsigned int level)
  ABSTRACT_ELEM_ITERATORS(active_local_subdomain_,subdomain_id_type sid)
  ABSTRACT_ELEM_ITERATORS(active_local_subdomain_set_,std::set<subdomain_id_type> ss)

  // Backwards compatibility
  virtual SimpleRange<element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) = 0;
  virtual SimpleRange<const_element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) const = 0;
  virtual SimpleRange<element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) = 0;
  virtual SimpleRange<const_element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) const = 0;
  virtual SimpleRange<element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) = 0;
  virtual SimpleRange<const_element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const = 0;

  // Discouraged from use - these iterators use outdated
  // pre-GhostingFunctor definitions and should be renamed if not
  // deprecated
  ABSTRACT_ELEM_ITERATORS(semilocal_,)         // active && Elem::is_semilocal()
  ABSTRACT_ELEM_ITERATORS(ghost_,)             // active && Elem::is_semilocal() && not local discouraged
  ABSTRACT_ELEM_ITERATORS(active_semilocal_,)

  // solution can be evaluated, with the given DoF map, for the given
  // variable number, or for all variables by default
  ABSTRACT_ELEM_ITERATORS(evaluable_,const DofMap & dof_map LIBMESH_COMMA unsigned int var_num = libMesh::invalid_uint)

  // solution can be evaluated for all variables of all given DoF maps
  ABSTRACT_ELEM_ITERATORS(multi_evaluable_,std::vector<const DofMap *> dof_maps)

#ifdef LIBMESH_ENABLE_AMR
  ABSTRACT_ELEM_ITERATORS(flagged_,unsigned char rflag)  // Elem::refinement_flag() == rflag

  // Elem::refinement_flag() == rflag && Elem::processor_id() == pid
  ABSTRACT_ELEM_ITERATORS(flagged_pid_,unsigned char rflag LIBMESH_COMMA processor_id_type pid)
#endif

  /*
   * node_iterator accessors
   *
   * The basic nodes_begin() and nodes_end() iterators iterate
   * over all nodes in a mesh, returning node pointers or const
   * node pointers when dereferenced (depending on whether the mesh
   * reference was const). range-for loops can be written using
   * node_ptr_range()
   *
   * Filtered versions of these iterators, which skip over all
   * nodes not matching some predicate, are also available, by
   * adding a prefix to the methods above.  E.g. local_ (in a form
   * like local_nodes_begin() or local_node_ptr_range()) will
   * iterate only over nodes whose processor_id() is the current
   * processor.
   *
   * All valid prefixes and their corresponding arguments can be found
   * in the macro invocations below.
   */
  ABSTRACT_NODE_ITERATORS(,)          // nodes_begin(), node_ptr_range(): all nodes
  ABSTRACT_NODE_ITERATORS(active_,)   // Node::active() == true; i.e. Node::id() != invalid_id
  ABSTRACT_NODE_ITERATORS(local_,)    // Node::processor_id() == this processor
  ABSTRACT_NODE_ITERATORS(bnd_,)      // BoundaryInfo::n_boundary_ids(node) > 0
  ABSTRACT_NODE_ITERATORS(pid_,processor_id_type pid)  // Node::processor_id() == pid
  ABSTRACT_NODE_ITERATORS(bid_,boundary_id_type bid)   // BoundaryInfo::has_boundary_id(node, bid)

  // solution can be evaluated, with the given DoF map, for the given
  // variable number, or for all variables by default
  ABSTRACT_NODE_ITERATORS(evaluable_,const DofMap & dof_map LIBMESH_COMMA unsigned int var_num = libMesh::invalid_uint)

  // solution can be evaluated for all variables of all given DoF maps
  ABSTRACT_NODE_ITERATORS(multi_evaluable_,std::vector<const DofMap *> dof_maps)

  /**
   * \returns A writable reference to the whole subdomain name map
   */
  std::map<subdomain_id_type, std::string> & set_subdomain_name_map ()
  { return _block_id_to_name; }
  const std::map<subdomain_id_type, std::string> & get_subdomain_name_map () const
  { return _block_id_to_name; }

  typedef std::vector<std::pair<std::pair<const Elem *, unsigned int>, Real>> constraint_rows_mapped_type;
  typedef std::map<const Node *, constraint_rows_mapped_type> constraint_rows_type;

  /**
   * Constraint rows accessors
   */
  constraint_rows_type & get_constraint_rows()
  { return _constraint_rows; }

  const constraint_rows_type & get_constraint_rows() const
  { return _constraint_rows; }

  dof_id_type n_constraint_rows() const;

  /**
   * Copy the constraints from the other mesh to this mesh
   */
  void copy_constraint_rows(const MeshBase & other_mesh);

  /**
   * Copy the constraints from the given matrix to this mesh.  The
   * \p constraint_operator should be an mxn matrix, where
   * m == this->n_nodes() and the operator indexing matches the
   * current node indexing.  This may require users to disable mesh
   * renumbering in between loading a mesh file and loading a
   * constraint matrix which matches it.
   *
   * If any "constraint" rows in the matrix are unit vectors, the node
   * corresponding to that row index will be left unconstrained, and
   * will be used to constrain any other nodes which have a non-zero
   * in the column index of that unit vector.
   *
   * For each matrix column index which does not correspond to an
   * existing node, a new NodeElem will be added to the mesh on which
   * to store the new unconstrained degree(s) of freedom.
   *
   * If \p precondition_constraint_operator is true, then the values
   * of those new unconstrained degrees of freedom may be scaled to
   * improve the conditioning of typical PDE matrices integrated on
   * constrained mesh elements.
   *
   * \p T for the constraint_operator in this function should be \p
   * Real or \p Number ... and the data should be \p Real - we just
   * allow complex \p T for the sake of subclasses which have to be
   * configured and compiled with only one runtime option.
   */
  template <typename T>
  void copy_constraint_rows(const SparseMatrix<T> & constraint_operator,
                            bool precondition_constraint_operator = false);

  /**
   * Prints (from processor 0) all mesh constraint rows.  If \p
   * print_nonlocal is true, then each constraint is printed once for
   * each processor that knows about it, which may be useful for \p
   * DistributedMesh debugging.
   */
  void print_constraint_rows(std::ostream & os=libMesh::out,
                             bool print_nonlocal=false) const;

  /**
   * Gets a string reporting all mesh constraint rows local to
   * this processor.  If \p print_nonlocal is true, then nonlocal
   * constraints which are locally known are included.
   */
  std::string get_local_constraints(bool print_nonlocal=false) const;

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * \deprecated This method has ben replaced by \p cache_elem_data which
   * caches data in addition to elem dimensions (e.g. elem subdomain ids)
   * Search the mesh and cache the different dimensions of the elements
   * present in the mesh.  This is done in prepare_for_use(), but can
   * be done manually by other classes after major mesh modifications.
   */
  void cache_elem_dims();
#endif // LIBMESH_ENABLE_DEPRECATED

  /*
   * Search the mesh and cache data for the elements
   * present in the mesh.  This is done in prepare_for_use(), but can
   * be done manually by other classes after major mesh modifications.
   * Data cached includes:
   *   - elem dimensions
   *   - elem subdomains
   */
  void cache_elem_data();

  /**
   * libMesh often expects all processors to know about names of all
   * subdomain ids, but distributed mesh generators may only know
   * about part of a mesh when creating names.  This method can
   * synchronize the subdomain id to name map across processors,
   * assuming no conflicts exist.
   */
  void sync_subdomain_name_map();

  /**
   * Search the mesh for elements that have a neighboring element
   * of dim+1 and set that element as the interior parent
   */
  void detect_interior_parents();

  /**
   * \return A mesh that may own interior parents of elements in this
   * mesh.  In most cases this mesh includes its own interior parents,
   * but in cases where a separate "interior" mesh was used to create
   * this mesh as a distinct lower-dimensional boundary (or boundary
   * subset) mesh, the original mesh will be returned here.
   */
  const MeshBase & interior_mesh() const { return *_interior_mesh; }

  /**
   * \return A writeable reference to the interior mesh.
   */
  MeshBase & interior_mesh() { return *_interior_mesh; }

  /**
   * Sets the interior mesh.  For advanced use only.
   */
  void set_interior_mesh(MeshBase & int_mesh) { _interior_mesh = &int_mesh; }

  /**
   * \return The cached mesh subdomains. As long as the mesh is prepared, this
   * should contain all the subdomain ids across processors. Relies on the mesh
   * being prepared
   */
  const std::set<subdomain_id_type> & get_mesh_subdomains() const
  { libmesh_assert(this->is_prepared()); return _mesh_subdomains; }

#ifdef LIBMESH_ENABLE_PERIODIC
  /**
   * Register a pair of boundaries as disjoint neighbor boundary pairs.
   */
  void add_disjoint_neighbor_boundary_pairs(const boundary_id_type b1,
                                   const boundary_id_type b2,
                                   const RealVectorValue & translation);

  PeriodicBoundaries * get_disjoint_neighbor_boundary_pairs();

  const PeriodicBoundaries * get_disjoint_neighbor_boundary_pairs() const;

  void remove_disjoint_boundary_pair(const boundary_id_type b1,
                                     const boundary_id_type b2);
#endif

  /**
   * Flags indicating in what ways a mesh has been prepared for use.
   */
  struct Preparation {
    bool is_partitioned = false,
         has_synched_id_counts = false,
         has_neighbor_ptrs = false,
         has_cached_elem_data = false,
         has_interior_parent_ptrs = false,
         has_removed_remote_elements = false,
         has_removed_orphaned_nodes = false,
         has_boundary_id_sets = false,
         has_reinit_ghosting_functors = false;

    operator bool() const {
      return is_partitioned &&
             has_synched_id_counts &&
             has_neighbor_ptrs &&
             has_cached_elem_data &&
             has_interior_parent_ptrs &&
             has_removed_remote_elements &&
             has_removed_orphaned_nodes &&
             has_reinit_ghosting_functors &&
             has_boundary_id_sets;
    }

    Preparation & operator= (bool set_all) {
      is_partitioned = set_all;
      has_synched_id_counts = set_all;
      has_neighbor_ptrs = set_all;
      has_cached_elem_data = set_all;
      has_interior_parent_ptrs = set_all;
      has_removed_remote_elements = set_all;
      has_removed_orphaned_nodes = set_all;
      has_reinit_ghosting_functors = set_all;
      has_boundary_id_sets = set_all;

      return *this;
    }

    bool operator== (const Preparation & other) {
      return is_partitioned == other.is_partitioned &&
             has_synched_id_counts == other.has_synched_id_counts &&
             has_neighbor_ptrs == other.has_neighbor_ptrs &&
             has_cached_elem_data == other.has_cached_elem_data &&
             has_interior_parent_ptrs == other.has_interior_parent_ptrs &&
             has_removed_remote_elements == other.has_removed_remote_elements &&
             has_removed_orphaned_nodes == other.has_removed_orphaned_nodes &&
             has_reinit_ghosting_functors == other.has_reinit_ghosting_functors &&
             has_boundary_id_sets == other.has_boundary_id_sets;
    }
  };

protected:

#ifdef LIBMESH_ENABLE_PERIODIC
  /// @brief The disjoint neighbor boundary id pairs.
  std::unique_ptr<PeriodicBoundaries> _disjoint_neighbor_boundary_pairs;
#endif

  /**
   * This class holds the boundary information.  It can store nodes, edges,
   * and faces with a corresponding id that facilitates setting boundary
   * conditions.
   *
   * Direct access to this class is now officially deprecated and will
   * be removed in future libMesh versions.  Use the \p get_boundary_info()
   * accessor instead.
   */
  std::unique_ptr<BoundaryInfo> boundary_info;

  /**
   * Moves any superclass data (e.g. GhostingFunctors that might rely
   * on element and nodal data (which is managed by subclasses!)
   * being already moved first.
   *
   * Must be manually called in dofobject-managing subclass move
   * operators.
   */
  void post_dofobject_moves(MeshBase && other_mesh);

  /**
   * Helper class to copy cached data, to synchronize with a possibly
   * unprepared \p other_mesh
   */
  void copy_cached_data (const MeshBase & other_mesh);

  /**
   * Shim to allow operator == (&) to behave like a virtual function
   * without having to be one.
   */
  virtual bool subclass_locally_equals (const MeshBase & other_mesh) const = 0;

  /**
   * Tests for equality of all elements and nodes in the mesh.  Helper
   * function for subclass_equals() in unstructured mesh subclasses.
   */
  bool nodes_and_elements_equal(const MeshBase & other_mesh) const;

  /**
   * \returns A writable reference to the number of partitions.
   */
  unsigned int & set_n_partitions ()
  { return _n_parts; }

  /**
   * The number of partitions the mesh has.  This is set by
   * the partitioners, and may not be changed directly by
   * the user.
   *
   * \note The number of partitions \e need \e not equal
   * this->n_processors(), consider for example the case where you
   * simply want to partition a mesh on one processor and view the
   * result in GMV.
   */
  unsigned int _n_parts;

  /**
   * The default mapping type (typically Lagrange) between master and
   * physical space to assign to newly added elements.
   */
  ElemMappingType _default_mapping_type;

  /**
   * The default mapping data (unused with Lagrange, used for nodal
   * weight lookup index with rational bases) to assign to newly added
   * elements.
   */
  unsigned char _default_mapping_data;

  /**
   * Flags indicating in what ways \p this mesh has been prepared.
   */
  Preparation _preparation;

  /**
   * A \p PointLocator class for this mesh.
   * This will not actually be built unless needed. Further, since we want
   * our \p point_locator() method to be \p const (yet do the dynamic allocating)
   * this needs to be mutable.  Since the PointLocatorBase::build() member is used,
   * and it operates on a constant reference to the mesh, this is OK.
   */
  mutable std::unique_ptr<PointLocatorBase> _point_locator;

  /**
   * Do we count lower dimensional elements in point locator refinement?
   * This is relevant in tree-based point locators, for example.
   */
  bool _count_lower_dim_elems_in_point_locator;

  /**
   * A partitioner to use at each prepare_for_use().
   *
   * This will be built in the constructor of each derived class, but
   * can be replaced by the user through the partitioner() accessor.
   */
  std::unique_ptr<Partitioner> _partitioner;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * The next available unique id for assigning ids to DOF objects
   */
  unique_id_type _next_unique_id;
#endif

  /**
   * Defaulting to \p this, a pointer to the mesh used to generate
   * boundary elements on \p this.
   */
  MeshBase *_interior_mesh;

  /**
   * If this is true then no partitioning should be done with the
   * possible exception of orphaned nodes.
   */
  bool _skip_noncritical_partitioning;

  /**
   * If this is true then no partitioning should be done.
   */
  bool _skip_all_partitioning;

  /**
   * If this is true then renumbering will be kept to a minimum.
   *
   * This is set when prepare_for_use() is called.
   */
  bool _skip_renumber_nodes_and_elements;

  /**
   * If this is \p true then we will skip \p find_neighbors in \p prepare_for_use
   */
  bool _skip_find_neighbors;

  /**
   * If this is false then even on DistributedMesh remote elements
   * will not be deleted during mesh preparation.
   *
   * This is true by default.
   */
  bool _allow_remote_element_removal;

  /**
   * The Exodus reader (and potentially other readers in the future?)
   * now supports setting Node and Elem unique_ids based on values
   * from within the Exodus file itself, rather than generating them
   * automatically in LibMesh. In this case, the unique_ids will not
   * necessarily be unique across the set of all _DofObjects_,
   * although they should still be unique within the individual sets
   * of Elems and Nodes. The reader can therefore set this Mesh flag
   * (which defaults to false) to indicate we should be less strict
   * when checking the "uniqueness" of unique_ids.
   */
  bool _allow_node_and_elem_unique_id_overlap;

  /**
   * This structure maintains the mapping of named blocks
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<subdomain_id_type, std::string> _block_id_to_name;

  /**
   * We cache the dimension of the elements present in the mesh.
   * So, if we have a mesh with 1D and 2D elements, this structure
   * will contain 1 and 2.
   */
  std::set<unsigned char> _elem_dims;

  /**
   * We cache the (default) order of the geometric elements present in
   * the mesh.  E.g. if we have a mesh with TRI3 and TRI6 elements,
   * this structure will contain FIRST and SECOND.
   */
  std::set<Order> _elem_default_orders;

  /**
   * We cache the maximum nodal order supported by all the mesh's
   * elements (the minimum supported_nodal_order() of any element)
   */
  Order _supported_nodal_order;

  /**
   * We cache the subdomain ids of the elements present in the mesh.
   */
  std::set<subdomain_id_type> _mesh_subdomains;

  /**
   * Map from "element set code" to list of set ids to which that element
   * belongs (and vice-versa). Remarks:
   * 1.) The elemset code is a dof_id_type because (if used) it is
   * stored as an extra_integer (named "elemset_code") on all elements,
   * and extra_integers are of type dof_id_type. Elements which do not
   * belong to any set should be assigned an elemset code of DofObject::invalid_id.
   * 2.) Element sets can be thought of as a generalization of the concept
   * of a subdomain. Subdomains have the following restrictions:
   *   a.) A given element can only belong to a single subdomain
   *   b.) When using Exodus file input/output, subdomains are (unfortunately)
   *   tied to the concept of exodus element blocks, which consist of a single
   *   geometric element type, somewhat limiting their generality.
   * 3.) The user is responsible for filling in the values of this map
   * in a consistent manner, unless the elemsets are read in from an
   * Exodus file, in which case the elemset codes will be set up
   * automatically. The codes can basically be chosen arbitrarily,
   * with the one requirement that elements which belong to no sets
   * should have a set code of DofObject::invalid_id.
   * 4.) We also keep a list of all the elemset ids which have been added in
   * order to support O(1) performance behavior in n_elemsets() calls.
   */
  std::map<dof_id_type, const MeshBase::elemset_type *> _elemset_codes;
  std::map<MeshBase::elemset_type, dof_id_type> _elemset_codes_inverse_map;
  MeshBase::elemset_type _all_elemset_ids;

  /**
   * The "spatial dimension" of the Mesh.  See the documentation for
   * Mesh::spatial_dimension() for more information.
   */
  unsigned char _spatial_dimension;

  /**
   * The array of names for integer data associated with each element
   * in the mesh
   */
  std::vector<std::string> _elem_integer_names;

  /**
   * The array of default initialization values for integer data
   * associated with each element in the mesh
   */
  std::vector<dof_id_type> _elem_integer_default_values;

  /**
   * The array of names for integer data associated with each node
   * in the mesh
   */
  std::vector<std::string> _node_integer_names;

  /**
   * The array of default initialization values for integer data
   * associated with each node in the mesh
   */
  std::vector<dof_id_type> _node_integer_default_values;

  /**
   * Size extra-integer arrays of all elements in the mesh
   */
  void size_elem_extra_integers();

  /**
   * Size extra-integer arrays of all nodes in the mesh
   */
  void size_node_extra_integers();

  /**
   * Merge extra-integer arrays from an \p other mesh.  Returns two
   * mappings from index values in \p other to (possibly newly created)
   * index values with the same string name in \p this mesh, the first
   * for element integers and the second for node integers.
   */
  std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
    merge_extra_integer_names(const MeshBase & other);

  /**
   * The default geometric GhostingFunctor, used to implement standard
   * libMesh element ghosting behavior.  We use a base class pointer
   * here to avoid dragging in more header dependencies.
   */
  std::unique_ptr<GhostingFunctor> _default_ghosting;

  /**
   * The list of all GhostingFunctor objects to be used when
   * distributing a DistributedMesh.
   *
   * Basically unused by ReplicatedMesh for now, but belongs to
   * MeshBase because the cost is trivial.
   */
  std::vector<GhostingFunctor *> _ghosting_functors;

  /**
   * Hang on to references to any GhostingFunctor objects we were
   * passed in shared_ptr form
   */
  std::map<GhostingFunctor *, std::shared_ptr<GhostingFunctor> > _shared_functors;

  // Keep track of any constraint equations that are inherent to the
  // mesh, such as FE nodes whose Rational Bernstein values need to be
  // constrained in terms of values on spline control nodes.
  //
  // _constraint_rows[constrained_node][i].first.first is an
  // element (e.g. a NodeElem for a spline control node),
  // _constraint_rows[constrained_node][i].first.second is the
  // local node id of that element which is a constraining node,
  // _constraint_rows[constrained_node][i].second is that node's
  // constraint coefficient.
  constraint_rows_type _constraint_rows;

  /**
   * If nonzero, we will call PointLocatorBase::set_close_to_point_tol()
   * on any PointLocators that we create.
   */
  Real _point_locator_close_to_point_tol;

  /**
   * The partitioner class is a friend so that it can set
   * the number of partitions.
   */
  friend class Partitioner;

  /**
   * The MeshInput classes are friends so that they can set the number
   * of partitions.
   */
  friend class MeshInput<MeshBase>;

  /**
   * Make the \p BoundaryInfo class a friend so that
   * it can create and interact with \p BoundaryMesh.
   */
  friend class BoundaryInfo;

  /**
   * Make the \p MeshCommunication class a friend so that
   * it can directly broadcast *_integer_names
   */
  friend class MeshCommunication;


  /**
   * The original iterator classes weren't properly const-safe;
   * relying on their const-incorrectness is now deprecated.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
typedef variant_filter_iterator<MeshBase::Predicate, Elem *> elem_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate,
                                Elem * const,
                                Elem * const &,
                                Elem * const *> const_elem_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate, Node *> node_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate,
                                Node * const,
                                Node * const &,
                                Node * const *> const_node_filter_iter;
#else
typedef variant_filter_iterator<MeshBase::Predicate,
                                Elem * const,
                                Elem * const &,
                                Elem * const *,
                                const Elem * const,
                                const Elem * const &,
                                const Elem * const *> elem_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate,
                                const Elem * const,
                                const Elem * const &,
                                const Elem * const *> const_elem_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate,
                                Node * const,
                                Node * const &,
                                Node * const *,
                                const Node * const,
                                const Node * const &,
                                const Node * const *> node_filter_iter;

typedef variant_filter_iterator<MeshBase::Predicate,
                                const Node * const,
                                const Node * const &,
                                const Node * const *> const_node_filter_iter;
#endif // LIBMESH_ENABLE_DEPRECATED

};











/**
 * The definition of the element_iterator struct.
 */
struct
MeshBase::element_iterator : MeshBase::elem_filter_iter
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  element_iterator (const IterType & d,
                    const IterType & e,
                    const PredType & p ) :
    elem_filter_iter(d,e,p) {}
};




/**
 * The definition of the const_element_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_element_iterator : MeshBase::const_elem_filter_iter
{
  /**
   * Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor.
   */
  template <typename PredType, typename IterType>
  const_element_iterator (const IterType & d,
                          const IterType & e,
                          const PredType & p ) :
    const_elem_filter_iter(d,e,p)  {}

  /**
   * The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
   * variant_filter_iterator copy constructor.
   *
   * \note This one is \e not templated!
   */
  const_element_iterator (const MeshBase::element_iterator & rhs) :
    const_elem_filter_iter(rhs) {}
};







/**
 * The definition of the node_iterator struct.
 */
struct
MeshBase::node_iterator : MeshBase::node_filter_iter
{
  /**
   * Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor.
   */
  template <typename PredType, typename IterType>
  node_iterator (const IterType & d,
                 const IterType & e,
                 const PredType & p ) :
    node_filter_iter(d,e,p)  {}
};




/**
 * The definition of the const_node_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_node_iterator : MeshBase::const_node_filter_iter
{
  /**
   * Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor.
   */
  template <typename PredType, typename IterType>
  const_node_iterator (const IterType & d,
                       const IterType & e,
                       const PredType & p ) :
    const_node_filter_iter(d,e,p)  {}

  /**
   * The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
   * variant_filter_iterator copy constructor.
   *
   * \note This one is *not* templated!
   */
  const_node_iterator (const MeshBase::node_iterator & rhs) :
    const_node_filter_iter(rhs) {}
};


template <typename T>
inline
unsigned int MeshBase::add_elem_datum(const std::string & name,
                                      bool allocate_data,
                                      const T * default_value)
{
  const std::size_t old_size = _elem_integer_names.size();

  unsigned int n_more_integers = (sizeof(T)-1)/sizeof(dof_id_type);
  std::vector<dof_id_type> int_data(n_more_integers+1, DofObject::invalid_id);
  if (default_value)
    std::memcpy(int_data.data(), default_value, sizeof(T));

  unsigned int start_idx = this->add_elem_integer(name, false, int_data[0]);
  for (unsigned int i=0; i != n_more_integers; ++i)
    this->add_elem_integer(name+"__"+std::to_string(i), false, int_data[i+1]);

  if (allocate_data && old_size != _elem_integer_names.size())
    this->size_elem_extra_integers();

  return start_idx;
}


template <typename T>
inline
std::vector<unsigned int> MeshBase::add_elem_data(const std::vector<std::string> & names,
                                                  bool allocate_data,
                                                  const std::vector<T> * default_values)
{
  libmesh_assert(!default_values || default_values->size() == names.size());

  std::vector<unsigned int> returnval(names.size());

  const std::size_t old_size = _elem_integer_names.size();

  for (auto i : index_range(names))
    returnval[i] =
      this->add_elem_datum<T>(names[i], false,
                              default_values ?
                              (*default_values)[i] : nullptr);

  if (allocate_data && old_size != _elem_integer_names.size())
    this->size_elem_extra_integers();

  return returnval;
}


template <typename T>
inline
unsigned int MeshBase::add_node_datum(const std::string & name,
                                      bool allocate_data,
                                      const T * default_value)
{
  const std::size_t old_size = _node_integer_names.size();

  unsigned int n_more_integers = (sizeof(T)-1)/sizeof(dof_id_type);
  std::vector<dof_id_type> int_data(n_more_integers+1, DofObject::invalid_id);
  if (default_value)
    std::memcpy(int_data.data(), default_value, sizeof(T));

  unsigned int start_idx = this->add_node_integer(name, false, int_data[0]);
  for (unsigned int i=0; i != n_more_integers; ++i)
    this->add_node_integer(name+"__"+std::to_string(i), false, int_data[i+1]);

  if (allocate_data && old_size != _node_integer_names.size())
    this->size_node_extra_integers();

  return start_idx;
}


template <typename T>
inline
std::vector<unsigned int> MeshBase::add_node_data(const std::vector<std::string> & names,
                                                  bool allocate_data,
                                                  const std::vector<T> * default_values)
{
  libmesh_assert(!default_values || default_values->size() == names.size());

  std::vector<unsigned int> returnval(names.size());

  const std::size_t old_size = _node_integer_names.size();

  for (auto i : index_range(names))
    returnval[i] =
      this->add_node_datum<T>(names[i], false,
                              default_values ?
                              (*default_values)[i] : nullptr);

  if (allocate_data && old_size != _node_integer_names.size())
    this->size_node_extra_integers();

  return returnval;
}



} // namespace libMesh

#endif // LIBMESH_MESH_BASE_H
