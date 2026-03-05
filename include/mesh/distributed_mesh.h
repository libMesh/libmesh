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



#ifndef LIBMESH_DISTRIBUTED_MESH_H
#define LIBMESH_DISTRIBUTED_MESH_H

// Local Includes
#include "libmesh_config.h"

#if LIBMESH_MAPVECTOR_CHUNK_SIZE == 1
#  include "libmesh/mapvector.h"
#else
#  include "libmesh/chunked_mapvector.h"
#endif

#include "libmesh/unstructured_mesh.h"

// C++ Includes
#include <cstddef>
#include <memory>
#include <set>

namespace libMesh
{

// Forward declarations
class Elem;
class Node;


/**
 * The \p DistributedMesh class is derived from the \p MeshBase class,
 * and is intended to provide identical functionality to the user
 * but be distributed rather than replicated across distributed-memory
 * systems.
 *
 * \author Roy Stogner
 * \date 2007
 * \brief Mesh data structure which is distributed across all processors.
 */
class DistributedMesh : public UnstructuredMesh
{
public:

  template <typename Obj>
#if LIBMESH_MAPVECTOR_CHUNK_SIZE == 1
  using dofobject_container = mapvector<Obj *, dof_id_type>;
#else
  using dofobject_container = chunked_mapvector<Obj *, dof_id_type, LIBMESH_MAPVECTOR_CHUNK_SIZE>;
#endif

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  DistributedMesh (const Parallel::Communicator & comm_in,
                   unsigned char dim=1);

  /**
   * Copy-constructor.  This should be able to take a
   * replicated or distributed mesh.
   */
  DistributedMesh (const MeshBase & other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * distributed mesh.
   */
  DistributedMesh (const DistributedMesh & other_mesh);

  /**
   * Move-constructor deleted in MeshBase.
   */
  DistributedMesh(DistributedMesh &&) = delete;

  /**
   * Copy assignment is not allowed.
   */
  DistributedMesh & operator= (const DistributedMesh &) = delete;

  /**
   * Overloaded operator= will move contents of other_mesh to calling
   * DistributedMesh object.
   */
  DistributedMesh & operator= (DistributedMesh && other_mesh);

  /**
   * Shim to call the move assignment operator for this class
  */
  virtual MeshBase & assign(MeshBase && other_mesh) override;

  /**
   * Shim to allow operator == (&) to behave like a virtual function
   * without having to be one.
   */
  virtual std::string_view subclass_first_difference_from (const MeshBase & other_mesh_base) const override;

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual std::unique_ptr<MeshBase> clone () const override
  {
    auto returnval = std::make_unique<DistributedMesh>(*this);
#ifdef DEBUG
    libmesh_assert(*returnval == *this);
#endif
    return returnval;
  }

  /**
   * Destructor.
   */
  virtual ~DistributedMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear() override;

  /**
   * Clear internal Elem data.
   */
  virtual void clear_elems() override;

  /**
   * Redistribute elements between processors.  This gets called
   * automatically by the Partitioner, and is a no-op in the case of a
   * serialized mesh.
   */
  virtual void redistribute () override;

  /**
   * Recalculate cached data after elements and nodes have been
   * repartitioned.
   */
  virtual void update_post_partitioning () override;

  /**
   * \returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const override final
  { return _is_serial; }

  /**
   * \returns \p true if all elements and nodes of the mesh
   * exist on the processor 0, \p false otherwise
   */
  virtual bool is_serial_on_zero () const override final
  { return _is_serial || _is_serial_on_proc_0; }

  /**
   * Asserts that not all elements and nodes of the mesh necessarily
   * exist on the current processor.
   */
  virtual void set_distributed () override final
  { _is_serial = false;
    _is_serial_on_proc_0 = false; }

  /**
   * \returns \p true if new elements and nodes can and should be
   * created in synchronization on all processors, \p false otherwise
   */
  virtual bool is_replicated () const override final
  { return false; }

  /**
   * Verify id, processor_id, and if applicable unique_id consistency
   * of a parallel objects container.
   * Calls libmesh_assert() on each possible failure in that container.
   */
  template <typename T>
  void libmesh_assert_valid_parallel_object_ids(const dofobject_container<T> &) const;

  /**
   * Verify id and processor_id consistency of our elements and
   * nodes containers.
   * Calls libmesh_assert() on each possible failure.
   */
  virtual void libmesh_assert_valid_parallel_ids() const override;

  /**
   * Verify p_level consistency of our elements containers.
   * Calls libmesh_assert() on each possible failure.
   */
  void libmesh_assert_valid_parallel_p_levels() const;

  /**
   * Verify refinement_flag and p_refinement_flag consistency of our
   * elements containers.
   * Calls libmesh_assert() on each possible failure.
   */
  void libmesh_assert_valid_parallel_flags() const;

  /**
   * Renumber a parallel objects container.
   *
   * \returns The smallest globally unused id for that container.
   */
  template <typename T>
  dof_id_type renumber_dof_objects (dofobject_container<T> &);

  /**
   * Remove nullptr elements from arrays.
   */
  virtual void renumber_nodes_and_elements () override;

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor
   */
  virtual void allgather() override;

  /**
   * Gathers all elements and nodes of the mesh onto
   * processor zero
   */
  virtual void gather_to_zero() override;

  /**
   * Deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   */
  virtual void delete_remote_elements() override;

  /**
   * Inserts the element _and_ adds it to a list of elements that
   * should not get deleted or have their descendants deleted by
   * delete_remote_elements.  This is handy for inserting otherwise
   * off-processor elements that you want to keep track of on this
   * processor.
   */
  virtual void add_extra_ghost_elem(Elem * e);

  /**
   * Clears extra ghost elements.
   */
  virtual void clear_extra_ghost_elems() { _extra_ghost_elems.clear(); }

  /**
   * Clears specified extra ghost elements
   */
  virtual void clear_extra_ghost_elems(const std::set<Elem *> & extra_ghost_elems);

  /**
   * Const accessor to the ghosted elements
   */
  const std::set<Elem *> & extra_ghost_elems() const { return _extra_ghost_elems; }

  // Cached methods that can be called in serial
  virtual dof_id_type n_nodes () const override final { return _n_nodes; }
  virtual dof_id_type max_node_id () const override final { return _max_node_id; }
  virtual void reserve_nodes (const dof_id_type) override final {}
  virtual dof_id_type n_elem () const override final { return _n_elem; }
  virtual dof_id_type n_active_elem () const override final;
  virtual dof_id_type max_elem_id () const override final { return _max_elem_id; }
  virtual void reserve_elem (const dof_id_type) override final {}

  // Parallel only method to update the caches
  virtual void update_parallel_id_counts () override;

  // And more parallel only methods to test non-cached values
  virtual dof_id_type parallel_n_nodes () const override;
  dof_id_type parallel_max_node_id () const;
  virtual dof_id_type parallel_n_elem () const override;
  dof_id_type parallel_max_elem_id () const;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const override;
  virtual void set_next_unique_id(unique_id_type id) override;
#endif

  virtual const Point & point (const dof_id_type i) const override final;

  virtual const Node * node_ptr (const dof_id_type i) const override final;
  virtual Node * node_ptr (const dof_id_type i) override final;

  virtual const Node * query_node_ptr (const dof_id_type i) const override final;
  virtual Node * query_node_ptr (const dof_id_type i) override final;

  virtual const Elem * elem_ptr (const dof_id_type i) const override final;
  virtual Elem * elem_ptr (const dof_id_type i) override final;

  virtual const Elem * query_elem_ptr (const dof_id_type i) const override final;
  virtual Elem * query_elem_ptr (const dof_id_type i) override final;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node * add_point (const Point & p,
                            const dof_id_type id = DofObject::invalid_id,
                            const processor_id_type proc_id = DofObject::invalid_processor_id) override final;
  virtual Node * add_node (Node * n) override final;
  virtual Node * add_node (std::unique_ptr<Node> n) override final;

  /**
   * Takes ownership of node \p n on this partition of a distributed
   * mesh, by setting n.processor_id() to this->processor_id(), as
   * well as changing n.id() and moving it in the mesh's internal
   * container to give it a new authoritative id.
   */
  virtual void own_node (Node & n) override final;

  virtual void delete_node (Node * n) override final;
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id) override final;
  virtual Elem * add_elem (Elem * e) override final;
  virtual Elem * add_elem (std::unique_ptr<Elem> e) override final;
  virtual Elem * insert_elem (Elem * e) override final;
  virtual Elem * insert_elem (std::unique_ptr<Elem> e) override final;
  virtual void delete_elem (Elem * e) override final;
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id) override final;

  /**
   * There is no reason for a user to ever call this function.
   *
   * This function restores a previously broken element/node numbering such that
   * \p mesh.node_ref(n).id() == n.
   */
  virtual void fix_broken_node_and_element_numbering () override;

public:
  /**
   * Elem and Node iterator accessor functions.  See MeshBase for
   * documentation.
   */
  DECLARE_ELEM_ITERATORS(,,);
  DECLARE_ELEM_ITERATORS(active_,,);
  DECLARE_ELEM_ITERATORS(ancestor_,,)
  DECLARE_ELEM_ITERATORS(subactive_,,)
  DECLARE_ELEM_ITERATORS(local_,,)
  DECLARE_ELEM_ITERATORS(unpartitioned_,,)
  DECLARE_ELEM_ITERATORS(facelocal_,,)
  DECLARE_ELEM_ITERATORS(level_, unsigned int level, level)
  DECLARE_ELEM_ITERATORS(pid_, processor_id_type pid, pid)
  DECLARE_ELEM_ITERATORS(type_, ElemType type, type)

  DECLARE_ELEM_ITERATORS(active_subdomain_, subdomain_id_type sid, sid)
  DECLARE_ELEM_ITERATORS(active_subdomain_set_, std::set<subdomain_id_type> ss, ss)

  // Backwards compatibility
  virtual SimpleRange<element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) override final { return active_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<const_element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) const override final { return active_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) override final { return active_local_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<const_element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) const override final { return active_local_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) override final { return active_subdomain_set_element_ptr_range(ss); }
  virtual SimpleRange<const_element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const override final { return active_subdomain_set_element_ptr_range(ss); }

  DECLARE_ELEM_ITERATORS(not_active_,,);
  DECLARE_ELEM_ITERATORS(not_ancestor_,,);
  DECLARE_ELEM_ITERATORS(not_subactive_,,);
  DECLARE_ELEM_ITERATORS(not_local_,,);
  DECLARE_ELEM_ITERATORS(not_level_, unsigned int level, level)

  DECLARE_ELEM_ITERATORS(active_local_,,)
  DECLARE_ELEM_ITERATORS(active_not_local_,,)
  DECLARE_ELEM_ITERATORS(active_unpartitioned_,,)
  DECLARE_ELEM_ITERATORS(active_type_, ElemType type, type)
  DECLARE_ELEM_ITERATORS(active_pid_, processor_id_type pid, pid)
  DECLARE_ELEM_ITERATORS(local_level_, unsigned int level, level)
  DECLARE_ELEM_ITERATORS(local_not_level_, unsigned int level, level)
  DECLARE_ELEM_ITERATORS(active_local_subdomain_, subdomain_id_type sid, sid)
  DECLARE_ELEM_ITERATORS(active_local_subdomain_set_, std::set<subdomain_id_type> ss, ss)

  DECLARE_ELEM_ITERATORS(semilocal_,,)
  DECLARE_ELEM_ITERATORS(ghost_,,)
  DECLARE_ELEM_ITERATORS(active_semilocal_,,)

  DECLARE_ELEM_ITERATORS(evaluable_, const DofMap & dof_map LIBMESH_COMMA unsigned int var_num = libMesh::invalid_uint, dof_map LIBMESH_COMMA var_num)
  DECLARE_ELEM_ITERATORS(multi_evaluable_, std::vector<const DofMap *> dof_maps, dof_maps)

#ifdef LIBMESH_ENABLE_AMR
  DECLARE_ELEM_ITERATORS(flagged_, unsigned char rflag, rflag)

  // Elem::refinement_flag() == rflag && Elem::processor_id() == pid
  DECLARE_ELEM_ITERATORS(flagged_pid_, unsigned char rflag LIBMESH_COMMA processor_id_type pid, rflag LIBMESH_COMMA pid)
#endif

  DECLARE_NODE_ITERATORS(,,)
  DECLARE_NODE_ITERATORS(active_,,)
  DECLARE_NODE_ITERATORS(local_,,)
  DECLARE_NODE_ITERATORS(bnd_,,)
  DECLARE_NODE_ITERATORS(pid_, processor_id_type pid, pid)
  DECLARE_NODE_ITERATORS(bid_, boundary_id_type bid, bid)

  DECLARE_NODE_ITERATORS(evaluable_, const DofMap & dof_map LIBMESH_COMMA unsigned int var_num = libMesh::invalid_uint, dof_map LIBMESH_COMMA var_num)

  DECLARE_NODE_ITERATORS(multi_evaluable_, std::vector<const DofMap *> dof_maps, dof_maps)


protected:

  /**
   * Move node and elements from a DistributedMesh.
   */
  virtual void move_nodes_and_elements(MeshBase && other_mesh) override;

  /**
   * The vertices (spatial coordinates) of the mesh.
   */
  dofobject_container<Node> _nodes;

  /**
   * The elements in the mesh.
   */
  dofobject_container<Elem> _elements;

  /**
   * A boolean remembering whether we're serialized or not
   */
  bool _is_serial;

  /**
   * A boolean remembering whether we're serialized to proc 0 or not
   */
  bool _is_serial_on_proc_0;

  /**
   * A boolean remembering whether we've recently deleted top-level elements or
   * not.  If so, we'll need to be extra careful when deleting "unused" nodes
   * they used to have, because those nodes might still be used by ghost elements.
   */
  bool _deleted_coarse_elements;

  /**
   * Cached data from the last renumber_nodes_and_elements call
   */
  dof_id_type _n_nodes, _n_elem, _max_node_id, _max_elem_id;

  /**
   * Guaranteed globally unused IDs for use when adding new
   * nodes or elements.
   */
  dof_id_type _next_free_local_node_id,
    _next_free_local_elem_id;
  dof_id_type _next_free_unpartitioned_node_id,
    _next_free_unpartitioned_elem_id;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * The next available unique id for assigning ids to unpartitioned DOF objects
   */
  unique_id_type _next_unpartitioned_unique_id;
#endif

  /**
   * These are extra ghost elements that we want to make sure
   * not to delete when we call delete_remote_elements()
   */
  std::set<Elem *> _extra_ghost_elems;

private:

  /**
   * Typedefs for the container implementation.
   */
  typedef dofobject_container<Elem>::veclike_iterator             elem_iterator_imp;
  typedef dofobject_container<Elem>::const_veclike_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.
   */
  typedef dofobject_container<Node>::veclike_iterator             node_iterator_imp;
  typedef dofobject_container<Node>::const_veclike_iterator const_node_iterator_imp;
};


} // namespace libMesh

#endif // LIBMESH_DISTRIBUTED_MESH_H
