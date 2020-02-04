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



#ifndef LIBMESH_DISTRIBUTED_MESH_H
#define LIBMESH_DISTRIBUTED_MESH_H

// Local Includes
#include "libmesh/mapvector.h"
#include "libmesh/unstructured_mesh.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ Includes
#include <cstddef>
#include <set>

namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
template <typename> class NodeTempl;


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
template <typename RealType = Real>
class DistributedMeshTempl : public UnstructuredMeshTempl<RealType>
{
public:
  typedef DistributedMeshTempl<RealType> DistributedMesh;
  typedef UnstructuredMeshTempl<RealType> UnstructuredMesh;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef NodeTempl<RealType> Node;
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;

  using typename MeshBase::element_iterator;
  using typename MeshBase::const_element_iterator;
  using typename MeshBase::node_iterator;
  using typename MeshBase::const_node_iterator;

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  DistributedMeshTempl (const Parallel::Communicator & comm_in,
                   unsigned char dim=1);

  /**
   * Copy-constructor.  This should be able to take a
   * replicated or distributed mesh.
   */
  DistributedMeshTempl (const UnstructuredMesh & other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * distributed mesh.
   */
  DistributedMeshTempl (const DistributedMesh & other_mesh);

  /**
   * Move-constructor.
   */
  DistributedMeshTempl(DistributedMesh &&) = default;

  /**
   * Copy and move assignment are not allowed.
   */
  DistributedMesh & operator= (const DistributedMesh &) = delete;
  DistributedMesh & operator= (DistributedMesh &&) = delete;

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual std::unique_ptr<MeshBase> clone () const override
  { return libmesh_make_unique<DistributedMesh>(*this); }

  /**
   * Destructor.
   */
  virtual ~DistributedMeshTempl();

  /**
   * Clear all internal data.
   */
  virtual void clear() override;

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
  virtual bool is_serial () const override
  { return _is_serial; }

  /**
   * \returns \p true if all elements and nodes of the mesh
   * exist on the processor 0, \p false otherwise
   */
  virtual bool is_serial_on_zero () const override
  { return _is_serial || _is_serial_on_proc_0; }

  /**
   * Asserts that not all elements and nodes of the mesh necessarily
   * exist on the current processor.
   */
  virtual void set_distributed () override
  { _is_serial = false;
    _is_serial_on_proc_0 = false; }

  /**
   * \returns \p true if new elements and nodes can and should be
   * created in synchronization on all processors, \p false otherwise
   */
  virtual bool is_replicated () const override
  { return false; }

  /**
   * Verify id, processor_id, and if applicable unique_id consistency
   * of a parallel objects container.
   * Calls libmesh_assert() on each possible failure in that container.
   */
  template <typename T>
  void libmesh_assert_valid_parallel_object_ids(const mapvector<T *,dof_id_type> &) const;

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
  dof_id_type renumber_dof_objects (mapvector<T *,dof_id_type> &);

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
  virtual dof_id_type n_nodes () const override { return _n_nodes; }
  virtual dof_id_type max_node_id () const override { return _max_node_id; }
  virtual void reserve_nodes (const dof_id_type) override {}
  virtual dof_id_type n_elem () const override { return _n_elem; }
  virtual dof_id_type n_active_elem () const override;
  virtual dof_id_type max_elem_id () const override { return _max_elem_id; }
  virtual void reserve_elem (const dof_id_type) override {}

  // Parallel only method to update the caches
  virtual void update_parallel_id_counts () override;

  // And more parallel only methods to test non-cached values
  virtual dof_id_type parallel_n_nodes () const override;
  dof_id_type parallel_max_node_id () const;
  virtual dof_id_type parallel_n_elem () const override;
  dof_id_type parallel_max_elem_id () const;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const override;
#endif

  virtual const Point & point (const dof_id_type i) const override;

  virtual const Node * node_ptr (const dof_id_type i) const override;
  virtual Node * node_ptr (const dof_id_type i) override;

  virtual const Node * query_node_ptr (const dof_id_type i) const override;
  virtual Node * query_node_ptr (const dof_id_type i) override;

  virtual const Elem * elem_ptr (const dof_id_type i) const override;
  virtual Elem * elem_ptr (const dof_id_type i) override;

  virtual const Elem * query_elem_ptr (const dof_id_type i) const override;
  virtual Elem * query_elem_ptr (const dof_id_type i) override;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node * add_point (const Point & p,
                            const dof_id_type id = DofObject::invalid_id,
                            const processor_id_type proc_id = DofObject::invalid_processor_id) override;
  virtual Node * add_node (Node * n) override;

  /**
   * Calls add_node().
   */
  virtual Node * insert_node(Node * n) override;

  /**
   * Takes ownership of node \p n on this partition of a distributed
   * mesh, by setting n.processor_id() to this->processor_id(), as
   * well as changing n.id() and moving it in the mesh's internal
   * container to give it a new authoritative id.
   */
  virtual void own_node (Node & n) override;

  virtual void delete_node (Node * n) override;
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id) override;
  virtual Elem * add_elem (Elem * e) override;
  virtual Elem * insert_elem (Elem * e) override;
  virtual void delete_elem (Elem * e) override;
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id) override;

  /**
   * There is no reason for a user to ever call this function.
   *
   * This function restores a previously broken element/node numbering such that
   * \p mesh.node_ref(n).id() == n.
   */
  virtual void fix_broken_node_and_element_numbering () override;

public:
  /**
   * Elem iterator accessor functions.
   */
  virtual element_iterator elements_begin () override;
  virtual element_iterator elements_end () override;
  virtual const_element_iterator elements_begin() const override;
  virtual const_element_iterator elements_end() const override;
  virtual SimpleRange<element_iterator> element_ptr_range() override { return {elements_begin(), elements_end()}; }
  virtual SimpleRange<const_element_iterator> element_ptr_range() const override { return {elements_begin(), elements_end()}; }

  virtual element_iterator active_elements_begin () override;
  virtual element_iterator active_elements_end () override;
  virtual const_element_iterator active_elements_begin() const override;
  virtual const_element_iterator active_elements_end() const override;
  virtual SimpleRange<element_iterator> active_element_ptr_range() override { return {active_elements_begin(), active_elements_end()}; }
  virtual SimpleRange<const_element_iterator> active_element_ptr_range() const override  { return {active_elements_begin(), active_elements_end()}; }

  virtual element_iterator ancestor_elements_begin () override;
  virtual element_iterator ancestor_elements_end () override;
  virtual const_element_iterator ancestor_elements_begin() const override;
  virtual const_element_iterator ancestor_elements_end() const override;

  virtual element_iterator subactive_elements_begin () override;
  virtual element_iterator subactive_elements_end () override;
  virtual const_element_iterator subactive_elements_begin() const override;
  virtual const_element_iterator subactive_elements_end() const override;

  virtual element_iterator not_active_elements_begin () override;
  virtual element_iterator not_active_elements_end () override;
  virtual const_element_iterator not_active_elements_begin() const override;
  virtual const_element_iterator not_active_elements_end() const override;

  virtual element_iterator not_ancestor_elements_begin () override;
  virtual element_iterator not_ancestor_elements_end () override;
  virtual const_element_iterator not_ancestor_elements_begin() const override;
  virtual const_element_iterator not_ancestor_elements_end() const override;

  virtual element_iterator not_subactive_elements_begin () override;
  virtual element_iterator not_subactive_elements_end () override;
  virtual const_element_iterator not_subactive_elements_begin() const override;
  virtual const_element_iterator not_subactive_elements_end() const override;

  virtual element_iterator local_elements_begin () override;
  virtual element_iterator local_elements_end () override;
  virtual const_element_iterator local_elements_begin () const override;
  virtual const_element_iterator local_elements_end () const override;

  virtual element_iterator semilocal_elements_begin () override;
  virtual element_iterator semilocal_elements_end () override;
  virtual const_element_iterator semilocal_elements_begin () const override;
  virtual const_element_iterator semilocal_elements_end () const override;

  virtual element_iterator active_semilocal_elements_begin () override;
  virtual element_iterator active_semilocal_elements_end () override;
  virtual const_element_iterator active_semilocal_elements_begin () const override;
  virtual const_element_iterator active_semilocal_elements_end () const override;

  virtual element_iterator facelocal_elements_begin () override;
  virtual element_iterator facelocal_elements_end () override;
  virtual const_element_iterator facelocal_elements_begin () const override;
  virtual const_element_iterator facelocal_elements_end () const override;

  virtual element_iterator not_local_elements_begin () override;
  virtual element_iterator not_local_elements_end () override;
  virtual const_element_iterator not_local_elements_begin () const override;
  virtual const_element_iterator not_local_elements_end () const override;

  virtual element_iterator active_local_elements_begin () override;
  virtual element_iterator active_local_elements_end () override;
  virtual const_element_iterator active_local_elements_begin () const override;
  virtual const_element_iterator active_local_elements_end () const override;
  virtual SimpleRange<element_iterator> active_local_element_ptr_range() override { return {active_local_elements_begin(), active_local_elements_end()}; }
  virtual SimpleRange<const_element_iterator> active_local_element_ptr_range() const override  { return {active_local_elements_begin(), active_local_elements_end()}; }

  virtual element_iterator active_not_local_elements_begin () override;
  virtual element_iterator active_not_local_elements_end () override;
  virtual const_element_iterator active_not_local_elements_begin () const override;
  virtual const_element_iterator active_not_local_elements_end () const override;

  virtual element_iterator level_elements_begin (unsigned int level) override;
  virtual element_iterator level_elements_end (unsigned int level) override;
  virtual const_element_iterator level_elements_begin (unsigned int level) const override;
  virtual const_element_iterator level_elements_end (unsigned int level) const override;

  virtual element_iterator not_level_elements_begin (unsigned int level) override;
  virtual element_iterator not_level_elements_end (unsigned int level) override;
  virtual const_element_iterator not_level_elements_begin (unsigned int level) const override;
  virtual const_element_iterator not_level_elements_end (unsigned int level) const override;

  virtual element_iterator local_level_elements_begin (unsigned int level) override;
  virtual element_iterator local_level_elements_end (unsigned int level) override;
  virtual const_element_iterator local_level_elements_begin (unsigned int level) const override;
  virtual const_element_iterator local_level_elements_end (unsigned int level) const override;

  virtual element_iterator local_not_level_elements_begin (unsigned int level) override;
  virtual element_iterator local_not_level_elements_end (unsigned int level) override;
  virtual const_element_iterator local_not_level_elements_begin (unsigned int level) const override;
  virtual const_element_iterator local_not_level_elements_end (unsigned int level) const override;

  virtual element_iterator pid_elements_begin (processor_id_type proc_id) override;
  virtual element_iterator pid_elements_end (processor_id_type proc_id) override;
  virtual const_element_iterator pid_elements_begin (processor_id_type proc_id) const override;
  virtual const_element_iterator pid_elements_end (processor_id_type proc_id) const override;

  virtual element_iterator type_elements_begin (ElemType type) override;
  virtual element_iterator type_elements_end (ElemType type) override;
  virtual const_element_iterator type_elements_begin (ElemType type) const override;
  virtual const_element_iterator type_elements_end (ElemType type) const override;

  virtual element_iterator active_type_elements_begin (ElemType type) override;
  virtual element_iterator active_type_elements_end (ElemType type) override;
  virtual const_element_iterator active_type_elements_begin (ElemType type) const override;
  virtual const_element_iterator active_type_elements_end (ElemType type) const override;

  virtual element_iterator active_pid_elements_begin (processor_id_type proc_id) override;
  virtual element_iterator active_pid_elements_end (processor_id_type proc_id) override;
  virtual const_element_iterator active_pid_elements_begin (processor_id_type proc_id) const override;
  virtual const_element_iterator active_pid_elements_end (processor_id_type proc_id) const override;

  virtual element_iterator unpartitioned_elements_begin () override;
  virtual element_iterator unpartitioned_elements_end () override;
  virtual const_element_iterator unpartitioned_elements_begin () const override;
  virtual const_element_iterator unpartitioned_elements_end () const override;

  virtual element_iterator active_unpartitioned_elements_begin () override;
  virtual element_iterator active_unpartitioned_elements_end () override;
  virtual const_element_iterator active_unpartitioned_elements_begin () const override;
  virtual const_element_iterator active_unpartitioned_elements_end () const override;

  virtual element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) override;
  virtual element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) override;
  virtual const_element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) const override;
  virtual const_element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) const override;
  virtual SimpleRange<element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) override
  { return {active_local_subdomain_elements_begin(subdomain_id), active_local_subdomain_elements_end(subdomain_id)}; }
  virtual SimpleRange<const_element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) const override
  { return {active_local_subdomain_elements_begin(subdomain_id), active_local_subdomain_elements_end(subdomain_id)}; }

  virtual element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) override;
  virtual element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) override;
  virtual const_element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) const override;
  virtual const_element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) const override;
  virtual SimpleRange<element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) override
  { return {active_subdomain_elements_begin(subdomain_id), active_subdomain_elements_end(subdomain_id)}; }
  virtual SimpleRange<const_element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) const override
  { return {active_subdomain_elements_begin(subdomain_id), active_subdomain_elements_end(subdomain_id)}; }

  virtual element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) override;
  virtual element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) override;
  virtual const_element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) const override;
  virtual const_element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) const override;
  virtual SimpleRange<element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) override
  { return {active_subdomain_set_elements_begin(ss), active_subdomain_set_elements_end(ss)}; }
  virtual SimpleRange<const_element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const override
  { return {active_subdomain_set_elements_begin(ss), active_subdomain_set_elements_end(ss)}; }

  virtual element_iterator ghost_elements_begin () override;
  virtual element_iterator ghost_elements_end () override;
  virtual const_element_iterator ghost_elements_begin () const override;
  virtual const_element_iterator ghost_elements_end () const override;

  virtual element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) override;

  virtual element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) override;

  virtual const_element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) const override;

  virtual const_element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) const override;

#ifdef LIBMESH_ENABLE_AMR
  virtual element_iterator flagged_elements_begin (unsigned char rflag) override;
  virtual element_iterator flagged_elements_end (unsigned char rflag) override;
  virtual const_element_iterator flagged_elements_begin (unsigned char rflag) const override;
  virtual const_element_iterator flagged_elements_end (unsigned char rflag) const override;

  virtual element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                       processor_id_type pid) override;
  virtual element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                     processor_id_type pid) override;
  virtual const_element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                             processor_id_type pid) const override;
  virtual const_element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                           processor_id_type pid) const override;
#endif // LIBMESH_ENABLE_AMR

  /**
   * Node iterator accessor functions.
   */
  virtual node_iterator nodes_begin () override;
  virtual node_iterator nodes_end () override;
  virtual const_node_iterator nodes_begin () const override;
  virtual const_node_iterator nodes_end () const override;
  virtual SimpleRange<node_iterator> node_ptr_range() override { return {nodes_begin(), nodes_end()}; }
  virtual SimpleRange<const_node_iterator> node_ptr_range() const override { return {nodes_begin(), nodes_end()}; }

  virtual node_iterator active_nodes_begin () override;
  virtual node_iterator active_nodes_end () override;
  virtual const_node_iterator active_nodes_begin () const override;
  virtual const_node_iterator active_nodes_end () const override;

  virtual node_iterator local_nodes_begin () override;
  virtual node_iterator local_nodes_end () override;
  virtual const_node_iterator local_nodes_begin () const override;
  virtual const_node_iterator local_nodes_end () const override;
  virtual SimpleRange<node_iterator> local_node_ptr_range() override { return {local_nodes_begin(), local_nodes_end()}; }
  virtual SimpleRange<const_node_iterator> local_node_ptr_range() const override { return {local_nodes_begin(), local_nodes_end()}; }

  virtual node_iterator pid_nodes_begin (processor_id_type proc_id) override;
  virtual node_iterator pid_nodes_end (processor_id_type proc_id) override;
  virtual const_node_iterator pid_nodes_begin (processor_id_type proc_id) const override;
  virtual const_node_iterator pid_nodes_end (processor_id_type proc_id) const override;

  virtual node_iterator bid_nodes_begin (boundary_id_type bndry_id) override;
  virtual node_iterator bid_nodes_end (boundary_id_type bndry_id) override;
  virtual const_node_iterator bid_nodes_begin (boundary_id_type bndry_id) const override;
  virtual const_node_iterator bid_nodes_end (boundary_id_type bndry_id) const override;

  virtual node_iterator bnd_nodes_begin () override;
  virtual node_iterator bnd_nodes_end () override;
  virtual const_node_iterator bnd_nodes_begin () const override;
  virtual const_node_iterator bnd_nodes_end () const override;

  virtual node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) override;
  virtual node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) override;
  virtual const_node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) const override;
  virtual const_node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) const override;


protected:

  /**
   * The vertices (spatial coordinates) of the mesh.
   */
  mapvector<Node *, dof_id_type> _nodes;

  /**
   * The elements in the mesh.
   */
  mapvector<Elem *, dof_id_type> _elements;

  /**
   * A boolean remembering whether we're serialized or not
   */
  bool _is_serial;

  /**
   * A boolean remembering whether we're serialized to proc 0 or not
   */
  bool _is_serial_on_proc_0;

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
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Elem *>.
   */
  typedef typename mapvector<Elem *, dof_id_type>::veclike_iterator             elem_iterator_imp;
  typedef typename mapvector<Elem *, dof_id_type>::const_veclike_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node *>.
   */
  typedef typename mapvector<Node *, dof_id_type>::veclike_iterator             node_iterator_imp;
  typedef typename mapvector<Node *, dof_id_type>::const_veclike_iterator const_node_iterator_imp;
};

typedef DistributedMeshTempl<Real> DistributedMesh;

} // namespace libMesh

#endif // LIBMESH_DISTRIBUTED_MESH_H
