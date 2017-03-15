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



#ifndef LIBMESH_DISTRIBUTED_MESH_H
#define LIBMESH_DISTRIBUTED_MESH_H

// Local Includes
#include "libmesh/mapvector.h"
#include "libmesh/unstructured_mesh.h"

// C++ Includes
#include <cstddef>
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

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  DistributedMesh (const Parallel::Communicator & comm_in,
                   unsigned char dim=1);

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  DistributedMesh (unsigned char dim=1);
#endif


  /**
   * Copy-constructor.  This should be able to take a
   * replicated or distributed mesh.
   */
  DistributedMesh (const UnstructuredMesh & other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * distributed mesh.
   */
  DistributedMesh (const DistributedMesh & other_mesh);

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual UniquePtr<MeshBase> clone () const libmesh_override
  { return UniquePtr<MeshBase>(new DistributedMesh(*this)); }

  /**
   * Destructor.
   */
  virtual ~DistributedMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear() libmesh_override;

  /**
   * Redistribute elements between processors.  This gets called
   * automatically by the Partitioner, and is a no-op in the case of a
   * serialized mesh.
   */
  virtual void redistribute () libmesh_override;

  /**
   * Recalculate cached data after elements and nodes have been
   * repartitioned.
   */
  virtual void update_post_partitioning () libmesh_override;

  /**
   * @returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const libmesh_override
  { return _is_serial; }

  /**
   * @returns \p true if new elements and nodes can and should be
   * created in synchronization on all processors, \p false otherwise
   */
  virtual bool is_replicated () const libmesh_override
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
  virtual void libmesh_assert_valid_parallel_ids() const libmesh_override;

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
   * Renumber a parallel objects container
   * Returns the smallest globally unused id for that
   * container.
   */
  template <typename T>
  dof_id_type renumber_dof_objects (mapvector<T *,dof_id_type> &);

  /**
   * Remove NULL elements from arrays
   */
  virtual void renumber_nodes_and_elements () libmesh_override;

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor
   */
  virtual void allgather() libmesh_override;

  /**
   * Gathers all elements and nodes of the mesh onto
   * processor zero
   */
  virtual void gather_to_zero() libmesh_override;

  /**
   * Deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   */
  virtual void delete_remote_elements() libmesh_override;

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

  // Cached methods that can be called in serial
  virtual dof_id_type n_nodes () const libmesh_override { return _n_nodes; }
  virtual dof_id_type max_node_id () const libmesh_override { return _max_node_id; }
  virtual void reserve_nodes (const dof_id_type) libmesh_override {}
  virtual dof_id_type n_elem () const libmesh_override { return _n_elem; }
  virtual dof_id_type n_active_elem () const libmesh_override;
  virtual dof_id_type max_elem_id () const libmesh_override { return _max_elem_id; }
  virtual void reserve_elem (const dof_id_type) libmesh_override {}

  // Parallel only method to update the caches
  virtual void update_parallel_id_counts () libmesh_override;

  // And more parallel only methods to test non-cached values
  virtual dof_id_type parallel_n_nodes () const libmesh_override;
  dof_id_type parallel_max_node_id () const;
  virtual dof_id_type parallel_n_elem () const libmesh_override;
  dof_id_type parallel_max_elem_id () const;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const libmesh_override;
#endif

  virtual const Point & point (const dof_id_type i) const libmesh_override;

  virtual const Node * node_ptr (const dof_id_type i) const libmesh_override;
  virtual Node * node_ptr (const dof_id_type i) libmesh_override;

  virtual const Node * query_node_ptr (const dof_id_type i) const libmesh_override;
  virtual Node * query_node_ptr (const dof_id_type i) libmesh_override;

  virtual const Elem * elem_ptr (const dof_id_type i) const libmesh_override;
  virtual Elem * elem_ptr (const dof_id_type i) libmesh_override;

  virtual const Elem * query_elem_ptr (const dof_id_type i) const libmesh_override;
  virtual Elem * query_elem_ptr (const dof_id_type i) libmesh_override;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node * add_point (const Point & p,
                            const dof_id_type id = DofObject::invalid_id,
                            const processor_id_type proc_id = DofObject::invalid_processor_id) libmesh_override;
  virtual Node * add_node (Node * n) libmesh_override;

  /**
   * Calls add_node().
   */
  virtual Node * insert_node(Node * n) libmesh_override;

  virtual void delete_node (Node * n) libmesh_override;
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id) libmesh_override;
  virtual Elem * add_elem (Elem * e) libmesh_override;
  virtual Elem * insert_elem (Elem * e) libmesh_override;
  virtual void delete_elem (Elem * e) libmesh_override;
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id) libmesh_override;

  /**
   * There is no reason for a user to ever call this function.
   *
   * This function restores a previously broken element/node numbering such that
   * \p mesh.node_ref(n).id() == n.
   */
  virtual void fix_broken_node_and_element_numbering () libmesh_override;

public:
  /**
   * Elem iterator accessor functions.
   */
  virtual element_iterator elements_begin () libmesh_override;
  virtual element_iterator elements_end () libmesh_override;
  virtual const_element_iterator elements_begin() const libmesh_override;
  virtual const_element_iterator elements_end() const libmesh_override;

  virtual element_iterator active_elements_begin () libmesh_override;
  virtual element_iterator active_elements_end () libmesh_override;
  virtual const_element_iterator active_elements_begin() const libmesh_override;
  virtual const_element_iterator active_elements_end() const libmesh_override;

  virtual element_iterator ancestor_elements_begin () libmesh_override;
  virtual element_iterator ancestor_elements_end () libmesh_override;
  virtual const_element_iterator ancestor_elements_begin() const libmesh_override;
  virtual const_element_iterator ancestor_elements_end() const libmesh_override;

  virtual element_iterator subactive_elements_begin () libmesh_override;
  virtual element_iterator subactive_elements_end () libmesh_override;
  virtual const_element_iterator subactive_elements_begin() const libmesh_override;
  virtual const_element_iterator subactive_elements_end() const libmesh_override;

  virtual element_iterator not_active_elements_begin () libmesh_override;
  virtual element_iterator not_active_elements_end () libmesh_override;
  virtual const_element_iterator not_active_elements_begin() const libmesh_override;
  virtual const_element_iterator not_active_elements_end() const libmesh_override;

  virtual element_iterator not_ancestor_elements_begin () libmesh_override;
  virtual element_iterator not_ancestor_elements_end () libmesh_override;
  virtual const_element_iterator not_ancestor_elements_begin() const libmesh_override;
  virtual const_element_iterator not_ancestor_elements_end() const libmesh_override;

  virtual element_iterator not_subactive_elements_begin () libmesh_override;
  virtual element_iterator not_subactive_elements_end () libmesh_override;
  virtual const_element_iterator not_subactive_elements_begin() const libmesh_override;
  virtual const_element_iterator not_subactive_elements_end() const libmesh_override;

  virtual element_iterator local_elements_begin () libmesh_override;
  virtual element_iterator local_elements_end () libmesh_override;
  virtual const_element_iterator local_elements_begin () const libmesh_override;
  virtual const_element_iterator local_elements_end () const libmesh_override;

  virtual element_iterator semilocal_elements_begin () libmesh_override;
  virtual element_iterator semilocal_elements_end () libmesh_override;
  virtual const_element_iterator semilocal_elements_begin () const libmesh_override;
  virtual const_element_iterator semilocal_elements_end () const libmesh_override;

  virtual element_iterator active_semilocal_elements_begin () libmesh_override;
  virtual element_iterator active_semilocal_elements_end () libmesh_override;
  virtual const_element_iterator active_semilocal_elements_begin () const libmesh_override;
  virtual const_element_iterator active_semilocal_elements_end () const libmesh_override;

  virtual element_iterator facelocal_elements_begin () libmesh_override;
  virtual element_iterator facelocal_elements_end () libmesh_override;
  virtual const_element_iterator facelocal_elements_begin () const libmesh_override;
  virtual const_element_iterator facelocal_elements_end () const libmesh_override;

  virtual element_iterator not_local_elements_begin () libmesh_override;
  virtual element_iterator not_local_elements_end () libmesh_override;
  virtual const_element_iterator not_local_elements_begin () const libmesh_override;
  virtual const_element_iterator not_local_elements_end () const libmesh_override;

  virtual element_iterator active_local_elements_begin () libmesh_override;
  virtual element_iterator active_local_elements_end () libmesh_override;
  virtual const_element_iterator active_local_elements_begin () const libmesh_override;
  virtual const_element_iterator active_local_elements_end () const libmesh_override;

  virtual element_iterator active_not_local_elements_begin () libmesh_override;
  virtual element_iterator active_not_local_elements_end () libmesh_override;
  virtual const_element_iterator active_not_local_elements_begin () const libmesh_override;
  virtual const_element_iterator active_not_local_elements_end () const libmesh_override;

  virtual element_iterator level_elements_begin (unsigned int level) libmesh_override;
  virtual element_iterator level_elements_end (unsigned int level) libmesh_override;
  virtual const_element_iterator level_elements_begin (unsigned int level) const libmesh_override;
  virtual const_element_iterator level_elements_end (unsigned int level) const libmesh_override;

  virtual element_iterator not_level_elements_begin (unsigned int level) libmesh_override;
  virtual element_iterator not_level_elements_end (unsigned int level) libmesh_override;
  virtual const_element_iterator not_level_elements_begin (unsigned int level) const libmesh_override;
  virtual const_element_iterator not_level_elements_end (unsigned int level) const libmesh_override;

  virtual element_iterator local_level_elements_begin (unsigned int level) libmesh_override;
  virtual element_iterator local_level_elements_end (unsigned int level) libmesh_override;
  virtual const_element_iterator local_level_elements_begin (unsigned int level) const libmesh_override;
  virtual const_element_iterator local_level_elements_end (unsigned int level) const libmesh_override;

  virtual element_iterator local_not_level_elements_begin (unsigned int level) libmesh_override;
  virtual element_iterator local_not_level_elements_end (unsigned int level) libmesh_override;
  virtual const_element_iterator local_not_level_elements_begin (unsigned int level) const libmesh_override;
  virtual const_element_iterator local_not_level_elements_end (unsigned int level) const libmesh_override;

  virtual element_iterator pid_elements_begin (processor_id_type proc_id) libmesh_override;
  virtual element_iterator pid_elements_end (processor_id_type proc_id) libmesh_override;
  virtual const_element_iterator pid_elements_begin (processor_id_type proc_id) const libmesh_override;
  virtual const_element_iterator pid_elements_end (processor_id_type proc_id) const libmesh_override;

  virtual element_iterator type_elements_begin (ElemType type) libmesh_override;
  virtual element_iterator type_elements_end (ElemType type) libmesh_override;
  virtual const_element_iterator type_elements_begin (ElemType type) const libmesh_override;
  virtual const_element_iterator type_elements_end (ElemType type) const libmesh_override;

  virtual element_iterator active_type_elements_begin (ElemType type) libmesh_override;
  virtual element_iterator active_type_elements_end (ElemType type) libmesh_override;
  virtual const_element_iterator active_type_elements_begin (ElemType type) const libmesh_override;
  virtual const_element_iterator active_type_elements_end (ElemType type) const libmesh_override;

  virtual element_iterator active_pid_elements_begin (processor_id_type proc_id) libmesh_override;
  virtual element_iterator active_pid_elements_end (processor_id_type proc_id) libmesh_override;
  virtual const_element_iterator active_pid_elements_begin (processor_id_type proc_id) const libmesh_override;
  virtual const_element_iterator active_pid_elements_end (processor_id_type proc_id) const libmesh_override;

  virtual element_iterator unpartitioned_elements_begin () libmesh_override;
  virtual element_iterator unpartitioned_elements_end () libmesh_override;
  virtual const_element_iterator unpartitioned_elements_begin () const libmesh_override;
  virtual const_element_iterator unpartitioned_elements_end () const libmesh_override;

  virtual element_iterator active_unpartitioned_elements_begin () libmesh_override;
  virtual element_iterator active_unpartitioned_elements_end () libmesh_override;
  virtual const_element_iterator active_unpartitioned_elements_begin () const libmesh_override;
  virtual const_element_iterator active_unpartitioned_elements_end () const libmesh_override;

  virtual element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) libmesh_override;
  virtual element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) libmesh_override;
  virtual const_element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) const libmesh_override;
  virtual const_element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) const libmesh_override;

  virtual element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) libmesh_override;
  virtual element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) libmesh_override;
  virtual const_element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) const libmesh_override;
  virtual const_element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) const libmesh_override;

  virtual element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) libmesh_override;
  virtual element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) libmesh_override;
  virtual const_element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) const libmesh_override;
  virtual const_element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) const libmesh_override;

  virtual element_iterator ghost_elements_begin () libmesh_override;
  virtual element_iterator ghost_elements_end () libmesh_override;
  virtual const_element_iterator ghost_elements_begin () const libmesh_override;
  virtual const_element_iterator ghost_elements_end () const libmesh_override;

  virtual element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) libmesh_override;

  virtual element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) libmesh_override;

  virtual const_element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) const libmesh_override;

  virtual const_element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) const libmesh_override;

#ifdef LIBMESH_ENABLE_AMR
  virtual element_iterator flagged_elements_begin (unsigned char rflag) libmesh_override;
  virtual element_iterator flagged_elements_end (unsigned char rflag) libmesh_override;
  virtual const_element_iterator flagged_elements_begin (unsigned char rflag) const libmesh_override;
  virtual const_element_iterator flagged_elements_end (unsigned char rflag) const libmesh_override;

  virtual element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                       processor_id_type pid) libmesh_override;
  virtual element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                     processor_id_type pid) libmesh_override;
  virtual const_element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                             processor_id_type pid) const libmesh_override;
  virtual const_element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                           processor_id_type pid) const libmesh_override;
#endif // LIBMESH_ENABLE_AMR

  /**
   * Node iterator accessor functions.
   */
  virtual node_iterator nodes_begin () libmesh_override;
  virtual node_iterator nodes_end () libmesh_override;
  virtual const_node_iterator nodes_begin () const libmesh_override;
  virtual const_node_iterator nodes_end () const libmesh_override;

  virtual node_iterator active_nodes_begin () libmesh_override;
  virtual node_iterator active_nodes_end () libmesh_override;
  virtual const_node_iterator active_nodes_begin () const libmesh_override;
  virtual const_node_iterator active_nodes_end () const libmesh_override;

  virtual node_iterator local_nodes_begin () libmesh_override;
  virtual node_iterator local_nodes_end () libmesh_override;
  virtual const_node_iterator local_nodes_begin () const libmesh_override;
  virtual const_node_iterator local_nodes_end () const libmesh_override;

  virtual node_iterator pid_nodes_begin (processor_id_type proc_id) libmesh_override;
  virtual node_iterator pid_nodes_end (processor_id_type proc_id) libmesh_override;
  virtual const_node_iterator pid_nodes_begin (processor_id_type proc_id) const libmesh_override;
  virtual const_node_iterator pid_nodes_end (processor_id_type proc_id) const libmesh_override;

  virtual node_iterator bid_nodes_begin (boundary_id_type bndry_id) libmesh_override;
  virtual node_iterator bid_nodes_end (boundary_id_type bndry_id) libmesh_override;
  virtual const_node_iterator bid_nodes_begin (boundary_id_type bndry_id) const libmesh_override;
  virtual const_node_iterator bid_nodes_end (boundary_id_type bndry_id) const libmesh_override;

  virtual node_iterator bnd_nodes_begin () libmesh_override;
  virtual node_iterator bnd_nodes_end () libmesh_override;
  virtual const_node_iterator bnd_nodes_begin () const libmesh_override;
  virtual const_node_iterator bnd_nodes_end () const libmesh_override;

  virtual node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) libmesh_override;
  virtual node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) libmesh_override;
  virtual const_node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) const libmesh_override;
  virtual const_node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) const libmesh_override;


protected:

  /**
   * The verices (spatial coordinates) of the mesh.
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
  typedef mapvector<Elem *, dof_id_type>::veclike_iterator             elem_iterator_imp;
  typedef mapvector<Elem *, dof_id_type>::const_veclike_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node *>.
   */
  typedef mapvector<Node *, dof_id_type>::veclike_iterator             node_iterator_imp;
  typedef mapvector<Node *, dof_id_type>::const_veclike_iterator const_node_iterator_imp;
};


} // namespace libMesh

#endif // LIBMESH_DISTRIBUTED_MESH_H
