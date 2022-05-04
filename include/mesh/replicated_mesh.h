// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_REPLICATED_MESH_H
#define LIBMESH_REPLICATED_MESH_H

// Local Includes
#include "libmesh/unstructured_mesh.h"

// C++ Includes
#include <cstddef>
#include <memory>
#include <unordered_map>

namespace libMesh
{

/**
 * The \p ReplicatedMesh class is derived from the \p MeshBase class,
 * and is used to store identical copies of a full mesh data
 * structure on each processor.
 *
 * Most methods for this class are found in MeshBase, and most
 * implementation details are found in UnstructuredMesh.
 *
 * \author Roy Stogner
 * \date 2007
 * \brief Mesh data structure replicated on all processors.
 */
class ReplicatedMesh : public UnstructuredMesh
{
public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  ReplicatedMesh (const Parallel::Communicator & comm_in,
                  unsigned char dim=1);

  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  ReplicatedMesh (const UnstructuredMesh & other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * serial mesh.
   */
  ReplicatedMesh (const ReplicatedMesh & other_mesh);

  /**
   * Move-constructor deleted in MeshBase.
   */
  ReplicatedMesh(ReplicatedMesh &&) = delete;

  /**
   * Copy assignment is not allowed.
   */
  ReplicatedMesh & operator= (const ReplicatedMesh &) = delete;

  /**
   * Move assignment operator.
  */
  ReplicatedMesh & operator= (ReplicatedMesh && other_mesh);

  /**
   * Shim to call the move assignment operator for this class
  */
  virtual MeshBase & assign(MeshBase && other_mesh) override;

  /**
   * Move node and elements from a ReplicatedMesh.
   */
  virtual void move_nodes_and_elements(MeshBase && other_mesh) override;

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual std::unique_ptr<MeshBase> clone () const override
  { return std::make_unique<ReplicatedMesh>(*this); }

  /**
   * Destructor.
   */
  virtual ~ReplicatedMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear() override;

  /**
   * Clear internal Elem data.
   */
  virtual void clear_elems() override;

  /**
   * Remove nullptr elements from arrays
   */
  virtual void renumber_nodes_and_elements () override;

  virtual dof_id_type n_nodes () const override final
  { return _n_nodes; }

  virtual dof_id_type parallel_n_nodes () const override final
  { return _n_nodes; }

  virtual dof_id_type max_node_id () const override final
  { return cast_int<dof_id_type>(_nodes.size()); }

  virtual void reserve_nodes (const dof_id_type nn) override final
  { _nodes.reserve (nn); }

  virtual dof_id_type n_elem () const override final
  { return _n_elem; }

  virtual dof_id_type parallel_n_elem () const override final
  { return _n_elem; }

  virtual dof_id_type n_active_elem () const override final;

  virtual dof_id_type max_elem_id () const override final
  { return cast_int<dof_id_type>(_elements.size()); }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const override final;
  virtual void set_next_unique_id(unique_id_type id) override final;
#endif

  virtual void reserve_elem (const dof_id_type ne) override final
  { _elements.reserve (ne); }

  virtual void update_parallel_id_counts () override;

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
   * Insert \p Node \p n into the Mesh at a location consistent with
   * n->id(), allocating extra storage if necessary.  Throws an error if:
   * .) n==nullptr
   * .) n->id() == DofObject::invalid_id
   * .) A node already exists in position n->id().
   *
   * This function differs from the ReplicatedMesh::add_node() function,
   * which is only capable of appending nodes at the end of the nodes
   * storage.
   */
  virtual Node * insert_node(Node * n) override final;
  virtual Node * insert_node(std::unique_ptr<Node> n) override final;

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

  /**
   * Return IDs of representative elements of all disconnected subdomains.
   * Subdomains are considered connected only when they are sharing at least
   * one d-1 dimensional object (side in 2D, face in 3D), where d is
   * the mesh dimension.
   * The optional argument can be used for getting the subdomain IDs of all
   * elements with element IDs as the index.
   * This function cannot be called for a mesh with hanging nodes from
   * adaptive mesh refinement.
   */
  std::vector<dof_id_type> get_disconnected_subdomains(std::vector<subdomain_id_type> * subdomain_ids = nullptr) const;

  /**
   * Return all points on boundary.
   * The key of the returned unordered map is the ID of a representative
   * element of all disconnected subdomains. Subdomains are considered
   * connected only when they are sharing at least one d-1 dimensional object
   * (side in 2D), where d is the mesh dimension.
   * The size of the unordered map value is the number of disconnected
   * boundaries for a subdomain. Boundaries are considered
   * connected only when they are sharing a d-2 dimensional object.
   * This function currently only works for 2D meshes.
   * The points of each boundary are ordered to form an enclosure.
   */
  std::unordered_map<dof_id_type, std::vector<std::vector<Point>>> get_boundary_points() const;

public:
  /**
   * Elem iterator accessor functions.
   */
  virtual element_iterator elements_begin () override final;
  virtual element_iterator elements_end () override final;
  virtual const_element_iterator elements_begin() const override final;
  virtual const_element_iterator elements_end() const override final;
  virtual SimpleRange<element_iterator> element_ptr_range() override final { return {elements_begin(), elements_end()}; }
  virtual SimpleRange<const_element_iterator> element_ptr_range() const override final  { return {elements_begin(), elements_end()}; }

  virtual element_iterator active_elements_begin () override final;
  virtual element_iterator active_elements_end () override final;
  virtual const_element_iterator active_elements_begin() const override final;
  virtual const_element_iterator active_elements_end() const override final;
  virtual SimpleRange<element_iterator> active_element_ptr_range() override final { return {active_elements_begin(), active_elements_end()}; }
  virtual SimpleRange<const_element_iterator> active_element_ptr_range() const override final  { return {active_elements_begin(), active_elements_end()}; }

  virtual element_iterator ancestor_elements_begin () override final;
  virtual element_iterator ancestor_elements_end () override final;
  virtual const_element_iterator ancestor_elements_begin() const override final;
  virtual const_element_iterator ancestor_elements_end() const override final;

  virtual element_iterator subactive_elements_begin () override final;
  virtual element_iterator subactive_elements_end () override final;
  virtual const_element_iterator subactive_elements_begin() const override final;
  virtual const_element_iterator subactive_elements_end() const override final;

  virtual element_iterator not_active_elements_begin () override final;
  virtual element_iterator not_active_elements_end () override final;
  virtual const_element_iterator not_active_elements_begin() const override final;
  virtual const_element_iterator not_active_elements_end() const override final;

  virtual element_iterator not_ancestor_elements_begin () override final;
  virtual element_iterator not_ancestor_elements_end () override final;
  virtual const_element_iterator not_ancestor_elements_begin() const override final;
  virtual const_element_iterator not_ancestor_elements_end() const override final;

  virtual element_iterator not_subactive_elements_begin () override final;
  virtual element_iterator not_subactive_elements_end () override final;
  virtual const_element_iterator not_subactive_elements_begin() const override final;
  virtual const_element_iterator not_subactive_elements_end() const override final;

  virtual element_iterator local_elements_begin () override final;
  virtual element_iterator local_elements_end () override final;
  virtual const_element_iterator local_elements_begin () const override final;
  virtual const_element_iterator local_elements_end () const override final;

  virtual element_iterator semilocal_elements_begin () override final;
  virtual element_iterator semilocal_elements_end () override final;
  virtual const_element_iterator semilocal_elements_begin () const override final;
  virtual const_element_iterator semilocal_elements_end () const override final;

  virtual element_iterator active_semilocal_elements_begin () override final;
  virtual element_iterator active_semilocal_elements_end () override final;
  virtual const_element_iterator active_semilocal_elements_begin () const override final;
  virtual const_element_iterator active_semilocal_elements_end () const override final;

  virtual element_iterator facelocal_elements_begin () override final;
  virtual element_iterator facelocal_elements_end () override final;
  virtual const_element_iterator facelocal_elements_begin () const override final;
  virtual const_element_iterator facelocal_elements_end () const override final;

  virtual element_iterator not_local_elements_begin () override final;
  virtual element_iterator not_local_elements_end () override final;
  virtual const_element_iterator not_local_elements_begin () const override final;
  virtual const_element_iterator not_local_elements_end () const override final;

  virtual element_iterator active_local_elements_begin () override final;
  virtual element_iterator active_local_elements_end () override final;
  virtual const_element_iterator active_local_elements_begin () const override final;
  virtual const_element_iterator active_local_elements_end () const override final;
  virtual SimpleRange<element_iterator> active_local_element_ptr_range() override final { return {active_local_elements_begin(), active_local_elements_end()}; }
  virtual SimpleRange<const_element_iterator> active_local_element_ptr_range() const override final  { return {active_local_elements_begin(), active_local_elements_end()}; }

  virtual element_iterator active_not_local_elements_begin () override final;
  virtual element_iterator active_not_local_elements_end () override final;
  virtual const_element_iterator active_not_local_elements_begin () const override final;
  virtual const_element_iterator active_not_local_elements_end () const override final;

  virtual element_iterator level_elements_begin (unsigned int level) override final;
  virtual element_iterator level_elements_end (unsigned int level) override final;
  virtual const_element_iterator level_elements_begin (unsigned int level) const override final;
  virtual const_element_iterator level_elements_end (unsigned int level) const override final;

  virtual element_iterator not_level_elements_begin (unsigned int level) override final;
  virtual element_iterator not_level_elements_end (unsigned int level) override final;
  virtual const_element_iterator not_level_elements_begin (unsigned int level) const override final;
  virtual const_element_iterator not_level_elements_end (unsigned int level) const override final;

  virtual element_iterator local_level_elements_begin (unsigned int level) override final;
  virtual element_iterator local_level_elements_end (unsigned int level) override final;
  virtual const_element_iterator local_level_elements_begin (unsigned int level) const override final;
  virtual const_element_iterator local_level_elements_end (unsigned int level) const override final;

  virtual element_iterator local_not_level_elements_begin (unsigned int level) override final;
  virtual element_iterator local_not_level_elements_end (unsigned int level) override final;
  virtual const_element_iterator local_not_level_elements_begin (unsigned int level) const override final;
  virtual const_element_iterator local_not_level_elements_end (unsigned int level) const override final;

  virtual element_iterator pid_elements_begin (processor_id_type proc_id) override final;
  virtual element_iterator pid_elements_end (processor_id_type proc_id) override final;
  virtual const_element_iterator pid_elements_begin (processor_id_type proc_id) const override final;
  virtual const_element_iterator pid_elements_end (processor_id_type proc_id) const override final;

  virtual element_iterator type_elements_begin (ElemType type) override final;
  virtual element_iterator type_elements_end (ElemType type) override final;
  virtual const_element_iterator type_elements_begin (ElemType type) const override final;
  virtual const_element_iterator type_elements_end (ElemType type) const override final;

  virtual element_iterator active_type_elements_begin (ElemType type) override final;
  virtual element_iterator active_type_elements_end (ElemType type) override final;
  virtual const_element_iterator active_type_elements_begin (ElemType type) const override final;
  virtual const_element_iterator active_type_elements_end (ElemType type) const override final;

  virtual element_iterator active_pid_elements_begin (processor_id_type proc_id) override final;
  virtual element_iterator active_pid_elements_end (processor_id_type proc_id) override final;
  virtual const_element_iterator active_pid_elements_begin (processor_id_type proc_id) const override final;
  virtual const_element_iterator active_pid_elements_end (processor_id_type proc_id) const override final;

  virtual element_iterator unpartitioned_elements_begin () override final;
  virtual element_iterator unpartitioned_elements_end () override final;
  virtual const_element_iterator unpartitioned_elements_begin () const override final;
  virtual const_element_iterator unpartitioned_elements_end () const override final;

  virtual element_iterator active_unpartitioned_elements_begin () override final;
  virtual element_iterator active_unpartitioned_elements_end () override final;
  virtual const_element_iterator active_unpartitioned_elements_begin () const override final;
  virtual const_element_iterator active_unpartitioned_elements_end () const override final;

  virtual element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) override final;
  virtual element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) override final;
  virtual const_element_iterator active_local_subdomain_elements_begin (subdomain_id_type subdomain_id) const override final;
  virtual const_element_iterator active_local_subdomain_elements_end (subdomain_id_type subdomain_id) const override final;
  virtual SimpleRange<element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) override final
  { return {active_local_subdomain_elements_begin(subdomain_id), active_local_subdomain_elements_end(subdomain_id)}; }
  virtual SimpleRange<const_element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) const override final
  { return {active_local_subdomain_elements_begin(subdomain_id), active_local_subdomain_elements_end(subdomain_id)}; }

  virtual element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) override final;
  virtual element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) override final;
  virtual const_element_iterator active_subdomain_elements_begin (subdomain_id_type subdomain_id) const override final;
  virtual const_element_iterator active_subdomain_elements_end (subdomain_id_type subdomain_id) const override final;
  virtual SimpleRange<element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) override final
  { return {active_subdomain_elements_begin(subdomain_id), active_subdomain_elements_end(subdomain_id)}; }
  virtual SimpleRange<const_element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type subdomain_id) const override final
  { return {active_subdomain_elements_begin(subdomain_id), active_subdomain_elements_end(subdomain_id)}; }

  virtual element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) override final;
  virtual element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) override final;
  virtual const_element_iterator active_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) const override final;
  virtual const_element_iterator active_subdomain_set_elements_end (std::set<subdomain_id_type> ss) const override final;
  virtual SimpleRange<element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) override final
  { return {active_subdomain_set_elements_begin(ss), active_subdomain_set_elements_end(ss)}; }
  virtual SimpleRange<const_element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const override final
  { return {active_subdomain_set_elements_begin(ss), active_subdomain_set_elements_end(ss)}; }

  virtual element_iterator active_local_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) override final;
  virtual element_iterator active_local_subdomain_set_elements_end (std::set<subdomain_id_type> ss) override final;
  virtual const_element_iterator active_local_subdomain_set_elements_begin (std::set<subdomain_id_type> ss) const override final;
  virtual const_element_iterator active_local_subdomain_set_elements_end (std::set<subdomain_id_type> ss) const override final;
  virtual SimpleRange<element_iterator> active_local_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) override final
  { return {active_local_subdomain_set_elements_begin(ss), active_local_subdomain_set_elements_end(ss)}; }
  virtual SimpleRange<const_element_iterator> active_local_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const override final
  { return {active_local_subdomain_set_elements_begin(ss), active_local_subdomain_set_elements_end(ss)}; }

  virtual element_iterator ghost_elements_begin () override final;
  virtual element_iterator ghost_elements_end () override final;
  virtual const_element_iterator ghost_elements_begin () const override final;
  virtual const_element_iterator ghost_elements_end () const override final;

  virtual element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) override final;

  virtual element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) override final;

  virtual const_element_iterator
  evaluable_elements_begin (const DofMap & dof_map,
                            unsigned int var_num = libMesh::invalid_uint) const override final;

  virtual const_element_iterator
  evaluable_elements_end (const DofMap & dof_map,
                          unsigned int var_num = libMesh::invalid_uint) const override final;

  virtual element_iterator
  multi_evaluable_elements_begin (std::vector<const DofMap *> dof_maps) override final;

  virtual element_iterator
  multi_evaluable_elements_end (std::vector<const DofMap *> dof_maps) override final;

  virtual const_element_iterator
  multi_evaluable_elements_begin (std::vector<const DofMap *> dof_maps) const override final;

  virtual const_element_iterator
  multi_evaluable_elements_end (std::vector<const DofMap *> dof_maps) const override final;

#ifdef LIBMESH_ENABLE_AMR
  virtual element_iterator flagged_elements_begin (unsigned char rflag) override final;
  virtual element_iterator flagged_elements_end (unsigned char rflag) override final;
  virtual const_element_iterator flagged_elements_begin (unsigned char rflag) const override final;
  virtual const_element_iterator flagged_elements_end (unsigned char rflag) const override final;

  virtual element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                       processor_id_type pid) override final;
  virtual element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                     processor_id_type pid) override final;
  virtual const_element_iterator flagged_pid_elements_begin (unsigned char rflag,
                                                             processor_id_type pid) const override final;
  virtual const_element_iterator flagged_pid_elements_end (unsigned char rflag,
                                                           processor_id_type pid) const override final;
#endif // LIBMESH_ENABLE_AMR

  /**
   * Node iterator accessor functions.
   */
  virtual node_iterator nodes_begin () override final;
  virtual node_iterator nodes_end () override final;
  virtual const_node_iterator nodes_begin () const override final;
  virtual const_node_iterator nodes_end () const override final;
  virtual SimpleRange<node_iterator> node_ptr_range() override final { return {nodes_begin(), nodes_end()}; }
  virtual SimpleRange<const_node_iterator> node_ptr_range() const override final { return {nodes_begin(), nodes_end()}; }

  virtual node_iterator active_nodes_begin () override final;
  virtual node_iterator active_nodes_end () override final;
  virtual const_node_iterator active_nodes_begin () const override final;
  virtual const_node_iterator active_nodes_end () const override final;

  virtual node_iterator local_nodes_begin () override final;
  virtual node_iterator local_nodes_end () override final;
  virtual const_node_iterator local_nodes_begin () const override final;
  virtual const_node_iterator local_nodes_end () const override final;
  virtual SimpleRange<node_iterator> local_node_ptr_range() override final { return {local_nodes_begin(), local_nodes_end()}; }
  virtual SimpleRange<const_node_iterator> local_node_ptr_range() const override final { return {local_nodes_begin(), local_nodes_end()}; }

  virtual node_iterator pid_nodes_begin (processor_id_type proc_id) override final;
  virtual node_iterator pid_nodes_end (processor_id_type proc_id) override final;
  virtual const_node_iterator pid_nodes_begin (processor_id_type proc_id) const override final;
  virtual const_node_iterator pid_nodes_end (processor_id_type proc_id) const override final;

  virtual node_iterator bid_nodes_begin (boundary_id_type bndry_id) override final;
  virtual node_iterator bid_nodes_end (boundary_id_type bndry_id) override final;
  virtual const_node_iterator bid_nodes_begin (boundary_id_type bndry_id) const override final;
  virtual const_node_iterator bid_nodes_end (boundary_id_type bndry_id) const override final;

  virtual node_iterator bnd_nodes_begin () override final;
  virtual node_iterator bnd_nodes_end () override final;
  virtual const_node_iterator bnd_nodes_begin () const override final;
  virtual const_node_iterator bnd_nodes_end () const override final;

  virtual node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) override final;
  virtual node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) override final;
  virtual const_node_iterator
  evaluable_nodes_begin (const DofMap & dof_map,
                         unsigned int var_num = libMesh::invalid_uint) const override final;
  virtual const_node_iterator
  evaluable_nodes_end (const DofMap & dof_map,
                       unsigned int var_num = libMesh::invalid_uint) const override final;

  virtual node_iterator
  multi_evaluable_nodes_begin (std::vector<const DofMap *> dof_maps) override final;
  virtual node_iterator
  multi_evaluable_nodes_end (std::vector<const DofMap *> dof_maps) override final;
  virtual const_node_iterator
  multi_evaluable_nodes_begin (std::vector<const DofMap *> dof_maps) const override final;
  virtual const_node_iterator
  multi_evaluable_nodes_end (std::vector<const DofMap *> dof_maps) const override final;

protected:

  /**
   * The vertices (spatial coordinates) of the mesh.
   */
  std::vector<Node *> _nodes;

  dof_id_type _n_nodes;

  /**
   * The elements in the mesh.
   */
  std::vector<Elem *> _elements;

  dof_id_type _n_elem;

private:

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Elem *>.
   */
  typedef std::vector<Elem *>::iterator             elem_iterator_imp;
  typedef std::vector<Elem *>::const_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node *>.
   */
  typedef std::vector<Node *>::iterator             node_iterator_imp;
  typedef std::vector<Node *>::const_iterator const_node_iterator_imp;
};

} // namespace libMesh



#endif // LIBMESH_REPLICATED_MESH_H
