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
  ReplicatedMesh (const MeshBase & other_mesh);

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
   * Shim to allow operator == (&) to behave like a virtual function
   * without having to be one.
   */
  virtual std::string_view subclass_first_difference_from (const MeshBase & other_mesh_base) const override;

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual std::unique_ptr<MeshBase> clone () const override
  {
    auto returnval = std::make_unique<ReplicatedMesh>(*this);
#ifdef DEBUG
    libmesh_assert(*returnval == *this);
#endif
    return returnval;
  }

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
   * @copydoc MeshBase::renumber_nodes_and_elements()
   *
   * Note that regardless of whether \p _skip_renumber_nodes_and_elements is true, this method will
   * always remove nullptr nodes and elements from arrays and consequently will always renumber. So in
   * practice \p _skip_renumber_nodes_and_elements means skip reordering of non-null nodes and
   * elements while allowing renumbering
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

  // Backwards compatibility
  virtual SimpleRange<element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) override final { return active_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<const_element_iterator> active_subdomain_elements_ptr_range(subdomain_id_type sid) const override final { return active_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) override final { return active_local_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<const_element_iterator> active_local_subdomain_elements_ptr_range(subdomain_id_type sid) const override final { return active_local_subdomain_element_ptr_range(sid); }
  virtual SimpleRange<element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) override final { return active_subdomain_set_element_ptr_range(ss); }
  virtual SimpleRange<const_element_iterator> active_subdomain_set_elements_ptr_range(std::set<subdomain_id_type> ss) const override final { return active_subdomain_set_element_ptr_range(ss); }

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
