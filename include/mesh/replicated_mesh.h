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



#ifndef LIBMESH_REPLICATED_MESH_H
#define LIBMESH_REPLICATED_MESH_H

// Local Includes
#include "libmesh/unstructured_mesh.h"
#include "libmesh/boundary_info.h"

// C++ Includes
#include <cstddef>

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

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  ReplicatedMesh (unsigned char dim=1);
#endif


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
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual UniquePtr<MeshBase> clone () const libmesh_override
  { return UniquePtr<MeshBase>(new ReplicatedMesh(*this)); }

  /**
   * Destructor.
   */
  virtual ~ReplicatedMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear() libmesh_override;

  /**
   * Remove NULL elements from arrays
   */
  virtual void renumber_nodes_and_elements () libmesh_override;

  virtual dof_id_type n_nodes () const libmesh_override
  { return cast_int<dof_id_type>(_nodes.size()); }

  virtual dof_id_type parallel_n_nodes () const libmesh_override
  { return cast_int<dof_id_type>(_nodes.size()); }

  virtual dof_id_type max_node_id () const libmesh_override
  { return cast_int<dof_id_type>(_nodes.size()); }

  virtual void reserve_nodes (const dof_id_type nn) libmesh_override
  { _nodes.reserve (nn); }

  virtual dof_id_type n_elem () const libmesh_override
  { return cast_int<dof_id_type>(_elements.size()); }

  virtual dof_id_type parallel_n_elem () const libmesh_override
  { return cast_int<dof_id_type>(_elements.size()); }

  virtual dof_id_type n_active_elem () const libmesh_override;

  virtual dof_id_type max_elem_id () const libmesh_override
  { return cast_int<dof_id_type>(_elements.size()); }

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  virtual unique_id_type parallel_max_unique_id () const libmesh_override;
#endif

  virtual void reserve_elem (const dof_id_type ne) libmesh_override
  { _elements.reserve (ne); }

  virtual void update_parallel_id_counts () libmesh_override;

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
   * Insert \p Node \p n into the Mesh at a location consistent with
   * n->id(), allocating extra storage if necessary.  Throws an error if:
   * .) n==NULL
   * .) n->id() == DofObject::invalid_id
   * .) A node already exists in position n->id().
   *
   * This function differs from the ReplicatedMesh::add_node() function,
   * which is only capable of appending nodes at the end of the nodes
   * storage.
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

  /**
   * Stitch \p other_mesh to this mesh so that this mesh is the union of the two meshes.
   * \p this_mesh_boundary and \p other_mesh_boundary are used to specify a dim-1 dimensional
   * surface on which we seek to merge any "overlapping" nodes, where we use the parameter
   * \p tol as a relative tolerance (relative to the smallest edge length on the surfaces
   * being stitched) to determine whether or not nodes are overlapping.
   * If \p clear_stitched_boundary_ids==true, this function clears boundary_info IDs in this
   * mesh associated \p this_mesh_boundary and \p other_mesh_boundary.
   * If \p use_binary_search is true, we use an optimized "sort then binary search" algorithm
   * for finding matching nodes. Otherwise we use a N^2 algorithm (which can be more reliable
   * at dealing with slightly misaligned meshes).
   * If \p enforce_all_nodes_match_on_boundaries is true, we throw an error if the number of
   * nodes on the specified boundaries don't match the number of nodes that were merged.
   * This is a helpful error check in some cases.
   * If \p skip_find_neighbors is true, a faster stitching method is used, where the lists of
   * neighbors for each elements are copied as well and patched, without calling the time-consuming
   * find_neighbors() function.
   */
  void stitch_meshes (ReplicatedMesh & other_mesh,
                      boundary_id_type this_mesh_boundary,
                      boundary_id_type other_mesh_boundary,
                      Real tol=TOLERANCE,
                      bool clear_stitched_boundary_ids=false,
                      bool verbose=true,
                      bool use_binary_search=true,
                      bool enforce_all_nodes_match_on_boundaries=false);

  /**
   * Similar to stitch_meshes, except that we stitch two adjacent surfaces within this mesh.
   */
  void stitch_surfaces (boundary_id_type boundary_id_1,
                        boundary_id_type boundary_id_2,
                        Real tol=TOLERANCE,
                        bool clear_stitched_boundary_ids=false,
                        bool verbose=true,
                        bool use_binary_search=true,
                        bool enforce_all_nodes_match_on_boundaries=false);

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
  std::vector<Node *> _nodes;

  /**
   * The elements in the mesh.
   */
  std::vector<Elem *> _elements;

private:

  /**
   * Helper function for stitch_meshes and stitch_surfaces
   * that does the mesh stitching.
   */
  void stitching_helper (ReplicatedMesh * other_mesh,
                         boundary_id_type boundary_id_1,
                         boundary_id_type boundary_id_2,
                         Real tol,
                         bool clear_stitched_boundary_ids,
                         bool verbose,
                         bool use_binary_search,
                         bool enforce_all_nodes_match_on_boundaries,
                         bool skip_find_neighbors);

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
