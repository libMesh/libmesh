// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SERIAL_MESH_H
#define LIBMESH_SERIAL_MESH_H

// Local Includes -----------------------------------
#include "libmesh/unstructured_mesh.h"
#include "libmesh/boundary_info.h"

// C++ Includes   -----------------------------------
#include <cstddef>

namespace libMesh
{






/**
 * The \p SerialMesh class is derived from the \p MeshBase class,
 * and currently represents the default Mesh implementation.
 * Most methods for this class are found in MeshBase, and most
 * implementation details are found in UnstructuredMesh.
*/

// ------------------------------------------------------------
// Mesh class definition
class SerialMesh : public UnstructuredMesh
{
 public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  SerialMesh (const Parallel::Communicator &comm,
	      unsigned int dim=1);

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  SerialMesh (unsigned int dim=1);
#endif


  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  SerialMesh (const UnstructuredMesh& other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * serial mesh.
   */
  SerialMesh (const SerialMesh& other_mesh);

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual AutoPtr<MeshBase> clone () const
    { return AutoPtr<MeshBase>(new SerialMesh(*this)); }

  /**
   * Destructor.
   */
  virtual ~SerialMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear();

  /**
   * Remove NULL elements from arrays
   */
  virtual void renumber_nodes_and_elements ();

  virtual dof_id_type n_nodes () const
  { return libmesh_cast_int<dof_id_type>(_nodes.size()); }

  virtual dof_id_type parallel_n_nodes () const
  { return libmesh_cast_int<dof_id_type>(_nodes.size()); }

  virtual dof_id_type max_node_id () const
  { return libmesh_cast_int<dof_id_type>(_nodes.size()); }

  virtual void reserve_nodes (const dof_id_type nn)
  { _nodes.reserve (nn); }

  virtual dof_id_type n_elem () const
  { return libmesh_cast_int<dof_id_type>(_elements.size()); }

  virtual dof_id_type parallel_n_elem () const
  { return libmesh_cast_int<dof_id_type>(_elements.size()); }

  virtual dof_id_type n_active_elem () const;

  virtual dof_id_type max_elem_id ()  const
  { return libmesh_cast_int<dof_id_type>(_elements.size()); }

  virtual void reserve_elem (const dof_id_type ne) { _elements.reserve (ne); }

  // SerialMesh has no caches to update
  virtual void update_parallel_id_counts () {}

  virtual const Point& point (const dof_id_type i) const ;
  virtual const Node&  node  (const dof_id_type i) const ;
  virtual Node& node (const dof_id_type i) ;
  virtual const Node* node_ptr (const dof_id_type i) const ;
  virtual Node* node_ptr (const dof_id_type i) ;
  virtual const Node* query_node_ptr (const dof_id_type i) const ;
  virtual Node* query_node_ptr (const dof_id_type i) ;
  virtual const Elem* elem (const dof_id_type i) const ;
  virtual Elem* elem (const dof_id_type i) ;
  virtual const Elem* query_elem (const dof_id_type i) const ;
  virtual Elem* query_elem (const dof_id_type i) ;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node* add_point (const Point& p,
			   const dof_id_type id =
			     DofObject::invalid_id,
			   const processor_id_type proc_id =
			     DofObject::invalid_processor_id);
  virtual Node* add_node (Node* n) ;
  virtual void delete_node (Node* n) ;
  virtual void renumber_node (dof_id_type old_id, dof_id_type new_id);
  virtual Elem* add_elem (Elem* e) ;
  virtual Elem* insert_elem (Elem* e) ;
  virtual void delete_elem (Elem* e) ;
  virtual void renumber_elem (dof_id_type old_id, dof_id_type new_id);

    /**
     * There is no reason for a user to ever call this function.
     *
     * This function restores a previously broken element/node numbering such that
     * \p mesh.node(n)->id() == n.
     */
  virtual void fix_broken_node_and_element_numbering ();

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
   * If enforce_all_nodes_match_on_boundaries is true, we throw an error if the number of
   * nodes on the specified boundaries don't match the number of nodes that were merged.
   * This is a helpful error check in some cases.
   */
  void stitch_meshes (SerialMesh& other_mesh,
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
  element_iterator elements_begin ();
  element_iterator elements_end   ();

  element_iterator active_elements_begin ();
  element_iterator active_elements_end   ();

  element_iterator ancestor_elements_begin ();
  element_iterator ancestor_elements_end   ();

  element_iterator subactive_elements_begin ();
  element_iterator subactive_elements_end   ();

  element_iterator not_active_elements_begin ();
  element_iterator not_active_elements_end   ();

  element_iterator not_ancestor_elements_begin ();
  element_iterator not_ancestor_elements_end   ();

  element_iterator not_subactive_elements_begin ();
  element_iterator not_subactive_elements_end   ();

  element_iterator local_elements_begin ();
  element_iterator local_elements_end   ();

  element_iterator not_local_elements_begin ();
  element_iterator not_local_elements_end   ();

  element_iterator active_local_elements_begin ();
  element_iterator active_local_elements_end   ();

  element_iterator active_not_local_elements_begin ();
  element_iterator active_not_local_elements_end   ();

  element_iterator level_elements_begin (const unsigned int level);
  element_iterator level_elements_end   (const unsigned int level);

  element_iterator not_level_elements_begin (const unsigned int level);
  element_iterator not_level_elements_end   (const unsigned int level);

  element_iterator local_level_elements_begin (const unsigned int level);
  element_iterator local_level_elements_end   (const unsigned int level);

  element_iterator local_not_level_elements_begin (const unsigned int level);
  element_iterator local_not_level_elements_end   (const unsigned int level);

  element_iterator pid_elements_begin (const processor_id_type proc_id);
  element_iterator pid_elements_end   (const processor_id_type proc_id);

  element_iterator type_elements_begin (const ElemType type);
  element_iterator type_elements_end   (const ElemType type);

  element_iterator active_type_elements_begin (const ElemType type);
  element_iterator active_type_elements_end   (const ElemType type);

  element_iterator active_pid_elements_begin (const processor_id_type proc_id);
  element_iterator active_pid_elements_end   (const processor_id_type proc_id);

  element_iterator unpartitioned_elements_begin ();
  element_iterator unpartitioned_elements_end ();

  element_iterator active_local_subdomain_elements_begin (const subdomain_id_type subdomain_id);
  element_iterator active_local_subdomain_elements_end   (const subdomain_id_type subdomain_id);

  element_iterator active_subdomain_elements_begin (const subdomain_id_type subdomain_id);
  element_iterator active_subdomain_elements_end   (const subdomain_id_type subdomain_id);

  /**
   * const Elem iterator accessor functions.
   */
  const_element_iterator elements_begin() const;
  const_element_iterator elements_end()   const;

  const_element_iterator active_elements_begin() const;
  const_element_iterator active_elements_end()   const;

  const_element_iterator ancestor_elements_begin() const;
  const_element_iterator ancestor_elements_end()   const;

  const_element_iterator subactive_elements_begin() const;
  const_element_iterator subactive_elements_end()   const;

  const_element_iterator not_active_elements_begin() const;
  const_element_iterator not_active_elements_end()   const;

  const_element_iterator not_ancestor_elements_begin() const;
  const_element_iterator not_ancestor_elements_end()   const;

  const_element_iterator not_subactive_elements_begin() const;
  const_element_iterator not_subactive_elements_end()   const;

  const_element_iterator local_elements_begin () const;
  const_element_iterator local_elements_end   () const;

  const_element_iterator not_local_elements_begin () const;
  const_element_iterator not_local_elements_end   () const;

  const_element_iterator active_local_elements_begin () const;
  const_element_iterator active_local_elements_end   () const;

  const_element_iterator active_not_local_elements_begin () const;
  const_element_iterator active_not_local_elements_end   () const;

  const_element_iterator level_elements_begin (const unsigned int level) const;
  const_element_iterator level_elements_end   (const unsigned int level) const;

  const_element_iterator not_level_elements_begin (const unsigned int level) const;
  const_element_iterator not_level_elements_end   (const unsigned int level) const;

  const_element_iterator local_level_elements_begin (const unsigned int level) const;
  const_element_iterator local_level_elements_end   (const unsigned int level) const;

  const_element_iterator local_not_level_elements_begin (const unsigned int level) const;
  const_element_iterator local_not_level_elements_end   (const unsigned int level) const;

  const_element_iterator pid_elements_begin (const processor_id_type proc_id) const;
  const_element_iterator pid_elements_end   (const processor_id_type proc_id) const;

  const_element_iterator type_elements_begin (const ElemType type) const;
  const_element_iterator type_elements_end   (const ElemType type) const;

  const_element_iterator active_type_elements_begin (const ElemType type) const;
  const_element_iterator active_type_elements_end   (const ElemType type) const;

  const_element_iterator active_pid_elements_begin (const processor_id_type proc_id) const;
  const_element_iterator active_pid_elements_end   (const processor_id_type proc_id) const;

  const_element_iterator unpartitioned_elements_begin () const;
  const_element_iterator unpartitioned_elements_end () const;

  const_element_iterator active_local_subdomain_elements_begin (const subdomain_id_type subdomain_id) const;
  const_element_iterator active_local_subdomain_elements_end   (const subdomain_id_type subdomain_id) const;

  const_element_iterator active_subdomain_elements_begin (const subdomain_id_type subdomain_id) const;
  const_element_iterator active_subdomain_elements_end   (const subdomain_id_type subdomain_id) const;






  /**
   * non-const Node iterator accessor functions.
   */
  node_iterator nodes_begin();
  node_iterator nodes_end();

  node_iterator active_nodes_begin();
  node_iterator active_nodes_end();

  node_iterator local_nodes_begin  ();
  node_iterator local_nodes_end    ();

  node_iterator pid_nodes_begin (const processor_id_type proc_id);
  node_iterator pid_nodes_end   (const processor_id_type proc_id);

  /**
   * const Node iterator accessor functions.
   */
  const_node_iterator nodes_begin() const;
  const_node_iterator nodes_end()   const;

  const_node_iterator active_nodes_begin() const;
  const_node_iterator active_nodes_end()   const;

  const_node_iterator local_nodes_begin  () const;
  const_node_iterator local_nodes_end    () const;

  const_node_iterator pid_nodes_begin (const processor_id_type proc_id) const;
  const_node_iterator pid_nodes_end   (const processor_id_type proc_id) const;

protected:

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  /**
   * Assign globally unique IDs to all DOF objects (Elements and Nodes)
   * if the library has been configured with unique_id support.
   */
  virtual void assign_unique_ids();
#endif

  /**
   * The verices (spatial coordinates) of the mesh.
   */
  std::vector<Node*> _nodes;

  /**
   * The elements in the mesh.
   */
  std::vector<Elem*> _elements;

private:

  /**
   * Helper function for stitch_meshes and stitch_surfaces
   * that does the mesh stitching.
   */
  void stitching_helper (SerialMesh* other_mesh,
                         boundary_id_type boundary_id_1,
                         boundary_id_type boundary_id_2,
                         Real tol,
                         bool clear_stitched_boundary_ids,
                         bool verbose,
                         bool use_binary_search,
                         bool enforce_all_nodes_match_on_boundaries);

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Elem*>.
   */
  typedef std::vector<Elem*>::iterator             elem_iterator_imp;
  typedef std::vector<Elem*>::const_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node*>.
   */
  typedef std::vector<Node*>::iterator             node_iterator_imp;
  typedef std::vector<Node*>::const_iterator const_node_iterator_imp;
};



} // namespace libMesh



#endif // LIBMESH_SERIAL_MESH_H
