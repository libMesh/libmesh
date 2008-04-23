// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __parallel_mesh_h__
#define __parallel_mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "elem.h"
#include "mapvector.h"
#include "node.h"
#include "unstructured_mesh.h"



// Prototype this function so we can delcare it a friend.
class ParallelMesh;
namespace MeshTools {
  namespace Private {
    void fix_broken_node_and_element_numbering (ParallelMesh &);
  }
}


/**
 * The \p ParallelMesh class is derived from the \p MeshBase class,
 * and is intended to provide identical functionality to the user
 * but be fully parallelized in memory.
 * By "is intended" I mean that it doesn't work that way yet.  Don't
 * use this class unless you're developing or debugging it.
*/

// ------------------------------------------------------------
// UnstructuredMesh class definition
class ParallelMesh : public UnstructuredMesh
{
 public:

  /**
   * Constructor.  Requires the dimension and optionally
   * a processor id.  Note that \p proc_id should always
   * be provided for multiprocessor applications.
   */
  ParallelMesh (unsigned int d);

  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  ParallelMesh (const UnstructuredMesh& other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * parallel mesh.
   */
  ParallelMesh (const ParallelMesh& other_mesh);

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual AutoPtr<MeshBase> clone () const
    { return AutoPtr<MeshBase>(new ParallelMesh(*this)); }

  /**
   * Destructor.
   */
  virtual ~ParallelMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear();

  /**
   * @returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const
    { return _is_serial; }

#ifdef DEBUG
  /**
   * Verify id and processor_id consistency of a parallel
   * objects container.
   * Calls libmesh_assert() on each possible failure in that container.
   */
  template <typename T>
  void libmesh_assert_valid_parallel_object_ids(const mapvector<T*>&) const;

  /**
   * Verify id and processor_id consistency of our elements and
   * nodes containers.
   * Calls libmesh_assert() on each possible failure.
   */
  void libmesh_assert_valid_parallel_ids() const;

  /**
   * Verify refinement_flag and p_refinement_flag consistency of our
   * elements containers.
   * Calls libmesh_assert() on each possible failure.
   */
  void libmesh_assert_valid_parallel_flags() const;
#endif

  /**
   * Renumber a parallel objects container
   * Returns the smallest globally unused id for that
   * container.
   */
  template <typename T>
  unsigned int renumber_dof_objects (mapvector<T*>&);

  /**
   * Remove NULL elements from arrays
   */
  virtual void renumber_nodes_and_elements ();

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor
   */
  virtual void allgather();

  /**
   * Deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   *
   */
  virtual void delete_remote_elements();

  // Cached methods that can be called in serial
  virtual unsigned int n_nodes () const { return _n_nodes; }
  virtual unsigned int max_node_id () const { return _max_node_id; }
  virtual void reserve_nodes (const unsigned int) { }
  virtual unsigned int n_elem () const { return _n_elem; }
  virtual unsigned int max_elem_id () const { return _max_elem_id; }
  virtual void reserve_elem (const unsigned int) { }

  // Parallel only method to update the caches
  virtual void update_parallel_id_counts ();

  // And more parallel only methods to test non-cached values
  unsigned int parallel_n_nodes () const;
  unsigned int parallel_max_node_id () const;
  unsigned int parallel_n_elem () const;
  unsigned int parallel_max_elem_id () const;

  /**
   * For meshes that don't store points/elems, these functions may be an issue!
   */
  virtual const Point& point (const unsigned int i) const ;
  virtual const Node&  node  (const unsigned int i) const ;
  virtual Node& node (const unsigned int i) ;
  virtual const Node* node_ptr (const unsigned int i) const ;
  virtual Node* & node_ptr (const unsigned int i) ;
  virtual Elem* elem (const unsigned int i) const ;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node* add_point (const Point& p,
			   const unsigned int id =
			     DofObject::invalid_id,
			   const unsigned int proc_id =
			     DofObject::invalid_processor_id);
  virtual Node* add_node (Node* n) ;
  virtual Node* insert_node (Node* n);
  virtual void delete_node (Node* n) ;
  virtual void renumber_node (unsigned int old_id, unsigned int new_id);
  virtual Elem* add_elem (Elem* e) ;
  virtual Elem* insert_elem (Elem* e) ;
  virtual void delete_elem (Elem* e) ;
  virtual void renumber_elem (unsigned int old_id, unsigned int new_id);

  friend void MeshTools::Private::fix_broken_node_and_element_numbering (ParallelMesh &);
  
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

  element_iterator pid_elements_begin (const unsigned int proc_id);
  element_iterator pid_elements_end   (const unsigned int proc_id);

  element_iterator type_elements_begin (const ElemType type);
  element_iterator type_elements_end   (const ElemType type);

  element_iterator active_type_elements_begin (const ElemType type);
  element_iterator active_type_elements_end   (const ElemType type);

  element_iterator active_pid_elements_begin (const unsigned int proc_id);
  element_iterator active_pid_elements_end   (const unsigned int proc_id);

  element_iterator unpartitioned_elements_begin ();
  element_iterator unpartitioned_elements_end ();

  
  
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

  const_element_iterator pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator pid_elements_end   (const unsigned int proc_id) const;

  const_element_iterator type_elements_begin (const ElemType type) const;
  const_element_iterator type_elements_end   (const ElemType type) const;

  const_element_iterator active_type_elements_begin (const ElemType type) const;
  const_element_iterator active_type_elements_end   (const ElemType type) const;

  const_element_iterator active_pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator active_pid_elements_end   (const unsigned int proc_id) const;

  const_element_iterator unpartitioned_elements_begin () const;
  const_element_iterator unpartitioned_elements_end () const;
  
  
  
  
  
  
  
  /**
   * non-const Node iterator accessor functions.
   */
  node_iterator nodes_begin();
  node_iterator nodes_end();
  
  node_iterator active_nodes_begin();
  node_iterator active_nodes_end();

  node_iterator local_nodes_begin  ();
  node_iterator local_nodes_end    ();
  
  node_iterator pid_nodes_begin (const unsigned int proc_id);
  node_iterator pid_nodes_end   (const unsigned int proc_id);

  /**
   * const Node iterator accessor functions.
   */
  const_node_iterator nodes_begin() const;
  const_node_iterator nodes_end()   const;

  const_node_iterator active_nodes_begin() const;
  const_node_iterator active_nodes_end()   const;

  const_node_iterator local_nodes_begin  () const;
  const_node_iterator local_nodes_end    () const;

  const_node_iterator pid_nodes_begin (const unsigned int proc_id) const;
  const_node_iterator pid_nodes_end   (const unsigned int proc_id) const;
 
protected:

  /**
   * The verices (spatial coordinates) of the mesh.
   */
  mapvector<Node*> _nodes;
  
  /**
   * The elements in the mesh.
   */
  mapvector<Elem*> _elements;

  /**
   * A boolean remembering whether we're serialized or not
   */
  bool _is_serial;

  /**
   * Cached data from the last renumber_nodes_and_elements call
   */
  unsigned int _n_nodes, _n_elem, _max_node_id, _max_elem_id;

  /**
   * Guaranteed globally unused IDs for use when adding new
   * nodes or elements.
   */
  unsigned int _next_free_local_node_id,
	       _next_free_local_elem_id;
  unsigned int _next_free_unpartitioned_node_id,
	       _next_free_unpartitioned_elem_id;

private:
  
  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Elem*>.
   */
  typedef mapvector<Elem*>::veclike_iterator             elem_iterator_imp;
  typedef mapvector<Elem*>::const_veclike_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node*>.
   */
  typedef mapvector<Node*>::veclike_iterator             node_iterator_imp;
  typedef mapvector<Node*>::const_veclike_iterator const_node_iterator_imp;
};


#endif
