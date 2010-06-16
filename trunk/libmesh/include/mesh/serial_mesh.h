// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __serial_mesh_h__
#define __serial_mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "unstructured_mesh.h"






// Prototype this function so we can delcare it a friend.
class SerialMesh;
namespace MeshTools {
  namespace Private {
    void fix_broken_node_and_element_numbering (SerialMesh &);
  }
}


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
  SerialMesh (unsigned int dim=1);

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

  virtual unsigned int n_nodes () const { return _nodes.size(); }
  virtual unsigned int max_node_id () const { return _nodes.size(); }
  virtual void reserve_nodes (const unsigned int nn) { _nodes.reserve (nn); }
  virtual unsigned int n_elem ()  const { return _elements.size(); }
  virtual unsigned int n_active_elem () const;
  virtual unsigned int max_elem_id ()  const { return _elements.size(); }
  virtual void reserve_elem (const unsigned int ne) { _elements.reserve (ne); }

  // SerialMesh has no caches to update
  virtual void update_parallel_id_counts () {}

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
  virtual void delete_node (Node* n) ;
  virtual void renumber_node (unsigned int old_id, unsigned int new_id);
  virtual Elem* add_elem (Elem* e) ;
  virtual Elem* insert_elem (Elem* e) ;
  virtual void delete_elem (Elem* e) ;
  virtual void renumber_elem (unsigned int old_id, unsigned int new_id);

  friend void MeshTools::Private::fix_broken_node_and_element_numbering (SerialMesh &);
  
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
  
  element_iterator active_local_subdomain_elements_begin (const unsigned int subdomain_id);
  element_iterator active_local_subdomain_elements_end   (const unsigned int subdomain_id);
  
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
  
  const_element_iterator active_local_subdomain_elements_begin (const unsigned int subdomain_id) const;
  const_element_iterator active_local_subdomain_elements_end   (const unsigned int subdomain_id) const;
  
  
  
  
  
  
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
  std::vector<Node*> _nodes;
  
  /**
   * The elements in the mesh.
   */
  std::vector<Elem*> _elements;

private:
  
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





#endif
