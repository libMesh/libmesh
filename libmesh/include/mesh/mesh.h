// $Id: mesh.h,v 1.15 2005-06-11 05:11:30 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __mesh_h__
#define __mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "mesh_base.h"

// Forward Declarations
class MeshData;


/**
 * The \p Mesh class is derived from the \p MeshBase class.
 * The user will typically want to instantiate and use the
 * Mesh class in her applications.
 * In order to use the adaptive mesh refinment capabilities
 * of the library, first instantiate a MeshRefinement object
 * with a reference to this class.  Then call the appropriate
 * refinement functions from that object.  To interact with the 
 * boundary, instantiate a BoundaryMesh with a reference to
 * this class, and then use that object's functionality.
*/

// ------------------------------------------------------------
// Mesh class definition
class Mesh : public MeshBase
{
 public:

  /**
   * Constructor.  Requires the dimension and optionally
   * a processor id.  Note that \p proc_id should always
   * be provided for multiprocessor applications.
   */
  Mesh (unsigned int d);

  /**
   * Copy-constructor.
   */
  Mesh (const Mesh& other_mesh);  

  /**
   * Destructor.
   */
  ~Mesh();

  /**
   * Reads the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.  This is now the only
   * way to read a mesh.  The \p Mesh then initializes its data
   * structures and is ready for use.
   *
   * In order to read the UNV and TetGen file types, you must
   * also pass a separate pointer to the MeshData object you will
   * use with this mesh, since these read methods expect it.
   */
  void read (const std::string& name,
	     MeshData* mesh_data=NULL);
  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   *
   * In order to write the UNV and TetGen file types, you must
   * also pass a separate pointer to the MeshData object you have been
   * using with this mesh, since these write methods expect it.
   */
  void write (const std::string& name,
	      MeshData* mesh_data=NULL);
  
  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write (const std::string& name,
	      const std::vector<Number>& values,
	      const std::vector<std::string>& variable_names);

  /**
   * Clear all internal data.
   */
  void clear();

  /**
   * Converts a (conforming, non-refined) mesh with linear 
   * elements into a mesh with second-order elements.  For 
   * example, a mesh consisting of \p Tet4 will be converted
   * to a mesh with \p Tet10 etc.  Note that for some elements
   * like \p Hex8 there exist @e two higher order equivalents,
   * \p Hex20 and \p Hex27.  When \p full_ordered is \p true
   * (default), then \p Hex27 is built.  Otherwise, \p Hex20
   * is built.  The same holds obviously for \p Quad4, \p Prism6
   * ...
   */
  void all_second_order (const bool full_ordered=true);
  
  /**
   * Generates a new mesh containing all the elements which
   * are assigned to processor \p pid.  This mesh is written
   * to the pid_mesh reference which you must create and pass
   * to the function.
   */
  void create_pid_mesh (Mesh& pid_mesh,
			const unsigned int pid) const;
  
  /**
   * Constructs a mesh called "new_mesh" from the current mesh by
   * iterating over the elements between it and it_end and adding
   * them to the new mesh.
   */
  void create_submesh (Mesh& new_mesh,
		       const_element_iterator& it,
		       const const_element_iterator& it_end) const;
  


  /**
   * Virtual Functions Inherited From MeshBase which must be
   * redefined.
   * TODO[JWP]: Remember to comment later!!!!
   */
  virtual unsigned int n_nodes () const { return _nodes.size(); }
  virtual void reserve_nodes (const unsigned int nn) { _nodes.reserve (nn); }
  virtual unsigned int n_elem ()  const { return _elements.size(); }
  virtual void reserve_elem (const unsigned int ne) { _elements.reserve (ne); }

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
  virtual Node* add_point (const Point& n) ;
  virtual void delete_node (Node* n) ;
  virtual Elem* add_elem (Elem* e) ;
  virtual void delete_elem (Elem* e) ;

  /**
   * Other functions from MeshBase requiring re-definition.
   */
  virtual void find_neighbors ();
  virtual void renumber_nodes_and_elements ();

  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract ();


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



public:
  /**
   * Elem iterator accessor functions.
   */
  element_iterator elements_begin ();
  element_iterator elements_end   ();

  element_iterator active_elements_begin ();
  element_iterator active_elements_end   ();

  element_iterator not_active_elements_begin ();
  element_iterator not_active_elements_end   ();

  element_iterator local_elements_begin ();
  element_iterator local_elements_end   ();

  element_iterator active_local_elements_begin ();
  element_iterator active_local_elements_end   ();

  element_iterator level_elements_begin (const unsigned int level);
  element_iterator level_elements_end   (const unsigned int level);

  element_iterator not_level_elements_begin (const unsigned int level);
  element_iterator not_level_elements_end   (const unsigned int level);

  element_iterator pid_elements_begin (const unsigned int proc_id);
  element_iterator pid_elements_end   (const unsigned int proc_id);

  element_iterator type_elements_begin (const ElemType type);
  element_iterator type_elements_end   (const ElemType type);

  element_iterator active_type_elements_begin (const ElemType type);
  element_iterator active_type_elements_end   (const ElemType type);

  element_iterator active_pid_elements_begin (const unsigned int proc_id);
  element_iterator active_pid_elements_end   (const unsigned int proc_id);

  
  
  /**
   * const Elem iterator accessor functions.
   */
  const_element_iterator elements_begin() const;
  const_element_iterator elements_end()   const;
  
  const_element_iterator active_elements_begin() const;
  const_element_iterator active_elements_end()   const;
  
  const_element_iterator not_active_elements_begin() const;
  const_element_iterator not_active_elements_end()   const;

  const_element_iterator local_elements_begin () const;
  const_element_iterator local_elements_end   () const;

  const_element_iterator active_local_elements_begin () const;
  const_element_iterator active_local_elements_end   () const;

  const_element_iterator level_elements_begin (const unsigned int level) const;
  const_element_iterator level_elements_end   (const unsigned int level) const;

  const_element_iterator not_level_elements_begin (const unsigned int level) const;
  const_element_iterator not_level_elements_end   (const unsigned int level) const;

  const_element_iterator pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator pid_elements_end   (const unsigned int proc_id) const;

  const_element_iterator type_elements_begin (const ElemType type) const;
  const_element_iterator type_elements_end   (const ElemType type) const;

  const_element_iterator active_type_elements_begin (const ElemType type) const;
  const_element_iterator active_type_elements_end   (const ElemType type) const;

  const_element_iterator active_pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator active_pid_elements_end   (const unsigned int proc_id) const;
  
  
  
  
  
  
  
  /**
   * non-const Node iterator accessor functions.
   */
  node_iterator nodes_begin();
  node_iterator nodes_end();
  
  node_iterator active_nodes_begin();
  node_iterator active_nodes_end();
  
  
  /**
   * const Node iterator accessor functions.
   */
  const_node_iterator nodes_begin() const;
  const_node_iterator nodes_end()   const;

  const_node_iterator active_nodes_begin() const;
  const_node_iterator active_nodes_end()   const;

  
};





#endif
