// $Id: mesh_base.h,v 1.53 2007-06-04 16:32:44 jwpeterson Exp $

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



#ifndef __mesh_base_h__
#define __mesh_base_h__



// C++ Includes   -----------------------------------
#include <string>

// forward declarations
class Elem;
class Node;
class Point;
class Partitioner;
class BoundaryInfo;

// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "enum_elem_type.h"
#include "variant_filter_iterator.h"
#include "multi_predicates.h"
#include "auto_ptr.h"




/**
 * This is the \p MeshBase class. This class provides all the data necessary
 * to describe a geometric entity.  It allows for the description of a
 * \p dim dimensional object that lives in \p DIM-dimensional space.
 * \par
 * A mesh is made of nodes and elements, and this class provides data
 * structures to store and access both.  A mesh may be partitioned into a
 * number of subdomains, and this class provides that functionality.
 * Furthermore, this class provides functions for reading and writing a
 * mesh to disk in various formats.
 *
 * \author  Benjamin S. Kirk
 * \date    $Date: 2007-06-04 16:32:44 $
 * \version $Revision: 1.53 $
 */


// ------------------------------------------------------------
// MeshBase class definition
class MeshBase
{
public:

  /**
   * Constructor.  Requires \p d, the dimension of the mesh.
   */
  MeshBase (unsigned int d);
  
  /**
   * Copy-constructor.
   */
  MeshBase (const MeshBase& other_mesh);

  /**
   * Destructor.
   */
  virtual ~MeshBase ();

  /**
   * This class holds the boundary information.  It can store nodes, edges,
   * and faces with a corresponding id that facilitates setting boundary
   * conditions.
   */
  AutoPtr<BoundaryInfo> boundary_info;
  
  /**
   * Deletes all the data that are currently stored.
   */
  virtual void clear ();
  
  /**
   * @returns \p true if the mesh has been prepared via a call
   * to \p prepare_for_use, \p false otherwise.
   */
  bool is_prepared () const
  { return _is_prepared; }
  
  /**
   * Returns the logical dimension of the mesh.
   */
  unsigned int mesh_dimension () const
  { return static_cast<unsigned int>(_dim); }
  
  /**
   * Returns the spatial dimension of the mesh.  Note that this is
   * defined at compile time in the header \p libmesh_common.h.
   */
  unsigned int spatial_dimension () const
  { return static_cast<unsigned int>(DIM); }
  
  /**
   * Returns the number of nodes in the mesh. This function and others must
   * be defined in derived classes since the MeshBase class has no specific
   * storage for nodes or elements.
   */
  virtual unsigned int n_nodes () const = 0; 

  /**
   * Reserves space for a known number of nodes.
   * Note that this method may or may not do anything, depending
   * on the actual \p Mesh implementation.  If you know the number
   * of nodes you will add and call this method before repeatedly
   * calling \p add_point() the implementation will be more efficient.
   */
  virtual void reserve_nodes (const unsigned int nn) = 0; 
  
  /**
   * Returns the number of elements in the mesh.
   */
  virtual unsigned int n_elem () const = 0; 

  /**
   * Reserves space for a known number of elements.
   * Note that this method may or may not do anything, depending
   * on the actual \p Mesh implementation.  If you know the number
   * of elements you will add and call this method before repeatedly
   * calling \p add_point() the implementation will be more efficient.
   */
  virtual void reserve_elem (const unsigned int ne) = 0; 

  /**
   * Returns the number of active elements in the mesh.  Implemented
   * in terms of active_element_iterators.
   */
  unsigned int n_active_elem () const;

  /**
   * Returns the number of elements on processor \p proc.
   */
  unsigned int n_elem_on_proc (const unsigned int proc) const;

  /**
   * Returns the number of elements on the local processor.
   */
  unsigned int n_local_elem () const
  { return this->n_elem_on_proc (libMesh::processor_id()); }

  /**
   * Returns the number of active elements on processor \p proc.
   */
  unsigned int n_active_elem_on_proc (const unsigned int proc) const;

  /**
   * Returns the number of active elements on the local processor.
   */
  unsigned int n_active_local_elem () const
  { return this->n_active_elem_on_proc (libMesh::processor_id()); }
  
  /**
   * This function returns the number of elements that will be written
   * out in the Tecplot format.  For example, a 9-noded quadrilateral will
   * be broken into 4 linear sub-elements for plotting purposes.  Thus, for
   * a mesh of 2 \p QUAD9 elements  \p n_tecplot_elem() will return 8.
   * Implemented in terms of element_iterators.
   */
  unsigned int n_sub_elem () const;

  /**
   * Same, but only counts active elements.
   */
  unsigned int n_active_sub_elem () const;
  
  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ point.
   */  
  virtual const Point& point (const unsigned int i) const = 0;

  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ node.
   */  
  virtual const Node& node (const unsigned int i) const = 0;
  
  /**
   * Return a reference to the \f$ i^{th} \f$ node.
   */  
  virtual Node& node (const unsigned int i) = 0;
  
  /**
   * Return a pointer to the \f$ i^{th} \f$ node.
   */  
  virtual const Node* node_ptr (const unsigned int i) const = 0;

  /**
   * Return a pointer to the \f$ i^{th} \f$ node.
   */  
  virtual Node* & node_ptr (const unsigned int i) = 0;

  /**
   * Return a pointer to the \f$ i^{th} \f$ element.
   */
  virtual Elem* elem (const unsigned int i) const = 0;

  /**
   * Add \p Node \p n to the vertex array.  The node will be appended to the
   * end of the vertex array. 
   */
  virtual Node* add_point (const Point& n) = 0;

  /**
   * Removes the Node n from the mesh.
   */
  virtual void delete_node (Node* n) = 0;
  
  /**
   * Add elem \p e to the end of the element array.
   */
  virtual Elem* add_elem (Elem* e) = 0;

  /**
   * Removes element \p e from the mesh. Note that calling this
   * method may produce isolated nodes, i.e. nodes not connected
   * to any element.
   */
  virtual void delete_elem (Elem* e) = 0;
		      
  /**
   * Locate element face (edge in 2D) neighbors.  This is done with the help
   * of a \p std::map that functions like a hash table.  When this function is
   * called only elements with \p NULL neighbor pointers are considered, so
   * the first call should take the longest.  Subsequent calls will only
   * consider new elements and the elements that lie on the boundary.
   * After this routine is called all the elements with a \p NULL neighbor
   * pointer are guaranteed to be on the boundary.  Thus this routine is
   * useful for automatically determining the boundaries of the domain.
   */
  virtual void find_neighbors () = 0;
  
  /**
   * After partitoning a mesh it is useful to renumber the nodes and elements
   * so that they lie in contiguous blocks on the processors.  This method
   * does just that.
   */
  virtual void renumber_nodes_and_elements () = 0;    

#ifdef ENABLE_AMR
  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract () = 0;
#endif

  /**
   * Prepare a newly created (or read) mesh for use.
   * This involves 3 steps:
   *  1.) call \p find_neighbors()
   *  2.) call \p partition()
   *  3.) call \p renumber_nodes_and_elements() 
   *
   * The read_xda_file boolean flag is true when prepare_for_use
   * is called from Mesh::read after reading an xda file.  It prevents
   * the renumbering of nodes and elements.  In general, leave this at
   * the default value of false.
   */
  void prepare_for_use (const bool skip_renumber_nodes_and_elements=false);
  
  /**
   * Call the default partitioner (currently \p metis_partition()).
   */
  void partition (const unsigned int n_parts=libMesh::n_processors());
  
  /**
   * Returns the number of subdomains in the global mesh. Note that it is
   * convenient to have one subdomain on each processor on parallel machines,
   * however this is not required. Multiple subdomains can exist on the same
   * processor.
   */
  unsigned int n_subdomains () const
  { return _n_sbd; }

  /**
   * Returns the number of partitions which have been defined via
   * a call to either mesh.partition() or by building a Partitioner
   * object and calling partition.  Note that the partitioner objects
   * are responsible for setting this value.
   */
  unsigned int n_partitions () const
  { return _n_parts; }
  
  /**
   * @returns the number of processors used in the
   * current simulation.
   */
  unsigned int n_processors () const
  { return libMesh::n_processors(); }

  /**
   * @returns the subdomain id for this processor.
   */
  unsigned int processor_id () const
  { return libMesh::processor_id(); }

  /**
   * @returns a string containing relevant information
   * about the mesh.
   */
  std::string get_info () const;

  /**
   * Prints relevant information about the mesh.
   */
  void print_info (std::ostream& os=std::cout) const;

  /**
   * Equivalent to calling print_info() above, but now you can write:
   * Mesh mesh;
   * std::cout << mesh << std::endl;
   */
  friend std::ostream& operator << (std::ostream& os, const MeshBase& m);

  /**
   * We need an empty, generic class to act as a predicate for this
   * and derived mesh classes.
   */
  typedef Predicates::multi_predicate Predicate;

  /**
   * structs for the element_iterator's.
   * Note that these iterators were designed so that derived mesh classes could use the
   * _same_ base class iterators interchangeably.  Their definition comes later in the
   * header file.
   */
  struct element_iterator;
  struct const_element_iterator;

  /**
   * structs for the node_iterator's.
   * Note that these iterators were designed so that derived mesh classes could use the
   * _same_ base class iterators interchangeably.  Their definition comes later in the
   * header file.
   */
  struct node_iterator;
  struct const_node_iterator;




public:


  
  /**
   * Elem iterator accessor functions.  These must be defined in
   * Concrete base classes.
   */
  virtual element_iterator elements_begin               () = 0;
  virtual element_iterator elements_end                 () = 0;
  virtual element_iterator active_elements_begin        () = 0;
  virtual element_iterator active_elements_end          () = 0;
  virtual element_iterator subactive_elements_begin     () = 0;
  virtual element_iterator subactive_elements_end       () = 0;
  virtual element_iterator not_active_elements_begin    () = 0;
  virtual element_iterator not_active_elements_end      () = 0;
  virtual element_iterator not_subactive_elements_begin () = 0;
  virtual element_iterator not_subactive_elements_end   () = 0;
  virtual element_iterator local_elements_begin         () = 0;
  virtual element_iterator local_elements_end           () = 0;
  virtual element_iterator active_local_elements_begin  () = 0;
  virtual element_iterator active_local_elements_end    () = 0;
  virtual element_iterator level_elements_begin         (const unsigned int level  ) = 0;
  virtual element_iterator level_elements_end           (const unsigned int level  ) = 0;
  virtual element_iterator not_level_elements_begin     (const unsigned int level  ) = 0;
  virtual element_iterator not_level_elements_end       (const unsigned int level  ) = 0;
  virtual element_iterator pid_elements_begin           (const unsigned int proc_id) = 0;
  virtual element_iterator pid_elements_end             (const unsigned int proc_id) = 0;
  virtual element_iterator type_elements_begin          (const ElemType type       ) = 0;
  virtual element_iterator type_elements_end            (const ElemType type       ) = 0;
  virtual element_iterator active_type_elements_begin   (const ElemType type       ) = 0;
  virtual element_iterator active_type_elements_end     (const ElemType type       ) = 0;
  virtual element_iterator active_pid_elements_begin    (const unsigned int proc_id) = 0;
  virtual element_iterator active_pid_elements_end      (const unsigned int proc_id) = 0;

  
  
  /**
   * const Elem iterator accessor functions.
   */
  virtual const_element_iterator elements_begin               () const = 0;
  virtual const_element_iterator elements_end                 () const = 0;
  virtual const_element_iterator active_elements_begin        () const = 0;
  virtual const_element_iterator active_elements_end          () const = 0;
  virtual const_element_iterator subactive_elements_begin     () const = 0;
  virtual const_element_iterator subactive_elements_end       () const = 0;
  virtual const_element_iterator not_active_elements_begin    () const = 0;
  virtual const_element_iterator not_active_elements_end      () const = 0;
  virtual const_element_iterator not_subactive_elements_begin () const = 0;
  virtual const_element_iterator not_subactive_elements_end   () const = 0;
  virtual const_element_iterator local_elements_begin         () const = 0;
  virtual const_element_iterator local_elements_end           () const = 0;
  virtual const_element_iterator active_local_elements_begin  () const = 0;
  virtual const_element_iterator active_local_elements_end    () const = 0;
  virtual const_element_iterator level_elements_begin         (const unsigned int level)   const = 0;
  virtual const_element_iterator level_elements_end           (const unsigned int level)   const = 0;
  virtual const_element_iterator not_level_elements_begin     (const unsigned int level)   const = 0;
  virtual const_element_iterator not_level_elements_end       (const unsigned int level)   const = 0;
  virtual const_element_iterator pid_elements_begin           (const unsigned int proc_id) const = 0;
  virtual const_element_iterator pid_elements_end             (const unsigned int proc_id) const = 0;
  virtual const_element_iterator type_elements_begin          (const ElemType type)        const = 0;
  virtual const_element_iterator type_elements_end            (const ElemType type)        const = 0;
  virtual const_element_iterator active_type_elements_begin   (const ElemType type)        const = 0;
  virtual const_element_iterator active_type_elements_end     (const ElemType type)        const = 0;
  virtual const_element_iterator active_pid_elements_begin    (const unsigned int proc_id) const = 0;
  virtual const_element_iterator active_pid_elements_end      (const unsigned int proc_id) const = 0;

  
  /**
   * non-const Node iterator accessor functions.
   */
  virtual node_iterator nodes_begin        () = 0;
  virtual node_iterator nodes_end          () = 0;
  virtual node_iterator active_nodes_begin () = 0;
  virtual node_iterator active_nodes_end   () = 0;
  virtual node_iterator local_nodes_begin  () = 0;
  virtual node_iterator local_nodes_end    () = 0;
  virtual node_iterator pid_nodes_begin    (const unsigned int proc_id) = 0;
  virtual node_iterator pid_nodes_end      (const unsigned int proc_id) = 0;


  /**
   * const Node iterator accessor functions.
   */
  virtual const_node_iterator nodes_begin        () const = 0;
  virtual const_node_iterator nodes_end          () const = 0;
  virtual const_node_iterator active_nodes_begin () const = 0;
  virtual const_node_iterator active_nodes_end   () const = 0;
  virtual const_node_iterator local_nodes_begin  () const = 0;
  virtual const_node_iterator local_nodes_end    () const = 0;
  virtual const_node_iterator pid_nodes_begin    (const unsigned int proc_id) const = 0;
  virtual const_node_iterator pid_nodes_end      (const unsigned int proc_id) const = 0;


  

  
  

  
  
protected:



  
  
  /**
   * Returns a writeable reference to the number of subdomains.
   */
  unsigned int& set_n_subdomains ()
  { return _n_sbd; }

  /**
   * Returns a writeable reference to the number of partitions.
   */
  unsigned int& set_n_partitions ()
  { return _n_parts; }
  
  /**
   * The number of subdomains the mesh has.
   * **NOTE** Not to be confused with the number of paritions!
   * The definition of subdomain can be anything the user wants,
   * e.g. a solid region bounded by a liquid region could be
   * referred to as subdomains 1 and 2, but those subdomains
   * could be partitioned over many processors.
   */
  unsigned int _n_sbd;

  /**
   * The number of partitions the mesh has.  This is set by
   * the partitioners, and may not be changed directly by
   * the user.
   * **NOTE** The number of partitions *need not* equal
   * libMesh::n_processors(), consider for example the case
   * where you simply want to partition a mesh on one
   * processor and view the result in GMV.
   */
  unsigned int _n_parts;
  
  /**
   * The logical dimension of the mesh.
   */     
  const unsigned int _dim;

  /**
   * Flag indicating if the mesh has been prepared for use.
   */
  bool _is_prepared;
  

  /**
   * The partitioner class is a friend so that it can set
   * the number of partitions.
   */
  friend class Partitioner;

  /**
   * Make the \p BoundaryInfo class a friend so that
   * it can create and interact with \p BoundaryMesh.
   */
  friend class BoundaryInfo;

};











/**
 * The definition of the element_iterator struct.
 */
struct
MeshBase::element_iterator :
variant_filter_iterator<MeshBase::Predicate,
			Elem*>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  element_iterator (const IterType& d,
		    const IterType& e,
		    const PredType& p ) :
    variant_filter_iterator<MeshBase::Predicate,
			    Elem*>(d,e,p) {}
};




/**
 * The definition of the const_element_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_element_iterator :
variant_filter_iterator<MeshBase::Predicate,
			Elem* const,
			Elem* const&,
			Elem* const*>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  const_element_iterator (const IterType& d,
			  const IterType& e,
			  const PredType& p ) :
    variant_filter_iterator<MeshBase::Predicate,
			    Elem* const,
			    Elem* const&,
			    Elem* const*>(d,e,p)  {}


  // The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
  // variant_filter_iterator copy constructor.  Note that this one is *not* templated!
  const_element_iterator (const MeshBase::element_iterator& rhs) :
    variant_filter_iterator<Predicate,
			    Elem* const,
			    Elem* const&,
			    Elem* const*>(rhs)
  {
    // std::cout << "Called element_iterator conversion-to-const ctor." << std::endl;
  }
};







/**
 * The definition of the node_iterator struct.
 */
struct
MeshBase::node_iterator :
variant_filter_iterator<MeshBase::Predicate,
			Node*>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  node_iterator (const IterType& d,
		 const IterType& e,
		 const PredType& p ) :
    variant_filter_iterator<MeshBase::Predicate,
			    Node*>(d,e,p) {}
};




/**
 * The definition of the const_node_iterator struct.  It is similar to the regular
 * iterator above, but also provides an additional conversion-to-const ctor.
 */
struct
MeshBase::const_node_iterator :
variant_filter_iterator<MeshBase::Predicate,
			Node* const,
			Node* const &,
			Node* const *>
{
  // Templated forwarding ctor -- forwards to appropriate variant_filter_iterator ctor
  template <typename PredType, typename IterType>
  const_node_iterator (const IterType& d,
		       const IterType& e,
		       const PredType& p ) :
    variant_filter_iterator<MeshBase::Predicate,
			    Node* const,
			    Node* const &,
			    Node* const *>(d,e,p)  {}


  // The conversion-to-const ctor.  Takes a regular iterator and calls the appropriate
  // variant_filter_iterator copy constructor.  Note that this one is *not* templated!
  const_node_iterator (const MeshBase::node_iterator& rhs) :
    variant_filter_iterator<Predicate,
			    Node* const,
			    Node* const &,
			    Node* const *>(rhs)
  {
    // std::cout << "Called node_iterator conversion-to-const ctor." << std::endl;
  }
};


#endif
