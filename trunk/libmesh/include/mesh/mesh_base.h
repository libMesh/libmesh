// $Id: mesh_base.h,v 1.33 2004-11-12 00:42:40 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm>
#include <vector>
#include <set>
#include <string>


// forward declarations
class Elem;
class EquationSystems;


// Local Includes -----------------------------------
#include "libmesh_common.h"
#include "boundary_info.h"
#include "node.h"
#include "enum_elem_type.h"
#include "sphere.h"
#include "enum_order.h"
//#include "elem_iterators.h"
//#include "node_iterators.h"
#include "partitioner.h"

#include "variant_filter_iterator.h"
#include "multi_predicates.h"
#include "elem.h"

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
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.33 $
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
   * Deletes all the data that are currently stored.
   */
  virtual void clear ();
  
  /**
   * This class holds the boundary information.  It can store nodes, edges,
   * and faces with a corresponding id that facilitates setting boundary
   * conditions.
   */
  BoundaryInfo boundary_info;
  
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
   * Returns the number of nodes in the mesh.
   */
  unsigned int n_nodes () const { return _nodes.size(); }

  /**
   * Returns the number of elements in the mesh.
   */
  unsigned int n_elem ()  const { return _elements.size(); }

  /**
   * Returns the number of active elements in the mesh.
   */
  unsigned int n_active_elem () const;

  /**
   * Return a vector of all
   * element types for the mesh.
   */
  std::vector<ElemType> elem_types () const;
  
  /**
   * Return the number of elements of type \p type.
   */
  unsigned int n_elem_of_type (const ElemType type) const;

  /**
   * Return the number of active elements of type \p type.
   */
  unsigned int n_active_elem_of_type (const ElemType type) const;

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
   */
  unsigned int n_sub_elem () const;

  /**
   * Same, but only counts active elements.
   */
  unsigned int n_active_sub_elem () const;

  /**
   * This function returns the sum over all the elemenents of the number
   * of nodes per element.  This can be useful for partitioning hybrid meshes.
   * A feasible load balancing scheme is to keep the weight per processor as
   * uniform as possible.
   */
  unsigned int total_weight () const;
  
  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ point.
   */  
  const Point& point (const unsigned int i) const;

  /**
   * Return a constant reference (for reading only) to the
   * \f$ i^{th} \f$ node.
   */  
  const Node& node (const unsigned int i) const;
  
  /**
   * Return a reference to the \f$ i^{th} \f$ node.
   */  
  Node& node (const unsigned int i);
  
  /**
   * Return a pointer to the \f$ i^{th} \f$ node.
   */  
  const Node* node_ptr (const unsigned int i) const;

  /**
   * Return a pointer to the \f$ i^{th} \f$ node.
   */  
  Node* & node_ptr (const unsigned int i);

  /**
   * Return a constant reference to the \p nodes vector holding the nodes.
   */
  // const std::vector<Node*> & get_nodes () const { return _nodes; }

  /**
   * Add \p Node \p n to the vertex array.  The node will be appended to the
   * end of the vertex array. 
   */
  Node* add_point (const Point& n);

  /**
   * Return a pointer to the \f$ i^{th} \f$ element.
   */
  Elem* elem (const unsigned int i) const;

  /**
   * Return a reference to the \p cells vector holding the elements.
   */
  // const std::vector<Elem*> & get_elem () const { return _elements; }

  /**
   * Add elem \p e to the end of the element array.
   */
  Elem* add_elem (Elem* e);

  /**
   * Removes element \p e from the mesh. Note that calling this
   * method may produce isolated nodes, i.e. nodes not connected
   * to any element.
   */
  void delete_elem (Elem* e);
    

#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * convenient typedef for origin coordinates; so that the
   * \p build_inf_elem() methods know whether the respective
   * coordinate was given or not
   */
  typedef std::pair<bool, double> InfElemOriginValue;

  /**
   * Build infinite elements atop a volume-based mesh,
   * determine origin automatically.  Also returns the
   * origin as a \p const \p Point to make it more obvious that
   * the origin should not change after the infinite elements 
   * have been built.  When symmetry planes are present, use 
   * the version with optional symmetry switches.  
   * The flag \p be_verbose enables some diagnostic output.
   */
  const Point build_inf_elem (const bool be_verbose = false);

  /**
   * @returns the origin of the infinite elements.
   * Builds infinite elements atop a volume-based mesh.
   * Finds all faces on the outer boundary and build infinite elements
   * on them.  Using the \p InfElemOriginValue the user can
   * prescribe only selected origin coordinates.  The remaining
   * coordinates are computed from the center of the bounding box
   * of the mesh.
   *
   * During the search for faces on which infinite elements are built, 
   * @e interior faces that are not on symmetry planes are found, too.  
   * When an (optional) pointer to \p inner_boundary_nodes is provided, 
   * then this vector will be filled with the nodes that lie on the
   * inner boundary.
   *
   * Faces which lie in at least one symmetry plane are skipped.
   * The three optional booleans \p x_sym, \p y_sym,
   * \p z_sym indicate symmetry planes (through the origin, obviously)
   * perpendicular to the \p x, \p y and \p z direction, 
   * respectively.  
   * The flag \p be_verbose enables some diagnostic output.
   */
  const Point build_inf_elem (const InfElemOriginValue& origin_x,
			      const InfElemOriginValue& origin_y,
			      const InfElemOriginValue& origin_z,
			      const bool x_sym = false,
			      const bool y_sym = false,
			      const bool z_sym = false,
			      const bool be_verbose = false,
			      std::vector<const Node*>* inner_boundary_nodes = NULL);

protected:

  /**
   * Build infinite elements atop a volume-based mesh.
   * Actual implementation.
   */
  void build_inf_elem (const Point& origin,
		       const bool x_sym = false,
		       const bool y_sym = false,
		       const bool z_sym = false,
		       const bool be_verbose = false,
		       std::set<std::pair<unsigned int,
		                          unsigned int> >* inner_faces = NULL);

public:

		      
#endif

		      
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
  void find_neighbors ();
  
  /**
   * After calling this function the input vector \p nodes_to_elem_map
   * will contain the node to element connectivity.  That is to say
   * \p nodes_to_elem_map[i][j] is the global number of \f$ j^{th} \f$
   * element connected to node \p i.
   */
  void build_nodes_to_elem_map (std::vector<std::vector<unsigned int> >&
				nodes_to_elem_map) const;
  
  /**
   * The same, except element pointers are returned instead of indices.
   */
  void build_nodes_to_elem_map (std::vector<std::vector<const Elem*> >&
				nodes_to_elem_map) const;


  /**
   * Calling this function on a 2D mesh will convert all the elements
   * to triangles.  \p QUAD4s will be converted to \p TRI3s, \p QUAD8s
   * and \p QUAD9s will be converted to \p TRI6s. 
   */
  void all_tri ();
    
  /**
   * Call the default partitioner (currently \p metis_partition()).
   */
  virtual void partition (const unsigned int n_parts=libMesh::n_processors());

  /**
   * After partitoning a mesh it is useful to renumber the nodes and elements
   * so that they lie in contiguous blocks on the processors.  This method
   * does just that.
   */
  void renumber_nodes_and_elements ();    

  /**
   * Randomly perturb the nodal locations.  This function will
   * move each node \p factor fraction of its minimum neighboring
   * node separation distance.  Nodes on the boundary are not moved
   * by default, however they may be by setting the flag
   * \p perturb_boundary true.
   */
  void distort (const Real factor, const bool perturb_boundary=false);
  
  /**
   * Translates the mesh.  The grid points are translated in the
   * \p x direction by \p xt, in the \p y direction by \p yt,
   * etc...
   */
  void translate (const Real xt=0., const Real yt=0., const Real zt=0.); 

  /**
   * Rotates the mesh.  The grid points are rotated about the 
   * \p x axis by \p xr , about the \p y axis by \p yr,
   * etc...  
   */
  void rotate (const Real xr, const Real yr=0., const Real zr=0.); 

  /**
   * Scales the mesh.  The grid points are scaled in the
   * \p x direction by \p xs, in the \p y direction by \p ys,
   * etc...  If only \p xs is specified then the scaling is
   * assumed uniform in all directions.
   */
  void scale (const Real xs, const Real ys=0., const Real zs=0.);
    
  /**
   * @returns two points defining a cartesian box that bounds the
   * mesh.  The first entry in the pair is the mininum, the second 
   * is the maximim.
   */
  std::pair<Point, Point> bounding_box () const;

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere bounding_sphere () const;
  
  /**
   * @returns two points defining a cartesian box that bounds the
   * elements belonging to processor pid.  If no processor id is specified
   * the bounding box for the whole mesh is returned.
   */
  std::pair<Point, Point> 
  processor_bounding_box (const unsigned int pid = libMesh::invalid_uint) const;

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere 
  processor_bounding_sphere (const unsigned int pid = libMesh::invalid_uint) const;

  /**
   * @returns two points defining a Cartesian box that bounds the
   * elements belonging to subdomain sid.  If no subdomain id is specified
   * the bounding box for the whole mesh is returned.
   */
  std::pair<Point, Point> 
  subdomain_bounding_box (const unsigned int sid = libMesh::invalid_uint) const;

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere 
  subdomain_bounding_sphere (const unsigned int pid = libMesh::invalid_uint) const;


//   /**
//    * Builds the connectivity graph. The matrix \p conn is such that the
//    * valence of each node is on the diagonal and there is a -1 for each
//    * node connected to the node.
//    */
//   void build_L_graph (PetscMatrix<Number>& conn) const;
  
//   /**
//    * Builds the connectivity graph. The matrix \p conn is such that the
//    * valence of each node is on the diagonal and there is a -1 for each
//    * node connected to the node.
//    */
//   void build_script_L_graph (PetscMatrix<Number>& conn) const;
  
  /**
   * Returns the number of subdomains in the global mesh. Note that it is
   * convenient to have one subdomain on each processor on parallel machines,
   * however this is not required. Multiple subdomains can exist on the same
   * processor.
   *
   */
  unsigned int n_subdomains () const { return _n_sbd; }

  /**
   * Returns the number of partitions which have been defined via
   * a call to either mesh.partition() or by building a Partitioner
   * object and calling partition.  Note that the partitioner objects
   * are responsible for setting this value.
   */
  unsigned int n_partitions () const { return _n_parts; }
  
  /**
   * @returns the number of processors used in the
   * current simulation.
   */
  unsigned int n_processors () const { return libMesh::n_processors(); }

  /**
   * @returns the subdomain id for this processor.
   */
  unsigned int processor_id () const { return libMesh::processor_id(); }

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
   * Convenient typedefs for the element_iterators returned by the accessor functions.
   * Note that these iterators were designed so that derived mesh classes could use the
   * _same_ typedefs.
   */
  typedef variant_filter_iterator<Elem*      , Predicate>       element_iterator;
  typedef variant_filter_iterator<Elem* const, Predicate> const_element_iterator;

  /**
   * Convenient typedefs for the node_iterators returned by the accessor functions.
   * Note that these iterators were designed so that derived classes could use the
   * _same_ typedefs.
   */
  typedef variant_filter_iterator<Node*      , Predicate>       node_iterator;
  typedef variant_filter_iterator<Node* const, Predicate> const_node_iterator;

private:
  // Typedefs for the container implementation.  In this case,
  // it's just a std::vector<Elem*>.
  typedef std::vector<Elem*>::iterator             elem_iterator_imp;
  typedef std::vector<Elem*>::const_iterator const_elem_iterator_imp;

  // Typedefs for the container implementation.  In this case,
  // it's just a std::vector<Node*>.
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

  element_iterator local_elements_begin ();
  element_iterator local_elements_end   ();

  element_iterator active_local_elements_begin ();
  element_iterator active_local_elements_end   ();

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

  const_element_iterator local_elements_begin () const;
  const_element_iterator local_elements_end   () const;

  const_element_iterator active_local_elements_begin () const;
  const_element_iterator active_local_elements_end   () const;

  const_element_iterator not_level_elements_begin (const unsigned int level) const;
  const_element_iterator not_level_elements_end   (const unsigned int level)   const;

  const_element_iterator pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator pid_elements_end   (const unsigned int proc_id)   const;

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


  

  
  /**
   * Fills the vector "on_boundary" with flags that tell whether each node
   * is on the domain boundary (true)) or not (false).
   */
  void find_boundary_nodes(std::vector<bool>& on_boundary) const;
  
  /**
   * Prepare a newly created (or read) mesh for use.
   * This involves 3 steps:
   *  1.) call \p find_neighbors()
   *  2.) call \p partition()
   *  3.) call \p renumber_nodes_and_elements() 
   */
  virtual void prepare_for_use ();
  
  
protected:



  
  
  /**
   * Returns a writeable reference to the number of subdomains.
   */
  unsigned int& set_n_subdomains () { return _n_sbd; }

  /**
   * Returns a writeable reference to the number of partitions.
   */
  unsigned int& set_n_partitions () { return _n_parts; }
  
  /**
   * The verices (spatial coordinates) of the mesh.
   */
  std::vector<Node*> _nodes;

  /**
   * The elements in the mesh.
   */
  std::vector<Elem*> _elements;

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
   * Make the \p BoundaryInfo class a friend so that
   * it can create and interact with \p BoundaryMesh.
   */
  friend class BoundaryInfo;

  /**
   * The partitioner class is a friend so that it can set
   * the number of partitions.
   */
  friend class Partitioner;

};








// ------------------------------------------------------------
// MeshBase inline methods
inline
Elem* MeshBase::elem (const unsigned int i) const
{
  assert (i < this->n_elem());
  assert (_elements[i] != NULL);
  
  return _elements[i];
}



inline
const Point& MeshBase::point (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  

  return (*_nodes[i]);
}



inline
const Node& MeshBase::node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return (*_nodes[i]);
}



inline
Node& MeshBase::node (const unsigned int i)
{
  if (i >= this->n_nodes())
    {
      std::cout << " i=" << i
		<< ", n_nodes()=" << this->n_nodes()
		<< std::endl;
      error();
    }
  
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);

  return (*_nodes[i]);
}



inline
const Node* MeshBase::node_ptr (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() != Node::invalid_id);  
  
  return _nodes[i];
}



inline
Node* & MeshBase::node_ptr (const unsigned int i)
{
  assert (i < this->n_nodes());

  return _nodes[i];
}


inline
void MeshBase::delete_elem(Elem* e)
{
  assert (e != NULL);

  std::vector<Elem*>::iterator pos = std::find (_elements.begin(),
						_elements.end(),
						e);

  // Huh? Element not in the vector?
  assert (pos != _elements.end());

  // delete the element
  delete e;
  
  // explicitly NULL the pointer
  e    = NULL;
  *pos = NULL;
  
  //_elements.erase(pos);

  return;
}





#endif
