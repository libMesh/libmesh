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



#ifndef __mesh_tools_h__
#define __mesh_tools_h__



// C++ Includes   -----------------------------------
#include <vector>
#include <set>

// Local Includes -----------------------------------
#include "libmesh.h"
#include "enum_elem_type.h"
#include "mesh_base.h"
#include "point.h" // some compilers want the full definition - I think so they can do 
                   // return-value-optimization for BoundingBox'es - BSK
// forward declarations
class SerialMesh;
class ParallelMesh;
class Sphere;
class Elem;

/**
 * Utility functions for operations on a \p Mesh object.  Here is where
 * useful functions for interfacing with a \p Mesh should be defined.
 * In general this namespace should be used to prevent the \p Mesh class
 * from becoming too cluttered.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision$
 */


// ------------------------------------------------------------
// MeshTools namespace
namespace MeshTools
{
  /**
   * Defines a Cartesian bounding box by the two
   * corner extremum.
   */
  typedef std::pair<Point, Point> BoundingBox;
  
  /**
   * This function returns the sum over all the elemenents of the number
   * of nodes per element.  This can be useful for partitioning hybrid meshes.
   * A feasible load balancing scheme is to keep the weight per processor as
   * uniform as possible.
   */
  unsigned int total_weight (const MeshBase& mesh);
  
  /**
   * This function returns the sum over all the elemenents on processor \p pid
   * of nodes per element.  This can be useful for partitioning hybrid meshes.
   * A feasible load balancing scheme is to keep the weight per processor as
   * uniform as possible.
   */
  unsigned int weight (const MeshBase& mesh, const unsigned int pid=libMesh::processor_id());
  
  /**
   * After calling this function the input vector \p nodes_to_elem_map
   * will contain the node to element connectivity.  That is to say
   * \p nodes_to_elem_map[i][j] is the global number of \f$ j^{th} \f$
   * element connected to node \p i.
   */
  void build_nodes_to_elem_map (const MeshBase& mesh,
				std::vector<std::vector<unsigned int> >& nodes_to_elem_map);
  
  /**
   * The same, except element pointers are returned instead of indices.
   */
  void build_nodes_to_elem_map (const MeshBase& mesh,
				std::vector<std::vector<const Elem*> >&	nodes_to_elem_map);


//   /**
//    * Calling this function on a 2D mesh will convert all the elements
//    * to triangles.  \p QUAD4s will be converted to \p TRI3s, \p QUAD8s
//    * and \p QUAD9s will be converted to \p TRI6s. 
//    */
//   void all_tri (MeshBase& mesh);

  /**
   * Fills the vector "on_boundary" with flags that tell whether each node
   * is on the domain boundary (true)) or not (false).
   */
  void find_boundary_nodes (const MeshBase& mesh,
			    std::vector<bool>& on_boundary);

  /**
   * @returns two points defining a cartesian box that bounds the
   * mesh.  The first entry in the pair is the mininum, the second 
   * is the maximim.
   */
  BoundingBox
  bounding_box (const MeshBase& mesh);

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere
  bounding_sphere (const MeshBase& mesh);
  
  /**
   * @returns two points defining a cartesian box that bounds the
   * elements belonging to processor pid.  If no processor id is specified
   * the bounding box for the whole mesh is returned.
   */
  BoundingBox
  processor_bounding_box (const MeshBase& mesh,
			  const unsigned int pid = libMesh::invalid_uint);

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere 
  processor_bounding_sphere (const MeshBase& mesh,
			     const unsigned int pid = libMesh::invalid_uint);

  /**
   * @returns two points defining a Cartesian box that bounds the
   * elements belonging to subdomain sid.  If no subdomain id is specified
   * the bounding box for the whole mesh is returned.
   */
  std::pair<Point, Point> 
  subdomain_bounding_box (const MeshBase& mesh,
			  const unsigned int sid = libMesh::invalid_uint);

  /**
   * Same, but returns a sphere instead of a box.
   */
  Sphere 
  subdomain_bounding_sphere (const MeshBase& mesh,
			     const unsigned int pid = libMesh::invalid_uint);


  /**
   * Return a vector of all element types for the mesh.  Implemented
   * in terms of element_iterators.
   */
  void elem_types (const MeshBase& mesh,
		   std::vector<ElemType>& et);
  
  /**
   * Return the number of elements of type \p type.  Implemented
   * in terms of type_element_iterators.
   */
  unsigned int n_elem_of_type (const MeshBase& mesh,
			       const ElemType type);

  /**
   * Return the number of active elements of type \p type.
   * Implemented in terms of active_type_element_iterators.
   */
  unsigned int n_active_elem_of_type (const MeshBase& mesh,
				      const ElemType type);

  /**
   * Return the number of elements of type \p type at the specified
   * refinement level.
   *
   * TODO: Replace all of the n_xxx_elem() functions like this with
   * a single function which takes a range of iterators and returns the
   * std::distance between them.
   */
  unsigned int n_non_subactive_elem_of_type_at_level(const MeshBase& mesh,
                                                     const ElemType type,
                                                     const unsigned int level);

  /**
   * Return the number of levels of refinement in the mesh.
   * Implemented by looping over all the local elements and finding the
   * maximum level, then summing in parallel.
   */
  unsigned int n_levels(const MeshBase& mesh);

  /**
   * Return the number of levels of refinement in the local mesh.
   * Implemented by looping over all the local elements and finding the
   * maximum level.
   */
  unsigned int n_local_levels(const MeshBase& mesh);
 
  /**
   * Return the number of levels of refinement in the active mesh.
   * Implemented by looping over all the active local elements and finding
   * the maximum level, then summing in parallel.
   */
  unsigned int n_active_levels(const MeshBase& mesh);

  /**
   * Return the number of levels of refinement in the active local mesh.
   * Implemented by looping over all the active local elements and finding
   * the maximum level.
   */
  unsigned int n_active_local_levels(const MeshBase& mesh);
 
  /**
   * Builds a set of node IDs for nodes which belong to non-subactive
   * elements.  Non-subactive elements are those which are either active
   * or inactive.  This is useful for determining which nodes should be
   * written to a data file, and is used by the XDA mesh writing methods.
   */
  void get_not_subactive_node_ids(const MeshBase& mesh, 
                                  std::set<unsigned int>& not_subactive_node_ids);

  /**
   * Count up the number of elements of a specific type (as defined by an iterator range).
   */
   unsigned int n_elem (MeshBase::const_element_iterator& begin,
                        MeshBase::const_element_iterator& end);

  /**
    * Given a mesh and a node in the mesh, the vector will be filled with
    * every node directly attached to the given one.
    */
   void find_nodal_neighbors(const MeshBase& mesh, const Node& n, 
                             std::vector<std::vector<const Elem*> >& nodes_to_elem_map, 
                             std::vector<const Node*>& neighbors);
   
   /**
    * Given a mesh hanging_nodes will be filled with an associative array keyed off the
    * global id of all the hanging nodes in the mesh.  It will hold an array of the 
    * parents of the node (meaning the two nodes to either side of it that make up
    * the side the hanging node is on.
    */
   void find_hanging_nodes_and_parents(const MeshBase& mesh, std::map<unsigned int, std::vector<unsigned int> >& hanging_nodes);

  // There is no reason for users to call functions in the MeshTools::Private namespace.
  namespace Private {
    /**
     * There is no reason for a user to ever call this function.
     *
     * This function determines partition-agnostic global indices for all nodes and elements 
     * in the mesh.  Note that after this function is called the mesh will likely be in an
     * inconsistent state, i.e. \p mesh.nodes(i)->id() != i in the nodes container.
     * Direct node/element access via the \p mesh.node(n) or \p mesh.elem(e) functions will
     * likely fail. The original numbering can (and should) be restored with a subsequent call to
     * \p fix_node_and_element_numbering().
     *
     */
    void globally_renumber_nodes_and_elements (MeshBase &);

  
    /**
     * There is no reason for a user to ever call this function.
     *
     * This function restores a previously broken element/node numbering such that
     * \p mesh.node(n)->id() == n. 
     */
    void fix_broken_node_and_element_numbering (SerialMesh &);
    

    /**
     * There is no reason for a user to ever call this function.
     *
     * This function restores a previously broken element/node numbering such that
     * \p mesh.node(n)->id() == n. 
     */
    void fix_broken_node_and_element_numbering (ParallelMesh &);
  } // end namespace Private
  
} // end namespace MeshTools


#endif // #define __mesh_tools_h__
