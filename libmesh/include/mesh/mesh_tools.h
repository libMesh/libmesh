// $Id: mesh_tools.h,v 1.2 2005-02-22 22:17:33 jwpeterson Exp $

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



#ifndef __mesh_tools_h__
#define __mesh_tools_h__



// C++ Includes   -----------------------------------
#include <vector>

// Local Includes -----------------------------------
#include "libmesh.h"

// forward declarations
class MeshBase;
class Sphere;
class Point;
class Elem;



/**
 * Utility functions for operations on a \p Mesh object.  Here is where
 * useful functions for interfacing with a \p Mesh should be defined.
 * In general this namespace should be used to prevent the \p Mesh class
 * from becoming too cluttered.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision: 1.2 $
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


} // end namespace MeshTools


#endif // #define __mesh_tools_h__
