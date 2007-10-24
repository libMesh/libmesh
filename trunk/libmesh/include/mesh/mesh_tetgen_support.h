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


#ifndef __mesh_tetgen_support_h__
#define __mesh_tetgen_support_h__

#include "libmesh_config.h"
#ifdef HAVE_TETGEN


// C++ includes
#include <vector>
#include <map>
#include <string>

// Local includes

// TetGen include file
#include "tetgen.h"


// Forward Declarations
class UnstructuredMesh;
class Node;

/**
 * Abstract class \p TetGenWrapper provides an interface for basic
 * access to TetGen data structures and methods.
 */
class TetGenWrapper
{
 public:

  /**
   * Constructor.
   */
  TetGenWrapper ();

  /**
   * Destructor.  Empty.
   */
  ~TetGenWrapper ();
  
  /**
   * Method set TetGen commandline switches
   */
  void set_switches(const std::string& s);

  /**
   * Method starts triangulization.
   */
  void run_tetgen();

  /**
   * Method returns number of tetrahedra in TetGen output.
   */
  int  get_numberoftetrahedra();

  /**
   * Method returns number of triangle surface elts. in TetGen output.
   */
  int  get_numberoftrifaces();

  /**
   * Method sets number of nodes in TetGen input.
   */
  void set_numberofpoints(const int i);

  /**
   * Method returns number of nodes in TetGen output.
   */
  int  get_numberofpoints();

  /**
   * Method sets number of facets in TetGen input.
   */
  void set_numberoffacets(const int i);

  /**
   * Method sets number of holes in TetGen input.
   */
  void set_numberofholes(const int i);

  /**
   * Method sets number of regions in TetGen input.
   */
  void set_numberofregions(const int i);

  /**
   * Method allocates memory, sets number of nodes in TetGen input.
   */
  void set_pointlist(const int numofpoints);

  /**
   * Method allocates memory, sets number of facets, holes in TetGen input.
   */
  void set_facetlist(const int numoffacets, const int numofholes);

  /**
   * Method allocates memory, sets number of regions in TetGen input.
   */
  void set_regionlist(const int numofregions);

  /**
   * Method sets coordinates of point i in TetGen input.
   */
  void set_node(const int i, const REAL x, const REAL y, const REAL z);

  /**
   * Method returns coordinates of point i in TetGen output.
   */
  void get_output_node(const int i, REAL& x, REAL& y, REAL& z);

  /**
   * Method returns index of jth node from element i in TetGen output.
   */
  int  get_element_node(const int i, const int j);

  /**
   * Method returns index of jth node from surface triangle i in TetGen output.
   */
  int  get_triface_node(const int i, const int j);

  /**
   * Method returns attribute of element i in TetGen output.
   */
  REAL get_element_attribute(const int i);
       
  /**
   * Method sets coordinates of hole i in TetGen input.
   */
  void set_hole(const int i, const REAL x, const REAL y, const REAL z);

  /**
   * Method sets number of polygons for facet i in TetGen input.
   */
  void set_facet_numberofpolygons(const int i, const int num);

  /**
   * Method sets number of holes for facet i in TetGen input.
   */
  void set_facet_numberofholes(const int i, const int num);

  /**
   * Method allocates memory, sets number of polygons for facet i
   * in TetGen input.
   */
  void set_facet_polygonlist(const int i, const int numofpolygons);

  /**
   * Method sets number of vertices for polygon j, facet i in TetGen input.
   */
  void set_polygon_numberofvertices(const int i, const int j, const int num);

  /**
   * Method allocates memory, sets number of vertices for polygon j,
   * facet i in TetGen input.
   */
  void set_polygon_vertexlist(const int i, const int j, const int numofvertices);

  /**
   * Method sets index of ith facet, jth polygon, kth vertex in
   * TetGen input.
   */
  void set_vertex(const int i, const int j, const int k, const int nodeindex);

  /**
   * Method sets coordinates, attribute and volume constraint for
   * region i in TetGen input.  Note that coordinates and attributes
   * will only be considered if the corresponding switches are
   * enabled.  See TetGen documentation for more details.
   */
  void set_region(const int i, const REAL x, const REAL y, const REAL z,
		  const REAL attribute, const REAL vol_constraint);

  /**
   * TetGen input structure.
   */
  tetgenio   tetgen_data;
  
  /**
   * TetGen output structure.
   */
  tetgenio*  tetgen_output;

  /**
   * TetGen mesh structure (from the TetGen library).
   */
  tetgenmesh      tetgen_mesh;
  
  /**
   * TetGen control class (from the TetGen library).
   */
  tetgenbehavior  tetgen_be; 

};



/**
 * Class \p TetGenMeshInterface provides an interface for
 * tetrahedrization of meshes using the TetGen library.  For
 * information about TetGen cf.
 * <a href="http://tetgen.berlios.de/">TetGen home page</a>.
 */
class TetGenMeshInterface
{
public:

  /**
   * Constructor. Takes a reference to the mesh.
   */
  TetGenMeshInterface (UnstructuredMesh& mesh);

  /**
   * Empty destructor.
   */
  ~TetGenMeshInterface() {};

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. 
   */
  void triangulate_pointset ();

  /** 
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Stores only 2D hull surface elements.
   */
  void pointset_convexhull ();

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Boundary constraints are taken from 
   * elements array.
   */
  void triangulate_conformingDelaunayMesh (const double quality_constraint,
					   const double volume_constraint);

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Boundary constraints are taken from 
   * elements array. Include carve-out functionality.
   */
  void triangulate_conformingDelaunayMesh_carvehole (const std::vector< Node *>& holes,
						     const double quality_constraint,
						     const double volume_constraint);

  /**
   * Help function. Returns _nodes index for *Node. 
   */
  int get_node_index (const Node* inode);


protected:

  /**
   * Local reference to the mesh we are working with.
   */
  UnstructuredMesh& _mesh;

};

#endif // HAVE_TETGEN

#endif 
