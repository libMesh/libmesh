// $Id: mesh_tetgen_support.h,v 1.12 2004-11-08 00:11:03 jwpeterson Exp $
 
// The libMesh Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh.h"

#ifdef TETGEN_13
#include "tetgen.h"
#endif // TETGEN_13


// Forward Declarations


/**
 * Abstract class \p TetGenWrapper provides an interface for basic
 * access to TetGen data structures and methods.  This class is
 * intended to simplify adaption of libmesh to new TetGen versions
 * (e.g. change of class nomenclature between tetgen1.2c and
 * tetgen1.3).
 */
class TetGenWrapper
{
 public:
  /**
   * Method set TetGen commandline switches
   */
  virtual void set_switches(std::string s) = 0;

  /**
   * Method starts triangulization.
   */
  virtual void run_tetgen() = 0;

  /**
   * Method returns number of tetrahedra in TetGen output.
   */
  virtual int  get_numberoftetrahedra() = 0;

  /**
   * Method returns number of triangle surface elts. in TetGen output.
   */
  virtual int  get_numberoftrifaces() = 0;

  /**
   * Method sets number of nodes in TetGen input.
   */
  virtual void set_numberofpoints(int i) = 0;

  /**
   * Method returns number of nodes in TetGen output.
   */
  virtual int  get_numberofpoints() = 0;

  /**
   * Method sets number of facets in TetGen input.
   */
  virtual void set_numberoffacets(int i) = 0;

  /**
   * Method sets number of holes in TetGen input.
   */
  virtual void set_numberofholes(int i) = 0;

  /**
   * Method allocates memory, sets number of nodes in TetGen input.
   */
  virtual void set_pointlist(int numofpoints) = 0;

  /**
   * Method allocates memory, sets number of facets, holes in TetGen input.
   */
  virtual void set_facetlist(int numoffacets, int numofholes) = 0;

  /**
   * Method sets coordinates of point i in TetGen input.
   */
  virtual void set_node(int i, REAL x, REAL y, REAL z) = 0;

  /**
   * Method returns coordinates of point i in TetGen output.
   */
  virtual void get_output_node(int i, REAL& x, REAL& y, REAL& z) = 0;

  /**
   * Method returns index of jth node from element i in TetGen output.
   */
  virtual int  get_element_node(int i, int j) = 0;

  /**
   * Method returns index of jth node from surface triangle i in TetGen output.
   */
  virtual int  get_triface_node(int i, int j) = 0;

  /**
   * Method sets coordinates of hole i in TetGen input. */
  virtual void set_hole(int i, REAL x, REAL y, REAL z) = 0;

  /**
   * Method sets number of polygons for facet i in TetGen input.
   */
  virtual void set_facet_numberofpolygons(int i, int num) = 0;

  /**
   * Method sets number of holes for facet i in TetGen input.
   */
  virtual void set_facet_numberofholes(int i, int num) = 0;

  /**
   * Method allocates memory, sets number of polygons for facet i
   * in TetGen input.
   */
  virtual void set_facet_polygonlist(int i, int numofpolygons) = 0;

  /**
   * Method sets number of vertices for polygon j, facet i in TetGen input.
   */
  virtual void set_polygon_numberofvertices(int i, int j, int num) = 0;

  /**
   * Method allocates memory, sets number of vertices for polygon j,
   * facet i in TetGen input.
   */
  virtual void set_polygon_vertexlist(int i, int j, int numofvertices) = 0;

  /**
   * Method sets index of ith facet, jth polygon, kth vertex in
   * TetGen input.
   */
  virtual void set_vertex(int i, int j, int k, int nodeindex) = 0;
};


/**
 * Semi-abstract implementation of TetGenWrapper class 
 * for TetGen 1.x data structures and methods.
 */
class TetGen1_wrapper : public TetGenWrapper
{
 public:
  
  /**
   * TetGen input structure.
   */
  tetgenio   tetgen_data;
  
  /**
   * TetGen output structure.
   */
  tetgenio*  tetgen_output;

  /**
   * Constructor class
   */
  TetGen1_wrapper();

  /**
   * Destructor class
   */
  virtual ~TetGen1_wrapper();

  void set_node(int i, REAL x, REAL y, REAL z);
  void set_hole(int i, REAL x, REAL y, REAL z);
  void set_numberofpoints(int i);
  int  get_numberofpoints();

  void get_output_node(int i, REAL& x, REAL& y, REAL& z);
  int  get_numberoftetrahedra();
  int  get_numberoftrifaces();
  int  get_element_node(int i, int j);
  int  get_triface_node(int i, int j);
};


#ifdef TETGEN_13

/**
 * Implementation of TetGenWrapper class provides an interface
 * for basic access to TetGen 1.3 data structures and methods.
 */
class TetGen13_wrapper : public TetGen1_wrapper
{
 public:
  /**
   * TetGen mesh structure.
   */
  tetgenmesh      tetgen_mesh;
  
  /**
   * TetGen control class.
   */
  tetgenbehavior  tetgen_be; 

  /**
   * String holding commandline switches for TetGen.
   */
  char*           switches;

  /**
   * Constructor class
   */
  TetGen13_wrapper();

  /**
   * Destructor class
   */
  ~TetGen13_wrapper();

  void set_switches(std::string s);
  void run_tetgen();
  void set_pointlist(int numofpoints);

  void set_numberoffacets(int i);
  void set_numberofholes(int i);
  void set_facetlist(int numoffacets, int numofholes);

  void set_facet_numberofpolygons(int i, int num);
  void set_facet_numberofholes(int i, int num);
  void set_facet_polygonlist(int i, int numofpolygons);
  void set_polygon_numberofvertices(int i, int j, int num);
  void set_polygon_vertexlist(int i, int j, int numofvertices);
  void set_vertex(int i, int j, int k, int nodeindex);
};

#endif // TETGEN_13


/**
 * Class \p TetGenMeshInterface provides an interface for
 * tetrahedrization of meshes using the TetGen library.  For
 * information about TetGen see: cf.
 * <a href="http://tetgen.berlios.de/">TetGen home page</a>.
 */
class TetGenMeshInterface
{
public:

  /**
   * Constructor. Takes a reference to the mesh.
   */
  TetGenMeshInterface (Mesh& mesh);


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
  void triangulate_conformingDelaunayMesh (double quality_constraint, double volume_constraint);

  /**
   * Method invokes TetGen library to compute a Delaunay tetrahedrization
   * from the nodes point set. Boundary constraints are taken from 
   * elements array. Include carve-out functionality.
   */
  void triangulate_conformingDelaunayMesh_carvehole
         (std::vector< Node *>& holes, double quality_constraint, double volume_constraint);

  /**
   * Help function. Returns _nodes index for *Node. 
   */
  int get_node_index (Node* inode);

  /**
   * Destructor.
   */
  ~TetGenMeshInterface();

protected:

  //-------------------------------------------------------------
  // local data
  
  /**
   * vector holding the nodes.
   */
  std::vector<Node*>& _nodes;    

  /**
   * vector holding the elements.
   */
  std::vector<Elem*>& _elements; 

  /**
   * stores new positions of nodes. Used when reading.
   */
  std::map<unsigned int,unsigned int> _assign_nodes; 

  /**
   * total number of nodes. Primarily used when reading.
   */
  unsigned int _num_nodes;

  /**
   * total number of elements. Primarily used when reading.
   */
  unsigned int _num_elements;

};

#endif // HAVE_TETGEN

#endif 
