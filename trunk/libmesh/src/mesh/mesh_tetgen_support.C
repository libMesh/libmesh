// $Id: mesh_tetgen_support.C,v 1.15 2004-11-12 00:42:42 jwpeterson Exp $
 
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

#include "libmesh_config.h"
#ifdef HAVE_TETGEN



// C++ includes
#include <sstream>


// Local includes
#include "mesh_tetgen_support.h"
#include "mesh_data.h"
#include "cell_tet4.h"
#include "face_tri3.h"
#include "mesh.h"

#ifdef TETGEN_13
#include "tetgen.h"
#endif // TETGEN_13


//----------------------------------------------------------------------
// TetGenMeshInterface class members
TetGenMeshInterface::TetGenMeshInterface (Mesh& mesh) :
  _mesh         (mesh)
{}



TetGenMeshInterface::~TetGenMeshInterface() {}







// =============================================================================
TetGen1_wrapper::TetGen1_wrapper()
{
  tetgen_output = new tetgenio;
}




TetGen1_wrapper::~TetGen1_wrapper()
{
  delete tetgen_output;
}




void TetGen1_wrapper::set_node(int i, REAL x, REAL y, REAL z)
{
  int index = i*3;
  tetgen_data.pointlist[index++] = x;
  tetgen_data.pointlist[index++] = y;
  tetgen_data.pointlist[index++] = z;
}




void TetGen1_wrapper::set_hole(int i, REAL x, REAL y, REAL z)
{
  int index = i*3;
  tetgen_data.holelist[index++] = x;
  tetgen_data.holelist[index++] = y;
  tetgen_data.holelist[index++] = z;
}




void TetGen1_wrapper::set_numberofpoints(int i)
{ tetgen_data.numberofpoints = i; }




void TetGen1_wrapper::get_output_node(int i, REAL& x, REAL& y, REAL& z)
{
  x = tetgen_output->pointlist[3*i];
  y = tetgen_output->pointlist[3*i+1];
  z = tetgen_output->pointlist[3*i+2];
}




int TetGen1_wrapper::get_numberoftetrahedra()
{ return tetgen_output->numberoftetrahedra; }




int TetGen1_wrapper::get_numberoftrifaces()
{ return tetgen_output->numberoftrifaces; }




int TetGen1_wrapper::get_numberofpoints()
{ return tetgen_output->numberofpoints; }




int TetGen1_wrapper::get_element_node(const int i, const int j)
{ return tetgen_output->tetrahedronlist[i*4+j]; }




int TetGen1_wrapper::get_triface_node(const int i, const int j)
{ return tetgen_output->trifacelist[i*3+j]; }




// =============================================================================

#ifdef TETGEN_13

TetGen13_wrapper::TetGen13_wrapper()
{
  this->tetgen_data.mesh_dim                = 3;
  this->tetgen_data.numberofpointattributes = 0;
  this->tetgen_data.firstnumber             = 0;
}




TetGen13_wrapper::~TetGen13_wrapper() {}




void TetGen13_wrapper::set_pointlist(const int numofpoints)
{
  this->set_numberofpoints(numofpoints);
  this->tetgen_data.pointlist = new REAL[tetgen_data.numberofpoints * 3];
}




void TetGen13_wrapper::set_switches(std::string s)
{
  // this allocates memory!!!!!  You must free it at some point!!!!
  switches = strdup(s.c_str()); 
  if (!tetgen_be.parse_commandline(switches)) 
    std::cout << "TetGen replies: Wrong switches!" << std::endl;
}




void TetGen13_wrapper::run_tetgen()
{
  // Call tetrahedralize from the TetGen library.
  tetrahedralize(&tetgen_be, &tetgen_data, tetgen_output);
}




void TetGen13_wrapper::set_numberoffacets(const int i)
{ this->tetgen_data.numberoffacets = i; }




void TetGen13_wrapper::set_numberofholes(const int i)
{ this->tetgen_data.numberofholes = i; }




void TetGen13_wrapper::set_facetlist(int numoffacets, int numofholes)
{
  set_numberoffacets(numoffacets);
  set_numberofholes(numofholes);
  this->tetgen_data.facetlist = new tetgenio::facet[this->tetgen_data.numberoffacets];
  for (int i=0; i<numoffacets; i++)
    this->tetgen_data.init(&(this->tetgen_data.facetlist[i]));
  this->tetgen_data.holelist = new REAL[this->tetgen_data.numberofholes * 3];
}




void TetGen13_wrapper::set_facet_numberofpolygons(int i, int num)
{ this->tetgen_data.facetlist[i].numberofpolygons = num; }




void TetGen13_wrapper::set_facet_numberofholes(int i, int num)
{ this->tetgen_data.facetlist[i].numberofholes = num; }




void TetGen13_wrapper::set_facet_polygonlist(int i, int numofpolygons)
{
  set_facet_numberofpolygons(i, numofpolygons);
  set_facet_numberofholes(i, 0);
  this->tetgen_data.facetlist[i].polygonlist = new tetgenio::polygon[numofpolygons];
  for (int j=0; j<this->tetgen_data.facetlist[i].numberofpolygons; j++)
    this->tetgen_data.init(&(this->tetgen_data.facetlist[i].polygonlist[j]));
}




void TetGen13_wrapper::set_polygon_numberofvertices(int i, int j, int num)
{ this->tetgen_data.facetlist[i].polygonlist[j].numberofvertices = num; }




void TetGen13_wrapper::set_polygon_vertexlist(int i, int j, int numofvertices)
{
  set_polygon_numberofvertices(i, j, numofvertices);
  this->tetgen_data.facetlist[i].polygonlist[j].vertexlist = new int[numofvertices];
}




void TetGen13_wrapper::set_vertex(int i, int j, int k, int nodeindex)
{
  this->tetgen_data.facetlist[i].polygonlist[j].vertexlist[k] = nodeindex;
}




// class type TetGen_access is cast to TetGen_Wrapper class for current TetGen version:
typedef TetGen13_wrapper TetGen_access;











// =============================================================================
void TetGenMeshInterface::triangulate_pointset ()   
{
  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 

  // fill input structure with point set data:
  tetgen_wrapper.set_pointlist( this->_mesh.n_nodes() );

  {
    int index = 0;
    MeshBase::node_iterator it  = this->_mesh.nodes_begin();
    const MeshBase::node_iterator end = this->_mesh.nodes_end();
    for ( ; it != end; ++it) 
      tetgen_wrapper.set_node(index++, (**it)(0), (**it)(1), (**it)(2));
  }
  
  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("-Q"); // TetGen switches: triangulation, Quiet mode
  tetgen_wrapper.run_tetgen();

  // save elements to mesh structure, nodes will not be changed:
  const unsigned int num_elements   = tetgen_wrapper.get_numberoftetrahedra();

  // Vector that temporarily holds the node labels defining element.
  unsigned long int node_labels[4];      

  for (unsigned int i=0; i<num_elements; ++i)
    {
      Elem* elem = new Tet4;

      // Get the nodes associated with this element
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);

      // Associate the nodes with this element
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
	elem->set_node(j) = this->_mesh.node_ptr( node_labels[j] );

      // Finally, add this element to the mesh.
      this->_mesh.add_elem(elem);
    }
} 










void TetGenMeshInterface::pointset_convexhull ()   
{
  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 

  // fill input structure with point set data:
  tetgen_wrapper.set_pointlist(this->_mesh.n_nodes());
  int index = 0;

  MeshBase::node_iterator it        = this->_mesh.nodes_begin();
  const MeshBase::node_iterator end = this->_mesh.nodes_end();
  for ( ; it != end; ++it) 
    tetgen_wrapper.set_node(index++, (**it)(0), (**it)(1), (**it)(2));

  
  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("-Q"); // TetGen switches: triangulation, Convex hull, Quiet mode
  tetgen_wrapper.run_tetgen();
  unsigned int num_elements   = tetgen_wrapper.get_numberoftrifaces();

  // Delete *all* old elements
  {
    MeshBase::element_iterator       it  = this->_mesh.elements_begin();
    const MeshBase::element_iterator end = this->_mesh.elements_end();
    for ( ; it != end; ++it) 
      this->_mesh.delete_elem (*it);
  }
  

  // Add the 2D elements which comprise the convex hull back to the mesh.
  // Vector that temporarily holds the node labels defining element.
  unsigned long int node_labels[3];
  
  for (unsigned int i=0; i<num_elements; ++i)
    {
      Elem* elem = new Tri3;

      // Get node labels associated with this element
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
	node_labels[j] = tetgen_wrapper.get_triface_node(i,j);

      // Associate nodes with this element
      for (unsigned int j=0; j<elem->n_nodes(); j++)
	elem->set_node(j) = this->_mesh.node_ptr( node_labels[j] );

      // Finally, add this element to the mesh.
      this->_mesh.add_elem(elem);
    }
} 











int TetGenMeshInterface::get_node_index (Node* inode)
{
  // This is NO elegant solution! (A linear search of the nodes.)
  unsigned int node_id;
  node_id = inode->id();

  for (unsigned int i=0; i<this->_mesh.n_nodes(); ++i)
    if (this->_mesh.node(i).id() == node_id)
      return i;

  std::cerr << "Error! Node not found in the mesh!" << std::endl;
  error();
  
  return 0;
}











void TetGenMeshInterface::triangulate_conformingDelaunayMesh (const double quality_constraint,
							      const double volume_constraint)   
{
  // >>> assert usage of TRI3 hull elements (no other elements! => mesh.all_tri() )

  // thorough check of hull integrity:
  // loop over hull with breadth-first search:
  std::set< Elem *>    visited;
  std::vector< Elem *> current;

  // Initialize the current vector with element 0
  current.push_back(this->_mesh.elem(0));
  
  while (!current.empty())
    {
      Elem* elem = current.back();

      // Attempt to insert this element into the visited set.
      visited.insert(elem);

      // Remove the last element from the vector.
      current.pop_back();

      // Loop over the element's neighbors
      for (unsigned int i=0; i<elem->n_neighbors(); ++i)
	{
	  // Attempt to find this element's neighbor in the visited set.
	  // If not found, insert this neighbor into the current vector.
	  if (visited.find(elem->neighbor(i)) == visited.end())
	    current.push_back(elem->neighbor(i));
	} 
    }
  
  if (visited.size() != this->_mesh.n_elem())
    {
      std::cerr << "triangulate: hull not connected: element(s) not reached by others.\n";
      error();
    } 

  // start triangulation method with empty holes list:
  std::vector< Node *> noholes;
  triangulate_conformingDelaunayMesh_carvehole(noholes, quality_constraint, volume_constraint);
}













void TetGenMeshInterface::triangulate_conformingDelaunayMesh_carvehole  (const std::vector< Node *>& holes,
									 const double quality_constraint,
									 const double volume_constraint)
{
  // Check mesh for validity: Must be composed of all TRI3 elements,
  // also it must be convex so each element must have non-NULL neighbors.
  {
    MeshBase::element_iterator it        = this->_mesh.elements_begin();
    const MeshBase::element_iterator end = this->_mesh.elements_end();

    for (; it != end ; ++it)
      {
	Elem* elem = *it;

	// Check for proper element type
	if (elem->type() != TRI3)
	  {
	    std::cerr << "ERROR: Some of the elements in the original mesh were not TRI3!" << std::endl;
	    error();
	  }

	for (unsigned int i=0; i<elem->n_neighbors(); ++i)
	  {
	    if (elem->neighbor(i) == NULL)
	      {
		std::cerr << "ERROR: Non-convex hull, cannot be tetrahedralized." << std::endl;
		error();
	      }
	  }
      }
  }
  
  // >>> fill input structure with point set data:
  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 

  tetgen_wrapper.set_pointlist(this->_mesh.n_nodes());
  {
    int index = 0;
    MeshBase::node_iterator it        = this->_mesh.nodes_begin();
    const MeshBase::node_iterator end = this->_mesh.nodes_end();
    for ( ; it != end; ++it) 
      tetgen_wrapper.set_node(index++, (**it)(0), (**it)(1), (**it)(2));
  }

  // >>> fill input structure "tetgenio" with facet data:
  int facet_num = this->_mesh.n_elem();

  // allocate memory in "tetgenio" structure:
  tetgen_wrapper.set_facetlist(facet_num, holes.size());


  {
    int insertnum = 0;
    MeshBase::element_iterator it        = this->_mesh.elements_begin();
    const MeshBase::element_iterator end = this->_mesh.elements_end();
    for (; it != end ; ++it)
      {
	tetgen_wrapper.set_facet_polygonlist(insertnum, 1);
	tetgen_wrapper.set_polygon_vertexlist(insertnum, 0, 3);

	Elem* elem = *it;

	for (unsigned int j=0; j<elem->n_nodes(); ++j)
	  tetgen_wrapper.set_vertex(insertnum, // facet number
				    0,         // polygon (always 0)
				    j,         // vertex in tetgen input
				    this->get_node_index(elem->get_node(j)));
	insertnum++;
      }
  }



  // fill hole list (if there are holes):
  if (holes.size() > 0)
    {
      std::vector< Node *>::const_iterator ihole;
      int hole_index = 0;
      for (ihole=holes.begin(); ihole!=holes.end(); ++ihole)
	tetgen_wrapper.set_hole(hole_index++, (**ihole)(0), (**ihole)(1), (**ihole)(2));
    }

  
  // >>> run TetGen triangulation method:
  // assemble switches:
  std::ostringstream oss; // string holding switches
  oss << "-pQ";

  if (quality_constraint != 0)
    oss << "q" << quality_constraint;
  
  if (volume_constraint != 0)
    oss << "a" << volume_constraint;
  
  std::string params = oss.str();

  tetgen_wrapper.set_switches(params); // TetGen switches: Piecewise linear complex, Quiet mode
  tetgen_wrapper.run_tetgen();

  // => nodes:
  unsigned int old_nodesnum = this->_mesh.n_nodes();
  REAL x=0., y=0., z=0.;
  const unsigned int num_nodes = tetgen_wrapper.get_numberofpoints();

  // Add additional nodes to the nodes vector.
  for (unsigned int i=old_nodesnum; i<=num_nodes; i++)
    {
      tetgen_wrapper.get_output_node(i, x,y,z);
      this->_mesh.add_point ( Point(x,y,z) );
    }
  
  // => tetrahedra:
  const unsigned int num_elements   = tetgen_wrapper.get_numberoftetrahedra();

  // Vector that temporarily holds the node labels defining element.
  unsigned long int node_labels[4];      

  for (unsigned int i=0; i<num_elements; i++)
    {
      //_elements[firstnumber+i] = new Tet4;     // TetGen only supports Tet4 elements.
      Elem* elem = new Tet4;

      // Fill up the the node_labels vector
      for (unsigned int j=0; j<elem->n_nodes(); j++) 
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);
      
      // Associate nodes with this element
      for (unsigned int j=0; j<elem->n_nodes(); j++) 
	elem->set_node(j) = this->_mesh.node_ptr( node_labels[j] );

      // Finally, add this element to the mesh
      this->_mesh.add_elem(elem);
    }
}

#endif

#endif // #ifdef HAVE_TETGEN
