// $Id: mesh_tetgen_support.C,v 1.11 2004-10-28 19:09:27 benkirk Exp $
 
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
TetGenMeshInterface::TetGenMeshInterface (Mesh& mesh, MeshData& data) :
  _nodes        (mesh._nodes),
  _elements     (mesh._elements),
  _num_nodes    (0),
  _num_elements (0),
  _mesh_data    (data)
{
}



TetGenMeshInterface::~TetGenMeshInterface()
{
}


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

int TetGen1_wrapper::get_element_node(int i, int j)
{ return tetgen_output->tetrahedronlist[i*4+j]; }

int TetGen1_wrapper::get_triface_node(int i, int j)
{ return tetgen_output->trifacelist[i*3+j]; }


// =============================================================================

#ifdef TETGEN_13

TetGen13_wrapper::TetGen13_wrapper()
{
  tetgen_data.mesh_dim                = 3;
  tetgen_data.numberofpointattributes = 0;
  tetgen_data.firstnumber             = 0;
}

TetGen13_wrapper::~TetGen13_wrapper()
{ }
  
void TetGen13_wrapper::set_pointlist(int numofpoints)
{
  set_numberofpoints(numofpoints);
  tetgen_data.pointlist = new REAL[tetgen_data.numberofpoints * 3];
}
  
void TetGen13_wrapper::set_switches(std::string s)
{ 
  switches = strdup(s.c_str()); 
  if (!tetgen_be.parse_commandline(switches)) {
    std::cout << "TetGen replies: Wrong switches!" << std::endl;
  } // if
}

void TetGen13_wrapper::run_tetgen()
{
  tetrahedralize(&tetgen_be, &tetgen_data, tetgen_output);
}

void TetGen13_wrapper::set_numberoffacets(int i)
{ tetgen_data.numberoffacets = i; }

void TetGen13_wrapper::set_numberofholes(int i)
{ tetgen_data.numberofholes = i; }

void TetGen13_wrapper::set_facetlist(int numoffacets, int numofholes)
{
  set_numberoffacets(numoffacets);
  set_numberofholes(numofholes);
  tetgen_data.facetlist = new tetgenio::facet[tetgen_data.numberoffacets];
  for (int i=0; i<numoffacets; i++)
    tetgen_data.init(&(tetgen_data.facetlist[i]));
  tetgen_data.holelist = new REAL[tetgen_data.numberofholes * 3];
}

void TetGen13_wrapper::set_facet_numberofpolygons(int i, int num)
{ tetgen_data.facetlist[i].numberofpolygons = num; }

void TetGen13_wrapper::set_facet_numberofholes(int i, int num)
{ tetgen_data.facetlist[i].numberofholes = num; }

void TetGen13_wrapper::set_facet_polygonlist(int i, int numofpolygons)
{
  set_facet_numberofpolygons(i, numofpolygons);
  set_facet_numberofholes(i, 0);
  tetgen_data.facetlist[i].polygonlist = new tetgenio::polygon[numofpolygons];
  for (int j=0; j<tetgen_data.facetlist[i].numberofpolygons; j++)
    tetgen_data.init(&(tetgen_data.facetlist[i].polygonlist[j]));
}

void TetGen13_wrapper::set_polygon_numberofvertices(int i, int j, int num)
{ tetgen_data.facetlist[i].polygonlist[j].numberofvertices = num; }

void TetGen13_wrapper::set_polygon_vertexlist(int i, int j, int numofvertices)
{
  set_polygon_numberofvertices(i, j, numofvertices);
  tetgen_data.facetlist[i].polygonlist[j].vertexlist = new int[numofvertices];
}

void TetGen13_wrapper::set_vertex(int i, int j, int k, int nodeindex)
{
  tetgen_data.facetlist[i].polygonlist[j].vertexlist[k] = nodeindex;
}

// class type TetGen_access is cast to TetGen_Wrapper class for current TetGen version:
typedef TetGen13_wrapper TetGen_access;


// =============================================================================


void TetGenMeshInterface::triangulate_pointset ()   
{
  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 

  _num_elements = 0; // counts additional elements 

  // fill input structure with point set data:
  tetgen_wrapper.set_pointlist(_nodes.size());
  int index = 0;
  std::vector< Node *>::iterator i;
  for (i=_nodes.begin(); i!=_nodes.end(); ++i) {
    tetgen_wrapper.set_node(index++, (**i)(0), (**i)(1), (**i)(2));
  } 

  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("-Q"); // TetGen switches: triangulation, Quiet mode
  tetgen_wrapper.run_tetgen();

  // save elements to mesh structure, nodes will not be changed:
  int n_nodes     = 4;                                  // (Tet4 elements)
  _num_elements   = tetgen_wrapper.get_numberoftetrahedra();
  int firstnumber = _elements.size()+1;                 // append position
  // Reserve space in the appropriate vector to avoid unnecessary allocations.
  _elements.resize (_num_elements+firstnumber);

  // Vector that assigns element nodes to their correct position:
  static const unsigned int assign_elm_nodes[] = { 0, 1, 2, 3};
  unsigned long int node_labels[4];      // Vector that temporarily holds the node labels defining element.
  for (unsigned int i=0; i<_num_elements; i++)
    {
      _elements[firstnumber+i] = new Tet4;     // TetGen only supports Tet4 elements.
      for (int j=0; j<n_nodes; j++)
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);
      // nodes are being stored in element
      for (int j=0; j<n_nodes; j++)
	_elements[firstnumber+i]->set_node(assign_elm_nodes[j]) = _nodes[node_labels[j]];
    } // for
} // triangulate_pointset


void TetGenMeshInterface::pointset_convexhull ()   
{
  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 

  _num_elements = 0; // counts additional elements 

  // fill input structure with point set data:
  tetgen_wrapper.set_pointlist(_nodes.size());
  int index = 0;
  std::vector< Node *>::iterator i;
  for (i=_nodes.begin(); i!=_nodes.end(); ++i) {
    tetgen_wrapper.set_node(index++, (**i)(0), (**i)(1), (**i)(2));
  } 

  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("-Q"); // TetGen switches: triangulation, Convex hull, Quiet mode
  tetgen_wrapper.run_tetgen();

  // save SURFACE elements to mesh structure, nodes will not be changed:
  int n_nodes     = 3;                                  // (Tri3 elements)
  _num_elements   = tetgen_wrapper.get_numberoftrifaces();
  //  int firstnumber = _elements.size()+1;                 // append position
  int firstnumber = 0;                                      // trivial append position

  // Delete old elements:
  std::vector< Elem *>::iterator j;
  for (j=_elements.begin(); j!=_elements.end(); ++j) {
    delete (*j);
  } 

  // Reserve space in the appropriate vector to avoid unnecessary allocations.
  _elements.resize (_num_elements+firstnumber);

  // Vector that assigns element nodes to their correct position:
  static const unsigned int assign_elm_nodes[] = { 0, 1, 2};
  unsigned long int node_labels[3];      // Vector that temporarily holds the node labels defining element.
  for (unsigned int i=0; i<_num_elements; i++)
    {
      _elements[firstnumber+i] = new Tri3;
      for (int j=0; j<n_nodes; j++)
	node_labels[j] = tetgen_wrapper.get_triface_node(i,j);
      // nodes are being stored in element
      for (int j=0; j<n_nodes; j++)
	_elements[firstnumber+i]->set_node(assign_elm_nodes[j]) = _nodes[node_labels[j]];
    } // for
} // pointset_convexhull




int TetGenMeshInterface::get_node_index (Node* inode)
{
  // This is NO elegant solution!
  unsigned int node_id;
  int result = 0;
  node_id = inode->id();

  for (unsigned int i=0; i<_nodes.size(); i++)
    if (_nodes[i]->id()==node_id) result=i;
  return result;
}


void TetGenMeshInterface::triangulate_conformingDelaunayMesh (double quality_constraint, double volume_constraint)   
{
  // >>> assert usage of TRI3 hull elements (no other elements! => mesh.all_tri() )

  // thorough check of hull integrity:
  // loop over hull with breadth-first search:
  std::set< Elem *>    visited;
  std::vector< Elem *> current;
  current.push_back(_elements[0]);
  while (!current.empty()) {
    Elem* actual = current.back();
    visited.insert(actual);
    current.pop_back();
    for (unsigned int i=0; i<actual->n_neighbors(); ++i) {
      int new_n = 0; // number of non-visited neighbours
      if (visited.find(actual->neighbor(i))==visited.end()) {
	new_n++;
	current.push_back(actual->neighbor(i));
      } // if
    } // for
  } // while
  if (visited.size() != _elements.size()) {
    std::cout << "triangulate: hull not connected: element(s) not reached by others.\n";
    error();
  } // if

  // start triangulation method with empty holes list:
  std::vector< Node *> noholes;
  triangulate_conformingDelaunayMesh_carvehole(noholes, quality_constraint, volume_constraint);
} // triangulate_conformingDelaunayMesh


void TetGenMeshInterface::triangulate_conformingDelaunayMesh_carvehole
      (std::vector< Node *>& holes, double quality_constraint, double volume_constraint)
{
  // >>> assert usage of TRI3 hull elements (no other elements! => mesh.all_tri() )
  // >>> intersection of inner and outer hull structures, connectedness etc. is not verified!

  // count number of facet elements:
  // (by the way check number of element neighbours - should be 3!)
  unsigned int tri3_num  = 0;
  unsigned int min_nb    = 4;
  unsigned int total_num = _elements.size();  
  for (unsigned int i=0; i<total_num; i++) {
    if (_elements[i]->type() == TRI3) { 
      tri3_num++;
      if (_elements[i]->n_neighbors() < min_nb) 
	min_nb = _elements[i]->n_neighbors();
    } // if
  } // for
  if (total_num > tri3_num) error();
  if (min_nb < 3) { 
    // hull cannot be closed
    std::cout << "triangulate_carvehole: hull not connected: element has " << min_nb << " neighbours.";
    error();
  }
  
  // >>> fill input structure with point set data:

  // class tetgen_wrapper allows library access on a basic level:
  TetGen_access tetgen_wrapper; 
  _num_elements = 0; // counts additional elements 

  tetgen_wrapper.set_pointlist(_nodes.size());
  int index = 0;
  std::vector< Node *>::iterator i;
  for (i=_nodes.begin(); i!=_nodes.end(); ++i) {
    tetgen_wrapper.set_node(index++, (**i)(0), (**i)(1), (**i)(2));
  } 

  // >>> fill input structure "tetgenio" with facet data:
  int facet_num = tri3_num; 

  // allocate memory in "tetgenio" structure:
  tetgen_wrapper.set_facetlist(facet_num, holes.size());

  // fill facetlist: each facet consists of only one polygon:
  int insertnum = 0;
  for (unsigned int i=0; i<total_num; i++) {
    tetgen_wrapper.set_facet_polygonlist(insertnum, 1);
    tetgen_wrapper.set_polygon_vertexlist(insertnum, 0, 3);
    for (int j=0; j<3; j++)
      tetgen_wrapper.set_vertex(insertnum,0,j, get_node_index(_elements[i]->get_node(j)));
    insertnum++;
  } // for

  // fill hole list (if there are holes):
  if (holes.size()>0) {
    std::vector< Node *>::iterator ihole;
    int hole_index = 0;
    for (ihole=holes.begin(); ihole!=holes.end(); ++ihole) {
      tetgen_wrapper.set_hole(hole_index++, (**ihole)(0), (**ihole)(1), (**ihole)(2));
    } // for
  } // if

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

  // >>> save elements to mesh structure, 
  // old nodes should not be changed, new nodes will be appended:

  // => nodes:
  // Reserve space in the _nodes vector to avoid unnecessary allocations.
  unsigned int old_nodesnum = _nodes.size();
  REAL x,y,z;
  _num_nodes       = tetgen_wrapper.get_numberofpoints();
  _nodes.resize (_num_nodes+1);
  // Add additional nodes to the nodes vector.
  for (unsigned int i=old_nodesnum; i<=_num_nodes; i++) {
    tetgen_wrapper.get_output_node(i, x,y,z);
    _nodes[i] = Node::build(x,y,z,i);
  }
  // => tetrahedra:
  int n_nodes     = 4;                                  // (Tet4 elements)
  _num_elements   = tetgen_wrapper.get_numberoftetrahedra();

  int firstnumber = _elements.size();                 // append position
  // Reserve space in the _elements vector to avoid unnecessary allocations.
  _elements.resize (_num_elements+firstnumber);

  // Vector that assigns element nodes to their correct position:
  static const unsigned int assign_elm_nodes[] = { 0, 1, 2, 3};
  unsigned long int node_labels[4];      // Vector that temporarily holds the node labels defining element.
  for (unsigned int i=0; i<_num_elements; i++)
    {
      _elements[firstnumber+i] = new Tet4;     // TetGen only supports Tet4 elements.
      for (int j=0; j<n_nodes; j++) 
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);
      
      // nodes are being stored in element
      for (int j=0; j<n_nodes; j++) 
	_elements[firstnumber+i]->set_node(assign_elm_nodes[j]) = _nodes[node_labels[j]];
    } // for
} // triangulate_conformingDelaunayMesh_carvehole

#endif

#endif
