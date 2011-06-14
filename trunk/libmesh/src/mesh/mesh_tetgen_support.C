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

#include "libmesh_config.h"
#ifdef LIBMESH_HAVE_TETGEN


// C++ includes
#include <sstream>
#include <set>

// Local includes
#include "cell_tet4.h"
#include "face_tri3.h"
#include "unstructured_mesh.h"
#include "mesh_tetgen_support.h"
#include "utility.h" // binary_find

namespace libMesh
{


//----------------------------------------------------------------------
// TetGenMeshInterface functions
TetGenWrapper::TetGenWrapper()
{
  tetgen_output = new tetgenio;

  this->tetgen_data.mesh_dim                = 3;
  this->tetgen_data.numberofpointattributes = 0;
  this->tetgen_data.firstnumber             = 0;
}



TetGenWrapper::~TetGenWrapper()
{
  delete tetgen_output;
}



void TetGenWrapper::set_node(unsigned i, REAL x, REAL y, REAL z)
{
  unsigned index = i*3;
  tetgen_data.pointlist[index++] = x;
  tetgen_data.pointlist[index++] = y;
  tetgen_data.pointlist[index++] = z;
}



void TetGenWrapper::set_hole(unsigned i, REAL x, REAL y, REAL z)
{
  unsigned index = i*3;
  tetgen_data.holelist[index++] = x;
  tetgen_data.holelist[index++] = y;
  tetgen_data.holelist[index++] = z;
}



void TetGenWrapper::set_numberofpoints(int i)
{
  // This is an int in tetgen, so use an int here even though it should be unsigned
  tetgen_data.numberofpoints = i;
}



void TetGenWrapper::get_output_node(unsigned i, REAL& x, REAL& y, REAL& z)
{
  // Bounds checking...
  if (i >= static_cast<unsigned>(tetgen_output->numberofpoints))
    {
      std::cerr << "Error, requested point " 
		<< i 
		<< ", but there are only " 
		<< tetgen_output->numberofpoints
		<< " points available."
		<< std::endl;
      libmesh_error();
    }

  x = tetgen_output->pointlist[3*i];
  y = tetgen_output->pointlist[3*i+1];
  z = tetgen_output->pointlist[3*i+2];
}



int TetGenWrapper::get_numberoftetrahedra()
{
  return tetgen_output->numberoftetrahedra;
}



int TetGenWrapper::get_numberoftrifaces()
{
  return tetgen_output->numberoftrifaces;
}



int TetGenWrapper::get_numberofpoints()
{
  return tetgen_output->numberofpoints;
}



int TetGenWrapper::get_element_node(unsigned i, unsigned j)
{
  return tetgen_output->tetrahedronlist[i*4+j];
}



int TetGenWrapper::get_triface_node(unsigned i, unsigned j)
{
  return tetgen_output->trifacelist[i*3+j];
}



REAL TetGenWrapper::get_element_attribute(unsigned i)
{
  libmesh_assert(tetgen_output->numberoftetrahedronattributes>0);
  return tetgen_output->tetrahedronattributelist[tetgen_output->numberoftetrahedronattributes*i];
}



void TetGenWrapper::allocate_pointlist(int numofpoints)
{
  // This is stored as an int in tetgen, so we store it that way as well.
  this->set_numberofpoints(numofpoints);

  // Is there previously-allocated memory here?
  if (this->tetgen_data.pointlist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }

  // We allocate memory here, the tetgenio destructor will delete it.
  this->tetgen_data.pointlist = new REAL[tetgen_data.numberofpoints * 3];
}



void TetGenWrapper::set_switches(const std::string& s)
{
  // Copy the string to a temporary buffer for passing to the C API
  char buffer[256];
  libmesh_assert (s.size() < sizeof(buffer)-1);
  buffer[ s.copy( buffer , sizeof( buffer ) - 1 ) ] = '\0' ;
  
  if (!tetgen_be.parse_commandline(buffer)) 
    libMesh::out << "TetGen replies: Wrong switches!" << std::endl;
}



void TetGenWrapper::run_tetgen()
{
  // Call tetrahedralize from the TetGen library.
  tetrahedralize(&tetgen_be, &tetgen_data, tetgen_output);
}



void TetGenWrapper::set_numberoffacets(int i)
{ 
  // This is stored as an int in TetGen
  this->tetgen_data.numberoffacets = i;
}



void TetGenWrapper::set_numberofholes(int i)
{ 
  // This is stored as an int in TetGen
  this->tetgen_data.numberofholes = i;
}



void TetGenWrapper::set_numberofregions(int i)
{
  // This is stored as an int in TetGen
  this->tetgen_data.numberofregions = i;
}



void TetGenWrapper::allocate_facetlist(int numoffacets, int numofholes)
{
  // These are both stored as ints in TetGen
  this->set_numberoffacets(numoffacets);
  this->set_numberofholes(numofholes);
  
  // Is there previously-allocated memory here?
  if (this->tetgen_data.facetlist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }

  // We allocate memory here, the tetgenio destructor cleans it up.
  this->tetgen_data.facetlist = new tetgenio::facet[this->tetgen_data.numberoffacets];

  for (int i=0; i<numoffacets; i++)
    this->tetgen_data.init(&(this->tetgen_data.facetlist[i]));


  // Is there previously-allocated memory here?
  if (this->tetgen_data.holelist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }
  
  this->tetgen_data.holelist = new REAL[this->tetgen_data.numberofholes * 3];
}



void TetGenWrapper::allocate_regionlist(int numofregions)
{
  this->set_numberofregions(numofregions);

  // Is there previously-allocated memory here?
  if (this->tetgen_data.regionlist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }
  
  // We allocate memory here, the tetgenio destructor cleans it up.
  this->tetgen_data.regionlist = new REAL[this->tetgen_data.numberofregions * 5];
}



void TetGenWrapper::set_facet_numberofpolygons(unsigned i, int num)
{
  // numberofpolygons is stored as an int in TetGen
  this->tetgen_data.facetlist[i].numberofpolygons = num;
}



void TetGenWrapper::set_facet_numberofholes(unsigned i, int num)
{
  // numberofholes is stored as an int in TetGen
  this->tetgen_data.facetlist[i].numberofholes = num;
}




void TetGenWrapper::allocate_facet_polygonlist(unsigned i, int numofpolygons)
{
  this->set_facet_numberofpolygons(i, numofpolygons);
  this->set_facet_numberofholes(i, 0);

  // Is there previously-allocated memory here?
  if (this->tetgen_data.facetlist[i].polygonlist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }

  // We allocate memory here, the tetgenio destructor cleans it up.
  this->tetgen_data.facetlist[i].polygonlist = new tetgenio::polygon[numofpolygons];

  for (int j=0; j<this->tetgen_data.facetlist[i].numberofpolygons; j++)
    this->tetgen_data.init(&(this->tetgen_data.facetlist[i].polygonlist[j]));
}



void TetGenWrapper::set_polygon_numberofvertices(unsigned i, unsigned j, int num)
{
  // numberofvertices is stored as an int in TetGen
  this->tetgen_data.facetlist[i].polygonlist[j].numberofvertices = num;
}



void TetGenWrapper::allocate_polygon_vertexlist(unsigned i, unsigned j, int numofvertices)
{
  this->set_polygon_numberofvertices(i, j, numofvertices);

  // Is there previously-allocated memory here?
  if (this->tetgen_data.facetlist[i].polygonlist[j].vertexlist != NULL)
    {
      libMesh::err << "Cannot allocate on top of previously allocated memory!" << std::endl;
      libmesh_error();
    }


  // We allocate memory here, the tetgenio destructor cleans it up.
  this->tetgen_data.facetlist[i].polygonlist[j].vertexlist = new int[numofvertices];
}




void TetGenWrapper::set_vertex(unsigned i, unsigned j, unsigned k, int nodeindex)
{
  // vertexlist entries are stored as ints in TetGen
  this->tetgen_data.facetlist[i].polygonlist[j].vertexlist[k] = nodeindex;
}



void TetGenWrapper::set_region(unsigned i, REAL x, REAL y, REAL z,
			       REAL attribute, REAL vol_constraint)
{
  unsigned index = i*5;
  tetgen_data.regionlist[index++] = x;
  tetgen_data.regionlist[index++] = y;
  tetgen_data.regionlist[index++] = z;
  tetgen_data.regionlist[index++] = attribute;
  tetgen_data.regionlist[index++] = vol_constraint;
}








//----------------------------------------------------------------------
// TetGenMeshInterface class members
TetGenMeshInterface::TetGenMeshInterface (UnstructuredMesh& mesh) :
  _mesh(mesh)
{
}



void TetGenMeshInterface::triangulate_pointset ()   
{
  // class tetgen_wrapper allows library access on a basic level:
  TetGenWrapper tetgen_wrapper; 

  // fill input structure with point set data:
  this->fill_pointlist(tetgen_wrapper);

  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("Q"); // TetGen switches: triangulation, Quiet mode
  tetgen_wrapper.run_tetgen();

  // save elements to mesh structure, nodes will not be changed:
  const unsigned int num_elements   = tetgen_wrapper.get_numberoftetrahedra();

  // Vector that temporarily holds the node labels defining element.
  unsigned int node_labels[4];      

  for (unsigned int i=0; i<num_elements; ++i)
    {
      Elem* elem = new Tet4;

      // Get the nodes associated with this element
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);

      // Associate the nodes with this element
      this->assign_nodes_to_elem(node_labels, elem);

      // Finally, add this element to the mesh.
      this->_mesh.add_elem(elem);
    }
} 



void TetGenMeshInterface::pointset_convexhull ()   
{
  // class tetgen_wrapper allows library access on a basic level
  TetGenWrapper tetgen_wrapper; 

  // Copy Mesh's node points into TetGen data structure
  this->fill_pointlist(tetgen_wrapper);

  // run TetGen triangulation method:
  tetgen_wrapper.set_switches("Q"); // TetGen switches: triangulation, Convex hull, Quiet mode
  tetgen_wrapper.run_tetgen();
  unsigned int num_elements   = tetgen_wrapper.get_numberoftrifaces();

  // Delete *all* old elements.  Yes, we legally delete elements while
  // iterating over them because no entries from the underlying container
  // are actually erased.
  {
    MeshBase::element_iterator       it  = this->_mesh.elements_begin();
    const MeshBase::element_iterator end = this->_mesh.elements_end();
    for ( ; it != end; ++it) 
      this->_mesh.delete_elem (*it);
  }
  

  // Add the 2D elements which comprise the convex hull back to the mesh.
  // Vector that temporarily holds the node labels defining element.
  unsigned int node_labels[3];
  
  for (unsigned int i=0; i<num_elements; ++i)
    {
      Elem* elem = new Tri3;

      // Get node labels associated with this element
      for (unsigned int j=0; j<elem->n_nodes(); ++j)
	node_labels[j] = tetgen_wrapper.get_triface_node(i,j);

      this->assign_nodes_to_elem(node_labels, elem);

      // Finally, add this element to the mesh.
      this->_mesh.add_elem(elem);
    }
} 





void TetGenMeshInterface::triangulate_conformingDelaunayMesh (double quality_constraint,
							      double volume_constraint)   
{
  // start triangulation method with empty holes list:
  std::vector< Node *> noholes;
  triangulate_conformingDelaunayMesh_carvehole(noholes, quality_constraint, volume_constraint);
}



void TetGenMeshInterface::triangulate_conformingDelaunayMesh_carvehole  (const std::vector< Node *>& holes,
									 double quality_constraint,
									 double volume_constraint)
{
  // Before calling this function, the Mesh must contain a convex hull
  // of TRI3 elements which define the boundary.
  unsigned hull_integrity_check = check_hull_integrity();

  // Possibly die if hull integrity check failed
  this->process_hull_integrity_result(hull_integrity_check);

  // class tetgen_wrapper allows library access on a basic level
  TetGenWrapper tetgen_wrapper; 

  // Copy Mesh's node points into TetGen data structure
  this->fill_pointlist(tetgen_wrapper);

  // >>> fill input structure "tetgenio" with facet data:
  int facet_num = this->_mesh.n_elem();

  // allocate memory in "tetgenio" structure:
  tetgen_wrapper.allocate_facetlist(facet_num, holes.size());


  // Set up tetgen data structures with existing facet information
  // from the convex hull.
  {
    int insertnum = 0;
    MeshBase::element_iterator it        = this->_mesh.elements_begin();
    const MeshBase::element_iterator end = this->_mesh.elements_end();
    for (; it != end ; ++it)
      {
	tetgen_wrapper.allocate_facet_polygonlist(insertnum, 1);
	tetgen_wrapper.allocate_polygon_vertexlist(insertnum, 0, 3);

	Elem* elem = *it;

	for (unsigned int j=0; j<elem->n_nodes(); ++j)
	  {
	    // We need to get the sequential index of elem->get_node(j), but
	    // it should already be stored in _sequential_to_libmesh_node_map...
	    unsigned libmesh_node_id = elem->node(j);

	    // The libmesh node IDs may not be sequential, but can we assume
	    // they are at least in order???  We will do so here.
	    std::vector<unsigned>::iterator it = 
	    Utility::binary_find(_sequential_to_libmesh_node_map.begin(),
				 _sequential_to_libmesh_node_map.end(),
				 libmesh_node_id);

	    // Check to see if not found: this could also indicate the sequential
	    // node map is not sorted...
	    if (it == _sequential_to_libmesh_node_map.end())
	      {
		libMesh::err << "Global node " << libmesh_node_id << " not found in sequential node map!"  << std::endl;
		libmesh_error();
	      }
	    
	    std::vector<unsigned>::iterator::difference_type 
	       sequential_index = std::distance(_sequential_to_libmesh_node_map.begin(), it);
	    
	    // Debugging:
//	    std::cout << "libmesh_node_id=" << libmesh_node_id
//		      << ", sequential_index=" << sequential_index
//		      << std::endl;
	    
	    tetgen_wrapper.set_vertex(insertnum, // facet number
				      0,         // polygon (always 0)
				      j,         // local vertex index in tetgen input
				      sequential_index);
	  }

	// Go to next facet in polygonlist
	insertnum++;
      }
  }



  // fill hole list (if there are holes):
  if (holes.size() > 0)
    {
      std::vector< Node *>::const_iterator ihole;
      unsigned hole_index = 0;
      for (ihole=holes.begin(); ihole!=holes.end(); ++ihole)
	tetgen_wrapper.set_hole(hole_index++, (**ihole)(0), (**ihole)(1), (**ihole)(2));
    }

  
  // >>> run TetGen triangulation method:
  // assemble switches:
  std::ostringstream oss; // string holding switches
  oss << "pQ";

  if (quality_constraint != 0)
    oss << "q" << std::fixed << quality_constraint;
  
  if (volume_constraint != 0)
    oss << "a" << std::fixed << volume_constraint;
  
  std::string params = oss.str();

  tetgen_wrapper.set_switches(params); // TetGen switches: Piecewise linear complex, Quiet mode
  tetgen_wrapper.run_tetgen();

  // => nodes:
  unsigned int old_nodesnum = this->_mesh.n_nodes();
  REAL x=0., y=0., z=0.;
  const unsigned int num_nodes = tetgen_wrapper.get_numberofpoints();

  // Debugging:
  // std::cout << "Original mesh had " << old_nodesnum << " nodes." << std::endl;
  // std::cout << "Reserving space for " << num_nodes << " total nodes." << std::endl;

  // Reserve space for additional nodes in the node map
  _sequential_to_libmesh_node_map.reserve(num_nodes);

  // Add additional nodes to the Mesh.
  // Original code had i<=num_nodes here (Note: the indexing is:
  // foo[3*i], [3*i+1], [3*i+2]) But according to the TetGen docs, "In
  // all cases, the first item in any array is stored starting at
  // index [0]."
  for (unsigned int i=old_nodesnum; i<num_nodes; i++) 
    {
      // Fill in x, y, z values
      tetgen_wrapper.get_output_node(i, x,y,z);

      // Catch the node returned by add_point()... this will tell us the ID
      // assigned by the Mesh.
      Node* new_node = this->_mesh.add_point ( Point(x,y,z) );

      // Store this new ID in our sequential-to-libmesh node mapping array
      _sequential_to_libmesh_node_map.push_back( new_node->id() );
    }
  
  // Debugging:
//  std::copy(_sequential_to_libmesh_node_map.begin(), 
//	    _sequential_to_libmesh_node_map.end(), 
//	    std::ostream_iterator<unsigned>(std::cout, " "));
//  std::cout << std::endl;
  

  // => tetrahedra:
  const unsigned int num_elements = tetgen_wrapper.get_numberoftetrahedra();

  // Vector that temporarily holds the node labels defining element connectivity.
  unsigned int node_labels[4];      

  for (unsigned int i=0; i<num_elements; i++)
    {
      // TetGen only supports Tet4 elements.
      Elem* elem = new Tet4;

      // Fill up the the node_labels vector
      for (unsigned int j=0; j<elem->n_nodes(); j++) 
	node_labels[j] = tetgen_wrapper.get_element_node(i,j);
      
      // Associate nodes with this element
      this->assign_nodes_to_elem(node_labels, elem);

      // Finally, add this element to the mesh
      this->_mesh.add_elem(elem);
    }

  // Delete original convex hull elements.  Is there ever a case where
  // we should not do this?
  this->delete_2D_hull_elements();
}





void TetGenMeshInterface::fill_pointlist(TetGenWrapper& wrapper)
{
  // fill input structure with point set data:
  wrapper.allocate_pointlist( this->_mesh.n_nodes() );

  // Make enough space to store a mapping between the implied sequential
  // node numbering used in tetgen and libmesh's (possibly) non-sequential
  // numbering scheme.
  _sequential_to_libmesh_node_map.clear();
  _sequential_to_libmesh_node_map.resize( this->_mesh.n_nodes() );

  {
    unsigned index = 0;
    MeshBase::node_iterator it  = this->_mesh.nodes_begin();
    const MeshBase::node_iterator end = this->_mesh.nodes_end();
    for ( ; it != end; ++it) 
      {
	_sequential_to_libmesh_node_map[index] = (*it)->id();
	wrapper.set_node(index++, (**it)(0), (**it)(1), (**it)(2));
      }
  }
}





void TetGenMeshInterface::assign_nodes_to_elem(unsigned* node_labels, Elem* elem)
{
  for (unsigned int j=0; j<elem->n_nodes(); ++j)
    {
      // Get the mapped node index to ask the Mesh for
      unsigned mapped_node_id = _sequential_to_libmesh_node_map[ node_labels[j] ];

      // Parallel mesh can return NULL pointers, this is bad...
      Node* current_node = this->_mesh.node_ptr( mapped_node_id ); 

      if (current_node == NULL)
	{
	  std::cerr << "Error! Mesh returned NULL node pointer!" << std::endl;
	  libmesh_error();
	}

      elem->set_node(j) = current_node;
    }
}





unsigned TetGenMeshInterface::check_hull_integrity()
{
  MeshBase::element_iterator it        = this->_mesh.elements_begin();
  const MeshBase::element_iterator end = this->_mesh.elements_end();

  for (; it != end ; ++it)
    {
      Elem* elem = *it;

      // Check for proper element type
      if (elem->type() != TRI3)
	{
	  //libMesh::err << "ERROR: Some of the elements in the original mesh were not TRI3!" << std::endl;
	  //libmesh_error();
	  return 1;
	}

      for (unsigned int i=0; i<elem->n_neighbors(); ++i)
	{
	  if (elem->neighbor(i) == NULL)
	    {
	      // libMesh::err << "ERROR: Non-convex hull, cannot be tetrahedralized." << std::endl;
	      // libmesh_error();
	      return 2;
	    }
	}
    }

  // If we made it here, return success!
  return 0;
}





void TetGenMeshInterface::process_hull_integrity_result(unsigned result)
{
  if (result != 0)
    {
      libMesh::err << "Error! Conforming Delaunay mesh tetrahedralization requires a convex hull." << std::endl;
      libMesh::err << "Consider calling TetGenMeshInterface::pointset_convexhull() followed " << std::endl;
      libMesh::err << "by Mesh::find_neighbors() first." << std::endl;
      libmesh_error();
    }
}




void TetGenMeshInterface::delete_2D_hull_elements()
{
  MeshBase::element_iterator it        = this->_mesh.elements_begin();
  const MeshBase::element_iterator end = this->_mesh.elements_end();

  for (; it != end ; ++it)
    {
      Elem* elem = *it;

      // Check for proper element type. Yes, we legally delete elements while
      // iterating over them because no entries from the underlying container
      // are actually erased.
      if (elem->type() == TRI3)
	_mesh.delete_elem(elem);
    }
}
 
 

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_TETGEN
