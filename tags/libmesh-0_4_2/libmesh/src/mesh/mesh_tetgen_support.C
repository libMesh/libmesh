// $Id: mesh_tetgen_support.C,v 1.2 2004-01-03 15:37:43 benkirk Exp $

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


// C++ includes
#include <stdio.h>

// Local includes
#include "mesh_tetgen_support.h"
#include "mesh_data.h"
#include "cell_tet4.h"



//----------------------------------------------------------------------
// MeshBase methods
void MeshBase::read_tetgen(const std::string& name)
{
  std::string name_node, name_ele, dummy;
	
  // tetgen only works in 3D
  assert (this->mesh_dimension() == 3);

  // Check name for *.node or *.ele extension.
  // Set std::istream for node_stream and ele_stream.
  //
  if (name.rfind(".node") < name.size()) 
    {
      name_node = name;
      dummy     = name;
      int position = dummy.rfind(".node");
      name_ele     = dummy.replace(position, 5, ".ele");
    }
  else if (name.rfind(".ele") < name.size()) 
    {
      name_ele = name;
      dummy    = name;
      int position = dummy.rfind(".ele");
      name_node    = dummy.replace(position, 4, ".node");
    }
  else
    {
      std::cerr << "ERROR: Unrecognized file name: "
		<< name << std::endl;
      error();
    }


  
  // Set the streams from which to read in
  std::ifstream node_stream (name_node.c_str());
  std::ifstream ele_stream  (name_ele.c_str());

  if ( !node_stream.good() || !ele_stream.good() )
    {
      std::cerr << "ERROR: One or both Input file(s) not good." << std::endl
		<< "Error checking files "
		<< name_node << " and "
		<< name_ele  << std::endl;
      error();
    }

  // Skip the comment lines at the beginning
  skip_comment_lines (node_stream, '#');
  skip_comment_lines (ele_stream, '#');
  
  TetGenMeshInterface tetgen_interface (_nodes,
					_elements,
					data);

  tetgen_interface.read (node_stream, ele_stream);

  node_stream.close();
  ele_stream.close();
}



//----------------------------------------------------------------------
// TetGenMeshInterface class members
TetGenMeshInterface::TetGenMeshInterface (std::vector<Node*>& nodes,
					  std::vector<Elem*>& elements,
					  MeshData& md) :
  _nodes        (nodes),
  _elements     (elements),
  _num_nodes    (0),
  _num_elements (0),
  _mesh_data    (md)
{
}



TetGenMeshInterface::~TetGenMeshInterface()
{
}




//----------------------------------------------------------------------
// Read in the mesh from node_stream and ele_stream.
void TetGenMeshInterface::read (std::istream& node_stream,
				std::istream& ele_stream)
{
  _num_nodes    = 0;
  _num_elements = 0;


  // Read all the datasets.
  this->node_in    (node_stream);
  this->element_in (ele_stream);

  // Tell the MeshData object that we are finished 
  // reading data.
  _mesh_data.close_foreign_id_maps ();

  // some more clean-up
  _assign_nodes.clear();
}



//----------------------------------------------------------------------
// Function to read in the node table.
void TetGenMeshInterface::node_in (std::istream& node_stream)
{
  // Check input buffer
  assert (node_stream.good());

  unsigned int dimension=0, nAttributes=0, BoundaryMarkers=0;

  node_stream >> _num_nodes       // Read the number of nodes from the stream
	      >> dimension        // Read the dimension from the stream
	      >> nAttributes      // Read the number of attributes from stream
	      >> BoundaryMarkers; // Read if or not boundary markers are included in *.node (0 or 1)

  // Read the nodal coordinates from the node_stream (*.node file).
  unsigned int node_lab=0;
  Real x=0., y=0., z=0.;
  Real dummy;

  // Reserve space in the _nodes vector to avoid unnecessary allocations.
  _nodes.resize (_num_nodes);
    
  for (unsigned int i=0; i<_num_nodes; i++)
    {
      node_stream >> node_lab  // node number
		  >> x         // x-coordinate value
		  >> y         // y-coordinate value
		  >> z;        // z-coordinate value

      // For the number of attributes read all into dummy.
      for (unsigned int j=0; j<nAttributes; j++)
	node_stream >> dummy;
      
      // Read boundary marker if BoundaryMarker=1.
      if (BoundaryMarkers == 1)
	node_stream >> dummy;

      // Store the new position of the node under its label.
      _assign_nodes.insert (std::make_pair(node_lab,i));

      // Add node to the nodes vector.
      _nodes[i] = Node::build(x,y,z,i);

      // Tell the MeshData object the foreign node id.
      _mesh_data.add_foreign_node_id (_nodes[i], node_lab);
    }
}

/*
 *----------------------------------------------------------------------
 * Function to read in the element table.
 */
void TetGenMeshInterface::element_in (std::istream& ele_stream)
{
  // Check input buffer
  assert (ele_stream.good());

  // Read the elements from the ele_stream (*.ele file). 
  unsigned int element_lab, n_nodes, nAttri=0;
  Real dummy=0.0;
  unsigned long int node_labels[4];    // Vector that temporarily holds the node labels defining element.
		
  ele_stream >> _num_elements // Read the number of tetrahedrons from the stream.
	     >> n_nodes       // Read the number of nodes per tetrahedron from the stream (defaults to 4).
	     >> nAttri;       // Read the number of attributes from stream.
    
  // Reserve space in the appropriate vector to avoid unnecessary allocations.
  _elements.resize (_num_elements);

  // Vector that assigns element nodes to their correct position.
  // TetGen is normaly 0-based
  // (right now this is strictly not necessary since it is the identity map,
  //  but in the future TetGen could change their numbering scheme.)
  static const unsigned int assign_elm_nodes[] = { 0, 1, 2, 3};
			
  for (unsigned int i=0; i<_num_elements; i++)
    {
      ele_stream >> element_lab;
			
      // TetGen only supports Tet4 elements.
      _elements[i] = new Tet4;
      
      // Read node labels
      for (unsigned int j=0; j<n_nodes; j++)
	ele_stream >> node_labels[j];
			
      // Read attributes from the stream.
      for (unsigned int j=0; j<nAttri; j++)
	ele_stream >> dummy;
			
      // nodes are being stored in element
      for (unsigned int j=0; j<n_nodes; j++)
	_elements[i]->set_node(assign_elm_nodes[j]) = _nodes[_assign_nodes[node_labels[j]]];
      
      // tell the MeshData object the foreign element id
      _mesh_data.add_foreign_elem_id (_elements[i], element_lab);
    }
}
