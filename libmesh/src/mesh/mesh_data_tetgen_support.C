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



// C++ includes
#include <fstream>

// Local includes
#include "mesh_data.h"




//---------------------------------------------------------------------
// MeshDate TetGen support function
void MeshData::read_tetgen (const std::string& name)
{
  std::string name_node, name_ele, dummy;
  std::string desc = name;
  
  
  // Check name for *.node or *.ele extension.
  // Set std::istream for node_stream and ele_stream.
  if (name.rfind(".node") < name.size()) 
    {
      name_node = name;
      dummy     = name;
      int position = dummy.rfind(".node");
      name_ele     = dummy.replace(position, 5, ".ele");
      desc.erase(position);
    }
  else if (name.rfind(".ele") < name.size()) 
    {
      name_ele = name;
      dummy    = name;
      int position = dummy.rfind(".ele");
      name_node    = dummy.replace(position, 4, ".node");
      desc.erase(position);
    }
  else
    {
      std::cerr << "ERROR: Unrecognized file name: "
		<< name << std::endl;
      error();
    }
  
  // Set the streams from which to read in.
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

  
  // Set the descriptive name.
  // TetGen won't give a name, so we use the filename.
  this->_data_descriptor = desc;

  
  //--------------------------------------------------
  // Read in the data associated with the nodes.
  {
    unsigned int n_node=0, f_n_id=0, nAttri=0, BoundMark=0;
    Real dummy=0.0;
    std::vector<Number> AttriValue;

    // Read the parameters from the node_stream.
    node_stream >> n_node     // Read the number of nodes
		>> dummy      // Read the dimension
		>> nAttri     // Read the number of attributes
		>> BoundMark; // (0 or 1) boundary markers are in the stream or not.

    // Resize the values vector.
    AttriValue.resize(nAttri);

    for (unsigned int i=0; i<n_node; i++)
      {
	node_stream >> f_n_id;

	
	// Read the nodal coordinates for this node into dummy,
	// since we don't need them.
	for (unsigned int j=0; j<3; j++)
	  node_stream >> dummy;

	// Read the attributes from the stream.
	for (unsigned int j=0; j<nAttri; j++)
	  node_stream >> AttriValue[j];

	// Read boundary marker if BoundaryMarker=1.
	if (BoundMark == 1)
    	  node_stream >> dummy;

	// For the foreign node id locate the Node*.
     	const Node* node = foreign_id_to_node(f_n_id);
				
	// Insert this node and the values in our _node_data.
	_node_data.insert (std::make_pair(node, AttriValue));
      }
  }

  
  //--------------------------------------------------
  // Read in the data associated with the elements.
  {
    unsigned int n_elem, f_e_id, n_nodes, nAttri=0;
    Real dummy=0.0;
    std::vector<Number> AttriValue;

    // Read the parameters from the ele_stream.
    ele_stream >> n_elem   // Read the number of tetrahedrons
	       >> n_nodes  // Read the points per tetrahedron
	       >> nAttri;  // Read the number of attributes

    // Resize the values vector.
    AttriValue.resize(nAttri);

    for (unsigned int i=0; i<n_elem; i++)
      {
	ele_stream >> f_e_id;


	// For the number of nodes for this element read them into dummy,
	// since we don't need them.
	for (unsigned int n=0; n<n_nodes; n++)
	  {
	    ele_stream >> dummy;
	  }

	// Read the attributes from the stream.
	for (unsigned int j=0; j<nAttri; j++)
	  {
	    ele_stream >> AttriValue[j];
	  }

	// For the foreign elem id locate the Elem*.
     	const Elem* elem = foreign_id_to_elem(f_e_id);
				
	// Insert this elem and the values in our _elem_data.
	_elem_data.insert (std::make_pair(elem, AttriValue));
      }
  }
  
  //--------------------------------------------------
  // Finished reading.  Now ready for use.
  this->_node_data_closed = true;
  this->_elem_data_closed = true;

  node_stream.close();
  ele_stream.close();
}
