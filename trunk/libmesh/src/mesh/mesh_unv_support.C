// $Id: mesh_unv_support.C,v 1.10 2003-05-14 11:54:37 ddreyer Exp $

// The Next Great Finite Element Library.
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
# include <stdio.h>


// Local includes
#include "mesh_unv_support.h"
#include "boundary_data.h"
#include "face_quad4.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad8.h"
#include "cell_tet4.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_tet10.h"
#include "cell_prism6.h"

//-----------------------------------------------------------------------------
// Mesh methods

void Mesh::read_unv(const std::string& name){
  std::ifstream file (name.c_str());
  read_unv(file);
  file.close();
  return;
}



void Mesh::read_unv(std::istream& in){
  UnvInterface i(in,_nodes,_elements,boundary_info);
}


//-----------------------------------------------------------------------------
// UnvInterface Methods

UnvInterface::UnvInterface(std::istream& in,
			   std::vector<Node*>& nodes,
			   std::vector<Elem*>& elements,
			   BoundaryInfo& boundary_info):
  _phys_file     (in),
  _nodes         (nodes),
  _elements      (elements),
  _boundary_info (boundary_info)
{
  if ( !_phys_file.good() )
  {
    std::cout << "Input file not good!" << std::endl;
    error();
  }

  _label_dataset_nodes = "2411";
  _label_dataset_elms  = "2412";
  _num_nodes = 0;
  _num_elements = 0;
  _need_D_to_e = true;


  // use temporary file name as buffer
  _temporary_file_name = tmpnam (NULL);
  std::cout << "Using temporary file: " << _temporary_file_name << std::endl;

  _temporary_file.open(_temporary_file_name, std::fstream::out);
  if  ( !_temporary_file.good() )
  {
    std::cout << "Error opening temporary file." << std::endl;
    error();
  }
   
  init();


  // writing operations in _temporary_file are finished, now we close the file
  // and open it for reading
  _temporary_file.close();
  _temporary_file.open(_temporary_file_name, std::fstream::in);

  if  ( !_temporary_file.good() )
  {
    std::cout << "Error re-opening temporary file for reading." << std::endl;
    error();
  }
   
  node_in();
  element_in();
}



UnvInterface::~UnvInterface()
{

  _temporary_file.close();

  if ( remove(_temporary_file_name) == -1 )
  {
    std::cout << "Error deleting temporary file." << std::endl;
    error();
  }

}




void UnvInterface::init()
{
  std::string olds,news;

  while (true)            // work through the file until break (caused by eof)
    {
      _phys_file >> olds >> news; // read two arguments from the stream

      // a "-1" followed by a number means the beginning of a dataset
      // stop combing at the end of the file
      while( ((olds != "-1") || (news == "-1") )
	     &&
	     !_phys_file.eof() )
	{
	  olds = news;    // go on reading
	  _phys_file >> news;
	}

      if(_phys_file.eof()) // end of file is reached
	{ break; }

      scan_dataset(news); // scan the dataset for important information
    }

}



void UnvInterface::scan_dataset(std::string ds_num){

  // some datasets need special treatment

  // dataset containing the nodes
  if (ds_num == _label_dataset_nodes)
    {
      //Write beginning of dataset to virtual file
      _temporary_file << "    -1\n"
		     << "  " << ds_num << "\n";

      // store the position of the dataset in the virtual file
      _ds_position[ds_num]=_temporary_file.tellp();

      // if _num_nodes is not 0 the dataset has already been scanned
      if (_num_nodes != 0)
	{
	  std::cerr << "Error: UnvInterface::scan_dataset():\n"
		    << "Trying to scan nodes twice!" << std::endl;
	  error();
	  return;
	}

      // Read from file and store readings in data
      // Then write the data to the virtual file, check, if Reals have
      // to be converted
      std::string data;                  // Sets of data to be read from file
      int i;                             // used for counting
      while (true)                       // read data until break
	{
	  _phys_file >> data;             // read the node label
	   	
	  if (data == "-1")
	  {     
	    // end of dataset is reached
	    _temporary_file << "    -1\n";
	    break;
	  }

	  _temporary_file << "\t" << data;

	  for(i=0;i<3;i++)
	  {
	    _phys_file >> data;
	    _temporary_file << "\t" << data;
	  }
	
	  _temporary_file << "\n"
		   << "  ";

	  if(_need_D_to_e == true)
	    {
	      for(i=0;i<3;i++)
	      {
		_phys_file >> data;
		_temporary_file << D_to_e(data) << "\t";
	      }
	    }
	  else
	    {
	      for(i=0;i<3;i++)
	      {
		_phys_file >> data;
		_temporary_file << data << "\t";
	      }
	    }

	  _temporary_file << "\n";

	  _num_nodes++;                   // count nodes
	}
    }

  // dataset containing the elements
  else if (ds_num == _label_dataset_elms)
    {
      //Write beginning of dataset to virtual file
      _temporary_file << "    -1\n"
		     << "  " << ds_num << "\n";

      // store the position of the dataset in the virtual file
      _ds_position[ds_num]=_temporary_file.tellp();

      // if _num_elements is not 0 the dataset has already been scanned
      if (_num_elements != 0)
	{
	  std::cerr << "Error: UnvInterface::scan_dataset():\n"
		    << "Trying to scan elements twice!" << std::endl;
	  error();
	  return;
	}

      std::string data;                  // Sets of data to be read from file
      unsigned int i;                    // used for counting

      while (true)                       // scan data until break
	{
	  _phys_file >> data;             // read element label
	  if (data == "-1") 
	  { 
	    // end of dataset is reached
	    _temporary_file << "    -1\n";
	    break;
	  }

	  _temporary_file << "\t" << data;
	
	  // Data of the node, including the number of nodes, which come next
	  // For further information read comments of element_in()
	  for(i=0;i<4;i++)
	  {
	    _phys_file >> data;
	    _temporary_file << "\t" << data;
	  }

	  _phys_file >> i;                 // number of nodes
	  _temporary_file << "\t" << i << "\n";
	  while (i > 0)                   // ignore the nodes
	    {
	      _phys_file >> data;
	      _temporary_file << "\t" << data;
	      i--;
	    }

	  _temporary_file << "\n";
	  _num_elements++;                 // count elements
	}
    }


  // datasets that are not of special interest are ignored
  else
    {
      std::string data;
      do {
	_phys_file >> data;                // read the beginning of every line
	_phys_file.ignore(256,'\n');}      // ignore the rest
      while (data != "-1");               // look for delimiter
    }
}



void UnvInterface::set_stream_pointer(std::string ds_num)
{
  // An error message is displayed if the specified
  // dataset does not exist
  if (static_cast<int>(_ds_position[ds_num]) == 0)
    {
      std::cerr << "Error: UnvInterface::set_stream_pointer(std::string ds_num):\n"
		<< "Dataset #" << ds_num 
		<< " not found in file." << std::endl;
      error();
      return;
    }

  // Move file pointer
  _temporary_file.seekg(_ds_position[ds_num],std::ios::beg);
}





void UnvInterface::node_in()
{
  // Variables needed for reading
  unsigned int node_lab;// label of the node
  unsigned long int exp_coord_sys_num,  // export coordinate system number(not supported yet)
                    disp_coord_sys_num, // displacement coordinate system number(not supp. yet)
                    color;              // color(not supported yet)
  Real x,y,z;           // coordinates of the node

  // allocate the correct amount of memory for the vector
  // that holds the nodes
  _nodes.resize(_num_nodes);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(_label_dataset_nodes);

  // read from virtual file
  for(unsigned int i=0;i<_num_nodes;i++)
    {
      _temporary_file >> node_lab                // read the node label
		     >> exp_coord_sys_num       // (not supported yet)
		     >> disp_coord_sys_num      // (not supported yet)
		     >> color                   // (not supported yet)
		     >> x                       // read x-coordinate
		     >> y                       // read y-coordinate
		     >> z;                      // read z-coordinate

      // store the new position of the node under its label
      _assign_nodes[node_lab]=i;

      // add node to the nodes vector.
      _nodes[i] = Node::build(x,y,z,i);

      // when there is a BoundaryData object, add the foreign id
      if (_boundary_info.has_boundary_data())
	  _boundary_info.get_boundary_data().add_foreign_node_id(_nodes[i],
								 node_lab);

    }


  // when there is a BoundaryData object, tell it we are finished with nodes
  if (_boundary_info.has_boundary_data())
      _boundary_info.get_boundary_data().close_node_map();
}



void UnvInterface::element_in()
{
  // Variables needed for reading from the file
  unsigned int element_lab,            // element label (not supported yet)
               n_nodes;                // number of nodes on element
  unsigned long int fe_descriptor_id,  // FE descriptor id
                    phys_prop_tab_num, // physical property table number (not supported yet)
                    mat_prop_tab_num,  // material property table number (not supported yet)
                    color;             // color (not supported yet)


  // vector that temporarily holds the node labels defining element
  // (this was named "nodes" but that clashes with UnvInterface::nodes
  unsigned long int node_labels[21];


  // vector that assigns element nodes to their correct position
  // for example:
  // 44:plane stress      | QUAD4
  // linear quadrilateral |
  // position in UNV-file | position in libmesh
  // assign_elm_node[1]   = 0
  // assign_elm_node[2]   = 3
  // assign_elm_node[3]   = 2
  // assign_elm_node[4]   = 1
  unsigned int assign_elm_nodes[21];

  // UNV is 1-based, we leave the 0th element of the vectors unused in order
  // to prevent confusion, this way we can store elements with up to 20 nodes

  // allocate the correct amount of memory for the vector
  // that holds the elements
  _elements.resize(_num_elements);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(_label_dataset_elms);

  // read from the virtual file
  for (unsigned int i=0;i<_num_elements;i++)
    {
      _temporary_file >> element_lab             // read element label
		     >> fe_descriptor_id        // read FE descriptor id
		     >> phys_prop_tab_num       // (not supported yet)
		     >> mat_prop_tab_num        // (not supported yet)
		     >> color                   // (not supported yet)
		     >> n_nodes;                // read number of nodes on element

      for (unsigned int j=1;j<=n_nodes;j++)
	_temporary_file >> node_labels[j];             // read node labels

      // UNV-FE descriptor id has to be recognized and translated to mesh
      switch (fe_descriptor_id)
	{
	
	case 41: // Plane Stress Linear Triangle
	case 91: // Thin Shell   Linear Triangle
	  _elements[i]=new Tri3;  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=2;
	  assign_elm_nodes[3]=1;
	  break;
	
	case 42: // Plane Stress Quadratic Triangle
	case 92: // Thin Shell   Quadratic Triangle
	  _elements[i]=new Tri6;  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=5;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=4;
	  assign_elm_nodes[5]=1;
	  assign_elm_nodes[6]=3;
	  break;
	
	case 43: // Plane Stress Cubic Triangle
	  std::cerr << "Error: UnvInterface::element_in():\n"
		    << "UNV-element type 43: Plane Stress Cubic Triangle"
		    << " not supported." << std::endl;
	  error();
	
	  break;
	
	case 44: // Plane Stress Linear Quadrilateral
	case 94: // Thin Shell   Linear Quadrilateral
	  _elements[i]=new Quad4; // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=3;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=1;
	  break;
	
	case 45: // Plane Stress Quadratic Quadrilateral
	case 95: // Thin Shell   Quadratic Quadrilateral
	  _elements[i]=new Quad8; // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=7;
	  assign_elm_nodes[3]=3;
	  assign_elm_nodes[4]=6;
	  assign_elm_nodes[5]=2;
	  assign_elm_nodes[6]=5;
	  assign_elm_nodes[7]=1;
	  assign_elm_nodes[8]=4;
	  break;
	
	case 46: // Plane Stress Cubic Quadrilateral
	  std::cerr << "Error: UnvInterface::element_in():\n"
		    << "UNV-element type 46: Plane Stress Cubic Quadrilateral"
		    << " not supported." << std::endl;
	  error();
	  break;
	
	
	case 111: // Solid Linear Tetrahedron
	  _elements[i]=new Tet4;  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=1;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=3;
	  break;


	case 112: // Solid Linear Prism
	  _elements[i]=new Prism6;  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=1;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=3;
	  assign_elm_nodes[5]=4;
	  assign_elm_nodes[6]=5;
	  break;
	
	case 115: // Solid Linear Brick
	  _elements[i]=new Hex8;  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=4;
	  assign_elm_nodes[3]=5;
	  assign_elm_nodes[4]=1;
	  assign_elm_nodes[5]=3;
	  assign_elm_nodes[6]=7;
	  assign_elm_nodes[7]=6;
	  assign_elm_nodes[8]=2;
	  break;
	
	
	case 116: // Solid Quadratic Brick
	  _elements[i]=new Hex20; // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=12;
	  assign_elm_nodes[3]=4;
	  assign_elm_nodes[4]=16;
	  assign_elm_nodes[5]=5;
	  assign_elm_nodes[6]=13;
	  assign_elm_nodes[7]=1;
	  assign_elm_nodes[8]=8;
	  assign_elm_nodes[9]=11;
	  assign_elm_nodes[10]=19;
	  assign_elm_nodes[11]=17;
	  assign_elm_nodes[12]=9;
	  assign_elm_nodes[13]=3;
	  assign_elm_nodes[14]=15;
	  assign_elm_nodes[15]=7;
	  assign_elm_nodes[16]=18;
	  assign_elm_nodes[17]=6;
	  assign_elm_nodes[18]=14;
	  assign_elm_nodes[19]=2;
	  assign_elm_nodes[20]=10;
	  break;
	
	
	case 117: // Solid Cubic Brick
	  std::cerr << "Error: UnvInterface::element_in():\n"
		    << "UNV-element type 117: Solid Cubic Brick"
		    << " not supported." << std::endl;
	  error();
	
	  break;

	
	case 118: // Solid Parabolic Tetrahedron
	  _elements[i]=new Tet10; // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=4;
	  assign_elm_nodes[3]=1;
	  assign_elm_nodes[4]=5;
	  assign_elm_nodes[5]=2;
	  assign_elm_nodes[6]=6;
	  assign_elm_nodes[7]=7;
	  assign_elm_nodes[8]=8;
	  assign_elm_nodes[9]=9;
	  assign_elm_nodes[10]=3;
	  break;
	
	
	default: // Unrecognized element type
	  std::cerr << "Error: UnvInterface::element_in():\n"
		    << "UNV-element type " << fe_descriptor_id
		    << " not supported." << std::endl;
	  error();
	
	  break;
	}

      // nodes are being stored in element
      for (unsigned int j=1;j<=n_nodes;j++)
	_elements[i]->set_node(assign_elm_nodes[j]) = _nodes[_assign_nodes[node_labels[j]]];


      // when there is a BoundaryData object, add the foreign elem id
      if (_boundary_info.has_boundary_data())
	  _boundary_info.get_boundary_data().add_foreign_elem_id(_elements[i],
								 element_lab);


    }

  // when there is a BoundaryData object, tell it we are finished with nodes
  if (_boundary_info.has_boundary_data())
      _boundary_info.get_boundary_data().close_elem_map();
}



// Method that converts strings containing float point numbers in D-notation
// to e-notation.
// e.g. 3.141592654D+00 --> 3.141592654e+00
// in order to make it readable for C++.
std::string& UnvInterface::D_to_e(std::string& number)
{
  // position of the "D" in the string
  unsigned int position;

  // find "D" in string, start looking at 6th element, to improve speed.
  // We dont expect a "D" earlier
  position = number.find("D",6);

  if(position!=std::string::npos)     // npos means no position
  {
    // replace "D" in string
    number.replace(position,1,"e"); 
  }
  else
  {
    // we assume that if this one number is written correctly, all numbers are
    _need_D_to_e = false;
  }

  return number;

}


