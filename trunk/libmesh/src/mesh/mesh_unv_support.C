// $Id: mesh_unv_support.C,v 1.2 2003-01-20 16:31:41 jwpeterson Exp $

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



//-----------------------------------------------------------------------------
// Mesh methods

void Mesh::read_unv(const std::string& name){
  std::ifstream file (name.c_str());
  read_unv(file);
  file.close();
  return;
};



void Mesh::read_unv(std::istream& in){
  UnvInterface i(in,_nodes,_elements);
};


//-----------------------------------------------------------------------------
// UnvInterface Methods

UnvInterface::UnvInterface(std::istream& _in,
			   std::vector<Node*>& _nodes,
			   std::vector<Elem*>& _elements):
  phys_file(_in),
  nodes(_nodes),
  elements(_elements)
{
  if ( !phys_file.good() )
  {
    std::cout << "Input file not good!" << std::endl;
    error();
  };

  label_dataset_nodes = "2411";
  label_dataset_elms  = "2412";
  num_nodes = 0;
  num_elements = 0;
  need_D_to_e = true;


  // use temporary file name as buffer
  temporary_file_name = tmpnam (NULL);
  std::cout << "Using temporary file: " << temporary_file_name << std::endl;

  temporary_file.open(temporary_file_name, std::fstream::out);
  if  ( !temporary_file.good() )
  {
    std::cout << "Error opening temporary file." << std::endl;
    error();
  };
   

  init();


  // writing operations in temporary_file are finished, now we close the file
  // and open it for reading
  temporary_file.close();
  temporary_file.open(temporary_file_name, std::fstream::in);

  if  ( !temporary_file.good() )
  {
    std::cout << "Error re-opening temporary file for reading." << std::endl;
    error();
  };
   

  node_in();
  element_in();
};



UnvInterface::~UnvInterface()
{

  temporary_file.close();

  if ( remove(temporary_file_name) == -1 )
  {
    std::cout << "Error deleting temporary file." << std::endl;
    error();
  };

};




void UnvInterface::init()
{
  std::string olds,news;

  while (true)            // work through the file until break (caused by eof)
    {
      phys_file >> olds >> news; // read two arguments from the stream

      // a "-1" followed by a number means the beginning of a dataset
      // stop combing at the end of the file
      while( ((olds != "-1") || (news == "-1") )
	     &&
	     !phys_file.eof() )
	{
	  olds = news;    // go on reading
	  phys_file >> news;
	};

      if(phys_file.eof()) // end of file is reached
	{ break; }

      scan_dataset(news); // scan the dataset for important information
    };

};



void UnvInterface::scan_dataset(std::string ds_num){

  // some datasets need special treatment

  // dataset containing the nodes
  if (ds_num == label_dataset_nodes)
    {
      //Write beginning of dataset to virtual file
      temporary_file << "    -1\n"
		     << "  " << ds_num << "\n";

      // store the position of the dataset in the virtual file
      ds_position[ds_num]=temporary_file.tellp();

      // if num_nodes is not 0 the dataset has already been scanned
      if (num_nodes != 0)
	{
	  std::cerr << "Error: UnvInterface::scan_dataset():\n"
		    << "Trying to scan nodes twice!" << std::endl;
	  error();
	  return;
	}

      // Read from file and store readings in data
      // Then write the data to the virtual file, check, if reals have
      // to be converted
      std::string data;                  // Sets of data to be read from file
      int i;                             // used for counting
      while (true)                       // read data until break
	{
	  phys_file >> data;             // read the node label
	   	
	  if (data == "-1")
	  {     
	    // end of dataset is reached
	    temporary_file << "    -1\n";
	    break;
	  }

	  temporary_file << "\t" << data;

	  for(i=0;i<3;i++)
	  {
	    phys_file >> data;
	    temporary_file << "\t" << data;
	  }
	
	  temporary_file << "\n"
		   << "  ";

	  if(need_D_to_e == true)
	    {
	      for(i=0;i<3;i++)
	      {
		phys_file >> data;
		temporary_file << D_to_e(data) << "\t";
	      }
	    }
	  else
	    {
	      for(i=0;i<3;i++)
	      {
		phys_file >> data;
		temporary_file << data << "\t";
	      }
	    }

	  temporary_file << "\n";

	  num_nodes++;                   // count nodes
	};
    }

  // dataset containing the elements
  else if (ds_num == label_dataset_elms)
    {
      //Write beginning of dataset to virtual file
      temporary_file << "    -1\n"
		     << "  " << ds_num << "\n";

      // store the position of the dataset in the virtual file
      ds_position[ds_num]=temporary_file.tellp();

      // if num_elements is not 0 the dataset has already been scanned
      if (num_elements != 0)
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
	  phys_file >> data;             // read element label
	  if (data == "-1") 
	  { 
	    // end of dataset is reached
	    temporary_file << "    -1\n";
	    break;
	  }

	  temporary_file << "\t" << data;
	
	  // Data of the node, including the number of nodes, which come next
	  // For further information read comments of element_in()
	  for(i=0;i<4;i++)
	  {
	    phys_file >> data;
	    temporary_file << "\t" << data;
	  }

	  phys_file >> i;                 // number of nodes
	  temporary_file << "\t" << i << "\n";
	  while (i > 0)                   // ignore the nodes
	    {
	      phys_file >> data;
	      temporary_file << "\t" << data;
	      i--;
	    }

	  temporary_file << "\n";
	  num_elements++;                 // count elements
	};
    }

  // datasets that are not of special interest are ignored
  else
    {
      std::string data;
      do {
	phys_file >> data;                // read the beginning of every line
	phys_file.ignore(256,'\n');}      // ignore the rest
      while (data != "-1");               // look for delimiter
    }
};



void UnvInterface::set_stream_pointer(std::string ds_num)
{
  // An error message is displayed if the specified
  // dataset does not exist
  if (static_cast<int>(ds_position[ds_num]) == 0)
    {
      std::cerr
	<< "Error: UnvInterface::set_stream_pointer(std::string ds_num):\n"
	<< "Dataset #" << ds_num << " not found in file." << std::endl;
      error();
      return;
    }

  // Move file pointer
  temporary_file.seekg(ds_position[ds_num],std::ios::beg);
};



// Method reads nodes from the virtual file and stores them in
// vector<Node> nodes in the order they come in.
// The original node labels are being stored in the
// map assign_nodes in order to assign the elements to
// the right nodes later.
void UnvInterface::node_in()
{
  // Variables needed for reading
  unsigned int node_lab;// label of the node
  unsigned long int exp_coord_sys_num,  // export coordinate system number(not supported yet)
                    disp_coord_sys_num, // displacement coordinate system number(not supp. yet)
                    color;              // color(not supported yet)
  real x,y,z;           // coordinates of the node

  // allocate the correct amount of memory for the vector
  // that holds the nodes
  nodes.resize(num_nodes);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(label_dataset_nodes);

  // read from virtual file
  for(unsigned int i=0;i<num_nodes;i++)
    {
      temporary_file >> node_lab                // read the node label
		     >> exp_coord_sys_num       // (not supported yet)
		     >> disp_coord_sys_num      // (not supported yet)
		     >> color                   // (not supported yet)
		     >> x                       // read x-coordinate
		     >> y                       // read y-coordinate
		     >> z;                      // read z-coordinate

      // store the new position of the node under its label
      assign_nodes[node_lab]=i;

      // add node to the nodes vector.
      nodes[i] = Node::build(x,y,z,i);
    };
};



// Method reads elements and stores them in
// vector<Elem*> elements in the same order as they
// come in. That means labels are ignored!
// Changes in node numbers are taken into account.
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
  elements.resize(num_elements);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(label_dataset_elms);

  // read from the virtual file
  for (unsigned int i=0;i<num_elements;i++)
    {
      temporary_file >> element_lab             // read element label
		     >> fe_descriptor_id        // read FE descriptor id
		     >> phys_prop_tab_num       // (not supported yet)
		     >> mat_prop_tab_num        // (not supported yet)
		     >> color                   // (not supported yet)
		     >> n_nodes;                // read number of nodes on element

      for (unsigned int j=1;j<=n_nodes;j++)
	temporary_file >> node_labels[j];             // read node labels

      // UNV-FE descriptor id has to be recognized and translated to mesh
      switch (fe_descriptor_id)
	{
	
	case 41: // Plane Stress Linear Triangle
	case 91: // Thin Shell   Linear Triangle
	  elements[i]=Elem::build(TRI3);  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=2;
	  assign_elm_nodes[3]=1;
	  break;
	
	case 42: // Plane Stress Quadratic Triangle
	case 92: // Thin Shell   Quadratic Triangle
	  elements[i]=Elem::build(TRI6);  // add new element
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
	  elements[i]=Elem::build(QUAD4); // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=3;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=1;
	  break;
	
	case 45: // Plane Stress Quadratic Quadrilateral
	case 95: // Thin Shell   Quadratic Quadrilateral
	  elements[i]=Elem::build(QUAD8); // add new element
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
	  elements[i]=Elem::build(TET4);  // add new element
	  assign_elm_nodes[1]=0;
	  assign_elm_nodes[2]=1;
	  assign_elm_nodes[3]=2;
	  assign_elm_nodes[4]=3;
	  break;

	
	case 115: // Solid Linear Brick
	  elements[i]=Elem::build(HEX8);  // add new element
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
	  elements[i]=Elem::build(HEX20); // add new element
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
	  elements[i]=Elem::build(TET10); // add new element
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
	elements[i]->set_node(assign_elm_nodes[j]) = nodes[assign_nodes[node_labels[j]]];

    }
};



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
    need_D_to_e = false;
  };

  return number;

};


