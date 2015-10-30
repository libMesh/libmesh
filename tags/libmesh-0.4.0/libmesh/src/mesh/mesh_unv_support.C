// $Id: mesh_unv_support.C,v 1.14 2003-05-20 22:43:10 benkirk Exp $

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
#include <stdio.h>
#include <iomanip>


// Local includes
#include "mesh_unv_support.h"
#include "mesh_data.h"
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
void Mesh::read_unv(const std::string& name)
{
  std::ifstream file (name.c_str());

  UnvMeshInterface unv_interface (_nodes,
				  _elements,
				  data);

  unv_interface.read (file);

  file.close();
}



void Mesh::write_unv(const std::string& name)
{
  std::ofstream file (name.c_str());

  UnvMeshInterface unv_interface (_nodes,
				  _elements,
				  data);

  unv_interface.write (file);

  file.close();
}





//-----------------------------------------------------------------------------
// UnvMeshInterface class static members
const std::string  UnvMeshInterface::_label_dataset_nodes       = "2411";
const std::string  UnvMeshInterface::_label_dataset_elements    = "2412";






//-----------------------------------------------------------------------------
// UnvMeshInterface class members
UnvMeshInterface::UnvMeshInterface (std::vector<Node*>& nodes,
				    std::vector<Elem*>& elements,
				    MeshData& md) :
  _nodes         (nodes),
  _elements      (elements),
  _mesh_data     (md)
{
  /*
   * Initialize these to dummy values
   */
  _num_nodes    = 0;
  _num_elements = 0;
  _need_D_to_e  = true;

  _assign_nodes.clear();
  _ds_position.clear();
}








UnvMeshInterface::~UnvMeshInterface()
{
}








void UnvMeshInterface::read (std::istream& in_stream,
			     bool verbose)
{
  if ( !in_stream.good() )
    {
      std::cerr << "ERROR: Input file not good." 
		<< std::endl;
      error();
    }

  _num_nodes    = 0;
  _num_elements = 0;
  _need_D_to_e  = true;

    
  /*
   * When reading, the data is stored in a temporary
   * file.  This is the tempfile.
   */
  std::string temp_buffer_name = tmpnam (NULL);
  if (verbose)
    {  
      std::cout << "Using temporary file: " 
		<< temp_buffer_name 
		<< std::endl;
    }

  std::fstream temp_buffer;
  temp_buffer.open (temp_buffer_name.c_str(), 
		    std::fstream::out);
  if  (!temp_buffer.good())
    {
      std::cerr << "ERROR: Could not open temporary file in write mode." 
		<< std::endl;
      error();
    }
   

  /*
   * locate the beginning of data sets,
   * store interesting datasets (through
   * \p buffer_interesting_datasets) in 
   * the tempfile, extract data useful for
   * faster re-reading of the tempfile
   */
  {
    std::string olds, news;

    while (true)
    {
      in_stream >> olds >> news;

      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while( ((olds != "-1") || (news == "-1") ) && !in_stream.eof() )
	{
	  olds = news;
	  in_stream >> news;
	}

      if(in_stream.eof())
	break;

      /*
       * if beginning of dataset, buffer it in
       * temp_buffer, if desired
       */
      if (news == _label_dataset_nodes)
        {
	  buffer_nodes (in_stream,
			temp_buffer);
        }

      else if (news == _label_dataset_elements)
        {
	  buffer_elements (in_stream,
			   temp_buffer);
        }

      else
        {
	  /*
	   * other datasets are ignored
	   */
        }

    }
  }


  /*
   * writing to temp_buffer is finished, 
   * close the file and open it for reading
   */
  {
    temp_buffer.close ();
    temp_buffer.open  (temp_buffer_name.c_str(), 
			    std::fstream::in);

    if  (!temp_buffer.good())
      {
        std::cout << "ERROR:  Could not open temporary file in read mode." 
		  << std::endl;
	error();
      }
  }


  /*
   * read all the datasets in the order
   * as given by the user, but only once!
   * Datasets that are already read are
   * added to the set \p already_read.
   */
  node_in    (temp_buffer);
  element_in (temp_buffer);


  /*
   * close the tempfile and try to delete it
   */
  {
    temp_buffer.close();

    if ( remove(temp_buffer_name.c_str()) == -1 )
      {
        std::cerr << "ERROR: Cannot delete temporary file." 
		  << std::endl;
	error();
      }
  }


  /*
   * tell the MeshData object that we are finished 
   * reading data
   */
  _mesh_data.close_foreign_id_maps ();


  /* 
   * some more clean-up
   */
  _assign_nodes.clear();
  _ds_position.clear();
}






void UnvMeshInterface::write (std::ostream& out_stream)
{
  if ( !out_stream.good() )
    {
      std::cerr << "ERROR: Input file not good." 
		<< std::endl;
      error();
    }


  /*
   * already know these data, so initialize
   * them.  Does not hurt.
   */
  _num_nodes    = _nodes.size();
  _num_elements = _elements.size();
  _need_D_to_e  = false;


  /*
   * write the nodes,  then the elements
   */
  node_out         (out_stream);
  element_out      (out_stream);
}





void UnvMeshInterface::buffer_nodes (std::istream& physical_file,
				     std::fstream& temp_file)
{
  std::string data;

  //Write beginning of dataset to virtual file
  temp_file << "    -1\n"
	    << "  " 
	    << _label_dataset_nodes 
	    << "\n";

  // store the position of the dataset in the virtual file
  _ds_position[_label_dataset_nodes]=temp_file.tellp();

  // if _num_nodes is not 0 the dataset has already been scanned
  if (_num_nodes != 0)
    {
      std::cerr << "Error: Trying to scan nodes twice!" 
		<< std::endl;
      error();
    }

  /* 
   * Read from file and buffer in temp_file,
   * check if Reals have to be converted
   */
  while (true)                       // read data until break
    {
      physical_file >> data;             // read the node label
	   	
      if (data == "-1")
        {     
	  // end of dataset is reached
	  temp_file << "    -1\n";
	  break;
	}

      temp_file << "\t" << data;

      for(unsigned int i=0;i<3;i++)
        {
	  physical_file >> data;
	  temp_file << "\t" << data;
	}
	
      temp_file << "\n  ";

      if(_need_D_to_e)
        {
	  for(unsigned int i=0;i<3;i++)
	    {  
	      physical_file >> data;
	      temp_file << D_to_e(data) << "\t";
	    }
	}
      else
        {
	  for(unsigned int i=0;i<3;i++)
	    {
	      physical_file >> data;
	      temp_file << data << "\t";
	    }
	}

      temp_file << "\n";

      _num_nodes++;                   // count nodes
    }
}






void UnvMeshInterface::buffer_elements (std::istream& physical_file,
					std::fstream& temp_file)
{
  std::string data;

  /*
   * buffer the element dataset in tempfile
   */
  temp_file << "    -1\n  " 
	    << _label_dataset_elements
	    << "\n";

  /*
   * remember the position of the dataset 
   * in the virtual file
   */
  _ds_position[_label_dataset_elements]=temp_file.tellp();

  if (_num_elements != 0)
    {
      std::cerr << "Error: Trying to scan elements twice!" 
		<< std::endl;
      error();
    }

  while (true)
    {
      physical_file >> data;             // read element label
      if (data == "-1") 
        { 
	  // end of dataset is reached
	  temp_file << "    -1\n";
	  break;
	}

      temp_file << "\t" << data;
	
      /*
       * Nodes & related data, including the 
       * number of nodes, just copy them
       */
      for(unsigned int cnt=0; cnt<4; cnt++)
        {
	  physical_file >> data;
	  temp_file << "\t" << data;
	}

      unsigned int n_nodes = 0;
      physical_file >> n_nodes;
      assert (n_nodes != 0);
      temp_file << "\t" << n_nodes << "\n";
      for (unsigned int cnt=0; cnt<n_nodes; cnt++)
        {
	  physical_file >> data;
	  temp_file << "\t" << data;
	}
      
      temp_file << "\n";
      _num_elements++;                 // count elements
    }
}








void UnvMeshInterface::set_stream_pointer(std::fstream& temp_file,
					  const std::string& ds_num)
{
  // An error message is displayed if the specified
  // dataset does not exist
  if (static_cast<int>(_ds_position[ds_num]) == 0)
    {
      std::cerr << "Error: Dataset #" 
		<< ds_num 
		<< " not found in file." 
		<< std::endl;
      error();
      return;
    }

  // Move file pointer
  temp_file.seekg(_ds_position[ds_num],std::ios::beg);
}





void UnvMeshInterface::node_in (std::fstream& temp_file)
{
  unsigned int node_lab;                // label of the node
  unsigned long int exp_coord_sys_num,  // export coordinate system number(not supported yet)
                    disp_coord_sys_num, // displacement coordinate system number(not supp. yet)
                    color;              // color(not supported yet)
  Real x,y,z;                           // coordinates of the node

  // allocate the correct amount of memory for the vector
  // that holds the nodes
  _nodes.resize(_num_nodes);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(temp_file,
		     _label_dataset_nodes);

  // read from virtual file
  for(unsigned int i=0;i<_num_nodes;i++)
    {
      temp_file >> node_lab                // read the node label
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

      // tell the MeshData object the foreign node id
      _mesh_data.add_foreign_node_id (_nodes[i],
				      node_lab);
    }

}



void UnvMeshInterface::element_in (std::fstream& temp_file)
{
  // Variables needed for reading from the file
  unsigned int element_lab,            // element label (not supported yet)
               n_nodes;                // number of nodes on element
  unsigned long int fe_descriptor_id,  // FE descriptor id
                    phys_prop_tab_num, // physical property table number (not supported yet)
                    mat_prop_tab_num,  // material property table number (not supported yet)
                    color;             // color (not supported yet)


  // vector that temporarily holds the node labels defining element
  // (this was named "nodes" but that clashes with UnvMeshInterface::nodes
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
  //
  // UNV is 1-based, we leave the 0th element of the vectors unused in order
  // to prevent confusion, this way we can store elements with up to 20 nodes
  unsigned int assign_elm_nodes[21];

  // allocate the correct amount of memory for the vector
  // that holds the elements
  _elements.resize(_num_elements);

  // put the file-pointer at the beginning of the dataset
  set_stream_pointer(temp_file,
		     _label_dataset_elements);

  // read from the virtual file
  for (unsigned int i=0;i<_num_elements;i++)
    {
      temp_file >> element_lab             // read element label
		>> fe_descriptor_id        // read FE descriptor id
		>> phys_prop_tab_num       // (not supported yet)
		>> mat_prop_tab_num        // (not supported yet)
		>> color                   // (not supported yet)
		>> n_nodes;                // read number of nodes on element

      for (unsigned int j=1;j<=n_nodes;j++)
	temp_file >> node_labels[j];             // read node labels

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
	  std::cerr << "ERROR: UNV-element type 43: Plane Stress Cubic Triangle"
		    << " not supported." 
		    << std::endl;
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
	  std::cerr << "ERROR: UNV-element type 46: Plane Stress Cubic Quadrilateral"
		    << " not supported." 
		    << std::endl;
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
	  std::cerr << "Error: UNV-element type 117: Solid Cubic Brick"
		    << " not supported." 
		    << std::endl;
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
	  std::cerr << "ERROR: UNV-element type " 
		    << fe_descriptor_id
		    << " not supported." 
		    << std::endl;
	  error();
	
	  break;
	}

      // nodes are being stored in element
      for (unsigned int j=1;j<=n_nodes;j++)
	_elements[i]->set_node(assign_elm_nodes[j]) = _nodes[_assign_nodes[node_labels[j]]];


      // tell the MeshData object the foreign element id
      _mesh_data.add_foreign_elem_id (_elements[i],
				      element_lab);

    }
}






void UnvMeshInterface::node_out (std::ostream& out_file)
{
  /*
   * we need the MeshData, otherwise we do not
   * know the foreign node id
   */
  if (!_mesh_data.active())
    {
      std::cerr << "ERROR: Need to have an active MeshData for writing UNV files."
		<< std::endl;
      error();
    }

  /*
   * Write beginning of dataset
   */
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_nodes
	   << "\n";


  unsigned int exp_coord_sys_dummy  = 0; // export coordinate sys. (not supported yet)
  unsigned int disp_coord_sys_dummy = 0; // displacement coordinate sys. (not supp. yet)
  unsigned int color_dummy          = 0; // color(not supported yet)


  for(unsigned int i=0;i<_nodes.size();i++)
    {
      const Node* current_node = _nodes[i];
      char buf[78];
      sprintf(buf, "%10d%10d%10d%10d\n", 
	      _mesh_data.node_to_foreign_id(current_node),
	      exp_coord_sys_dummy,
	      disp_coord_sys_dummy,
	      color_dummy);
      out_file << buf;

      // the coordinates
      if (DIM == 3)
	  sprintf(buf, "%25.16E%25.16E%25.16E\n", 
		  (*current_node)(0),
		  (*current_node)(1),
		  (*current_node)(2));
      else if (DIM == 2)
	  sprintf(buf, "%25.16E%25.16E\n", 
		  (*current_node)(0),
		  (*current_node)(1));
      else
	  sprintf(buf, "%25.16E\n", 
		  (*current_node)(0));

      out_file << buf;
    }


  /*
   * Write end of dataset
   */
  out_file << "    -1\n";
}






void UnvMeshInterface::element_out(std::ostream& out_file)
{
  /*
   * we need the MeshData, otherwise we do not
   * know the foreign node id
   */
  if (!_mesh_data.active())
    {
      std::cerr << "ERROR: Need to have an active MeshData for writing UNV files."
		<< std::endl;
      error();
    }


  /*
   * Write beginning of dataset
   */
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_elements
	   << "\n";

  unsigned int elem_n_nodes;                 /* number of nodes on element */
  unsigned long int fe_descriptor_id;        /* FE descriptor id */
  unsigned long int phys_prop_tab_dummy = 2; /* physical property (not supported yet) */
  unsigned long int mat_prop_tab_dummy = 1;  /* material property (not supported yet) */
  unsigned long int color_dummy = 7;         /* color (not supported yet) */


  /*
   * vector that assigns element nodes to their correct position
   * currently only elements with up to 20 nodes
   *
   * Example:
   * QUAD4               | 44:plane stress
   *                     | linear quad
   * position in libMesh | UNV numbering
   * (note: 0-based)     | (note: 1-based)
   *     
   * assign_elm_node[0]  = 1
   * assign_elm_node[3]  = 2
   * assign_elm_node[2]  = 3
   * assign_elm_node[1]  = 4
   */
  unsigned int assign_elm_nodes[20];


  for (unsigned int el_cnt=0;el_cnt<_elements.size();el_cnt++)
    {
      const Elem* elem = _elements[el_cnt];

      elem_n_nodes = elem->n_nodes();


      switch (elem->type())
	{
	
	case TRI3:
	  fe_descriptor_id = 41; // Plane Stress Linear Triangle
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[2] = 2;
	  assign_elm_nodes[1] = 3;
	  break;

	case TRI6:
	  fe_descriptor_id = 42; // Plane Stress Quadratic Triangle
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[5] = 2;
	  assign_elm_nodes[2] = 3;
	  assign_elm_nodes[4] = 4;
	  assign_elm_nodes[1] = 5;
	  assign_elm_nodes[3] = 6;
	  break;

	case QUAD4:
	  fe_descriptor_id = 44; // Plane Stress Linear Quadrilateral
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[3] = 2;
	  assign_elm_nodes[2] = 3;
	  assign_elm_nodes[1] = 4;
	  break;
	
	case QUAD8:
	  fe_descriptor_id = 45; // Plane Stress Quadratic Quadrilateral
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[7] = 2;
	  assign_elm_nodes[3] = 3;
	  assign_elm_nodes[6] = 4;
	  assign_elm_nodes[2] = 5;
	  assign_elm_nodes[5] = 6;
	  assign_elm_nodes[1] = 7;
	  assign_elm_nodes[4] = 8;
	  break;
	
	case TET4:
	  fe_descriptor_id = 111; // Solid Linear Tetrahedron
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[1] = 2;
	  assign_elm_nodes[2] = 3;
	  assign_elm_nodes[3] = 4;
	  break;

	case PRISM6:
	  fe_descriptor_id = 112; // Solid Linear Prism
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[1] = 2;
	  assign_elm_nodes[2] = 3;
	  assign_elm_nodes[3] = 4;
	  assign_elm_nodes[4] = 5;
	  assign_elm_nodes[5] = 6;
	  break;

	case HEX8:
	  fe_descriptor_id = 115; // Solid Linear Brick
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[4] = 2;
	  assign_elm_nodes[5] = 3;
	  assign_elm_nodes[1] = 4;
	  assign_elm_nodes[3] = 5;
	  assign_elm_nodes[7] = 6;
	  assign_elm_nodes[6] = 7;
	  assign_elm_nodes[2] = 8;
	  break;
	
	case HEX20:
	  fe_descriptor_id = 116; // Solid Quadratic Brick
	  assign_elm_nodes[ 0] = 1;
	  assign_elm_nodes[12] = 2;
	  assign_elm_nodes[ 4] = 3;
	  assign_elm_nodes[16] = 4;
	  assign_elm_nodes[ 5] = 5;
	  assign_elm_nodes[13] = 6;
	  assign_elm_nodes[ 1] = 7;
	  assign_elm_nodes[ 8] = 8;
	  assign_elm_nodes[11] = 9;
	  assign_elm_nodes[19] = 10;
	  assign_elm_nodes[17] = 11;
	  assign_elm_nodes[ 9] = 12;
	  assign_elm_nodes[ 3] = 13;
	  assign_elm_nodes[15] = 14;
	  assign_elm_nodes[ 7] = 15;
	  assign_elm_nodes[18] = 16;
	  assign_elm_nodes[ 6] = 17;
	  assign_elm_nodes[14] = 18;
	  assign_elm_nodes[ 2] = 19;
	  assign_elm_nodes[10] = 20;
	  break;
		
	case TET10:
	  fe_descriptor_id = 118; // Solid Parabolic Tetrahedron
	  assign_elm_nodes[0] = 1;
	  assign_elm_nodes[4] = 2;
	  assign_elm_nodes[1] = 3;
	  assign_elm_nodes[5] = 4;
	  assign_elm_nodes[2] = 5;
	  assign_elm_nodes[6] = 6;
	  assign_elm_nodes[7] = 7;
	  assign_elm_nodes[8] = 8;
	  assign_elm_nodes[9] = 9;
	  assign_elm_nodes[3] = 10;
	  break;
	
	
	default:
	  std::cerr << "ERROR: Element type = " 
		    << elem->type() 
		    << "not supported in "
		    << "UnvMeshInterface!"
		    << std::endl;
	  error();	
	  break;
	}


      out_file << std::setw(10) << _mesh_data.elem_to_foreign_id(elem)   /* element ID */
	       << std::setw(10) << fe_descriptor_id                      /* type of element */
	       << std::setw(10) << phys_prop_tab_dummy                   /* not supported */
	       << std::setw(10) << mat_prop_tab_dummy                    /* not supported */
	       << std::setw(10) << color_dummy                           /* not supported */
	       << std::setw(10) << elem_n_nodes << std::endl;            /* No. of nodes per element */


      for (unsigned int j=0;j<elem_n_nodes;j++)
	{
	  /*
	   * assign_elm_nodes[j]-th node: i.e., j loops over the
	   * libMesh numbering, and assign_elm_nodes[j] over the
	   * UNV numbering.
	   */
	  const Node* node_in_unv_order = elem->get_node(assign_elm_nodes[j]-1);

	  // write foreign label for this node
	  out_file << std::setw(10) << _mesh_data.node_to_foreign_id(node_in_unv_order);
	}

      out_file << std::endl;

    }


  /*
   * Write end of dataset
   */
  out_file << "    -1\n";
}











std::string& UnvMeshInterface::D_to_e(std::string& number)
{
  // position of the "D" in the string
  unsigned int position;

  // find "D" in string, start looking at 6th element, to improve speed.
  // We dont expect a "D" earlier
  position = number.find("D",6);

  if(position!=std::string::npos)     // npos means no position
      // replace "D" in string
      number.replace(position,1,"e"); 
  else
      // we assume that if this one number is written correctly, all numbers are
      _need_D_to_e = false;

  return number;

}
