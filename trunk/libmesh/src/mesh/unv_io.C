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
#include <iomanip>
#include <cstdio>   // for std::sprintf
#include <algorithm> // for std::sort
#include <fstream>

// Local includes
#include "libmesh_config.h"
#include "unv_io.h"
#include "mesh_data.h"
#include "mesh_base.h"
#include "face_quad4.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "cell_tet4.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_tet10.h"
#include "cell_prism6.h"

#ifdef HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif





//-----------------------------------------------------------------------------
// UNVIO class static members
const std::string UNVIO::_label_dataset_nodes    = "2411";
const std::string UNVIO::_label_dataset_elements = "2412";



// ------------------------------------------------------------
// UNVIO class members
void UNVIO::clear ()
{
  /*
   * Initialize these to dummy values
   */
  this->_n_nodes     = 0;
  this->_n_elements  = 0;
  this->_need_D_to_e = true;

  this->_assign_nodes.clear();
  this->_ds_position.clear();
}





void UNVIO::read (const std::string& file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef HAVE_GZSTREAM
      
      igzstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);
      
#else
      
      std::cerr << "ERROR:  You must have the zlib.h header "
		<< "files and libraries to read and write "
		<< "compressed streams."
		<< std::endl;
      libmesh_error();
      
#endif
      return;      
    }
  
  else
    {
      std::ifstream in_stream (file_name.c_str());
      this->read_implementation (in_stream);
      return;
    }
}


void UNVIO::read_implementation (std::istream& in_stream)
{
  // clear everything, so that
  // we can start from scratch
  this->clear ();


   // Note that we read this file
   // @e twice.  First time to
   // detect the number of nodes
   // and elements (and possible
   // conversion tasks like D_to_e)
   // and the order of datasets
   // (nodes first, then elements,
   // or the other way around),
   // and second to do the actual
   // read.
  std::vector<std::string> order_of_datasets;
  order_of_datasets.reserve(2);
  
  {
    // the first time we read the file,
    // merely to obtain overall info
    if ( !in_stream.good() )
      {
        std::cerr << "ERROR: Input file not good." 
		  << std::endl;
	libmesh_error();
      }

    
    // Count nodes and elements, then let 
    // other methods read the element and 
    // node data.  Also remember which
    // dataset comes first: nodes or elements
    if (this->verbose())
      std::cout << "  Counting nodes and elements" << std::endl;
    

//    bool reached_eof = false;
    bool found_node  = false;
    bool found_elem  = false;


    std::string olds, news;

    while (in_stream.good())
      {
	in_stream >> olds >> news;
	
	// a "-1" followed by a number means the beginning of a dataset
	// stop combing at the end of the file
	while ( ((olds != "-1") || (news == "-1") ) && !in_stream.eof() )
	  {	  
	    olds = news;
	    in_stream >> news;
	  }

//  	if (in_stream.eof())
//  	  {
//  	    reached_eof = true;
//  	    break;
//  	  }

	
	// if beginning of dataset, buffer it in
	// temp_buffer, if desired
	if (news == _label_dataset_nodes)
	  {
	    found_node = true;
	    order_of_datasets.push_back (_label_dataset_nodes);
	    this->count_nodes (in_stream);

	    // we can save some time scanning the file
	    // when we know we already have everything
	    // we want
	    if (found_elem)
	      break;
	  }

	else if (news == _label_dataset_elements)
	  {
	    found_elem = true;
	    order_of_datasets.push_back (_label_dataset_elements);
	    this->count_elements (in_stream);

	    // we can save some time scanning the file
	    // when we know we already have everything
	    // we want
	    if (found_node)
	      break;
	  }
      }


    // Here we should better have found
    // the datasets for nodes and elements,
    // otherwise the unv files is bad!
    if (!found_elem)
      {
	std::cerr << "ERROR: Could not find elements!" << std::endl;
	libmesh_error();
      }

    if (!found_node)
      {
	std::cerr << "ERROR: Could not find nodes!" << std::endl;
	libmesh_error();
      }


    // Don't close, just seek to the beginning
    in_stream.seekg(0, std::ios::beg);
    
    if (!in_stream.good() )
      {
        std::cerr << "ERROR: Cannot re-read input file." 
		  << std::endl;
	libmesh_error();
      }
  }





  // We finished scanning the file,
  // and our member data 
  // \p this->_n_nodes,
  // \p this->_n_elements,
  // \p this->_need_D_to_e
  // should be properly initialized.
  {
    // Read the datasets in the order that
    // we already know
    libmesh_assert (order_of_datasets.size()==2);

    for (unsigned int ds=0; ds < order_of_datasets.size(); ds++)
      {
	if (order_of_datasets[ds] == _label_dataset_nodes)
	  this->node_in    (in_stream);

	else if (order_of_datasets[ds] == _label_dataset_elements)
	  this->element_in (in_stream);

	else
	  libmesh_error();
      }


    // tell the MeshData object that we are finished 
    // reading data
    this->_mesh_data.close_foreign_id_maps ();

    if (this->verbose())
      std::cout << "  Finished." << std::endl << std::endl;
  }
  
  // save memory
  this->_assign_nodes.clear();
  this->_ds_position.clear();
}





void UNVIO::write (const std::string& file_name)
{
  if (file_name.rfind(".gz") < file_name.size())
    {
#ifdef HAVE_GZSTREAM
      
      ogzstream out_stream(file_name.c_str());
      this->write_implementation (out_stream);
      
#else
      
      std::cerr << "ERROR:  You must have the zlib.h header "
		<< "files and libraries to read and write "
		<< "compressed streams."
		<< std::endl;
      libmesh_error();
      
#endif
      
      return;      
    }
  
  else
    {
      std::ofstream out_stream (file_name.c_str());
      this->write_implementation (out_stream);
      return;
    }
}




void UNVIO::write_implementation (std::ostream& out_file)
{
  if ( !out_file.good() )
    {
      std::cerr << "ERROR: Output file not good." 
		<< std::endl;
      libmesh_error();
    }


  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // already know these data, so initialize
  // them.  Does not hurt.
  this->_n_nodes      = mesh.n_nodes();
  this->_n_elements   = mesh.n_elem();
  this->_need_D_to_e  = false;


 
  // we need the MeshData, otherwise we do not
  // know the foreign node id
  if (!this->_mesh_data.active())
    if (!this->_mesh_data.compatibility_mode())
      {
	std::cerr << std::endl
		  << "*************************************************************************" << std::endl
		  << "* WARNING: MeshData neither active nor in compatibility mode.           *" << std::endl
		  << "*          Enable compatibility mode for MeshData.  Use this Universal  *" << std::endl
		  << "*          file with caution: libMesh node and element ids are used.    *" << std::endl
		  << "*************************************************************************" << std::endl
		  << std::endl;
	this->_mesh_data.enable_compatibility_mode();
      }



  // write the nodes,  then the elements
  this->node_out    (out_file);
  this->element_out (out_file);
}





void UNVIO::count_nodes (std::istream& in_file)
{
  // if this->_n_nodes is not 0 the dataset 
  // has already been scanned
  if (this->_n_nodes != 0)
    {
      std::cerr << "Error: Trying to scan nodes twice!" 
		<< std::endl;
      libmesh_error();
    }


  // Read from file, count nodes,
  // check if floats have to be converted
  std::string data;

  in_file >> data; // read the first node label


  if (data == "-1")
    {
      std::cerr << "ERROR: Bad, already reached end of dataset before even starting to read nodes!"
		<< std::endl;
      libmesh_error();
    }

 
  // ignore the misc data for this node
  in_file.ignore(256,'\n');	      


  
  // Now we are there to verify whether we need
  // to convert from D to e or not
  in_file >> data;

  // When this "data" contains a "D", then
  // we have to convert each and every float...
  // But also assume when _this_ specific
  // line does not contain a "D", then the
  // other lines won't, too.
  {
// #ifdef __HP_aCC
//     // Use an "int" instead of unsigned int,
//     // otherwise HP aCC may crash!
//     const int position          = data.find("D",6);
// #else
//     const unsigned int position = data.find("D",6);
// #endif
    std::string::size_type position = data.find("D",6);

    if (position!=std::string::npos) // npos means no position
      {
	this->_need_D_to_e = true;
	
	if (this->verbose())
	  std::cout << "  Convert from \"D\" to \"e\"" << std::endl;
      }
    else
      this->_need_D_to_e = false;
  }

  // read the remaining two coordinates
  in_file >> data;
  in_file >> data;


  // this was our first node
  this->_n_nodes++;



  // proceed _counting_ the remaining
  // nodes.
  while (in_file.good())
    {
      // read the node label
      in_file >> data; 
	   	
      if (data == "-1")
	// end of dataset is reached
	break;
      
      // ignore the remaining data (coord_sys_label, color etc)
      in_file.ignore (256, '\n');
      // ignore the coordinates
      in_file.ignore (256, '\n');

      this->_n_nodes++;
    }


  if (in_file.eof())
    {
      std::cerr << "ERROR: File ended before end of node dataset!"
		<< std::endl;
      libmesh_error();
    }

  if (this->verbose())
    std::cout << "  Nodes   : " << this->_n_nodes << std::endl;
}






void UNVIO::count_elements (std::istream& in_file)
{
  if (this->_n_elements != 0)
    {
      std::cerr << "Error: Trying to scan elements twice!" 
		<< std::endl;
      libmesh_error();
    }


  // Simply read the element
  // dataset for the @e only
  // purpose to count nodes!
  
  std::string data;
  unsigned int fe_id;

  while (!in_file.eof())
    {
      // read element label
      in_file >> data;
      
      // end of dataset?
      if (data == "-1") 
	break;
	
      // read fe_id
      in_file >> fe_id;

      // Skip related data,
      // and node number list
      in_file.ignore (256,'\n');
      in_file.ignore (256,'\n');

      // For some elements the node numbers
      // are given more than one record

      // TET10 or QUAD9
      if (fe_id == 118 || fe_id == 300)
	  in_file.ignore (256,'\n');

      // HEX20
      if (fe_id == 116)
	{
	  in_file.ignore (256,'\n');
	  in_file.ignore (256,'\n');
	}
     
      this->_n_elements++;
    }


  if (in_file.eof())
    {
      std::cerr << "ERROR: File ended before end of element dataset!"
		<< std::endl;
      libmesh_error();
    }

  if (this->verbose())
    std::cout << "  Elements: " << this->_n_elements << std::endl;
}



void UNVIO::node_in (std::istream& in_file)
{
  if (this->verbose())
    std::cout << "  Reading nodes" << std::endl;

  // adjust the \p istream to our position
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_nodes);

  if (!ok)
    {
      std::cerr << "ERROR: Could not find node dataset!" << std::endl;
      libmesh_error();
    }

  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  unsigned int node_lab;           // label of the node
  unsigned int exp_coord_sys_num,  // export coordinate system number       (not supported yet)
               disp_coord_sys_num, // displacement coordinate system number (not supported yet)
               color;              // color                                 (not supported yet)

  // allocate the correct amount 
  // of memory for the node vector
  this->_assign_nodes.reserve (this->_n_nodes);


  // always 3 coordinates in the UNV file, no matter
  // which dimensionality libMesh is in
  //std::vector<Real> xyz (3);
  Point xyz;
      
  // depending on whether we have to convert each
  // coordinate (float), we offer two versions.
  // Note that \p count_nodes() already verified
  // whether this file uses "D" of "e"
  if (this->_need_D_to_e)
    {
      // ok, convert...
      std::string num_buf;

      for(unsigned int i=0; i<this->_n_nodes; i++)
        {
	  libmesh_assert (!in_file.eof());

	  in_file >> node_lab                // read the node label
		  >> exp_coord_sys_num       // (not supported yet)
		  >> disp_coord_sys_num      // (not supported yet)
		  >> color;                  // (not supported yet)

	  // take care of the
	  // floating-point data
	  for (unsigned int d=0; d<3; d++)
	    {
	      in_file >> num_buf;
	      xyz(d) = this->D_to_e (num_buf);
	    }

	  // set up the id map
	  this->_assign_nodes.push_back (node_lab);
	  
	  // add node to the Mesh & 
	  // tell the MeshData object the foreign node id
	  // (note that mesh.add_point() returns a pointer to the new node)
	  this->_mesh_data.add_foreign_node_id (mesh.add_point(xyz,i), node_lab);
	}
    }

  else
    {
      // very well, no need to convert anything,
      // just plain import.
      for (unsigned int i=0;i<this->_n_nodes;i++)
        {
	  libmesh_assert (!in_file.eof());

	  in_file >> node_lab                // read the node label
		  >> exp_coord_sys_num       // (not supported yet)
		  >> disp_coord_sys_num      // (not supported yet)
		  >> color                   // (not supported yet)
		  >> xyz(0)                  // read x-coordinate
		  >> xyz(1)                  // read y-coordinate
		  >> xyz(2);                 // read z-coordinate

	  // set up the id map
	  this->_assign_nodes.push_back (node_lab);
	  
	  // add node to the Mesh & 
	  // tell the MeshData object the foreign node id
	  // (note that mesh.add_point() returns a pointer to the new node)
	  this->_mesh_data.add_foreign_node_id (mesh.add_point(xyz,i), node_lab);
	}
    }

  // now we need to sort the _assign_nodes vector so we can
  // search it efficiently like a map
  std::sort (this->_assign_nodes.begin(),
	     this->_assign_nodes.end());
}





void UNVIO::element_in (std::istream& in_file)
{
  if (this->verbose())
    std::cout << "  Reading elements" << std::endl;

  MeshBase& mesh = MeshInput<MeshBase>::mesh();

  // adjust the \p istream to our
  // position
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_elements);

  if (!ok)
    {
      std::cerr << "ERROR: Could not find element dataset!" << std::endl;
      libmesh_error();
    }


  unsigned int      element_lab,       // element label (not supported yet)
                    n_nodes;           // number of nodes on element
  unsigned long int fe_descriptor_id,  // FE descriptor id
                    phys_prop_tab_num, // physical property table number (not supported yet)
                    mat_prop_tab_num,  // material property table number (not supported yet)
                    color;             // color (not supported yet)


  // vector that temporarily holds the node labels defining element
  std::vector<unsigned int> node_labels (21);


  // vector that assigns element nodes to their correct position
  // for example:
  // 44:plane stress      | QUAD4
  // linear quadrilateral |
  // position in UNV-file | position in libmesh
  // assign_elem_node[1]   = 0
  // assign_elem_node[2]   = 3
  // assign_elem_node[3]   = 2
  // assign_elem_node[4]   = 1
  //
  // UNV is 1-based, we leave the 0th element of the vectors unused in order
  // to prevent confusion, this way we can store elements with up to 20 nodes
  unsigned int assign_elem_nodes[21];
  

  // Get the beginning and end of the _assign_nodes vector
  // to eliminate repeated function calls
  const std::vector<unsigned int>::const_iterator it_begin =
    this->_assign_nodes.begin();
  
  const std::vector<unsigned int>::const_iterator it_end   =
    this->_assign_nodes.end();


      
  // read from the virtual file
  for (unsigned int i=0; i<this->_n_elements; i++)
    {
      in_file >> element_lab             // read element label
	      >> fe_descriptor_id        // read FE descriptor id
	      >> phys_prop_tab_num       // (not supported yet)
	      >> mat_prop_tab_num        // (not supported yet)
	      >> color                   // (not supported yet)
	      >> n_nodes;                // read number of nodes on element

      for (unsigned int j=1; j<=n_nodes; j++)
	in_file >> node_labels[j];       // read node labels

      Elem* elem = NULL;                 // element pointer
      
      switch (fe_descriptor_id)
	{
	
	case 41: // Plane Stress Linear Triangle
	case 91: // Thin Shell   Linear Triangle
	  {
	    elem = new Tri3;  // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=2;
	    assign_elem_nodes[3]=1;
	    break;
	  }
	
	case 42: // Plane Stress Quadratic Triangle
	case 92: // Thin Shell   Quadratic Triangle
	  {
	    elem = new Tri6;  // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=5;
	    assign_elem_nodes[3]=2;
	    assign_elem_nodes[4]=4;
	    assign_elem_nodes[5]=1;
	    assign_elem_nodes[6]=3;
	    break;
	  }
	
	case 43: // Plane Stress Cubic Triangle
	  {
	    std::cerr << "ERROR: UNV-element type 43: Plane Stress Cubic Triangle"
		      << " not supported." 
		      << std::endl;
	    libmesh_error();	    
	    break;
	  }
	
	case 44: // Plane Stress Linear Quadrilateral
	case 94: // Thin Shell   Linear Quadrilateral
	  {
	    elem = new Quad4; // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=3;
	    assign_elem_nodes[3]=2;
	    assign_elem_nodes[4]=1;
	    break;
	  }
	
	case 45: // Plane Stress Quadratic Quadrilateral
	case 95: // Thin Shell   Quadratic Quadrilateral
	  {
	    elem = new Quad8; // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=7;
	    assign_elem_nodes[3]=3;
	    assign_elem_nodes[4]=6;
	    assign_elem_nodes[5]=2;
	    assign_elem_nodes[6]=5;
	    assign_elem_nodes[7]=1;
	    assign_elem_nodes[8]=4;
	    break;
	  }
	
	case 300: // Thin Shell   Quadratic Quadrilateral (nine nodes)
	  {
	    elem = new Quad9; // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=7;
	    assign_elem_nodes[3]=3;
	    assign_elem_nodes[4]=6;
	    assign_elem_nodes[5]=2;
	    assign_elem_nodes[6]=5;
	    assign_elem_nodes[7]=1;
	    assign_elem_nodes[8]=4;
	    assign_elem_nodes[9]=8;
	    break;
	  }
	
	case 46: // Plane Stress Cubic Quadrilateral
	  {
	    std::cerr << "ERROR: UNV-element type 46: Plane Stress Cubic Quadrilateral"
		      << " not supported." 
		      << std::endl;
	    libmesh_error();
	    break;
	  }
	
	case 111: // Solid Linear Tetrahedron
	  {
	    elem = new Tet4;  // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=1;
	    assign_elem_nodes[3]=2;
	    assign_elem_nodes[4]=3;
	    break;
	  }

	case 112: // Solid Linear Prism
	  {
	    elem = new Prism6;  // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=1;
	    assign_elem_nodes[3]=2;
	    assign_elem_nodes[4]=3;
	    assign_elem_nodes[5]=4;
	    assign_elem_nodes[6]=5;
	    break;
	  }
	
	case 115: // Solid Linear Brick
	  {
	    elem = new Hex8;  // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=4;
	    assign_elem_nodes[3]=5;
	    assign_elem_nodes[4]=1;
	    assign_elem_nodes[5]=3;
	    assign_elem_nodes[6]=7;
	    assign_elem_nodes[7]=6;
	    assign_elem_nodes[8]=2;
	    break;
	  }
		
	case 116: // Solid Quadratic Brick
	  {
	    elem = new Hex20; // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=12;
	    assign_elem_nodes[3]=4;
	    assign_elem_nodes[4]=16;
	    assign_elem_nodes[5]=5;
	    assign_elem_nodes[6]=13;
	    assign_elem_nodes[7]=1;
	    assign_elem_nodes[8]=8;

	    assign_elem_nodes[9]=11;
	    assign_elem_nodes[10]=19;
	    assign_elem_nodes[11]=17;
	    assign_elem_nodes[12]=9;

	    assign_elem_nodes[13]=3;
	    assign_elem_nodes[14]=15;
	    assign_elem_nodes[15]=7;
	    assign_elem_nodes[16]=18;
	    assign_elem_nodes[17]=6;
	    assign_elem_nodes[18]=14;
	    assign_elem_nodes[19]=2;
	    assign_elem_nodes[20]=10;
	    break;
	  }	
	
	case 117: // Solid Cubic Brick
	  {
	    std::cerr << "Error: UNV-element type 117: Solid Cubic Brick"
		      << " not supported." 
		      << std::endl;
	    libmesh_error();	    
	    break;
	  }
	
	case 118: // Solid Quadratic Tetrahedron
	  {
	    elem = new Tet10; // create new element
	    
	    assign_elem_nodes[1]=0;
	    assign_elem_nodes[2]=4;
	    assign_elem_nodes[3]=1;
	    assign_elem_nodes[4]=5;
	    assign_elem_nodes[5]=2;
	    assign_elem_nodes[6]=6;
	    assign_elem_nodes[7]=7;
	    assign_elem_nodes[8]=8;
	    assign_elem_nodes[9]=9;
	    assign_elem_nodes[10]=3;
	    break;
	  }	
	
	default: // Unrecognized element type
	  {
	    std::cerr << "ERROR: UNV-element type " 
		      << fe_descriptor_id
		      << " not supported." 
		      << std::endl;
	    libmesh_error();	    
	    break;
	  }
	}

      // nodes are being stored in element
      for (unsigned int j=1; j<=n_nodes; j++)
	{
	  // Find the position of node_labels[j] in the _assign_nodes vector.
	  const std::pair<std::vector<unsigned int>::const_iterator,
	                  std::vector<unsigned int>::const_iterator>	    
	    it = std::equal_range (it_begin,
				   it_end,
				   node_labels[j]);

	  // it better be there, so libmesh_assert that it was found.
	  libmesh_assert (it.first  != it.second);
	  libmesh_assert (*(it.first) == node_labels[j]);

	  // Now, the distance between this UNV id and the beginning of
	  // the _assign_nodes vector will give us a unique id in the
	  // range [0,n_nodes) that we can use for defining a contiguous
	  // connectivity.
	  const unsigned int assigned_node = std::distance (it_begin,
							    it.first);

	  // Make sure we didn't get an out-of-bounds id
	  libmesh_assert (assigned_node < this->_n_nodes);
	  
	  elem->set_node(assign_elem_nodes[j]) =
	    mesh.node_ptr(assigned_node);
	}
      
      // add elem to the Mesh & 
      // tell the MeshData object the foreign elem id
      // (note that mesh.add_elem() returns a pointer to the new element)
      elem->set_id(i);
      this->_mesh_data.add_foreign_elem_id (mesh.add_elem(elem), element_lab);
    }
}






void UNVIO::node_out (std::ostream& out_file)
{
  
  libmesh_assert (this->_mesh_data.active() ||
	  this->_mesh_data.compatibility_mode());


  if (this->verbose())
    std::cout << "  Writing " << this->_n_nodes << " nodes" << std::endl;

  // Write beginning of dataset
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_nodes
	   << '\n';


  unsigned int exp_coord_sys_dummy  = 0; // export coordinate sys. (not supported yet)
  unsigned int disp_coord_sys_dummy = 0; // displacement coordinate sys. (not supp. yet)
  unsigned int color_dummy          = 0; // color(not supported yet)

  // A reference to the parent class's mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  MeshBase::const_node_iterator       nd  = mesh.nodes_begin();
  const MeshBase::const_node_iterator end = mesh.nodes_end();

  for (; nd != end; ++nd)
    {
      const Node* current_node = *nd;
      
      char buf[78];
      std::sprintf(buf, "%10d%10d%10d%10d\n", 
		   this->_mesh_data.node_to_foreign_id(current_node),
		   exp_coord_sys_dummy,
		   disp_coord_sys_dummy,
		   color_dummy);
      out_file << buf;

      // the coordinates
      if (mesh.spatial_dimension() == 3)
	std::sprintf(buf, "%25.16E%25.16E%25.16E\n", 
		     static_cast<double>((*current_node)(0)),
		     static_cast<double>((*current_node)(1)),
		     static_cast<double>((*current_node)(2)));
      else if (mesh.spatial_dimension() == 2)
	std::sprintf(buf, "%25.16E%25.16E\n", 
		     static_cast<double>((*current_node)(0)),
		     static_cast<double>((*current_node)(1)));
      else
	std::sprintf(buf, "%25.16E\n", 
		     static_cast<double>((*current_node)(0)));
      
      out_file << buf;
    }
  
  
  // Write end of dataset
  out_file << "    -1\n";
}






void UNVIO::element_out(std::ostream& out_file)
{
  libmesh_assert (this->_mesh_data.active() ||
	  this->_mesh_data.compatibility_mode());

  if (this->verbose())
    std::cout << "  Writing elements" << std::endl;
  
  // Write beginning of dataset
  out_file << "    -1\n"
	   << "  " 
	   << _label_dataset_elements
	   << "\n";

  unsigned long int fe_descriptor_id = 0;    // FE descriptor id 
  unsigned long int phys_prop_tab_dummy = 2; // physical property (not supported yet)
  unsigned long int mat_prop_tab_dummy = 1;  // material property (not supported yet)
  unsigned long int color_dummy = 7;         // color (not supported yet)


  // vector that assigns element nodes to their correct position
  // currently only elements with up to 20 nodes
  //  
  // Example:
  // QUAD4               | 44:plane stress
  //                     | linear quad
  // position in libMesh | UNV numbering
  // (note: 0-based)     | (note: 1-based)
  //     
  // assign_elem_node[0]  = 0
  // assign_elem_node[1]  = 3
  // assign_elem_node[2]  = 2
  // assign_elem_node[3]  = 1
  unsigned int assign_elem_nodes[20];

  unsigned int n_elem_written=0;

  // A reference to the parent class's mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();

  MeshBase::const_element_iterator it  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; it != end; ++it)
    {
      const Elem* elem = *it;

      elem->n_nodes();
      
      switch (elem->type())
	{
	
	case TRI3:
	  {
	    fe_descriptor_id = 41; // Plane Stress Linear Triangle
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 2;
	    assign_elem_nodes[2] = 1;
	    break;
	  }

	case TRI6:
	  {
	    fe_descriptor_id = 42; // Plane Stress Quadratic Triangle
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 5;
	    assign_elem_nodes[2] = 2;
	    assign_elem_nodes[3] = 4;
	    assign_elem_nodes[4] = 1;
	    assign_elem_nodes[5] = 3;
	    break;
	  }
	  
	case QUAD4:
	  {
	    fe_descriptor_id = 44; // Plane Stress Linear Quadrilateral
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 3;
	    assign_elem_nodes[2] = 2;
	    assign_elem_nodes[3] = 1;
	    break;
	  }
	
	case QUAD8:
	  {
	    fe_descriptor_id = 45; // Plane Stress Quadratic Quadrilateral
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 7;
	    assign_elem_nodes[2] = 3;
	    assign_elem_nodes[3] = 6;
	    assign_elem_nodes[4] = 2;
	    assign_elem_nodes[5] = 5;
	    assign_elem_nodes[6] = 1;
	    assign_elem_nodes[7] = 4;
	    break;
	  }
	
	case QUAD9:
	  {
	    fe_descriptor_id = 300; // Plane Stress Quadratic Quadrilateral
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 7;
	    assign_elem_nodes[2] = 3;
	    assign_elem_nodes[3] = 6;
	    assign_elem_nodes[4] = 2;
	    assign_elem_nodes[5] = 5;
	    assign_elem_nodes[6] = 1;
	    assign_elem_nodes[7] = 4;
	    assign_elem_nodes[8] = 8;
	    break;
	  }
	
	case TET4:
	  {
	    fe_descriptor_id = 111; // Solid Linear Tetrahedron
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 1;
	    assign_elem_nodes[2] = 2;
	    assign_elem_nodes[3] = 3;
	    break;
	  }

	case PRISM6:
	  {
	    fe_descriptor_id = 112; // Solid Linear Prism
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 1;
	    assign_elem_nodes[2] = 2;
	    assign_elem_nodes[3] = 3;
	    assign_elem_nodes[4] = 4;
	    assign_elem_nodes[5] = 5;
	    break;
	  }

	case HEX8:
	  {
	    fe_descriptor_id = 115; // Solid Linear Brick
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 4;
	    assign_elem_nodes[2] = 5;
	    assign_elem_nodes[3] = 1;
	    assign_elem_nodes[4] = 3;
	    assign_elem_nodes[5] = 7;
	    assign_elem_nodes[6] = 6;
	    assign_elem_nodes[7] = 2;
	    break;
	  }
	
	case HEX20:
	  {
	    fe_descriptor_id = 116; // Solid Quadratic Brick
	    assign_elem_nodes[ 0] = 0;
	    assign_elem_nodes[ 1] = 12;
	    assign_elem_nodes[ 2] = 4;
	    assign_elem_nodes[ 3] = 16;
	    assign_elem_nodes[ 4] = 5;
	    assign_elem_nodes[ 5] = 13;
	    assign_elem_nodes[ 6] = 1;
	    assign_elem_nodes[ 7] = 8;

	    assign_elem_nodes[ 8] = 11;
	    assign_elem_nodes[ 9] = 19;
	    assign_elem_nodes[10] = 17;
	    assign_elem_nodes[11] = 9;

	    assign_elem_nodes[12] = 3;
	    assign_elem_nodes[13] = 15;
	    assign_elem_nodes[14] = 7;
	    assign_elem_nodes[15] = 18;
	    assign_elem_nodes[16] = 6;
	    assign_elem_nodes[17] = 14;
	    assign_elem_nodes[18] = 2;
	    assign_elem_nodes[19] = 10;


	    break;
	  }
		
	case TET10:
	  {
	    fe_descriptor_id = 118; // Solid Quadratic Tetrahedron
	    assign_elem_nodes[0] = 0;
	    assign_elem_nodes[1] = 4;
	    assign_elem_nodes[2] = 1;
	    assign_elem_nodes[3] = 5;
	    assign_elem_nodes[4] = 2;
	    assign_elem_nodes[5] = 6;
	    assign_elem_nodes[6] = 7;
	    assign_elem_nodes[7] = 8;
	    assign_elem_nodes[8] = 9;
	    assign_elem_nodes[9] = 3;
	    break;
	  }
		
	default:
	  {
	    std::cerr << "ERROR: Element type = " 
		      << elem->type() 
		      << " not supported in "
		      << "UNVIO!"
		      << std::endl;
	    libmesh_error();	
	    break;
	  }
	}


      out_file << std::setw(10) << this->_mesh_data.elem_to_foreign_id(elem)  // element ID
	       << std::setw(10) << fe_descriptor_id                           // type of element
	       << std::setw(10) << phys_prop_tab_dummy                        // not supported 
	       << std::setw(10) << mat_prop_tab_dummy                         // not supported 
	       << std::setw(10) << color_dummy                                // not supported 
	       << std::setw(10) << elem->n_nodes()                            // No. of nodes per element
	       << '\n';

      for (unsigned int j=0; j<elem->n_nodes(); j++)
	{
	  // assign_elem_nodes[j]-th node: i.e., j loops over the
	  // libMesh numbering, and assign_elem_nodes[j] over the
	  // UNV numbering.
	  const Node* node_in_unv_order = elem->get_node(assign_elem_nodes[j]);

	  // new record after 8 id entries
	  if (j==8 || j==16)
	    out_file << '\n';

	  // write foreign label for this node
	  out_file << std::setw(10) << this->_mesh_data.node_to_foreign_id(node_in_unv_order);


	}

      out_file << '\n';

      n_elem_written++;
    }

  if (this->verbose())
    std::cout << "  Finished writing " << n_elem_written << " elements" << std::endl;

  // Write end of dataset
  out_file << "    -1\n";
}

