// $Id: mesh_unv_support.C,v 1.26 2003-09-12 03:28:56 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson

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
#include <algorithm>


// Local includes
#include "mesh_unv_support.h"
#include "mesh_data.h"
// #include "elem.h"
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
// MeshBase methods
void MeshBase::read_unv(const std::string& name)
{
  UnvMeshInterface unv_interface (_nodes,
				  _elements,
				  data,
				  true);

  unv_interface.read (name);
}



void MeshBase::write_unv(const std::string& name)
{
  UnvMeshInterface unv_interface (_nodes,
				  _elements,
				  data,
				  true);

  unv_interface.write (name);
}





//-----------------------------------------------------------------------------
// UnvMeshInterface class static members
const std::string  UnvMeshInterface::_label_dataset_nodes       = "2411";
const std::string  UnvMeshInterface::_label_dataset_elements    = "2412";






//-----------------------------------------------------------------------------
// UnvMeshInterface class members
UnvMeshInterface::UnvMeshInterface (std::vector<Node*>& nodes,
				    std::vector<Elem*>& elements,
				    MeshData& md,
				    const bool be_verbose) :
  _verbose       (be_verbose),
  _nodes         (nodes),
  _elements      (elements),
  _n_nodes       (0),
  _n_elements    (0),
  _mesh_data     (md)
{
  // read/write only works in 2D/3D
  assert (DIM > 1);
  this->clear();

  if (_verbose)
    std::cout << std::endl << " Universal file Interface:" << std::endl;
}




UnvMeshInterface::~UnvMeshInterface()
{
}




void UnvMeshInterface::clear ()
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





void UnvMeshInterface::read (const std::string& file_name)
{
  /*
   * clear everything, so that
   * we can start from scratch
   */
  this->clear ();


  /*
   * Note that we read this file
   * @e twice.  First time to
   * detect the number of nodes
   * and elements (and possible
   * conversion tasks like D_to_e)
   * and the order of datasets
   * (nodes first, then elements,
   * or the other way around),
   * and second to do the actual
   * read.
   */
  std::fstream in_stream;

  std::vector<std::string> order_of_datasets;
  order_of_datasets.reserve(2);

  {
    /*
     * the first time we read the file,
     * merely to obtain overall info
     */
    in_stream.open (file_name.c_str(), std::fstream::in);

    if ( !in_stream.good() )
      {
        std::cerr << "ERROR: Input file not good." 
		  << std::endl;
	error();
      }

    
    /*
     * Count nodes and elements, then let 
     * other methods read the element and 
     * node data.  Also remember which
     * dataset comes first: nodes or elements
     */
    if (_verbose)
      std::cout << "  Counting nodes and elements" << std::endl;
    

    bool reached_eof = false;
    bool found_node  = false;
    bool found_elem  = false;


    std::string olds, news;

    while (true)
    {
      in_stream >> olds >> news;

      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while ( ((olds != "-1") || (news == "-1") ) && !in_stream.eof() )
	{	  
	  olds = news;
	  in_stream >> news;
	}

      if (in_stream.eof())
        {
	  break;
	  reached_eof = true;
	}

      /*
       * if beginning of dataset, buffer it in
       * temp_buffer, if desired
       */
      if (news == _label_dataset_nodes)
        {
	  found_node = true;
	  order_of_datasets.push_back (_label_dataset_nodes);
	  this->count_nodes (in_stream);

	  /*
	   * we can save some time scanning the file
	   * when we know we already have everything
	   * we want
	   */
	  if (found_elem)
	    break;
        }

      else if (news == _label_dataset_elements)
        {
	  found_elem = true;
	  order_of_datasets.push_back (_label_dataset_elements);
	  this->count_elements (in_stream);

	  /*
	   * we can save some time scanning the file
	   * when we know we already have everything
	   * we want
	   */
	  if (found_node)
	    break;
        }
    }


    /*
     * Here we should better have found
     * the datasets for nodes and elements,
     * otherwise the unv files is bad!
     */
    if (!found_elem)
      {
	std::cerr << "ERROR: Could not find elements!" << std::endl;
	error();
      }

    if (!found_node)
      {
	std::cerr << "ERROR: Could not find nodes!" << std::endl;
	error();
      }


    in_stream.close();    
  }





  /*
   * We finished scanning the file,
   * and our member data 
   * \p this->_n_nodes,
   * \p this->_n_elements,
   * \p this->_need_D_to_e
   * should be properly initialized.
   */
  {
    in_stream.open (file_name.c_str(), std::fstream::in);


    if ( !in_stream.good() )
      {
        std::cerr << "ERROR: Cannot re-read input file." 
		  << std::endl;
	error();
      }


    /*
     * Read the datasets in the order that
     * we already know
     */
    assert (order_of_datasets.size()==2);

    for (unsigned int ds=0; ds < order_of_datasets.size(); ds++)
      {
	if (order_of_datasets[ds] == _label_dataset_nodes)
	  this->node_in    (in_stream);

	else if (order_of_datasets[ds] == _label_dataset_elements)
	  this->element_in (in_stream);

	else
	  error();
      }


    /*
     * tell the MeshData object that we are finished 
     * reading data
     */
    this->_mesh_data.close_foreign_id_maps ();

    in_stream.close();    

    if (_verbose)
      std::cout << "  Finished." << std::endl << std::endl;
  }



  /* 
   * save memory
   */
  this->_assign_nodes.clear();
  this->_ds_position.clear();
}






void UnvMeshInterface::write (const std::string& file_name)
{
  std::ofstream out_file (file_name.c_str());

  if ( !out_file.good() )
    {
      std::cerr << "ERROR: Output file not good." 
		<< std::endl;
      error();
    }


  /*
   * already know these data, so initialize
   * them.  Does not hurt.
   */
  this->_n_nodes      = this->_nodes.size();
  this->_n_elements   = this->_elements.size();
  this->_need_D_to_e  = false;


  /*
   * we need the MeshData, otherwise we do not
   * know the foreign node id
   */
  if (!_mesh_data.active())
    {
      if (!_mesh_data.compatibility_mode())
        {
	  std::cerr << std::endl
		    << "*************************************************************************" << std::endl
		    << "* WARNING: MeshData neither active nor in compatibility mode.           *" << std::endl
		    << "*          Enable compatibility mode for MeshData.  Use this Universal  *" << std::endl
		    << "*          file with caution: libMesh node and element ids are used.    *" << std::endl
		    << "*************************************************************************" << std::endl
		    << std::endl;
	  _mesh_data.enable_compatibility_mode();
	}
    }


  /*
   * write the nodes,  then the elements
   */
  this->node_out    (out_file);
  this->element_out (out_file);
}





void UnvMeshInterface::count_nodes (std::istream& in_file)
{
  /*
   * if this->_n_nodes is not 0 the dataset 
   * has already been scanned
   */
  if (this->_n_nodes != 0)
    {
      std::cerr << "Error: Trying to scan nodes twice!" 
		<< std::endl;
      error();
    }


  /* 
   * Read from file, count nodes,
   * check if floats have to be converted
   */
  std::string data;

  in_file >> data; // read the first node label


  if (data == "-1")
    {
      std::cerr << "ERROR: Bad, already reached end of dataset before even starting to read nodes!"
		<< std::endl;
      error();
    }

 
  // ignore the misc data for this node
  in_file.ignore(256,'\n');	      


  /*
   * Now we are there to verify whether we need
   * to convert from D to e or not
   */
  in_file >> data;

  /*
   * When this "data" contains a "D", then
   * we have to convert each and every float...
   * But also assume when _this_ specific
   * line does not contain a "D", then the
   * other lines won't, too.
   *
   */
  {
#ifdef __HP_aCC
    // Use an "int" instead of unsigned int,
    // otherwise HP aCC may crash!
    const int position          = data.find("D",6);
#else
    const unsigned int position = data.find("D",6);
#endif

    if(position!=std::string::npos)     // npos means no position
      {
	this->_need_D_to_e = true;
	
	if (_verbose)
	  std::cout << "  Convert from \"D\" to \"e\"" << std::endl;
      }
    else
      this->_need_D_to_e = false;
  }

  /*
   * read the remaining two coordinates
   */
  in_file >> data;
  in_file >> data;


  /*
   * this was our first node
   */
  this->_n_nodes++;



  /*
   * proceed _counting_ the remaining
   * nodes.
   */
  while (!in_file.eof())
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
      error();
    }

  if (_verbose)
    std::cout << "  Nodes   : " << this->_n_nodes << std::endl;
}






void UnvMeshInterface::count_elements (std::istream& in_file)
{
  if (this->_n_elements != 0)
    {
      std::cerr << "Error: Trying to scan elements twice!" 
		<< std::endl;
      error();
    }


  /*
   * Simply read the element
   * dataset for the @e only
   * purpose to count nodes!
   */

  std::string data;

  while (!in_file.eof())
    {
      // read element label
      in_file >> data;
      
      
      // end of dataset?
      if (data == "-1") 
	break;
	
      /*
       * Skip related data,
       * and node number list
       */
      in_file.ignore (256,'\n');
      in_file.ignore (256,'\n');
      
      this->_n_elements++;
    }


  if (in_file.eof())
    {
      std::cerr << "ERROR: File ended before end of element dataset!"
		<< std::endl;
      error();
    }

  if (_verbose)
    std::cout << "  Elements: " << this->_n_elements << std::endl;
}









void UnvMeshInterface::node_in (std::istream& in_file)
{
  if (_verbose)
    std::cout << "  Reading nodes" << std::endl;

  /*
   * adjust the \p istream to our
   * position
   */
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_nodes);

  if (!ok)
    {
      std::cerr << "ERROR: Could not find node dataset!" << std::endl;
      error();
    }

  unsigned int node_lab;           // label of the node
  unsigned int exp_coord_sys_num,  // export coordinate system number       (not supported yet)
               disp_coord_sys_num, // displacement coordinate system number (not supported yet)
               color;              // color                                 (not supported yet)

  /*
   * allocate the correct amount 
   * of memory for the node vector
   */
  this->_nodes.resize         (this->_n_nodes);
  this->_assign_nodes.reserve (this->_n_nodes);


  /**
   * always 3 coordinates in the UNV file, no matter
   * which dimensionality libMesh is in
   */
  std::vector<Real> xyz (3);
      
  /*
   * depending on whether we have to convert each
   * coordinate (float), we offer two versions.
   * Note that \p count_nodes() already verified
   * whether this file uses "D" of "e"
   */
  if (this->_need_D_to_e)
    {
      /*
       * ok, convert...
       */
      std::string num_buf;

      for(unsigned int i=0;i<this->_n_nodes;i++)
        {
	  assert (!in_file.eof());

	  in_file >> node_lab                // read the node label
		  >> exp_coord_sys_num       // (not supported yet)
		  >> disp_coord_sys_num      // (not supported yet)
		  >> color;                  // (not supported yet)

	  /*
	   * take care of the
	   * floating-point data
	   */
	  for (unsigned int d=0; d<3; d++)
	    {
	      in_file >> num_buf;
	      xyz[d] = this->D_to_e (num_buf);
	    }

	  // set up the id map
	  this->_assign_nodes.push_back (node_lab);
	  
	  // add node to the nodes vector.
	  this->_nodes[i] = Node::build (xyz[0],
					 xyz[1],
					 xyz[2],
					 i);

	  // tell the MeshData object the foreign node id
	  this->_mesh_data.add_foreign_node_id (this->_nodes[i], node_lab);
	}
    }

  else

    {
      /*
       * very well, no need to convert anything,
       * just plain import.
       */
      for (unsigned int i=0;i<this->_n_nodes;i++)
        {
	  assert (!in_file.eof());

	  in_file >> node_lab                // read the node label
		  >> exp_coord_sys_num       // (not supported yet)
		  >> disp_coord_sys_num      // (not supported yet)
		  >> color                   // (not supported yet)
		  >> xyz[0]                  // read x-coordinate
		  >> xyz[1]                  // read y-coordinate
		  >> xyz[2];                 // read z-coordinate

	  // set up the id map
	  this->_assign_nodes.push_back (node_lab);
	  
	  // add node to the nodes vector.
	  this->_nodes[i] = Node::build (xyz[0], xyz[1], xyz[2], i);

	  // tell the MeshData object the foreign node id
	  this->_mesh_data.add_foreign_node_id (this->_nodes[i], node_lab);
	}
    }

  // now we need to sort the _assign_nodes vector so we can
  // search it efficiently like a map
  std::sort (this->_assign_nodes.begin(),
	     this->_assign_nodes.end());
}





void UnvMeshInterface::element_in (std::istream& in_file)
{
  if (_verbose)
    std::cout << "  Reading elements" << std::endl;
  

  /*
   * adjust the \p istream to our
   * position
   */
  const bool ok = this->beginning_of_dataset(in_file, _label_dataset_elements);

  if (!ok)
    {
      std::cerr << "ERROR: Could not find element dataset!" << std::endl;
      error();
    }


  unsigned int element_lab,            // element label (not supported yet)
               n_nodes;                // number of nodes on element
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
  std::vector<unsigned int> assign_elem_nodes (21);


  // allocate the correct amount of memory for the vector
  // that holds the elements
  this->_elements.resize (this->_n_elements);
  

  // Get the beginning and end of the _assign_nodes vector
  // to eliminate repeated function calls
  const std::vector<unsigned int>::const_iterator it_begin =
    this->_assign_nodes.begin();
  
  const std::vector<unsigned int>::const_iterator it_end   =
    this->_assign_nodes.end();


      
  // read from the virtual file
  for (unsigned int i=0;i<this->_n_elements;i++)
    {
      in_file >> element_lab             // read element label
	      >> fe_descriptor_id        // read FE descriptor id
	      >> phys_prop_tab_num       // (not supported yet)
	      >> mat_prop_tab_num        // (not supported yet)
	      >> color                   // (not supported yet)
	      >> n_nodes;                // read number of nodes on element

      for (unsigned int j=1;j<=n_nodes;j++)
	in_file >> node_labels[j];             // read node labels


      switch (fe_descriptor_id)
	{
	
	case 41: // Plane Stress Linear Triangle
	case 91: // Thin Shell   Linear Triangle
	  _elements[i]=new Tri3;  // add new element
	  // _elements[i]= Elem::build(TRI3);
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=2;
	  assign_elem_nodes[3]=1;
	  break;
	
	case 42: // Plane Stress Quadratic Triangle
	case 92: // Thin Shell   Quadratic Triangle
	  // _elements[i]= Elem::build(TRI6);
	  _elements[i]=new Tri6;  // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=5;
	  assign_elem_nodes[3]=2;
	  assign_elem_nodes[4]=4;
	  assign_elem_nodes[5]=1;
	  assign_elem_nodes[6]=3;
	  break;
	
	case 43: // Plane Stress Cubic Triangle
	  std::cerr << "ERROR: UNV-element type 43: Plane Stress Cubic Triangle"
		    << " not supported." 
		    << std::endl;
	  error();
	
	  break;
	
	case 44: // Plane Stress Linear Quadrilateral
	case 94: // Thin Shell   Linear Quadrilateral
	  // _elements[i]= Elem::build(QUAD4);
	  _elements[i]=new Quad4; // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=3;
	  assign_elem_nodes[3]=2;
	  assign_elem_nodes[4]=1;
	  break;
	
	case 45: // Plane Stress Quadratic Quadrilateral
	case 95: // Thin Shell   Quadratic Quadrilateral
	  // _elements[i]= Elem::build(QUAD8);
	  _elements[i]=new Quad8; // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=7;
	  assign_elem_nodes[3]=3;
	  assign_elem_nodes[4]=6;
	  assign_elem_nodes[5]=2;
	  assign_elem_nodes[6]=5;
	  assign_elem_nodes[7]=1;
	  assign_elem_nodes[8]=4;
	  break;
	
	case 46: // Plane Stress Cubic Quadrilateral
	  std::cerr << "ERROR: UNV-element type 46: Plane Stress Cubic Quadrilateral"
		    << " not supported." 
		    << std::endl;
	  error();
	  break;
	
	
	case 111: // Solid Linear Tetrahedron
	  // _elements[i]= Elem::build(TET4);
	  _elements[i]=new Tet4;  // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=1;
	  assign_elem_nodes[3]=2;
	  assign_elem_nodes[4]=3;
	  break;


	case 112: // Solid Linear Prism
	  // _elements[i]= Elem::build(PRISM6);
	  _elements[i]=new Prism6;  // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=1;
	  assign_elem_nodes[3]=2;
	  assign_elem_nodes[4]=3;
	  assign_elem_nodes[5]=4;
	  assign_elem_nodes[6]=5;
	  break;
	
	case 115: // Solid Linear Brick
	  // _elements[i]= Elem::build(HEX8);
	  _elements[i]=new Hex8;  // add new element
	  assign_elem_nodes[1]=0;
	  assign_elem_nodes[2]=4;
	  assign_elem_nodes[3]=5;
	  assign_elem_nodes[4]=1;
	  assign_elem_nodes[5]=3;
	  assign_elem_nodes[6]=7;
	  assign_elem_nodes[7]=6;
	  assign_elem_nodes[8]=2;
	  break;
	
	
	case 116: // Solid Quadratic Brick
	  // _elements[i]= Elem::build(HEX20);
	  _elements[i]=new Hex20; // add new element
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
	
	
	case 117: // Solid Cubic Brick
	  std::cerr << "Error: UNV-element type 117: Solid Cubic Brick"
		    << " not supported." 
		    << std::endl;
	  error();
	
	  break;

	
	case 118: // Solid Quadratic Tetrahedron
	  // _elements[i]= Elem::build(TET10);
	  _elements[i]=new Tet10; // add new element
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
	{
	  // Find the position of node_labels[j] in the _assign_nodes vector.
	  const std::pair<std::vector<unsigned int>::const_iterator,
	                  std::vector<unsigned int>::const_iterator>	    
	    it = std::equal_range (it_begin,
				   it_end,
				   node_labels[j]);

	  // it better be there, so assert that it was found.
	  assert (it.first  != it.second);
	  assert (*(it.first) == node_labels[j]);

	  // Now, the distance between this UNV id and the beginning of
	  // the _assign_nodes vector will give us a unique id in the
	  // range [0,n_nodes) that we can use for defining a contiguous
	  // connectivity.
	  const unsigned int assigned_node = std::distance (it_begin,
							    it.first);

	  // Make sure we didn't get an out-of-bounds id
	  assert (assigned_node < this->_nodes.size());
	  
	  this->_elements[i]->set_node(assign_elem_nodes[j]) =
	    this->_nodes[assigned_node];
	}


      // tell the MeshData object the foreign element id
      this->_mesh_data.add_foreign_elem_id (this->_elements[i], element_lab);
    }
}






void UnvMeshInterface::node_out (std::ostream& out_file)
{
  assert (_mesh_data.active() || _mesh_data.compatibility_mode());


  if (_verbose)
    std::cout << "  Writing " << this->_nodes.size() << " nodes" << std::endl;

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


  for(unsigned int i=0;i<this->_nodes.size();i++)
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
  assert (_mesh_data.active() || _mesh_data.compatibility_mode());

  if (_verbose)
    std::cout << "  Writing elements" << std::endl;
  
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
   * assign_elem_node[0]  = 1
   * assign_elem_node[3]  = 2
   * assign_elem_node[2]  = 3
   * assign_elem_node[1]  = 4
   */
  unsigned int assign_elem_nodes[20];

  unsigned int n_elem_written=0;

  const_elem_iterator           elem_it  (std::make_pair(_elements.begin(), _elements.end()));
  const const_elem_iterator     elem_end (std::make_pair(_elements.end(),   _elements.end()));
  
  
  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;
      
      elem_n_nodes = elem->n_nodes();
      

      switch (elem->type())
	{
	
	case TRI3:
	  fe_descriptor_id = 41; // Plane Stress Linear Triangle
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[2] = 2;
	  assign_elem_nodes[1] = 3;
	  break;

	case TRI6:
	  fe_descriptor_id = 42; // Plane Stress Quadratic Triangle
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[5] = 2;
	  assign_elem_nodes[2] = 3;
	  assign_elem_nodes[4] = 4;
	  assign_elem_nodes[1] = 5;
	  assign_elem_nodes[3] = 6;
	  break;

	case QUAD4:
	  fe_descriptor_id = 44; // Plane Stress Linear Quadrilateral
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[3] = 2;
	  assign_elem_nodes[2] = 3;
	  assign_elem_nodes[1] = 4;
	  break;
	
	case QUAD8:
	  fe_descriptor_id = 45; // Plane Stress Quadratic Quadrilateral
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[7] = 2;
	  assign_elem_nodes[3] = 3;
	  assign_elem_nodes[6] = 4;
	  assign_elem_nodes[2] = 5;
	  assign_elem_nodes[5] = 6;
	  assign_elem_nodes[1] = 7;
	  assign_elem_nodes[4] = 8;
	  break;
	
	case TET4:
	  fe_descriptor_id = 111; // Solid Linear Tetrahedron
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[1] = 2;
	  assign_elem_nodes[2] = 3;
	  assign_elem_nodes[3] = 4;
	  break;

	case PRISM6:
	  fe_descriptor_id = 112; // Solid Linear Prism
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[1] = 2;
	  assign_elem_nodes[2] = 3;
	  assign_elem_nodes[3] = 4;
	  assign_elem_nodes[4] = 5;
	  assign_elem_nodes[5] = 6;
	  break;

	case HEX8:
	  fe_descriptor_id = 115; // Solid Linear Brick
	  assign_elem_nodes[0] = 1;
	  assign_elem_nodes[4] = 2;
	  assign_elem_nodes[5] = 3;
	  assign_elem_nodes[1] = 4;
	  assign_elem_nodes[3] = 5;
	  assign_elem_nodes[7] = 6;
	  assign_elem_nodes[6] = 7;
	  assign_elem_nodes[2] = 8;
	  break;
	
	case HEX20:
	  fe_descriptor_id = 116; // Solid Quadratic Brick
	  assign_elem_nodes[ 0] = 1;
	  assign_elem_nodes[ 1] = 13;
	  assign_elem_nodes[ 2] = 5;
	  assign_elem_nodes[ 3] = 17;
	  assign_elem_nodes[ 4] = 6;
	  assign_elem_nodes[ 5] = 14;
	  assign_elem_nodes[ 6] = 2;
	  assign_elem_nodes[ 7] = 9;

	  assign_elem_nodes[ 8] = 12;
	  assign_elem_nodes[ 9] = 20;
	  assign_elem_nodes[10] = 18;
	  assign_elem_nodes[11] = 10;

	  assign_elem_nodes[12] = 4;
	  assign_elem_nodes[13] = 16;
	  assign_elem_nodes[14] = 8;
	  assign_elem_nodes[15] = 19;
	  assign_elem_nodes[16] = 7;
	  assign_elem_nodes[17] = 15;
	  assign_elem_nodes[18] = 3;
	  assign_elem_nodes[19] = 11;
	  break;
		
	case TET10:
	  fe_descriptor_id = 118; // Solid Quadratic Tetrahedron

	  assign_elem_nodes[0] = 1; //0
	  assign_elem_nodes[1] = 5; //4
	  assign_elem_nodes[2] = 2; //1
	  assign_elem_nodes[3] = 6; //5
	  assign_elem_nodes[4] = 3; //2
	  assign_elem_nodes[5] = 7; //6
	  assign_elem_nodes[6] = 8; //7
	  assign_elem_nodes[7] = 9; //8
	  assign_elem_nodes[8] = 10; //9
	  assign_elem_nodes[9] = 4; //3
	  break;
	
	
	default:
	  std::cerr << "ERROR: Element type = " 
		    << elem->type() 
		    << " not supported in "
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
	   * assign_elem_nodes[j]-th node: i.e., j loops over the
	   * libMesh numbering, and assign_elem_nodes[j] over the
	   * UNV numbering.
	   */
	  const Node* node_in_unv_order = elem->get_node(assign_elem_nodes[j]-1);

	  // write foreign label for this node
	  out_file << std::setw(10) << _mesh_data.node_to_foreign_id(node_in_unv_order);
	}

      out_file << std::endl;

      n_elem_written++;
    }

  if (_verbose)
    std::cout << "  Finished writing " << n_elem_written << " elements" << std::endl;

  /*
   * Write end of dataset
   */
  out_file << "    -1\n";
}

