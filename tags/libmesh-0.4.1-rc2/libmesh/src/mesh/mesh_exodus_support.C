// $Id: mesh_exodus_support.C,v 1.9 2003-09-16 15:59:31 benkirk Exp $

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
#include <iomanip>

// Local includes
#include "mesh_base.h"
#include "boundary_mesh.h"
#include "boundary_info.h"
#include "mesh.h"


#ifdef HAVE_EXODUS_API
namespace exII {
extern "C" {
#include "exodusII.h"
}
}

#include "mesh_exodus_support.h"






/**
 * ExodusII Class static member variables
 */

// 2D node map definitions
const int ExodusII::ElementMaps::quad4_node_map[4] = {0, 1, 2, 3};
  
const int ExodusII::ElementMaps::quad8_node_map[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  
const int ExodusII::ElementMaps::quad9_node_map[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  
const int ExodusII::ElementMaps::tri3_node_map[3]  = {0, 1, 2};
  
const int ExodusII::ElementMaps::tri6_node_map[6]  = {0, 1, 2, 3, 4, 5};

// 2D edge map definitions
const int ExodusII::ElementMaps::tri_edge_map[3] = {0, 1, 2};

const int ExodusII::ElementMaps::quad_edge_map[4] = {0, 1, 2, 3};



// 3D node map definitions
const int ExodusII::ElementMaps::hex8_node_map[8]   = {0, 1, 2, 3, 4, 5, 6, 7};
  
const int ExodusII::ElementMaps::hex20_node_map[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
  
const int ExodusII::ElementMaps::hex27_node_map[27] = { 1,  5,  6,  2,  0,  4,  7,  3, 13, 17, 14,  9,  8, 16,
						       18, 10, 12, 19, 15, 11, 24, 25, 22, 26, 21, 23, 20};
  
const int ExodusII::ElementMaps::tet4_node_map[4]   = {0, 1, 2, 3};
  
const int ExodusII::ElementMaps::tet10_node_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  

// 3D face map definitions
const int ExodusII::ElementMaps::tet_face_map[4] =   {1, 2, 3, 0};

const int ExodusII::ElementMaps::hex_face_map[6] =   {1, 2, 3, 4, 0, 5};

const int ExodusII::ElementMaps::hex27_face_map[6] = {1, 0, 3, 5, 4, 2};



/**
 * ExodusII Class member functions
 */

ExodusII::~ExodusII()
{
  delete [] title;
  delete [] elem_type;
}



void ExodusII::check_err(const int err, const std::string msg)
{
  if (err < 0)
    {
      std::cout << msg << std::endl;
      error();
    }
}



void ExodusII::message(const std::string msg)
{
  if (verbose) std::cout << msg << std::endl;
}



void ExodusII::message(const std::string msg, int i)
{
  if (verbose) std::cout << msg << i << "." << std::endl;
}



void ExodusII::open(const char* filename)
{
  ex_id = exII::ex_open(filename,
			EX_READ,
			&comp_ws,
			&io_ws,
			&ex_version);
  
  check_err(ex_id, "Error opening ExodusII mesh file.");
  if (verbose) std::cout << "File opened successfully." << std::endl;
}



void ExodusII::read_header()
{
  ex_err = exII::ex_get_init(ex_id,
			     title,
			     &num_dim,
			     &num_nodes,
			     &num_elem,
			     &num_elem_blk,
			     &num_node_sets,
			     &num_side_sets);

  check_err(ex_err, "Error retrieving header info.");
  message("Exodus header info retrieved successfully.");
}




void ExodusII::print_header()
{
  if (verbose)
    std::cout << "Title: \t" << title << std::endl
	      << "Mesh Dimension: \t"   << num_dim << std::endl
	      << "Number of Nodes: \t" << num_nodes << std::endl
	      << "Number of elements: \t" << num_elem << std::endl
	      << "Number of elt blocks: \t" << num_elem_blk << std::endl
	      << "Number of node sets: \t" << num_node_sets << std::endl
	      << "Number of side sets: \t" << num_side_sets << std::endl;
}



void ExodusII::read_nodes()
{
  x.resize(num_nodes);
  y.resize(num_nodes); 
  z.resize(num_nodes); 

  ex_err = exII::ex_get_coord(ex_id,
			      static_cast<void*>(&x[0]),
			      static_cast<void*>(&y[0]),
			      static_cast<void*>(&z[0]));
  
  check_err(ex_err, "Error retrieving nodal data.");
  message("Nodal data retrieved successfully."); 
}



void ExodusII::print_nodes()
{
  for (int i=0; i<num_nodes; i++)
    {
      std::cout << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
    }
}



void ExodusII::read_block_info()
{
  block_ids.resize(num_elem_blk);
  ex_err = exII::ex_get_elem_blk_ids(ex_id,
				     &block_ids[0]); // Get all element block IDs.
  // Usually, there is only one
  // block since there is only
  // one type of element.
  // However, there could be more.

  check_err(ex_err, "Error getting block IDs.");
  message("All block IDs retrieved successfully."); 
}



void ExodusII::read_elem_in_block(int block)
{
  ex_err = exII::ex_get_elem_block(ex_id,
				   block_ids[block],
				   elem_type,
				   &num_elem_this_blk,
				   &num_nodes_per_elem,
				   &num_attr);
  check_err(ex_err, "Error getting block info.");
  message("Info retrieved successfully for block: ", block); 
  
  
  
  // Read in the connectivity of the elements of this block
  connect.resize(num_nodes_per_elem*num_elem_this_blk);
  ex_err = exII::ex_get_elem_conn(ex_id,
				  block_ids[block],
				  &connect[0]);
  
  check_err(ex_err, "Error reading block connectivity.");
  message("Connectivity retrieved successfully for block: ", block); 
}



void ExodusII::read_sideset_info()
{
  ss_ids.resize(num_side_sets);
  ex_err = exII::ex_get_side_set_ids(ex_id,
				     &ss_ids[0]);
  check_err(ex_err, "Error retrieving sideset information.");
  message("All sideset information retrieved successfully."); 

  // Resize appropriate data structures -- only do this once outside the loop
  num_sides_per_set.resize(num_side_sets);
  num_df_per_set.resize(num_side_sets);
  
  req_info = EX_INQ_SS_ELEM_LEN; // Inquire about the length of the
                                 // concatenated side sets element list
  ex_err = exII::ex_inquire(ex_id,
			    req_info,
			    &ret_int,
			    &ret_float,
			    &ret_char);
  check_err(ex_err, "Error inquiring about side set element list length.");

  //std::cout << "Value returned by ex_inquire was: " << ret_int << std::endl;
  num_elem_all_sidesets = ret_int;
  elem_list.resize (num_elem_all_sidesets);
  side_list.resize (num_elem_all_sidesets);
  id_list.resize   (num_elem_all_sidesets);
}


void ExodusII::read_sideset(int id, int offset)
{
  ex_err = exII::ex_get_side_set_param(ex_id,
				       ss_ids[id],
				       &num_sides_per_set[id],
				       &num_df_per_set[id]);
  check_err(ex_err, "Error retrieving sideset parameters.");
  message("Parameters retrieved successfully for sideset: ", id);

  ex_err = exII::ex_get_side_set(ex_id,
				 ss_ids[id],
				 &elem_list[offset],
				 &side_list[offset]);
  check_err(ex_err, "Error retrieving sideset data.");
  message("Data retrieved successfully for sideset: ", id);

  for (int i=0; i<num_sides_per_set[id]; i++)
    id_list[i+offset] = ss_ids[id];
}



void ExodusII::print_sideset_info()
{
  for (int i=0; i<num_elem_all_sidesets; i++)
    {
      std::cout << elem_list[i] << " " << side_list[i] << std::endl;
    }
}



void ExodusII::close()
{
  ex_err = exII::ex_close(ex_id);
  check_err(ex_err, "Error closing Exodus file.");
  message("Exodus file closed successfully."); 
}


/**
 * Conversion Class member functions
 *
 */
const ExodusII::Conversion ExodusII::ElementMaps::assign_conversion(const std::string type)
{
  if (type == "QUAD4") 
    return assign_conversion(QUAD4);

  else if (type == "QUAD8")
    return assign_conversion(QUAD8);

  else if (type == "QUAD9")
    return assign_conversion(QUAD9);

  else if (type == "TRI3")
    return assign_conversion(TRI3);

  else if (type == "TRI6")
    return assign_conversion(TRI6);

  else if (type == "HEX8")
    return assign_conversion(HEX8);

  else if (type == "HEX20")
    return assign_conversion(HEX20);

  else if (type == "HEX27")
    return assign_conversion(HEX27);

  else if (type == "TETRA4")
    return assign_conversion(TET4);

  else if (type == "TETRA10")
    return assign_conversion(TET10);

  else
    {
      std::cerr << "ERROR! Unrecognized element type: " << type << std::endl;
      error();
    }

  error();
  
  const Conversion conv(tri3_node_map, tri_edge_map, TRI3); // dummy
  return conv;  
}



const ExodusII::Conversion ExodusII::ElementMaps::assign_conversion(const ElemType type)
{
  switch (type)
    {

    case QUAD4:
      {
	const Conversion conv(quad4_node_map, quad_edge_map, QUAD4);
	return conv;
      }

    case QUAD8:
      {
	const Conversion conv(quad8_node_map, quad_edge_map, QUAD8);
	return conv;
      }
      
    case QUAD9:
      {
	const Conversion conv(quad9_node_map, quad_edge_map, QUAD9);
	return conv;
      }
      
    case TRI3:
      {
	const Conversion conv(tri3_node_map, tri_edge_map, TRI3);
	return conv;
      }
      
    case TRI6:
      {
	const Conversion conv(tri6_node_map, tri_edge_map, TRI6);
	return conv;
      }
      
    case HEX8:
      {
	const Conversion conv(hex8_node_map, hex_face_map, HEX8);
	return conv;
      }
      
    case HEX20:
      {
	const Conversion conv(hex20_node_map, hex_face_map, HEX20);
	return conv;
      }
      
    case HEX27:
      {
	const Conversion conv(hex27_node_map, hex27_face_map, HEX27);
	return conv;
      }
      
    case TET4:
      {
	const Conversion conv(tet4_node_map, tet_face_map, TET4);
	return conv;
      }
      
    case TET10:
      {
	const Conversion conv(tet10_node_map, tet_face_map, TET10);
	return conv;
      }
      
    default:
      error();
    }

  error();
  
  const Conversion conv(tri3_node_map, tri_edge_map, TRI3); // dummy
  return conv;  
}

#endif



void Mesh::read_exd(const std::string& name)
{
  // Generate error if API is not defined

  /**
   * Clear any existing mesh data
   */
  clear();

#ifndef HAVE_EXODUS_API
  {
    std::cerr << "ERROR, ExodusII API is not defined." << std::endl
	      << "Input " << name << " cannot be read" << std::endl;
    error();
  }

#else
  
  {    
    // Begin interfacing with the ExodusII data file.

    assert(_dim != 1);      // No support for 1D ExodusII meshes


    int verbose = 0;
    
#ifdef DEBUG    
    verbose = 1;
#endif
    
    ExodusII ex(verbose);     // Instantiate ExodusII interface: if ex(1), verbose=true 
    ExodusII::ElementMaps em; // Instantiate the ElementMaps interface
    
    ex.open(name.c_str());    // Open the exodus file, if possible

    ex.read_header();         // Get header information from exodus file

    ex.print_header();        // Print header information

    assert((unsigned int) ex.get_num_dim() == _dim);  // Be sure number of dimensions
                                                     // is equal to the number of 
                                                     // dimensions in the mesh supplied.

    ex.read_nodes();                      // Read nodes from the exodus file
    _nodes.resize(ex.get_num_nodes());  // Resize the nodes vector

    // Loop over the nodes, create Nodes.
    for (int i=0; i<ex.get_num_nodes(); i++)
      node_ptr(i) = Node::build(ex.get_x(i),
				ex.get_y(i),
				ex.get_z(i),
				i);

    
    // Get information about all the blocks
    ex.read_block_info();
   
    // Read in the element connectivity for each block.
    // Resize the elements vector
    _elements.resize(ex.get_num_elem());
    
    int nelem_last_block = 0;

    // Loop over all the blocks
    for (int i=0; i<ex.get_num_elem_blk(); i++)
      {
	// Read the information for block i
	ex.read_elem_in_block(i);

	// Set any relevant node/edge maps for this element
	const std::string type(ex.get_elem_type());
        const ExodusII::Conversion conv = em.assign_conversion(type); 
      
	// Loop over all the faces in this block
	int jmax = nelem_last_block+ex.get_num_elem_this_blk();
	for (int j=nelem_last_block; j<jmax; j++)
	  {
	    _elements[j] = Elem::build(conv.get_canonical_type());
	    _elements[j]->set_id (j);
	    
	    // Set all the nodes for this element
	    for (int k=0; k<ex.get_num_nodes_per_elem(); k++)
	      {
		int gi = (j-nelem_last_block)*ex.get_num_nodes_per_elem() + conv.get_node_map(k); // global index 
		int node_number      = ex.get_connect(gi);        // Global node number (1-based)
		elem(j)->set_node(k) = node_ptr((node_number-1)); // Set node number
		// Subtract 1 since
		// exodus is internally 1-based
	      }
	  }
	nelem_last_block += ex.get_num_elem_this_blk();
      }

    
    // Read in sideset information -- this is useful for applying boundary conditions
    {
      ex.read_sideset_info(); // Get basic information about ALL sidesets
      int offset=0;
      for (int i=0; i<ex.get_num_side_sets(); i++)
	{
	  offset += (i > 0 ? ex.get_num_sides_per_set(i-1) : 0); // Compute new offset
	  ex.read_sideset(i, offset);
	}
      
      
      //ex.print_sideset_info();
      
      const std::vector<int>& elem_list = ex.get_elem_list();
      const std::vector<int>& side_list = ex.get_side_list();
      const std::vector<int>& id_list   = ex.get_id_list();
      
      for (unsigned int e=0; e<elem_list.size(); e++)
	{
	  // Set any relevant node/edge maps for this element
	  const ExodusII::Conversion conv = em.assign_conversion(elem(elem_list[e]-1)->type());
	  
	  boundary_info.add_side(elem_list[e]-1,
				 conv.get_side_map(side_list[e]-1),
				 id_list[e]);
	}
    }

    
    ex.close();            // Close the exodus file, if possible
  }

#endif
}
