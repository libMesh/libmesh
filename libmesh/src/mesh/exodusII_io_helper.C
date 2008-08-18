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


#include "exodusII_io_helper.h"


#ifdef HAVE_EXODUS_API

#include <cstring>

#include "boundary_info.h"
#include "enum_elem_type.h"
#include "elem.h"
#include "system.h"
#include "numeric_vector.h"

// ------------------------------------------------------------
// ExodusII_IO_Helper::ElementMaps static data

// 2D node map definitions
const int ExodusII_IO_Helper::ElementMaps::quad4_node_map[4] = {0, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::quad8_node_map[8] = {0, 1, 2, 3, 4, 5, 6, 7};
const int ExodusII_IO_Helper::ElementMaps::quad9_node_map[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const int ExodusII_IO_Helper::ElementMaps::tri3_node_map[3]  = {0, 1, 2};
const int ExodusII_IO_Helper::ElementMaps::tri6_node_map[6]  = {0, 1, 2, 3, 4, 5};

// 2D edge map definitions
const int ExodusII_IO_Helper::ElementMaps::tri_edge_map[3] = {0, 1, 2};
const int ExodusII_IO_Helper::ElementMaps::quad_edge_map[4] = {0, 1, 2, 3};

// 3D node map definitions
const int ExodusII_IO_Helper::ElementMaps::hex8_node_map[8]   = {0, 1, 2, 3, 4, 5, 6, 7};
const int ExodusII_IO_Helper::ElementMaps::hex20_node_map[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
							10, 11, 12, 13, 14, 15, 16, 17, 18, 19};  

// Perhaps an older Hex27 node numbering?  This no longer works.
//const int ExodusII_IO_Helper::ElementMaps::hex27_node_map[27] = { 1,  5,  6,  2,  0,  4,  7,  3, 13, 17, 14,  9,  8, 16,
//							  18, 10, 12, 19, 15, 11, 24, 25, 22, 26, 21, 23, 20};

// After trial-and-error, we find a nearly identical mapping with Exodus,
// only 20 and 26 are transposed.
const int ExodusII_IO_Helper::ElementMaps::hex27_node_map[27] = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  // Mid-face nodes and centroid
  26, 21, 22, 23, 24, 25, 20};
//20  21  22  23  24  25  26 // LibMesh indices


const int ExodusII_IO_Helper::ElementMaps::tet4_node_map[4]   = {0, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::tet10_node_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

const int ExodusII_IO_Helper::ElementMaps::prism6_node_map[6]   = {0, 1, 2, 3, 4, 5};
const int ExodusII_IO_Helper::ElementMaps::pyramid5_node_map[5] = {0, 1, 2, 3, 4};
  
// 3D face map definitions
const int ExodusII_IO_Helper::ElementMaps::tet_face_map[4]     = {1, 2, 3, 0};
const int ExodusII_IO_Helper::ElementMaps::hex_face_map[6]     = {1, 2, 3, 4, 0, 5};
const int ExodusII_IO_Helper::ElementMaps::hex27_face_map[6]   = {1, 0, 3, 5, 4, 2};
const int ExodusII_IO_Helper::ElementMaps::prism_face_map[5]   = {-1,-1,-1,-1,-1}; // Not Implemented!
const int ExodusII_IO_Helper::ElementMaps::pyramid_face_map[5] = {-1,-1,-1,-1,-1}; // Not Implemented!


// ------------------------------------------------------------
// ExodusII_IO_Helper class members
ExodusII_IO_Helper::~ExodusII_IO_Helper()
{
  //delete [] title;
  //delete [] elem_type;
}

void ExodusII_IO_Helper::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;
}



void ExodusII_IO_Helper::check_err(const int err, const std::string msg)
{
  if (err < 0)
    {
      std::cout << msg << std::endl;
      libmesh_error();
    }
}



void ExodusII_IO_Helper::message(const std::string msg)
{
  if (_verbose) std::cout << msg << std::endl;
}



void ExodusII_IO_Helper::message(const std::string msg, int i)
{
  if (_verbose) std::cout << msg << i << "." << std::endl;
}



void ExodusII_IO_Helper::open(const char* filename)
{
  ex_id = exII::ex_open(filename,
			EX_READ,
			&comp_ws,
			&io_ws,
			&ex_version);
  
  check_err(ex_id, "Error opening ExodusII mesh file.");
  if (_verbose) std::cout << "File opened successfully." << std::endl;
}



void ExodusII_IO_Helper::read_header()
{
  ex_err = exII::ex_get_init(ex_id,
			     &title[0],
			     &num_dim,
			     &num_nodes,
			     &num_elem,
			     &num_elem_blk,
			     &num_node_sets,
			     &num_side_sets);

  check_err(ex_err, "Error retrieving header info.");

  num_time_steps = inquire(exII::EX_INQ_TIME, "Error retrieving time steps");

  exII::ex_get_var_param(ex_id, "n", &num_nodal_vars);

  message("Exodus header info retrieved successfully.");
}




void ExodusII_IO_Helper::print_header()
{
  if (_verbose)
    std::cout << "Title: \t" << &title[0] << std::endl
	      << "Mesh Dimension: \t"   << num_dim << std::endl
	      << "Number of Nodes: \t" << num_nodes << std::endl
	      << "Number of elements: \t" << num_elem << std::endl
	      << "Number of elt blocks: \t" << num_elem_blk << std::endl
	      << "Number of node sets: \t" << num_node_sets << std::endl
	      << "Number of side sets: \t" << num_side_sets << std::endl;
}



void ExodusII_IO_Helper::read_nodes()
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



void ExodusII_IO_Helper::print_nodes()
{
  for (int i=0; i<num_nodes; i++)
    {
      std::cout << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
    }
}



void ExodusII_IO_Helper::read_block_info()
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



int ExodusII_IO_Helper::get_block_id(int block)
{
  libmesh_assert(static_cast<unsigned int>(block) < block_ids.size());
    
  return block_ids[block];
}



void ExodusII_IO_Helper::read_elem_in_block(int block)
{
  ex_err = exII::ex_get_elem_block(ex_id,
				   block_ids[block],
				   &elem_type[0],
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



void ExodusII_IO_Helper::read_sideset_info()
{
  ss_ids.resize(num_side_sets);
  if (num_side_sets > 0)
    {
      ex_err = exII::ex_get_side_set_ids(ex_id,
					 &ss_ids[0]);
      check_err(ex_err, "Error retrieving sideset information.");
      message("All sideset information retrieved successfully."); 

      // Resize appropriate data structures -- only do this once outside the loop
      num_sides_per_set.resize(num_side_sets);
      num_df_per_set.resize(num_side_sets);

      // Inquire about the length of the
      // concatenated side sets element list
      req_info = exII::EX_INQ_SS_ELEM_LEN; 
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
}


void ExodusII_IO_Helper::read_sideset(int id, int offset)
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



void ExodusII_IO_Helper::print_sideset_info()
{
  for (int i=0; i<num_elem_all_sidesets; i++)
    {
      std::cout << elem_list[i] << " " << side_list[i] << std::endl;
    }
}



void ExodusII_IO_Helper::close()
{
  ex_err = exII::ex_close(ex_id);
  check_err(ex_err, "Error closing Exodus file.");
  message("Exodus file closed successfully."); 
}

int ExodusII_IO_Helper::inquire(int req_info, std::string error_msg)
{
  ex_err = exII::ex_inquire(ex_id,
			    req_info,
			    &ret_int,
			    &ret_float,
			    &ret_char);
    
  check_err(ex_err, error_msg);

  return ret_int;
}

const std::vector<double>& ExodusII_IO_Helper::get_time_steps()
{
  time_steps.resize(num_time_steps);
  exII::ex_get_all_times(ex_id, &time_steps[0]);
  return time_steps;
}



const std::vector<std::string>& ExodusII_IO_Helper::get_nodal_var_names()
{
  // Allocate enough space for our variable name strings.
  nodal_var_names.resize(num_nodal_vars);
    
  // Use the vvc and strings objects to emulate the behavior of
  // a char** object.
  vvc.resize(num_nodal_vars);
  strings.resize(num_nodal_vars);
  for (int i=0;i<num_nodal_vars;i++)
    vvc[i].resize(MAX_STR_LENGTH+1);

  for (int i=0;i<num_nodal_vars;i++)
    strings[i]=&(vvc[i][0]); // set pointer into vvc only *after* all resizing is complete
    
  exII::ex_get_var_names(ex_id,
			 "n",
			 num_nodal_vars,
			 &strings[0]//var_names
			 );

  if (_verbose)
    {
      std::cout << "Read the variable(s) from the file:" << std::endl;
      for (int i=0; i<num_nodal_vars; i++)
	std::cout << "strings[" << i << "]=" << strings[i] << std::endl;
    }

  
  // Copy the char buffers into strings.  
  for (int i=0;i<num_nodal_vars;i++)
    nodal_var_names[i]=strings[i]; // calls string::op=(const char*)

  
  return nodal_var_names;
}




const std::vector<double>& ExodusII_IO_Helper::get_nodal_var_values(std::string nodal_var_name, int time_step)
{
  nodal_var_values.resize(num_nodes);
    
  this->get_nodal_var_names();

  //See if we can find the variable we are looking for
  unsigned int var_index = 0;
  bool found = false;

  found = nodal_var_names[var_index] == nodal_var_name;
    
  while(!found && var_index < nodal_var_names.size())
    {
      var_index++;
      found = nodal_var_names[var_index] == nodal_var_name;
    }

  if(!found)
    {
      std::cerr << "Unable to locate variable named: " << nodal_var_name << std::endl;
      return nodal_var_values;
    }

  exII::ex_get_nodal_var(ex_id, time_step, var_index+1, num_nodes, &nodal_var_values[0]);
  
  return nodal_var_values;
}




// For Writing Solutions

void ExodusII_IO_Helper::create(std::string filename)
{
  //Store things as doubles
  comp_ws = 8;
  io_ws = 8;
    
  ex_id = exII::ex_create(filename.c_str(), EX_CLOBBER, &comp_ws, &io_ws);
    
  ex_id = exII::ex_open(filename.c_str(),
			EX_WRITE,
			&comp_ws,
			&io_ws,
			&ex_version);
  
  check_err(ex_id, "Error creating ExodusII mesh file.");
  if (_verbose) std::cout << "File created successfully." << std::endl;

  _created = true;
}



void ExodusII_IO_Helper::initialize(std::string str_title, const MeshBase & mesh)
{
  num_dim = mesh.mesh_dimension();
  num_nodes = mesh.n_nodes();
  num_elem = mesh.n_elem();
  num_elem_blk = 1;
  num_node_sets = 0;
  num_side_sets = 0;

  ex_err = exII::ex_put_init(ex_id,
			     str_title.c_str(),
			     num_dim,
			     num_nodes,
			     num_elem,
			     num_elem_blk,
			     num_node_sets,
			     num_side_sets);
    
  check_err(ex_err, "Error initializing new Exodus file.");
}

void ExodusII_IO_Helper::write_nodal_coordinates(const MeshBase & mesh)
{
  x.resize(num_nodes);
  y.resize(num_nodes);
  z.resize(num_nodes);

  for (/* unsigned */ int i=0; i<num_nodes; ++i)
    {
      x[i]=(*mesh.node_ptr(i))(0);
      y[i]=(*mesh.node_ptr(i))(1);
      z[i]=(*mesh.node_ptr(i))(2);
    }

  ex_err = exII::ex_put_coord(ex_id, &x[0], &y[0], &z[0]);

  check_err(ex_err, "Error writing coordinates to Exodus file.");
}



// FIXME: This apparently assumes only a single element block, i.e.
// a mesh having only a single type of element.
void ExodusII_IO_Helper::write_elements(const MeshBase & mesh)
{
#ifdef DEBUG
  std::vector<ElemType> et;
  MeshTools::Generation::elem_types(mesh, et);

  if (et.size() > 1)
    {
      std::cerr << "Warning! write_elements() currently assumes\n";
      std::cerr << "a Mesh with a single element type!\n";
      std::cerr << "It will probably fail spectacularly and\n";
      std::cerr << "unpredictably on hybrid meshes!" << std::endl;
    }
#endif
  
  ExodusII_IO_Helper::ElementMaps em;
  const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(0)->type());

  num_nodes_per_elem = mesh.elem(0)->n_nodes();
      
  ex_err = exII::ex_put_elem_block(ex_id, 1, conv.exodus_elem_type().c_str(), num_elem,num_nodes_per_elem, 0);
  check_err(ex_err, "Error writing element block.");

  connect.resize(num_elem*num_nodes_per_elem);

  for (int i=0; i<num_elem; i++)
    {
      Elem * elem = mesh.elem(i);

      for (int j=0; j<num_nodes_per_elem; j++)
	{
	  const unsigned int connect_index   = (i*num_nodes_per_elem)+j;
	  const unsigned int elem_node_index = conv.get_node_map(j);
	  if (_verbose)
	    {
	      std::cout << "Exodus node index: " << j
			<< "=LibMesh node index " << elem_node_index << std::endl;
	    }
	  connect[connect_index] = elem->node(elem_node_index)+1;
	}
    }

  ex_err = exII::ex_put_elem_conn(ex_id, 1, &connect[0]);
  check_err(ex_err, "Error writing element connectivities");
}



void ExodusII_IO_Helper::initialize_nodal_variables(std::vector<std::string> names)
{
  num_nodal_vars = names.size();

  ex_err = exII::ex_put_var_param(ex_id, "n", num_nodal_vars);
  check_err(ex_err, "Error setting number of nodal vars.");

  // Use the vvc and strings objects to emulate the behavior of
  // a char** object.
  vvc.resize(num_nodal_vars);
  strings.resize(num_nodal_vars);

  // For each string in names, allocate enough space in vvc and copy from
  // the C++ string into vvc for passing to the C interface.
  for (int i=0; i<num_nodal_vars; i++)
    {
      vvc[i].resize(names[i].size()+1);
      std::strcpy(&(vvc[i][0]), names[i].c_str());
    }

  for (int i=0; i<num_nodal_vars; i++)
    strings[i] = &(vvc[i][0]); // set pointer into vvc only *after* all resizing is complete

  if (_verbose)
    {
      std::cout << "Writing variable name(s) to file: " << std::endl;
      for (int i=0;i<num_nodal_vars;i++)
	std::cout << "strings[" << i << "]=" << strings[i] << std::endl;
    }
    
  ex_err = exII::ex_put_var_names(ex_id,
				  "n",
				  num_nodal_vars,
				  &strings[0]//var_names
				  );
    
  check_err(ex_err, "Error setting nodal variable names.");
}



void ExodusII_IO_Helper::write_timestep(int timestep, double time)
{
  ex_err = exII::ex_put_time(ex_id, timestep, &time);
  check_err(ex_err, "Error writing timestep.");
}



void ExodusII_IO_Helper::write_nodal_values(int var_id, const std::vector<Number> & values, int timestep)
{
  ex_err = exII::ex_put_nodal_var(ex_id, timestep, var_id, num_nodes, &values[0]);
  check_err(ex_err, "Error writing nodal values.");
}  



bool ExodusII_IO_Helper::created()
{
  return _created;
}



// ------------------------------------------------------------
// ExodusII_IO_Helper::Conversion class members
ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(const std::string type)
{
  if ((type == "QUAD4") || (type == "QUAD"))
    return assign_conversion(QUAD4);

  else if (type == "QUAD8")
    return assign_conversion(QUAD8);

  else if (type == "QUAD9")
    return assign_conversion(QUAD9);

  else if ((type == "TRI3") || (type == "TRIANGLE"))
    return assign_conversion(TRI3);

  else if (type == "TRI6")
    return assign_conversion(TRI6);

  else if ((type == "HEX8") || (type == "HEX") || (type=="hex8")) 
    return assign_conversion(HEX8);

  else if (type == "HEX20")
    return assign_conversion(HEX20);

  else if (type == "HEX27")
    return assign_conversion(HEX27);

  else if ((type == "TETRA4") || (type == "TETRA"))
    return assign_conversion(TET4);

  else if (type == "TETRA10")
    return assign_conversion(TET10);

  else if (type == "WEDGE")
    return assign_conversion(PRISM6);

  else if (type == "PYRAMID" || type == "PYRAMID5")
    return assign_conversion(PYRAMID5);

  else
    {
      std::cerr << "ERROR! Unrecognized element type: " << type << std::endl;
      libmesh_error();
    }

  libmesh_error();
  
  const Conversion conv(tri3_node_map, tri_edge_map, TRI3,"TRI3"); // dummy
  return conv;  
}



ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(const ElemType type)
{
  switch (type)
    {

    case QUAD4:
      {
	const Conversion conv(quad4_node_map, quad_edge_map, QUAD4, "QUAD4");
	return conv;
      }

    case QUAD8:
      {
	const Conversion conv(quad8_node_map, quad_edge_map, QUAD8, "QUAD8");
	return conv;
      }
      
    case QUAD9:
      {
	const Conversion conv(quad9_node_map, quad_edge_map, QUAD9, "QUAD9");
	return conv;
      }
      
    case TRI3:
      {
	const Conversion conv(tri3_node_map, tri_edge_map, TRI3, "TRI3");
	return conv;
      }
      
    case TRI6:
      {
	const Conversion conv(tri6_node_map, tri_edge_map, TRI6, "TRI6");
	return conv;
      }
      
    case HEX8:
      {
	const Conversion conv(hex8_node_map, hex_face_map, HEX8, "HEX8");
	return conv;
      }
      
    case HEX20:
      {
	const Conversion conv(hex20_node_map, hex_face_map, HEX20, "HEX20");
	return conv;
      }
      
    case HEX27:
      {
	const Conversion conv(hex27_node_map, hex27_face_map, HEX27, "HEX27");
	return conv;
      }
      
    case TET4:
      {
	const Conversion conv(tet4_node_map, tet_face_map, TET4, "TETRA4");
	return conv;
      }
      
    case TET10:
      {
	const Conversion conv(tet10_node_map, tet_face_map, TET10, "TETRA10");
	return conv;
      }

    case PRISM6:
      {
	const Conversion conv(prism6_node_map, prism_face_map, PRISM6, "WEDGE");
	return conv;
      }

    case PYRAMID5:
      {
	const Conversion conv(pyramid5_node_map, pyramid_face_map, PYRAMID5, "PYRAMID5");
	return conv;
      }
	
    default:
      libmesh_error();
    }
    
  libmesh_error();
    
  const Conversion conv(tri3_node_map, tri_edge_map, TRI3, "TRI3"); // dummy
  return conv;  
}



#endif // #ifdef HAVE_EXODUS_API
