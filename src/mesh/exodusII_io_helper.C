// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/exodusII_io_helper.h"


#ifdef LIBMESH_HAVE_EXODUS_API

#include <algorithm>
#include <functional>

#include "libmesh/boundary_info.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"

#ifdef DEBUG
#include "libmesh/mesh_tools.h"  // for elem_types warning
#endif

// This macro returns the length of the array a.  Don't
// try using it on empty arrays, since it accesses the
// zero'th element.
#define ARRAY_LENGTH(a) (sizeof((a))/sizeof((a)[0]))

// Anonymous namespace for file local data
namespace
{
  // Define equivalence classes of Cubit/Exodus element types that map to
  // libmesh ElemTypes
  std::map<std::string, libMeshEnums::ElemType> element_equivalence_map;

  // This function initializes the element_equivalence_map the first time it
  // is called, and returns early all other times.
  void init_element_equivalence_map()
  {
    if (element_equivalence_map.empty())
      {
        // EDGE2 equivalences
        element_equivalence_map["EDGE2"]  = EDGE2;
        element_equivalence_map["TRUSS"]  = EDGE2;
        element_equivalence_map["BEAM"]   = EDGE2;
        element_equivalence_map["BAR"]    = EDGE2;
        element_equivalence_map["TRUSS2"] = EDGE2;
        element_equivalence_map["BEAM2"]  = EDGE2;
        element_equivalence_map["BAR2"]   = EDGE2;

        // EDGE3 equivalences
        element_equivalence_map["EDGE3"]  = EDGE3;
        element_equivalence_map["TRUSS3"] = EDGE3;
        element_equivalence_map["BEAM3"]  = EDGE3;
        element_equivalence_map["BAR3"]   = EDGE3;

        // QUAD4 equivalences
        element_equivalence_map["QUAD"]   = QUAD4;
        element_equivalence_map["QUAD4"]  = QUAD4;
        // element_equivalence_map["SHELL"]  = QUAD4;
        // element_equivalence_map["SHELL4"] = QUAD4;

        // QUAD8 equivalences
        element_equivalence_map["QUAD8"]  = QUAD8;
        // element_equivalence_map["SHELL8"] = QUAD8;

        // QUAD9 equivalences
        element_equivalence_map["QUAD9"]  = QUAD9;
        // element_equivalence_map["SHELL9"] = QUAD9;

        // TRI3 equivalences
        element_equivalence_map["TRI"]       = TRI3;
        element_equivalence_map["TRI3"]      = TRI3;
        element_equivalence_map["TRIANGLE"]  = TRI3;
        // element_equivalence_map["TRISHELL"]  = TRI3;
        // element_equivalence_map["TRISHELL3"] = TRI3;

        // TRI6 equivalences
        element_equivalence_map["TRI6"]      = TRI6;
        // element_equivalence_map["TRISHELL6"] = TRI6;

        // HEX8 equivalences
        element_equivalence_map["HEX"]  = HEX8;
        element_equivalence_map["HEX8"] = HEX8;

        // HEX20 equivalences
        element_equivalence_map["HEX20"] = HEX20;

        // HEX27 equivalences
        element_equivalence_map["HEX27"] = HEX27;

        // TET4 equivalences
        element_equivalence_map["TETRA"]  = TET4;
        element_equivalence_map["TETRA4"] = TET4;

        // TET10 equivalences
        element_equivalence_map["TETRA10"] = TET10;

        // PRISM6 equivalences
        element_equivalence_map["WEDGE"] = PRISM6;

        // PRISM15 equivalences
        element_equivalence_map["WEDGE15"] = PRISM15;

        // PRISM18 equivalences
        element_equivalence_map["WEDGE18"] = PRISM18;

        // PYRAMID5 equivalences
        element_equivalence_map["PYRAMID"]  = PYRAMID5;
        element_equivalence_map["PYRAMID5"] = PYRAMID5;
      }
  }
}



namespace libMesh
{

// ------------------------------------------------------------
// ExodusII_IO_Helper::ElementMaps static data

// 1D node map definitions
const int ExodusII_IO_Helper::ElementMaps::edge2_node_map[2] = {0, 1};
const int ExodusII_IO_Helper::ElementMaps::edge3_node_map[3] = {0, 1, 2};

// 1D edge maps
// FIXME: This notion may or may not be defined in ExodusII
const int ExodusII_IO_Helper::ElementMaps::edge_edge_map[2] = {0, 1};
const int ExodusII_IO_Helper::ElementMaps::edge_inverse_edge_map[2] = {1, 2};

// 2D node map definitions
const int ExodusII_IO_Helper::ElementMaps::quad4_node_map[4] = {0, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::quad8_node_map[8] = {0, 1, 2, 3, 4, 5, 6, 7};
const int ExodusII_IO_Helper::ElementMaps::quad9_node_map[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const int ExodusII_IO_Helper::ElementMaps::tri3_node_map[3]  = {0, 1, 2};
const int ExodusII_IO_Helper::ElementMaps::tri6_node_map[6]  = {0, 1, 2, 3, 4, 5};

// 2D edge map definitions
const int ExodusII_IO_Helper::ElementMaps::tri_edge_map[3] = {0, 1, 2};
const int ExodusII_IO_Helper::ElementMaps::quad_edge_map[4] = {0, 1, 2, 3};

//These take a libMesh ID and turn it into an Exodus ID
const int ExodusII_IO_Helper::ElementMaps::tri_inverse_edge_map[3] = {1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::quad_inverse_edge_map[4] = {1, 2, 3, 4};

// 3D node map definitions
const int ExodusII_IO_Helper::ElementMaps::hex8_node_map[8]   = {0, 1, 2, 3, 4, 5, 6, 7};
const int ExodusII_IO_Helper::ElementMaps::hex20_node_map[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
							10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

// Perhaps an older Hex27 node numbering?  This no longer works.
//const int ExodusII_IO_Helper::ElementMaps::hex27_node_map[27] = { 1,  5,  6,  2,  0,  4,  7,  3, 13, 17, 14,  9,  8, 16,
//							  18, 10, 12, 19, 15, 11, 24, 25, 22, 26, 21, 23, 20};

//DRG: Using the newest exodus documentation available on sourceforge and using Table 2 to see
// where all of the nodes over 21 are supposed to go... we come up with:
const int ExodusII_IO_Helper::ElementMaps::hex27_node_map[27] = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  // Mid-face nodes and centroid
  21, 25, 24, 26, 23, 22, 20};
//20  21  22  23  24  25  26 // LibMesh indices

// The hex27 appears to be the only element without a 1:1 map between its
// node numbering and libmesh's.  Therefore when writing out hex27's we need
// to invert this map...
const int ExodusII_IO_Helper::ElementMaps::hex27_inverse_node_map[27] = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  // Mid-face nodes and centroid
  26, 20, 25, 24, 22, 21, 23};
//20  21  22  23  24  25  26


const int ExodusII_IO_Helper::ElementMaps::tet4_node_map[4]   = {0, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::tet10_node_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

const int ExodusII_IO_Helper::ElementMaps::prism6_node_map[6]   = {0, 1, 2, 3, 4, 5};
const int ExodusII_IO_Helper::ElementMaps::prism15_node_map[15]   = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,
								     10, 11, 12, 13, 14};
const int ExodusII_IO_Helper::ElementMaps::prism18_node_map[18]   = {0, 1, 2, 3, 4, 5, 6,  7,  8,  9,
								     10, 11, 12, 13, 14, 15, 16, 17};
const int ExodusII_IO_Helper::ElementMaps::pyramid5_node_map[5] = {0, 1, 2, 3, 4};

// 3D face map definitions
const int ExodusII_IO_Helper::ElementMaps::tet_face_map[4]     = {1, 2, 3, 0};

const int ExodusII_IO_Helper::ElementMaps::hex_face_map[6]     = {1, 2, 3, 4, 0, 5};
const int ExodusII_IO_Helper::ElementMaps::hex27_face_map[6]   = {1, 2, 3, 4, 0, 5};
//const int ExodusII_IO_Helper::ElementMaps::hex27_face_map[6]   = {1, 0, 3, 5, 4, 2};
const int ExodusII_IO_Helper::ElementMaps::prism_face_map[5]   = {1, 2, 3, 0, 4};
const int ExodusII_IO_Helper::ElementMaps::pyramid_face_map[5] = {-1,-1,-1,-1,-1}; // Not Implemented!

//These take a libMesh ID and turn it into an Exodus ID
const int ExodusII_IO_Helper::ElementMaps::tet_inverse_face_map[4]     = {4, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::hex_inverse_face_map[6]     = {5, 1, 2, 3, 4, 6};
const int ExodusII_IO_Helper::ElementMaps::hex27_inverse_face_map[6]   = {5, 1, 2, 3, 4, 6};
//const int ExodusII_IO_Helper::ElementMaps::hex27_inverse_face_map[6]   = {2, 1, 6, 3, 5, 4};
const int ExodusII_IO_Helper::ElementMaps::prism_inverse_face_map[5]   = {4, 1, 2, 3, 5};
const int ExodusII_IO_Helper::ElementMaps::pyramid_inverse_face_map[5] = {-1,-1,-1,-1,-1}; // Not Implemented!


// ------------------------------------------------------------
// ExodusII_IO_Helper class members

  ExodusII_IO_Helper::ExodusII_IO_Helper(const ParallelObject &parent,
                                         bool v,
                                         bool run_only_on_proc0) :
    ParallelObject(parent),
    ex_id(0),
    ex_err(0),
    num_dim(0),
    num_globals(0),
    num_nodes(0),
    num_elem(0),
    num_elem_blk(0),
    num_node_sets(0),
    num_side_sets(0),
    num_elem_this_blk(0),
    num_nodes_per_elem(0),
    num_attr(0),
    num_elem_all_sidesets(0),
    num_time_steps(0),
    num_nodal_vars(0),
    num_elem_vars(0),
    _created(false),
    _opened(false),
    _verbose(v),
    _run_only_on_proc0(run_only_on_proc0),
    _elem_vars_initialized(false),
    _global_vars_initialized(false),
    _use_mesh_dimension_instead_of_spatial_dimension(false)
  {
    title.resize(MAX_LINE_LENGTH+1);
    elem_type.resize(MAX_STR_LENGTH);
  }



ExodusII_IO_Helper::~ExodusII_IO_Helper()
{
}



void ExodusII_IO_Helper::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;
}



const char* ExodusII_IO_Helper::get_elem_type() const
{
  return &elem_type[0];
}



void ExodusII_IO_Helper::message(const std::string msg)
{
  if (_verbose) libMesh::out << msg << std::endl;
}



void ExodusII_IO_Helper::message(const std::string msg, int i)
{
  if (_verbose) libMesh::out << msg << i << "." << std::endl;
}



void ExodusII_IO_Helper::open(const char* filename)
{
  // Version of Exodus you are using
  float ex_version = 0.;

  // Word size in bytes of the floating point variables used in the
  // application program (0, 4, or 8)
  int comp_ws = sizeof(Real);

  // Word size in bytes of the floating point data as they are stored
  // in the ExodusII file
  int io_ws = 0;

  ex_id = exII::ex_open(filename,
			EX_READ,
			&comp_ws,
			&io_ws,
			&ex_version);

  std::string err_msg = std::string("Error opening ExodusII mesh file: ") + std::string(filename);
  EX_CHECK_ERR(ex_id, err_msg);
  if (_verbose) libMesh::out << "File opened successfully." << std::endl;

  _opened = true;
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

  EX_CHECK_ERR(ex_err, "Error retrieving header info.");

  num_time_steps = inquire(exII::EX_INQ_TIME, "Error retrieving time steps");

  ex_err = exII::ex_get_var_param(ex_id, "n", &num_nodal_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of nodal variables.");

  ex_err = exII::ex_get_var_param(ex_id, "e", &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  message("Exodus header info retrieved successfully.");
}




void ExodusII_IO_Helper::print_header()
{
  if (_verbose)
    libMesh::out << "Title: \t" << &title[0] << std::endl
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

  EX_CHECK_ERR(ex_err, "Error retrieving nodal data.");
  message("Nodal data retrieved successfully.");
}



void ExodusII_IO_Helper::read_node_num_map ()
{
  node_num_map.resize(num_nodes);

  ex_err = exII::ex_get_node_num_map (ex_id,
				      node_num_map.empty() ? NULL : &node_num_map[0]);

  EX_CHECK_ERR(ex_err, "Error retrieving nodal number map.");
  message("Nodal numbering map retrieved successfully.");

  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] node_num_map[i] = ";
      for (unsigned int i=0; i<static_cast<unsigned int>(std::min(10, num_nodes-1)); ++i)
        libMesh::out << node_num_map[i] << ", ";
      libMesh::out << "... " << node_num_map.back() << std::endl;
    }
}


void ExodusII_IO_Helper::print_nodes(std::ostream &out_stream)
{
  for (int i=0; i<num_nodes; i++)
    out_stream << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
}



void ExodusII_IO_Helper::read_block_info()
{
  block_ids.resize(num_elem_blk);
  // Get all element block IDs.
  ex_err = exII::ex_get_elem_blk_ids(ex_id,
				     block_ids.empty() ? NULL : &block_ids[0]);
  // Usually, there is only one
  // block since there is only
  // one type of element.
  // However, there could be more.

  EX_CHECK_ERR(ex_err, "Error getting block IDs.");
  message("All block IDs retrieved successfully.");

  char name_buffer[MAX_STR_LENGTH+1];
  for (int i=0; i<num_elem_blk; ++i)
  {
    ex_err = exII::ex_get_name(ex_id, exII::EX_ELEM_BLOCK,
                               block_ids[i], name_buffer);
    EX_CHECK_ERR(ex_err, "Error getting block name.");
    id_to_block_names[block_ids[i]] = name_buffer;
  }
  message("All block names retrieved successfully.");
}



int ExodusII_IO_Helper::get_block_id(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), block_ids.size());

  return block_ids[index];
}



std::string ExodusII_IO_Helper::get_block_name(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), block_ids.size());

  return id_to_block_names[block_ids[index]];
}



int ExodusII_IO_Helper::get_side_set_id(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), ss_ids.size());

  return ss_ids[index];
}



std::string ExodusII_IO_Helper::get_side_set_name(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), ss_ids.size());

  return id_to_ss_names[ss_ids[index]];
}



int ExodusII_IO_Helper::get_node_set_id(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), nodeset_ids.size());

  return nodeset_ids[index];
}



std::string ExodusII_IO_Helper::get_node_set_name(int index)
{
  libmesh_assert_less (static_cast<unsigned int>(index), nodeset_ids.size());

  return id_to_ns_names[nodeset_ids[index]];
}




void ExodusII_IO_Helper::read_elem_in_block(int block)
{
  libmesh_assert_less (static_cast<unsigned int>(block), block_ids.size());

  ex_err = exII::ex_get_elem_block(ex_id,
				   block_ids[block],
				   &elem_type[0],
				   &num_elem_this_blk,
				   &num_nodes_per_elem,
				   &num_attr);
  if (_verbose)
    libMesh::out << "Reading a block of " << num_elem_this_blk
	      << " " << &elem_type[0] << "(s)"
	      << " having " << num_nodes_per_elem
	      << " nodes per element." << std::endl;

  EX_CHECK_ERR(ex_err, "Error getting block info.");
  message("Info retrieved successfully for block: ", block);



  // Read in the connectivity of the elements of this block,
  // watching out for the case where we actually have no
  // elements in this block (possible with parallel files)
  connect.resize(num_nodes_per_elem*num_elem_this_blk);

  if (!connect.empty())
    {
      ex_err = exII::ex_get_elem_conn(ex_id,
				      block_ids[block],
				      &connect[0]);

      EX_CHECK_ERR(ex_err, "Error reading block connectivity.");
      message("Connectivity retrieved successfully for block: ", block);
    }
}




void ExodusII_IO_Helper::read_elem_num_map ()
{
  elem_num_map.resize(num_elem);

  ex_err = exII::ex_get_elem_num_map (ex_id,
				      elem_num_map.empty() ? NULL : &elem_num_map[0]);

  EX_CHECK_ERR(ex_err, "Error retrieving element number map.");
  message("Element numbering map retrieved successfully.");


  if (_verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] elem_num_map[i] = ";
      for (unsigned int i=0; i<static_cast<unsigned int>(std::min(10, num_elem-1)); ++i)
        libMesh::out << elem_num_map[i] << ", ";
      libMesh::out << "... " << elem_num_map.back() << std::endl;
    }
}



void ExodusII_IO_Helper::read_sideset_info()
{
  ss_ids.resize(num_side_sets);
  if (num_side_sets > 0)
    {
      ex_err = exII::ex_get_side_set_ids(ex_id,
					 &ss_ids[0]);
      EX_CHECK_ERR(ex_err, "Error retrieving sideset information.");
      message("All sideset information retrieved successfully.");

      // Resize appropriate data structures -- only do this once outside the loop
      num_sides_per_set.resize(num_side_sets);
      num_df_per_set.resize(num_side_sets);

      // Inquire about the length of the concatenated side sets element list
      num_elem_all_sidesets = inquire(exII::EX_INQ_SS_ELEM_LEN, "Error retrieving length of the concatenated side sets element list!");

      elem_list.resize (num_elem_all_sidesets);
      side_list.resize (num_elem_all_sidesets);
      id_list.resize   (num_elem_all_sidesets);
    }

  char name_buffer[MAX_STR_LENGTH+1];
  for (int i=0; i<num_side_sets; ++i)
  {
    ex_err = exII::ex_get_name(ex_id, exII::EX_SIDE_SET,
                               ss_ids[i], name_buffer);
    EX_CHECK_ERR(ex_err, "Error getting side set name.");
    id_to_ss_names[ss_ids[i]] = name_buffer;
  }
  message("All side set names retrieved successfully.");
}

void ExodusII_IO_Helper::read_nodeset_info()
{
  nodeset_ids.resize(num_node_sets);
  if (num_node_sets > 0)
    {
      ex_err = exII::ex_get_node_set_ids(ex_id,
					 &nodeset_ids[0]);
      EX_CHECK_ERR(ex_err, "Error retrieving nodeset information.");
      message("All nodeset information retrieved successfully.");

      // Resize appropriate data structures -- only do this once outnode the loop
      num_nodes_per_set.resize(num_node_sets);
      num_node_df_per_set.resize(num_node_sets);
    }

  char name_buffer[MAX_STR_LENGTH+1];
  for (int i=0; i<num_node_sets; ++i)
  {
    ex_err = exII::ex_get_name(ex_id, exII::EX_NODE_SET,
                               nodeset_ids[i], name_buffer);
    EX_CHECK_ERR(ex_err, "Error getting node set name.");
    id_to_ns_names[nodeset_ids[i]] = name_buffer;
  }
  message("All node set names retrieved successfully.");
}



void ExodusII_IO_Helper::read_sideset(int id, int offset)
{
  libmesh_assert_less (static_cast<unsigned int>(id), ss_ids.size());
  libmesh_assert_less (static_cast<unsigned int>(id), num_sides_per_set.size());
  libmesh_assert_less (static_cast<unsigned int>(id), num_df_per_set.size());
  libmesh_assert_less_equal (static_cast<unsigned int>(offset), elem_list.size());
  libmesh_assert_less_equal (static_cast<unsigned int>(offset), side_list.size());

  ex_err = exII::ex_get_side_set_param(ex_id,
				       ss_ids[id],
				       &num_sides_per_set[id],
				       &num_df_per_set[id]);
  EX_CHECK_ERR(ex_err, "Error retrieving sideset parameters.");
  message("Parameters retrieved successfully for sideset: ", id);


  // It's OK for offset==elem_list.size() as long as num_sides_per_set[id]==0
  // because in that case we don't actually read anything...
#ifdef DEBUG
  if (static_cast<unsigned int>(offset) == elem_list.size() ||
      static_cast<unsigned int>(offset) == side_list.size() )
    libmesh_assert_equal_to (num_sides_per_set[id], 0);
#endif


  // Don't call ex_get_side_set unless there are actually sides there to get.
  // Exodus prints an annoying warning in DEBUG mode otherwise...
  if (num_sides_per_set[id] > 0)
    {
      ex_err = exII::ex_get_side_set(ex_id,
				     ss_ids[id],
				     &elem_list[offset],
				     &side_list[offset]);
      EX_CHECK_ERR(ex_err, "Error retrieving sideset data.");
      message("Data retrieved successfully for sideset: ", id);

      for (int i=0; i<num_sides_per_set[id]; i++)
	id_list[i+offset] = ss_ids[id];
    }
}



void ExodusII_IO_Helper::read_nodeset(int id)
{
  libmesh_assert_less (static_cast<unsigned int>(id), nodeset_ids.size());
  libmesh_assert_less (static_cast<unsigned int>(id), num_nodes_per_set.size());
  libmesh_assert_less (static_cast<unsigned int>(id), num_node_df_per_set.size());

  ex_err = exII::ex_get_node_set_param(ex_id,
				       nodeset_ids[id],
				       &num_nodes_per_set[id],
				       &num_node_df_per_set[id]);
  EX_CHECK_ERR(ex_err, "Error retrieving nodeset parameters.");
  message("Parameters retrieved successfully for nodeset: ", id);

  node_list.resize(num_nodes_per_set[id]);

  // Don't call ex_get_node_set unless there are actually nodes there to get.
  // Exodus prints an annoying warning message in DEBUG mode otherwise...
  if (num_nodes_per_set[id] > 0)
    {
      ex_err = exII::ex_get_node_set(ex_id,
				     nodeset_ids[id],
				     &node_list[0]);

      EX_CHECK_ERR(ex_err, "Error retrieving nodeset data.");
      message("Data retrieved successfully for nodeset: ", id);
    }
}



void ExodusII_IO_Helper::close()
{
  // Always call close on processor 0.
  // If we're running on multiple processors, i.e. as one of several Nemesis files,
  // we call close on all processors...
  if ((this->processor_id() == 0) || (!_run_only_on_proc0))
    {
      ex_err = exII::ex_close(ex_id);
      EX_CHECK_ERR(ex_err, "Error closing Exodus file.");
      message("Exodus file closed successfully.");
    }
}



int ExodusII_IO_Helper::inquire(int req_info_in, std::string error_msg)
{
  int ret_int = 0;
  char ret_char = 0;
  float ret_float = 0.;

  ex_err = exII::ex_inquire(ex_id,
			    req_info_in,
			    &ret_int,
			    &ret_float,
			    &ret_char);

  EX_CHECK_ERR(ex_err, error_msg);

  return ret_int;
}



void ExodusII_IO_Helper::read_time_steps()
{
  if (num_time_steps > 0)
    {
      time_steps.resize(num_time_steps);
      ex_err = exII::ex_get_all_times(ex_id, &time_steps[0]);
      EX_CHECK_ERR(ex_err, "Error reading timesteps!");
    }
}



void ExodusII_IO_Helper::read_nodal_var_names()
{
  NamesData names_table(num_nodal_vars, MAX_STR_LENGTH);

  ex_err = exII::ex_get_var_names(ex_id,
                                  "n",
                                  num_nodal_vars,
                                  names_table.get_char_star_star()
                                  );
  EX_CHECK_ERR(ex_err, "Error reading nodal variable names!");


  if (_verbose)
    {
      libMesh::out << "Read the variable(s) from the file:" << std::endl;
      for (int i=0; i<num_nodal_vars; i++)
	libMesh::out << names_table.get_char_star(i) << std::endl;
    }

  // Allocate enough space for our variable name strings.
  nodal_var_names.resize(num_nodal_vars);

  // Copy the char buffers into strings.
  for (int i=0; i<num_nodal_vars; i++)
    nodal_var_names[i] = names_table.get_char_star(i); // calls string::op=(const char*)
}




void ExodusII_IO_Helper::read_nodal_var_values(std::string nodal_var_name, int time_step)
{
  // Read the nodal variable names from file, so we can see if we have the one we're looking for
  this->read_nodal_var_names();

  // See if we can find the variable we are looking for
  unsigned int var_index = 0;
  bool found = false;

  // Do a linear search for nodal_var_name in nodal_var_names
  for (; var_index<nodal_var_names.size(); ++var_index)
    {
      found = (nodal_var_names[var_index] == nodal_var_name);
      if (found)
        break;
    }

  if (!found)
    {
      libMesh::err << "Unable to locate variable named: " << nodal_var_name << std::endl;
      libMesh::err << "Available variables: " << std::endl;
      for (unsigned int i=0; i<nodal_var_names.size(); ++i)
        libMesh::err << nodal_var_names[i] << std::endl;

      libmesh_error();
    }

  // Allocate enough space to store the nodal variable values
  nodal_var_values.resize(num_nodes);

  // Call the Exodus API to read the nodal variable values
  ex_err = exII::ex_get_nodal_var(ex_id,
                                  time_step,
                                  var_index+1,
                                  num_nodes,
                                  &nodal_var_values[0]);
  EX_CHECK_ERR(ex_err, "Error reading nodal variable values!");
}



void ExodusII_IO_Helper::read_elemental_var_names()
{
  NamesData names_table(num_elem_vars, MAX_STR_LENGTH);

  ex_err = exII::ex_get_var_names(ex_id,
                                  "e",
                                  num_elem_vars,
                                  names_table.get_char_star_star()
                                  );
  EX_CHECK_ERR(ex_err, "Error reading elemental variable names!");

  if (_verbose)
    {
      libMesh::out << "Read the variable(s) from the file:" << std::endl;
      for (int i=0; i<num_elem_vars; i++)
        libMesh::out << names_table.get_char_star(i) << std::endl;
    }

  // Allocate enough space for our variable name strings.
  elem_var_names.resize(num_elem_vars);

  // Copy the char buffers into strings.
  for (int i=0; i<num_elem_vars; i++)
    elem_var_names[i] = names_table.get_char_star(i); // calls string::op=(const char*)
}



void ExodusII_IO_Helper::read_elemental_var_values(std::string elemental_var_name, int time_step)
{
  // CAUTION: this assumes that libMesh element numbering is identical to exodus block-by-block element numbering
  // There is no way how to get the whole elemental field from the exodus file, so we have to go block by block

  elem_var_values.resize(num_elem);

  this->read_elemental_var_names();

  // See if we can find the variable we are looking for
  unsigned int var_index = 0;
  bool found = false;

  // Do a linear search for nodal_var_name in nodal_var_names
  for (; var_index<elem_var_names.size(); ++var_index)
    {
      found = (elem_var_names[var_index] == elemental_var_name);
      if (found)
        break;
    }

  if (!found)
    {
      libMesh::err << "Unable to locate variable named: " << elemental_var_name << std::endl;
      libMesh::err << "Available variables: " << std::endl;
      for (unsigned int i=0; i<elem_var_names.size(); ++i)
        libMesh::err << elem_var_names[i] << std::endl;

      libmesh_error();
    }

  unsigned int ex_el_num = 0;
  for (unsigned int i=0; i<static_cast<unsigned int>(num_elem_blk); i++)
    {
      int n_blk_elems = 0;
      ex_err = exII::ex_get_elem_block(ex_id,
                                       block_ids[i],
                                       NULL,
                                       &n_blk_elems,
                                       NULL,
                                       NULL);
      EX_CHECK_ERR(ex_err, "Error getting number of elements in block.");

      std::vector<Real> block_elem_var_values(num_elem);
      ex_err = exII::ex_get_elem_var(ex_id,
                                     time_step,
                                     var_index+1,
                                     block_ids[i],
                                     n_blk_elems,
                                     &block_elem_var_values[0]);
      EX_CHECK_ERR(ex_err, "Error getting elemental values.");

      for (unsigned int j=0; j<static_cast<unsigned int>(n_blk_elems); j++)
        {
          elem_var_values[ex_el_num] = block_elem_var_values[j];
          ex_el_num++;
        }
    }
}


// For Writing Solutions

void ExodusII_IO_Helper::create(std::string filename)
{
  // If we're processor 0, always create the file.
  // If we running on all procs, e.g. as one of several Nemesis files, also
  // call create there.
  if ((this->processor_id() == 0) || (!_run_only_on_proc0))
    {
      // Fall back on double precision when necessary since ExodusII
      // doesn't seem to support long double
      int comp_ws = std::min(sizeof(Real), sizeof(double));
      int io_ws = std::min(sizeof(Real), sizeof(double));

      ex_id = exII::ex_create(filename.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

      EX_CHECK_ERR(ex_id, "Error creating ExodusII mesh file.");

      if (_verbose)
        libMesh::out << "File created successfully." << std::endl;
    }

  _created = true;
}




void ExodusII_IO_Helper::initialize_discontinuous(std::string str_title, const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_use_mesh_dimension_instead_of_spatial_dimension)
    num_dim = mesh.mesh_dimension();
  else
    num_dim = mesh.spatial_dimension();

  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for (; it!=end; ++it)
	num_nodes += (*it)->n_nodes();

  num_elem = mesh.n_elem();

  std::vector<boundary_id_type> unique_side_boundaries;
  std::vector<boundary_id_type> unique_node_boundaries;

  mesh.boundary_info->build_side_boundary_ids(unique_side_boundaries);
  mesh.boundary_info->build_node_boundary_ids(unique_node_boundaries);

  num_side_sets = unique_side_boundaries.size();
  num_node_sets = unique_node_boundaries.size();

  //loop through element and map between block and element vector
  std::map<subdomain_id_type, std::vector<unsigned int>  > subdomain_map;

  for (it=mesh.active_elements_begin(); it!=end; ++it)
    {
      const Elem * elem = *it;
      subdomain_id_type cur_subdomain = elem->subdomain_id();

      subdomain_map[cur_subdomain].push_back(elem->id());
    }
  num_elem_blk = subdomain_map.size();

  if (str_title.size() > MAX_LINE_LENGTH)
    {
      libMesh::err << "Warning, Exodus files cannot have titles longer than "
		   << MAX_LINE_LENGTH
		   << " characters.  Your title will be truncated."
		   << std::endl;
      str_title.resize(MAX_LINE_LENGTH);
    }

  ex_err = exII::ex_put_init(ex_id,
			     str_title.c_str(),
			     num_dim,
			     num_nodes,
			     num_elem,
			     num_elem_blk,
			     num_node_sets,
			     num_side_sets);

  EX_CHECK_ERR(ex_err, "Error initializing new Exodus file.");
}



void ExodusII_IO_Helper::initialize(std::string str_title, const MeshBase & mesh)
{
  // n_active_elem() is a parallel_only function
  unsigned int n_active_elem = mesh.n_active_elem();

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_use_mesh_dimension_instead_of_spatial_dimension)
    num_dim = mesh.mesh_dimension();
  else
    num_dim = mesh.spatial_dimension();

  num_nodes = mesh.n_nodes();
  num_elem = mesh.n_elem();

  std::vector<boundary_id_type> unique_side_boundaries;
  std::vector<boundary_id_type> unique_node_boundaries;

  mesh.boundary_info->build_side_boundary_ids(unique_side_boundaries);
  mesh.boundary_info->build_node_boundary_ids(unique_node_boundaries);

  num_side_sets = unique_side_boundaries.size();
  num_node_sets = unique_node_boundaries.size();

  //loop through element and map between block and element vector
  std::map<subdomain_id_type, std::vector<unsigned int>  > subdomain_map;

  MeshBase::const_element_iterator it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for (; it!=end; ++it)
  {
    const Elem * elem = *it;
    subdomain_id_type cur_subdomain = elem->subdomain_id();

    subdomain_map[cur_subdomain].push_back(elem->id());
  }
  num_elem_blk = subdomain_map.size();

  if (str_title.size() > MAX_LINE_LENGTH)
    {
      libMesh::err << "Warning, Exodus files cannot have titles longer than "
		   << MAX_LINE_LENGTH
		   << " characters.  Your title will be truncated."
		   << std::endl;
      str_title.resize(MAX_LINE_LENGTH);
    }

  ex_err = exII::ex_put_init(ex_id,
			     str_title.c_str(),
			     num_dim,
			     num_nodes,
			     n_active_elem,
			     num_elem_blk,
			     num_node_sets,
			     num_side_sets);

  EX_CHECK_ERR(ex_err, "Error initializing new Exodus file.");
}



void ExodusII_IO_Helper::write_nodal_coordinates(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  x.resize(num_nodes);
  y.resize(num_nodes);
  z.resize(num_nodes);

  // Use a node iterator instead of looping over i!
  {
    unsigned i = 0;
    MeshBase::const_node_iterator it = mesh.nodes_begin();
    const MeshBase::const_node_iterator end = mesh.nodes_end();
    for (; it!=end; ++it, ++i)
      {
	const Node* node = *it;

	x[i] = (*node)(0) + _coordinate_offset(0);

#if LIBMESH_DIM > 1
	y[i]=(*node)(1) + _coordinate_offset(1);
#else
	y[i]=0.;
#endif
#if LIBMESH_DIM > 2
	z[i]=(*node)(2) + _coordinate_offset(2);
#else
	z[i]=0.;
#endif
      }
  }

  ex_err = exII::ex_put_coord(ex_id,
			      x.empty() ? NULL : &x[0],
			      y.empty() ? NULL : &y[0],
			      z.empty() ? NULL : &z[0]);

  EX_CHECK_ERR(ex_err, "Error writing coordinates to Exodus file.");
}



void ExodusII_IO_Helper::write_nodal_coordinates_discontinuous(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  x.resize(num_nodes);
  y.resize(num_nodes);
  z.resize(num_nodes);

  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  unsigned int i = 0;
  for (; it!=end; ++it)
    for (unsigned int n=0; n<(*it)->n_nodes(); n++)
    {
      x[i]=(*it)->point(n)(0);
#if LIBMESH_DIM > 1
      y[i]=(*it)->point(n)(1);
#else
      y[i]=0.;
#endif
#if LIBMESH_DIM > 2
      z[i]=(*it)->point(n)(2);
#else
      z[i]=0.;
#endif
      i++;
    }

  ex_err = exII::ex_put_coord(ex_id,
			      x.empty() ? NULL : &x[0],
			      y.empty() ? NULL : &y[0],
			      z.empty() ? NULL : &z[0]);

  EX_CHECK_ERR(ex_err, "Error writing coordinates to Exodus file.");
}



void ExodusII_IO_Helper::write_elements(const MeshBase & mesh)
{
  // n_active_elem() is a parallel_only function
  unsigned int n_active_elem = mesh.n_active_elem();

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  std::map<unsigned int, std::vector<unsigned int>  > subdomain_map;

  MeshBase::const_element_iterator mesh_it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  //loop through element and map between block and element vector
  for (; mesh_it!=end; ++mesh_it)
  {
    const Elem * elem = *mesh_it;

    unsigned int cur_subdomain = elem->subdomain_id();

    subdomain_map[cur_subdomain].push_back(elem->id());
  }

  std::vector<int> elem_num_map_out;

  std::map<unsigned int, std::vector<unsigned int>  >::iterator it;

  // element map vector
  num_elem_blk = subdomain_map.size();
  block_ids.resize(num_elem_blk);
  std::vector<unsigned int> elem_map(n_active_elem);
  std::vector<unsigned int>::iterator curr_elem_map_end = elem_map.begin();

  // Note: It appears that there is a bug in exodusII::ex_put_name where
  // the index returned from the ex_id_lkup is erronously used.  For now
  // the work around is to use the alternative function ex_put_names, but
  // this function requires a char** datastructure.
  NamesData names_table(num_elem_blk, MAX_STR_LENGTH);

  unsigned int counter=0;
  for (it=subdomain_map.begin(); it!=subdomain_map.end(); it++)
    {
      block_ids[counter] = (*it).first;
      names_table.push_back_entry(mesh.subdomain_name((*it).first));

      std::vector<unsigned int> & tmp_vec = (*it).second;

      ExodusII_IO_Helper::ElementMaps em;

      //Use the first element in this block to get representative information.
      //Note that Exodus assumes all elements in a block are of the same type!
      //We are using that same assumption here!
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(tmp_vec[0])->type());
      num_nodes_per_elem = mesh.elem(tmp_vec[0])->n_nodes();

      ex_err = exII::ex_put_elem_block(ex_id, (*it).first, conv.exodus_elem_type().c_str(), tmp_vec.size(),num_nodes_per_elem,0);

      EX_CHECK_ERR(ex_err, "Error writing element block.");

      connect.resize(tmp_vec.size()*num_nodes_per_elem);

      for (unsigned int i=0; i<tmp_vec.size(); i++)
        {
          unsigned int elem_id = tmp_vec[i];
          elem_num_map_out.push_back(elem_id);
          libmesh_elem_num_to_exodus[elem_id] = elem_num_map_out.size();

          const Elem* elem = mesh.elem(elem_id);

          // Exodus/Nemesis want every block to have the same element type
          // libmesh_assert_equal_to (elem->type(), conv.get_canonical_type());

	  // But we can get away with writing e.g. HEX8 and INFHEX8 in
	  // the same block...
          libmesh_assert_equal_to (elem->n_nodes(), Elem::build(conv.get_canonical_type(), NULL)->n_nodes());

          for (unsigned int j=0; j<static_cast<unsigned int>(num_nodes_per_elem); j++)
            {
              const unsigned int connect_index   = (i*num_nodes_per_elem)+j;
              const unsigned int elem_node_index = conv.get_inverse_node_map(j); // inverse node map is for writing.
              if (_verbose)
                {
                  libMesh::out << "Exodus node index: " << j
                                << "=LibMesh node index " << elem_node_index << std::endl;
                }
              connect[connect_index] = elem->node(elem_node_index)+1;
            }
        }
    ex_err = exII::ex_put_elem_conn(ex_id, (*it).first, &connect[0]);
    EX_CHECK_ERR(ex_err, "Error writing element connectivities");

    // write out the element number map
    curr_elem_map_end = std::transform(tmp_vec.begin(), tmp_vec.end(), curr_elem_map_end,
                   std::bind2nd(std::plus<unsigned int>(), 1));  // Add one to each id for exodus!
    ex_err = exII::ex_put_elem_num_map(ex_id, (int *)&elem_map[0]);
    EX_CHECK_ERR(ex_err, "Error writing element map");

    counter++;
  }
//  ex_err = exII::ex_put_elem_num_map(ex_id, &elem_num_map_out[0]);
//  EX_CHECK_ERR(ex_err, "Error writing element connectivities");

  // Write out the block names
  if (num_elem_blk > 0)
    {
      ex_err = exII::ex_put_names(ex_id, exII::EX_ELEM_BLOCK, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing element names");
    }
}




void ExodusII_IO_Helper::write_elements_discontinuous(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  std::map<unsigned int, std::vector<unsigned int>  > subdomain_map;

  MeshBase::const_element_iterator mesh_it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  // loop through element and map between block and element vector
  for (; mesh_it!=end; ++mesh_it)
    {
      const Elem * elem = *mesh_it;

      // Only write out the active elements
      if (elem->active())
      {
        unsigned int cur_subdomain = elem->subdomain_id();

        subdomain_map[cur_subdomain].push_back(elem->id());
      }
    }

  std::vector<int> elem_num_map_out;

  std::map<unsigned int, std::vector<unsigned int>  >::iterator it;

  for (it=subdomain_map.begin(); it!=subdomain_map.end(); it++)
    {
      std::vector<unsigned int> & tmp_vec = (*it).second;

      ExodusII_IO_Helper::ElementMaps em;

      //Use the first element in this block to get representative information.
      //Note that Exodus assumes all elements in a block are of the same type!
      //We are using that same assumption here!
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(tmp_vec[0])->type());
      num_nodes_per_elem = mesh.elem(tmp_vec[0])->n_nodes();

      ex_err = exII::ex_put_elem_block(ex_id, (*it).first, conv.exodus_elem_type().c_str(), tmp_vec.size(),num_nodes_per_elem,0);

      EX_CHECK_ERR(ex_err, "Error writing element block.");

      connect.resize(tmp_vec.size()*num_nodes_per_elem);

      for (unsigned int i=0; i<tmp_vec.size(); i++)
        {
          unsigned int elem_id = tmp_vec[i];
          elem_num_map_out.push_back(elem_id);
          libmesh_elem_num_to_exodus[elem_id] = elem_num_map_out.size();

          for (unsigned int j=0; j<static_cast<unsigned int>(num_nodes_per_elem); j++)
            {
              const unsigned int connect_index   = (i*num_nodes_per_elem)+j;
              const unsigned int elem_node_index = conv.get_inverse_node_map(j); // Inverse node map is for writing
              if (_verbose)
                {
                  libMesh::out << "Exodus node index: " << j
                                << "=LibMesh node index " << elem_node_index << std::endl;
                }
              connect[connect_index] = i*num_nodes_per_elem+elem_node_index+1;
            }
        }
    ex_err = exII::ex_put_elem_conn(ex_id, (*it).first, &connect[0]);
    EX_CHECK_ERR(ex_err, "Error writing element connectivities");

    // Create temporary integer storage for the element number map
    std::vector<int> elem_map(tmp_vec.size());

    // Add one to each id for exodus!
    std::transform(tmp_vec.begin(),
                   tmp_vec.end(),
                   elem_map.begin(),
                   std::bind2nd(std::plus<int>(), 1));

    // And write to file
    ex_err = exII::ex_put_elem_num_map(ex_id, &elem_map[0]);
    EX_CHECK_ERR(ex_err, "Error writing element map");
  }

//  ex_err = exII::ex_put_elem_num_map(ex_id, &elem_num_map_out[0]);
//  EX_CHECK_ERR(ex_err, "Error writing element connectivities");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_sidesets(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  ExodusII_IO_Helper::ElementMaps em;

  std::vector< dof_id_type > el;
  std::vector< unsigned short int > sl;
  std::vector< boundary_id_type > il;

  mesh.boundary_info->build_side_list(el, sl, il);

  // Maps from sideset id to the element and sides
  std::map<int, std::vector<int> > elem;
  std::map<int, std::vector<int> > side;

  // Accumulate the vectors to pass into ex_put_side_set
  for (unsigned int i=0; i<el.size(); i++)
  {
    std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
    /**
     * We need to build up active elements if AMR is enabled and add
     * them to the exodus sidesets instead of the potentially inactive "parent" elements
     */
    mesh.elem(el[i])->active_family_tree_by_side(family, sl[i], false);
#else
    family.push_back(mesh.elem(el[i]));
#endif

    for (unsigned int j=0; j<family.size(); ++j)
    {
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(family[j]->id())->type());

      // Use the libmesh to exodus datastructure map to get the proper sideset IDs
      // The datastructure contains the "collapsed" contiguous ids
      elem[il[i]].push_back(libmesh_elem_num_to_exodus[family[j]->id()]);
      side[il[i]].push_back(conv.get_inverse_side_map(sl[i]));
    }
  }

  std::vector<boundary_id_type> side_boundary_ids;
  mesh.boundary_info->build_side_boundary_ids(side_boundary_ids);

  // Write out the sideset names, but only if there is something to write
  if (side_boundary_ids.size() > 0)
    {
      NamesData names_table(side_boundary_ids.size(), MAX_STR_LENGTH);

      for (unsigned int i=0; i<side_boundary_ids.size(); i++)
        {
          int ss_id = side_boundary_ids[i];

          int actual_id = ss_id;

          names_table.push_back_entry(mesh.boundary_info->sideset_name(ss_id));

          ex_err = exII::ex_put_side_set_param(ex_id, actual_id, elem[ss_id].size(), 0);
          EX_CHECK_ERR(ex_err, "Error writing sideset parameters");

          ex_err = exII::ex_put_side_set(ex_id, actual_id, &elem[ss_id][0], &side[ss_id][0]);
          EX_CHECK_ERR(ex_err, "Error writing sidesets");
        }

      ex_err = exII::ex_put_names(ex_id, exII::EX_SIDE_SET, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing sideset names");
    }
}



void ExodusII_IO_Helper::write_nodesets(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  std::vector< dof_id_type > nl;
  std::vector< boundary_id_type > il;

  mesh.boundary_info->build_node_list(nl, il);

  // Maps from nodeset id to the nodes
  std::map<boundary_id_type, std::vector<int> > node;

  // Accumulate the vectors to pass into ex_put_node_set
  for (unsigned int i=0; i<nl.size(); i++)
    node[il[i]].push_back(nl[i]+1);

  std::vector<boundary_id_type> node_boundary_ids;
  mesh.boundary_info->build_node_boundary_ids(node_boundary_ids);

  // Write out the nodeset names, but only if there is something to write
  if (node_boundary_ids.size() > 0)
    {
      NamesData names_table(node_boundary_ids.size(), MAX_STR_LENGTH);

      for (unsigned int i=0; i<node_boundary_ids.size(); i++)
        {
          int nodeset_id = node_boundary_ids[i];

          int actual_id = nodeset_id;

          names_table.push_back_entry(mesh.boundary_info->nodeset_name(nodeset_id));

          ex_err = exII::ex_put_node_set_param(ex_id, actual_id, node[nodeset_id].size(), 0);
          EX_CHECK_ERR(ex_err, "Error writing nodeset parameters");

          ex_err = exII::ex_put_node_set(ex_id, actual_id, &node[nodeset_id][0]);
          EX_CHECK_ERR(ex_err, "Error writing nodesets");
        }

      // Write out the nodeset names
      ex_err = exII::ex_put_names(ex_id, exII::EX_NODE_SET, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing nodeset names");
    }
}



void ExodusII_IO_Helper::initialize_element_variables(const MeshBase & /* mesh */,
                                                      std::vector<std::string> names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  num_elem_vars = names.size();

  if (num_elem_vars == 0)
    return;

  if (_elem_vars_initialized)
    return;

  // Set the flag so we can skip this stuff on subsequent calls to
  // initialize_element_variables()
  _elem_vars_initialized = true;

  ex_err = exII::ex_put_var_param(ex_id,
                                  "e",
                                  num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error setting number of element vars.");

  // Form the element variable truth table and send to Exodus.
  // This tells which variables are written to which blocks,
  // and can dramatically speed up writing element variables

  // We really should initialize all entries in the truth table to 0
  // and then loop over all subdomains, setting their entries to 1
  // if a given variable exists on that subdomain.  However,
  // we don't have that information, and the element variables
  // passed to us are padded with zeroes for the blocks where
  // they aren't defined.  To be consistent with that, fill
  // the truth table with ones.

  std::vector<int> truth_tab(num_elem_blk*num_elem_vars, 1);
  ex_err = exII::ex_put_elem_var_tab(ex_id,
                                     num_elem_blk,
                                     num_elem_vars,
                                     &truth_tab[0]);
  EX_CHECK_ERR(ex_err, "Error writing element truth table.");

  NamesData names_table(num_elem_vars, MAX_STR_LENGTH);

  // Store the input names in the format required by Exodus.
  for (int i=0; i<num_elem_vars; ++i)
    names_table.push_back_entry(names[i]);

  if (_verbose)
    {
      libMesh::out << "Writing variable name(s) to file: " << std::endl;
      for (int i=0; i<num_elem_vars; ++i)
	libMesh::out << names_table.get_char_star(i) << std::endl;
    }

  ex_err = exII::ex_put_var_names(ex_id,
				  "e",
				  num_elem_vars,
				  names_table.get_char_star_star()
				  );

  EX_CHECK_ERR(ex_err, "Error setting element variable names.");
}



void ExodusII_IO_Helper::initialize_nodal_variables(std::vector<std::string> names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  num_nodal_vars = names.size();

  ex_err = exII::ex_put_var_param(ex_id, "n", num_nodal_vars);
  EX_CHECK_ERR(ex_err, "Error setting number of nodal vars.");

  if (num_nodal_vars > 0)
    {
      NamesData names_table(num_nodal_vars, MAX_STR_LENGTH);

      for (int i=0; i<num_nodal_vars; i++)
        names_table.push_back_entry(names[i]);

      if (_verbose)
        {
          libMesh::out << "Writing variable name(s) to file: " << std::endl;
          for (int i=0; i<num_nodal_vars; i++)
            libMesh::out << names_table.get_char_star(i) << std::endl;
        }

      ex_err = exII::ex_put_var_names(ex_id,
                                      "n",
                                      num_nodal_vars,
                                      names_table.get_char_star_star()
                                      );

      EX_CHECK_ERR(ex_err, "Error setting nodal variable names.");
    }
}



void ExodusII_IO_Helper::initialize_global_variables(const std::vector<std::string> & names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_global_vars_initialized)
    return;

  _global_vars_initialized = true;

  num_globals = names.size();

  ex_err = exII::ex_put_var_param(ex_id, "g", num_globals);
  EX_CHECK_ERR(ex_err, "Error setting number of global vars.");

  if (num_globals > 0)
    {
      NamesData names_table(num_globals, MAX_STR_LENGTH);

      for (int i=0; i<num_globals; i++)
        names_table.push_back_entry(names[i]);

      if (_verbose)
        {
          libMesh::out << "Writing variable name(s) to file: " << std::endl;
          for (int i=0; i<num_globals; ++i)
            libMesh::out << names_table.get_char_star(i) << std::endl;
        }

      ex_err = exII::ex_put_var_names(ex_id,
                                      "g",
                                      num_globals,
                                      names_table.get_char_star_star()
                                      );

      EX_CHECK_ERR(ex_err, "Error setting global variable names.");
    }
}



void ExodusII_IO_Helper::write_timestep(int timestep, Real time)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  ex_err = exII::ex_put_time(ex_id, timestep, &time);
  EX_CHECK_ERR(ex_err, "Error writing timestep.");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_element_values(const MeshBase & mesh, const std::vector<Number> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Loop over the element blocks and write the data one block at a time
  std::map<unsigned int, std::vector<unsigned int> > subdomain_map;

  const unsigned int num_vars = values.size() / num_elem;

  MeshBase::const_element_iterator mesh_it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  // loop through element and map between block and element vector
  for (; mesh_it!=end; ++mesh_it)
    {
      const Elem * elem = *mesh_it;
      subdomain_map[elem->subdomain_id()].push_back(elem->id());
    }

  // For each variable, create a 'data' array which holds all the elemental variable
  // values *for a given block* on this processor, then write that data vector to file
  // before moving onto the next block.
  for (unsigned int i=0; i<num_vars; ++i)
    {
      // The size of the subdomain map is the number of blocks.
      std::map<unsigned int, std::vector<unsigned int> >::iterator it = subdomain_map.begin();

      for (unsigned int j=0; it!=subdomain_map.end(); ++it, ++j)
        {
          const std::vector<unsigned int> & elem_nums = (*it).second;
          const unsigned int num_elems_this_block = elem_nums.size();
          std::vector<Number> data(num_elems_this_block);

          for (unsigned int k=0; k<num_elems_this_block; ++k)
            data[k] = values[i*num_elem + elem_nums[k]];

          ex_err = exII::ex_put_elem_var(ex_id,
                                         timestep,
                                         i+1,
                                         this->get_block_id(j),
                                         num_elems_this_block,
                                         &data[0]);
          EX_CHECK_ERR(ex_err, "Error writing element values.");
        }
    }

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_nodal_values(int var_id, const std::vector<Number> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  ex_err = exII::ex_put_nodal_var(ex_id, timestep, var_id, num_nodes, &values[0]);
  EX_CHECK_ERR(ex_err, "Error writing nodal values.");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_information_records(const std::vector<std::string>& records)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  int num_records = records.size();

  if (num_records > 0)
    {
      NamesData info(num_records, MAX_LINE_LENGTH);

      // If an entry is longer than MAX_LINE_LENGTH characters it's not an error, we just
      // write the first MAX_LINE_LENGTH characters to the file.
      for (unsigned i=0; i<records.size(); ++i)
        info.push_back_entry(records[i]);

      ex_err = exII::ex_put_info(ex_id, num_records, info.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing global values.");

      ex_err = exII::ex_update(ex_id);
      EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
    }
}



void ExodusII_IO_Helper::write_global_values(const std::vector<Number> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  ex_err = exII::ex_put_glob_vars(ex_id, timestep, num_globals, &values[0]);
  EX_CHECK_ERR(ex_err, "Error writing global values.");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::use_mesh_dimension_instead_of_spatial_dimension(bool val)
{
  _use_mesh_dimension_instead_of_spatial_dimension = val;
}

void ExodusII_IO_Helper::set_coordinate_offset(Point p)
{
  _coordinate_offset = p;
}


bool ExodusII_IO_Helper::created()
{
  return _created;
}



bool ExodusII_IO_Helper::opened()
{
  return _opened;
}



// ------------------------------------------------------------
// ExodusII_IO_Helper::Conversion class members
ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(std::string type_str)
{
  init_element_equivalence_map();

  // Do only upper-case comparisons
  std::transform(type_str.begin(), type_str.end(), type_str.begin(), ::toupper);

  std::map<std::string, libMeshEnums::ElemType>::iterator it =
    element_equivalence_map.find(type_str);

  if (it != element_equivalence_map.end())
    return assign_conversion( it->second );
  else
    {
      libMesh::err << "ERROR! Unrecognized element type_str: " << type_str << std::endl;
      libmesh_error();
    }

  libmesh_error();

  // dummy return value, we won't get here
  return assign_conversion (EDGE2);
}



ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(const ElemType type)
{
  switch (type)
    {
    case EDGE2:
      {
	const Conversion conv(edge2_node_map,
			      ARRAY_LENGTH(edge2_node_map),
			      edge2_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(edge2_node_map),
			      edge_edge_map,
			      ARRAY_LENGTH(edge_edge_map),
			      edge_inverse_edge_map,
			      ARRAY_LENGTH(edge_inverse_edge_map),
			      EDGE2, "EDGE2");
	return conv;
      }
    case EDGE3:
      {
	const Conversion conv(edge3_node_map,
			      ARRAY_LENGTH(edge3_node_map),
			      edge3_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(edge3_node_map),
			      edge_edge_map,
			      ARRAY_LENGTH(edge_edge_map),
			      edge_inverse_edge_map,
			      ARRAY_LENGTH(edge_inverse_edge_map),
			      EDGE3, "EDGE3");
	return conv;
      }
    case QUAD4:
      {
	const Conversion conv(quad4_node_map,
			      ARRAY_LENGTH(quad4_node_map),
			      quad4_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(quad4_node_map),
			      quad_edge_map,
			      ARRAY_LENGTH(quad_edge_map),
			      quad_inverse_edge_map,
			      ARRAY_LENGTH(quad_inverse_edge_map),
			      QUAD4,
			      "QUAD4");
	return conv;
      }

    case QUAD8:
      {
	const Conversion conv(quad8_node_map,
			      ARRAY_LENGTH(quad8_node_map),
			      quad8_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(quad8_node_map),
			      quad_edge_map,
			      ARRAY_LENGTH(quad_edge_map),
			      quad_inverse_edge_map,
			      ARRAY_LENGTH(quad_inverse_edge_map),
			      QUAD8,
			      "QUAD8");
	return conv;
      }

    case QUAD9:
      {
	const Conversion conv(quad9_node_map,
			      ARRAY_LENGTH(quad9_node_map),
			      quad9_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(quad9_node_map),
			      quad_edge_map,
			      ARRAY_LENGTH(quad_edge_map),
			      quad_inverse_edge_map,
			      ARRAY_LENGTH(quad_inverse_edge_map),
			      QUAD9,
			      "QUAD9");
	return conv;
      }

    case TRI3:
      {
	const Conversion conv(tri3_node_map,
			      ARRAY_LENGTH(tri3_node_map),
			      tri3_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(tri3_node_map),
			      tri_edge_map,
			      ARRAY_LENGTH(tri_edge_map),
			      tri_inverse_edge_map,
			      ARRAY_LENGTH(tri_inverse_edge_map),
			      TRI3,
			      "TRI3");
	return conv;
      }

    case TRI6:
      {
	const Conversion conv(tri6_node_map,
			      ARRAY_LENGTH(tri6_node_map),
			      tri6_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(tri6_node_map),
			      tri_edge_map,
			      ARRAY_LENGTH(tri_edge_map),
			      tri_inverse_edge_map,
			      ARRAY_LENGTH(tri_inverse_edge_map),
			      TRI6,
			      "TRI6");
	return conv;
      }

    case HEX8:
      {
	const Conversion conv(hex8_node_map,
			      ARRAY_LENGTH(hex8_node_map),
			      hex8_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(hex8_node_map),
			      hex_face_map,
			      ARRAY_LENGTH(hex_face_map),
			      hex_inverse_face_map,
			      ARRAY_LENGTH(hex_inverse_face_map),
			      HEX8,
			      "HEX8");
	return conv;
      }

    case HEX20:
      {
	const Conversion conv(hex20_node_map,
			      ARRAY_LENGTH(hex20_node_map),
			      hex20_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(hex20_node_map),
			      hex_face_map,
			      ARRAY_LENGTH(hex_face_map),
			      hex_inverse_face_map,
			      ARRAY_LENGTH(hex_inverse_face_map),
			      HEX20,
			      "HEX20");
	return conv;
      }

    case HEX27:
      {
	const Conversion conv(hex27_node_map,
			      ARRAY_LENGTH(hex27_node_map),
			      hex27_inverse_node_map, // different inverse node map for Hex27!
			      ARRAY_LENGTH(hex27_inverse_node_map),
			      hex27_face_map,
			      ARRAY_LENGTH(hex27_face_map),
			      hex27_inverse_face_map,
			      ARRAY_LENGTH(hex27_inverse_face_map),
			      HEX27,
			      "HEX27");
	return conv;
      }

    case TET4:
      {
	const Conversion conv(tet4_node_map,
			      ARRAY_LENGTH(tet4_node_map),
			      tet4_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(tet4_node_map),
			      tet_face_map,
			      ARRAY_LENGTH(tet_face_map),
			      tet_inverse_face_map,
			      ARRAY_LENGTH(tet_inverse_face_map),
			      TET4,
			      "TETRA4");
	return conv;
      }

    case TET10:
      {
	const Conversion conv(tet10_node_map,
			      ARRAY_LENGTH(tet10_node_map),
			      tet10_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(tet10_node_map),
			      tet_face_map,
			      ARRAY_LENGTH(tet_face_map),
			      tet_inverse_face_map,
			      ARRAY_LENGTH(tet_inverse_face_map),
			      TET10,
			      "TETRA10");
	return conv;
      }

    case PRISM6:
      {
	const Conversion conv(prism6_node_map,
			      ARRAY_LENGTH(prism6_node_map),
			      prism6_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(prism6_node_map),
			      prism_face_map,
			      ARRAY_LENGTH(prism_face_map),
			      prism_inverse_face_map,
			      ARRAY_LENGTH(prism_inverse_face_map),
			      PRISM6,
			      "WEDGE");
	return conv;
      }

    case PRISM15:
      {
	const Conversion conv(prism15_node_map,
			      ARRAY_LENGTH(prism15_node_map),
			      prism15_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(prism15_node_map),
			      prism_face_map,
			      ARRAY_LENGTH(prism_face_map),
			      prism_inverse_face_map,
			      ARRAY_LENGTH(prism_inverse_face_map),
			      PRISM15,
			      "WEDGE15");
	return conv;
      }

    case PRISM18:
      {
	const Conversion conv(prism18_node_map,
			      ARRAY_LENGTH(prism18_node_map),
			      prism18_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(prism18_node_map),
			      prism_face_map,
			      ARRAY_LENGTH(prism_face_map),
			      prism_inverse_face_map,
			      ARRAY_LENGTH(prism_inverse_face_map),
			      PRISM18,
			      "WEDGE18");
	return conv;
      }

    case PYRAMID5:
      {
	const Conversion conv(pyramid5_node_map,
			      ARRAY_LENGTH(pyramid5_node_map),
			      pyramid5_node_map, // inverse node map same as forward node map
			      ARRAY_LENGTH(pyramid5_node_map),
			      pyramid_face_map,
			      ARRAY_LENGTH(pyramid_face_map),
			      pyramid_inverse_face_map,
			      ARRAY_LENGTH(pyramid_inverse_face_map),
			      PYRAMID5,
			      "PYRAMID5");
	return conv;
      }

    default:
      libmesh_error();
    }

  libmesh_error();

  // dummy return value, we will never get here
  const Conversion conv(tri3_node_map,
			ARRAY_LENGTH(tri3_node_map),
			tri3_node_map, // inverse node map same as forward node map
			ARRAY_LENGTH(tri3_node_map),
			tri_edge_map,
			ARRAY_LENGTH(tri_edge_map),
			tri_inverse_edge_map,
			ARRAY_LENGTH(tri_inverse_edge_map),
			TRI3,
			"TRI3");
  return conv;
}



ExodusII_IO_Helper::NamesData::NamesData(size_t n_strings, size_t string_length) :
    data_table(n_strings),
    data_table_pointers(n_strings),
    counter(0),
    table_size(n_strings)
{
  for (size_t i=0; i<n_strings; ++i)
    {
      data_table[i].resize(string_length + 1);

      // NULL-terminate these strings, just to be safe.
      data_table[i][0] = '\0';

      // Set pointer into the data_table
      data_table_pointers[i] = &(data_table[i][0]);
    }
}



void ExodusII_IO_Helper::NamesData::push_back_entry(const std::string & name)
{
  libmesh_assert_less (counter, table_size);

  // 1.) Copy the C++ string into the vector<char>...
  size_t num_copied = name.copy(&data_table[counter][0], data_table[counter].size()-1);

  // 2.) ...And null-terminate it.
  data_table[counter][num_copied] = '\0';

  // Go to next row
  ++counter;
}



char** ExodusII_IO_Helper::NamesData::get_char_star_star()
{
  return &data_table_pointers[0];
}



char* ExodusII_IO_Helper::NamesData::get_char_star(int i)
{
  if (static_cast<unsigned>(i) >= table_size)
    {
      libMesh::err << "Requested char* " << i << " but only have " << table_size << "!" << std::endl;
      libmesh_error();
    }
  else
    return &(data_table[i][0]);
}


} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_EXODUS_API
