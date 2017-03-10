// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <sstream>
#include <cstdlib> // std::strtol

#include "libmesh/boundary_info.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/elem.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/string_to_enum.h"

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
using namespace libMesh;

// Define equivalence classes of Cubit/Exodus element types that map to
// libmesh ElemTypes
std::map<std::string, ElemType> element_equivalence_map;

// This function initializes the element_equivalence_map the first time it
// is called, and returns early all other times.
void init_element_equivalence_map()
{
  if (element_equivalence_map.empty())
    {
      // We use an ExodusII SPHERE element to represent a NodeElem
      element_equivalence_map["SPHERE"] = NODEELEM;

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

      // QUADSHELL4 equivalences
      element_equivalence_map["SHELL"]  = QUADSHELL4;
      element_equivalence_map["SHELL4"] = QUADSHELL4;

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

      // TRISHELL3 equivalences
      element_equivalence_map["TRISHELL"]  = TRISHELL3;
      element_equivalence_map["TRISHELL3"] = TRISHELL3;

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

      // PYRAMID13 equivalences
      element_equivalence_map["PYRAMID13"] = PYRAMID13;

      // PYRAMID14 equivalences
      element_equivalence_map["PYRAMID14"] = PYRAMID14;
    }
}
}



namespace libMesh
{

// ------------------------------------------------------------
// ExodusII_IO_Helper::ElementMaps static data

// 0D node map definitions
const int ExodusII_IO_Helper::ElementMaps::nodeelem_node_map[1] = {0};

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
const int ExodusII_IO_Helper::ElementMaps::trishell3_edge_map[3] = {0, 1, 2};
const int ExodusII_IO_Helper::ElementMaps::quadshell4_edge_map[4] = {0, 1, 2, 3};

//These take a libMesh ID and turn it into an Exodus ID
const int ExodusII_IO_Helper::ElementMaps::tri_inverse_edge_map[3] = {1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::quad_inverse_edge_map[4] = {1, 2, 3, 4};
const int ExodusII_IO_Helper::ElementMaps::trishell3_inverse_edge_map[3] = {3, 4, 5};
const int ExodusII_IO_Helper::ElementMaps::quadshell4_inverse_edge_map[4] = {3, 4, 5, 6};

// 3D node map definitions
const int ExodusII_IO_Helper::ElementMaps::hex8_node_map[8]   = {0, 1, 2, 3, 4, 5, 6, 7};
const int ExodusII_IO_Helper::ElementMaps::hex20_node_map[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                                  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

// Perhaps an older Hex27 node numbering?  This no longer works.
//const int ExodusII_IO_Helper::ElementMaps::hex27_node_map[27] = { 1,  5,  6,  2,  0,  4,  7,  3, 13, 17, 14,  9,  8, 16,
//  18, 10, 12, 19, 15, 11, 24, 25, 22, 26, 21, 23, 20};

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
const int ExodusII_IO_Helper::ElementMaps::pyramid13_node_map[13] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
const int ExodusII_IO_Helper::ElementMaps::pyramid14_node_map[14] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

// Shell element face maps
const int ExodusII_IO_Helper::ElementMaps::trishell3_shellface_map[2] = {0, 1};
const int ExodusII_IO_Helper::ElementMaps::quadshell4_shellface_map[2] = {0, 1};

//These take a libMesh ID and turn it into an Exodus ID
const int ExodusII_IO_Helper::ElementMaps::trishell3_inverse_shellface_map[2] = {1, 2};
const int ExodusII_IO_Helper::ElementMaps::quadshell4_inverse_shellface_map[2] = {1, 2};

// 3D face map definitions
const int ExodusII_IO_Helper::ElementMaps::tet_face_map[4]     = {1, 2, 3, 0};

const int ExodusII_IO_Helper::ElementMaps::hex_face_map[6]     = {1, 2, 3, 4, 0, 5};
const int ExodusII_IO_Helper::ElementMaps::hex27_face_map[6]   = {1, 2, 3, 4, 0, 5};
//const int ExodusII_IO_Helper::ElementMaps::hex27_face_map[6]   = {1, 0, 3, 5, 4, 2};
const int ExodusII_IO_Helper::ElementMaps::prism_face_map[5]   = {1, 2, 3, 0, 4};

// Pyramids are not included in the ExodusII specification. The ordering below matches
// the sideset ordering that CUBIT generates.
const int ExodusII_IO_Helper::ElementMaps::pyramid_face_map[5] = {0, 1, 2, 3, 4};

//These take a libMesh ID and turn it into an Exodus ID
const int ExodusII_IO_Helper::ElementMaps::tet_inverse_face_map[4]     = {4, 1, 2, 3};
const int ExodusII_IO_Helper::ElementMaps::hex_inverse_face_map[6]     = {5, 1, 2, 3, 4, 6};
const int ExodusII_IO_Helper::ElementMaps::hex27_inverse_face_map[6]   = {5, 1, 2, 3, 4, 6};
//const int ExodusII_IO_Helper::ElementMaps::hex27_inverse_face_map[6]   = {2, 1, 6, 3, 5, 4};
const int ExodusII_IO_Helper::ElementMaps::prism_inverse_face_map[5]   = {4, 1, 2, 3, 5};

// Pyramids are not included in the ExodusII specification. The ordering below matches
// the sideset ordering that CUBIT generates.
const int ExodusII_IO_Helper::ElementMaps::pyramid_inverse_face_map[5] = {1, 2, 3, 4, 5};



// ExodusII_IO_Helper::Conversion static data
const int ExodusII_IO_Helper::Conversion::invalid_id = std::numeric_limits<int>::max();

// ------------------------------------------------------------
// ExodusII_IO_Helper class members

ExodusII_IO_Helper::ExodusII_IO_Helper(const ParallelObject & parent,
                                       bool v,
                                       bool run_only_on_proc0,
                                       bool single_precision) :
  ParallelObject(parent),
  ex_id(0),
  ex_err(0),
  num_dim(0),
  num_global_vars(0),
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
  verbose(v),
  opened_for_writing(false),
  opened_for_reading(false),
  _run_only_on_proc0(run_only_on_proc0),
  _elem_vars_initialized(false),
  _global_vars_initialized(false),
  _nodal_vars_initialized(false),
  _use_mesh_dimension_instead_of_spatial_dimension(false),
  _write_as_dimension(0),
  _single_precision(single_precision)
{
  title.resize(MAX_LINE_LENGTH+1);
  elem_type.resize(MAX_STR_LENGTH);
}



ExodusII_IO_Helper::~ExodusII_IO_Helper()
{
}



const char * ExodusII_IO_Helper::get_elem_type() const
{
  return &elem_type[0];
}



void ExodusII_IO_Helper::message(const std::string & msg)
{
  if (verbose) libMesh::out << msg << std::endl;
}



void ExodusII_IO_Helper::message(const std::string & msg, int i)
{
  if (verbose) libMesh::out << msg << i << "." << std::endl;
}



void ExodusII_IO_Helper::open(const char * filename, bool read_only)
{
  // Version of Exodus you are using
  float ex_version = 0.;

  // Word size in bytes of the floating point variables used in the
  // application program.  Exodus only supports 4-byte and 8-byte
  // floats.
  int comp_ws = std::min(sizeof(Real), std::size_t(8));

  // Word size in bytes of the floating point data as they are stored
  // in the ExodusII file.  "If this argument is 0, the word size of the
  // floating point data already stored in the file is returned"
  int io_ws = 0;

  ex_id = exII::ex_open(filename,
                        read_only ? EX_READ : EX_WRITE,
                        &comp_ws,
                        &io_ws,
                        &ex_version);

  std::string err_msg = std::string("Error opening ExodusII mesh file: ") + std::string(filename);
  EX_CHECK_ERR(ex_id, err_msg);
  if (verbose) libMesh::out << "File opened successfully." << std::endl;

  if (read_only)
    opened_for_reading = true;
  else
    opened_for_writing = true;

  current_filename = std::string(filename);
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

  this->read_num_time_steps();

  ex_err = exII::ex_get_var_param(ex_id, "n", &num_nodal_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of nodal variables.");

  ex_err = exII::ex_get_var_param(ex_id, "e", &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  ex_err = exII::ex_get_var_param(ex_id, "g", &num_global_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of global variables.");

  message("Exodus header info retrieved successfully.");
}




void ExodusII_IO_Helper::read_qa_records()
{
  // The QA records are four MAX_STR_LENGTH-byte character strings.
  int num_qa_rec =
    this->inquire(exII::EX_INQ_QA, "Error retrieving number of QA records");

  if (verbose)
    libMesh::out << "Found "
                 << num_qa_rec
                 << " QA record(s) in the Exodus file."
                 << std::endl;

  if (num_qa_rec > 0)
    {
      // How to dynamically allocate an array of fixed-size char * arrays in C++.
      // http://stackoverflow.com/questions/8529359/creating-a-dynamic-sized-array-of-fixed-sized-int-arrays-in-c
      typedef char * inner_array_t[4];
      inner_array_t * qa_record = new inner_array_t[num_qa_rec];

      for (int i=0; i<num_qa_rec; i++)
        for (int j=0; j<4; j++)
          qa_record[i][j] = new char[MAX_STR_LENGTH+1];

      ex_err = exII::ex_get_qa (ex_id, qa_record);
      EX_CHECK_ERR(ex_err, "Error reading the QA records.");

      // Print the QA records
      if (verbose)
        {
          for (int i=0; i<num_qa_rec; i++)
            {
              libMesh::out << "QA Record: " << i << std::endl;
              for (int j=0; j<4; j++)
                libMesh::out << qa_record[i][j] << std::endl;
            }
        }


      // Clean up dynamically-allocated memory
      for (int i=0; i<num_qa_rec; i++)
        for (int j=0; j<4; j++)
          delete [] qa_record[i][j];

      delete [] qa_record;
    }
}




void ExodusII_IO_Helper::print_header()
{
  if (verbose)
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
                              static_cast<void *>(&x[0]),
                              static_cast<void *>(&y[0]),
                              static_cast<void *>(&z[0]));

  EX_CHECK_ERR(ex_err, "Error retrieving nodal data.");
  message("Nodal data retrieved successfully.");
}



void ExodusII_IO_Helper::read_node_num_map ()
{
  node_num_map.resize(num_nodes);

  ex_err = exII::ex_get_node_num_map (ex_id,
                                      node_num_map.empty() ? libmesh_nullptr : &node_num_map[0]);

  EX_CHECK_ERR(ex_err, "Error retrieving nodal number map.");
  message("Nodal numbering map retrieved successfully.");

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] node_num_map[i] = ";
      for (unsigned int i=0; i<static_cast<unsigned int>(std::min(10, num_nodes-1)); ++i)
        libMesh::out << node_num_map[i] << ", ";
      libMesh::out << "... " << node_num_map.back() << std::endl;
    }
}


void ExodusII_IO_Helper::print_nodes(std::ostream & out_stream)
{
  for (int i=0; i<num_nodes; i++)
    out_stream << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
}



void ExodusII_IO_Helper::read_block_info()
{
  block_ids.resize(num_elem_blk);
  // Get all element block IDs.
  ex_err = exII::ex_get_elem_blk_ids(ex_id,
                                     block_ids.empty() ? libmesh_nullptr : &block_ids[0]);
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
  if (verbose)
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
                                      elem_num_map.empty() ? libmesh_nullptr : &elem_num_map[0]);

  EX_CHECK_ERR(ex_err, "Error retrieving element number map.");
  message("Element numbering map retrieved successfully.");


  if (verbose)
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
      // Don't close the file if it was never opened, this raises an Exodus error
      if (opened_for_writing || opened_for_reading)
        {
          ex_err = exII::ex_close(ex_id);
          EX_CHECK_ERR(ex_err, "Error closing Exodus file.");
          message("Exodus file closed successfully.");
        }
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
  // Make sure we have an up-to-date count of the number of time steps in the file.
  this->read_num_time_steps();

  if (num_time_steps > 0)
    {
      time_steps.resize(num_time_steps);
      ex_err = exII::ex_get_all_times(ex_id, &time_steps[0]);
      EX_CHECK_ERR(ex_err, "Error reading timesteps!");
    }
}



void ExodusII_IO_Helper::read_num_time_steps()
{
  num_time_steps =
    this->inquire(exII::EX_INQ_TIME, "Error retrieving number of time steps");
}



void ExodusII_IO_Helper::read_nodal_var_values(std::string nodal_var_name, int time_step)
{
  // Read the nodal variable names from file, so we can see if we have the one we're looking for
  this->read_var_names(NODAL);

  // See if we can find the variable we are looking for
  std::size_t var_index = 0;
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
      libMesh::err << "Available variables: " << std::endl;
      for (std::size_t i=0; i<nodal_var_names.size(); ++i)
        libMesh::err << nodal_var_names[i] << std::endl;

      libmesh_error_msg("Unable to locate variable named: " << nodal_var_name);
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



void ExodusII_IO_Helper::read_var_names(ExodusVarType type)
{
  switch (type)
    {
    case NODAL:
      this->read_var_names_impl("n", num_nodal_vars, nodal_var_names);
      break;
    case ELEMENTAL:
      this->read_var_names_impl("e", num_elem_vars, elem_var_names);
      break;
    case GLOBAL:
      this->read_var_names_impl("g", num_global_vars, global_var_names);
      break;
    default:
      libmesh_error_msg("Unrecognized ExodusVarType " << type);
    }
}



void ExodusII_IO_Helper::read_var_names_impl(const char * var_type,
                                             int & count,
                                             std::vector<std::string> & result)
{
  // First read and store the number of names we have
  ex_err = exII::ex_get_var_param(ex_id, var_type, &count);
  EX_CHECK_ERR(ex_err, "Error reading number of variables.");

  // Do nothing if no variables are detected
  if (count == 0)
    return;

  // Second read the actual names and convert them into a format we can use
  NamesData names_table(count, MAX_STR_LENGTH);

  ex_err = exII::ex_get_var_names(ex_id,
                                  var_type,
                                  count,
                                  names_table.get_char_star_star()
                                  );
  EX_CHECK_ERR(ex_err, "Error reading variable names!");

  if (verbose)
    {
      libMesh::out << "Read the variable(s) from the file:" << std::endl;
      for (int i=0; i<count; i++)
        libMesh::out << names_table.get_char_star(i) << std::endl;
    }

  // Allocate enough space for our variable name strings.
  result.resize(count);

  // Copy the char buffers into strings.
  for (int i=0; i<count; i++)
    result[i] = names_table.get_char_star(i); // calls string::op=(const char *)
}




void ExodusII_IO_Helper::write_var_names(ExodusVarType type, std::vector<std::string> & names)
{
  switch (type)
    {
    case NODAL:
      this->write_var_names_impl("n", num_nodal_vars, names);
      break;
    case ELEMENTAL:
      this->write_var_names_impl("e", num_elem_vars, names);
      break;
    case GLOBAL:
      this->write_var_names_impl("g", num_global_vars, names);
      break;
    default:
      libmesh_error_msg("Unrecognized ExodusVarType " << type);
    }
}



void ExodusII_IO_Helper::write_var_names_impl(const char * var_type, int & count, std::vector<std::string> & names)
{
  // Update the count variable so that it's available to other parts of the class.
  count = cast_int<int>(names.size());

  // Write that number of variables to the file.
  ex_err = exII::ex_put_var_param(ex_id, var_type, count);
  EX_CHECK_ERR(ex_err, "Error setting number of vars.");

  if (names.size() > 0)
    {
      NamesData names_table(names.size(), MAX_STR_LENGTH);

      // Store the input names in the format required by Exodus.
      for (std::size_t i=0; i<names.size(); ++i)
        names_table.push_back_entry(names[i]);

      if (verbose)
        {
          libMesh::out << "Writing variable name(s) to file: " << std::endl;
          for (std::size_t i=0; i<names.size(); ++i)
            libMesh::out << names_table.get_char_star(i) << std::endl;
        }

      ex_err = exII::ex_put_var_names(ex_id,
                                      var_type,
                                      cast_int<int>(names.size()),
                                      names_table.get_char_star_star()
                                      );

      EX_CHECK_ERR(ex_err, "Error writing variable names.");
    }
}



void ExodusII_IO_Helper::read_elemental_var_values(std::string elemental_var_name,
                                                   int time_step,
                                                   std::map<dof_id_type, Real> & elem_var_value_map)
{
  this->read_var_names(ELEMENTAL);

  // See if we can find the variable we are looking for
  std::size_t var_index = 0;
  bool found = false;

  // Do a linear search for nodal_var_name in nodal_var_names
  for (; var_index<elem_var_names.size(); ++var_index)
    if (elem_var_names[var_index] == elemental_var_name)
      {
        found = true;
        break;
      }

  if (!found)
    {
      libMesh::err << "Available variables: " << std::endl;
      for (std::size_t i=0; i<elem_var_names.size(); ++i)
        libMesh::err << elem_var_names[i] << std::endl;

      libmesh_error_msg("Unable to locate variable named: " << elemental_var_name);
    }

  // Sequential index which we can use to look up the element ID in the elem_num_map.
  unsigned ex_el_num = 0;

  for (unsigned i=0; i<static_cast<unsigned>(num_elem_blk); i++)
    {
      ex_err = exII::ex_get_elem_block(ex_id,
                                       block_ids[i],
                                       libmesh_nullptr,
                                       &num_elem_this_blk,
                                       libmesh_nullptr,
                                       libmesh_nullptr);
      EX_CHECK_ERR(ex_err, "Error getting number of elements in block.");

      std::vector<Real> block_elem_var_values(num_elem);
      ex_err = exII::ex_get_elem_var(ex_id,
                                     time_step,
                                     var_index+1,
                                     block_ids[i],
                                     num_elem_this_blk,
                                     &block_elem_var_values[0]);
      EX_CHECK_ERR(ex_err, "Error getting elemental values.");

      for (unsigned j=0; j<static_cast<unsigned>(num_elem_this_blk); j++)
        {
          // Use the elem_num_map to obtain the ID of this element in the Exodus file,
          // and remember to subtract 1 since libmesh is zero-based and Exodus is 1-based.
          unsigned mapped_elem_id = this->elem_num_map[ex_el_num] - 1;

          // Store the elemental value in the map.
          elem_var_value_map[mapped_elem_id] = block_elem_var_values[j];

          // Go to the next sequential element ID.
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
      int
        comp_ws = 0,
        io_ws = 0;

      if(_single_precision)
        {
          comp_ws = cast_int<int>(sizeof(float));
          io_ws = cast_int<int>(sizeof(float));
        }
      // Fall back on double precision when necessary since ExodusII
      // doesn't seem to support long double
      else
        {
          comp_ws = cast_int<int>
            (std::min(sizeof(Real), sizeof(double)));
          io_ws = cast_int<int>
            (std::min(sizeof(Real), sizeof(double)));
        }

      ex_id = exII::ex_create(filename.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

      EX_CHECK_ERR(ex_id, "Error creating ExodusII mesh file.");

      if (verbose)
        libMesh::out << "File created successfully." << std::endl;
    }

  opened_for_writing = true;
  current_filename = filename;
}



void ExodusII_IO_Helper::initialize(std::string str_title, const MeshBase & mesh, bool use_discontinuous)
{
  // n_active_elem() is a parallel_only function
  unsigned int n_active_elem = mesh.n_active_elem();

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // If _write_as_dimension is nonzero, use it to set num_dim in the Exodus file.
  if (_write_as_dimension)
    num_dim = _write_as_dimension;
  else if (_use_mesh_dimension_instead_of_spatial_dimension)
    num_dim = mesh.mesh_dimension();
  else
    num_dim = mesh.spatial_dimension();

  num_elem = mesh.n_elem();

  if (!use_discontinuous)
    {
      // Don't rely on mesh.n_nodes() here.  If nodes have been
      // deleted recently, it will be incorrect.
      num_nodes = cast_int<int>(std::distance(mesh.nodes_begin(),
                                              mesh.nodes_end()));
    }
  else
    {
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();
      for (; it!=end; ++it)
        num_nodes += (*it)->n_nodes();
    }

  std::vector<boundary_id_type> unique_side_boundaries;
  std::vector<boundary_id_type> unique_node_boundaries;

  mesh.get_boundary_info().build_side_boundary_ids(unique_side_boundaries);
  {
    // Add shell face boundaries to the list of side boundaries, since ExodusII
    // treats these the same way.
    std::vector<boundary_id_type> shellface_boundaries;
    mesh.get_boundary_info().build_shellface_boundary_ids(shellface_boundaries);
    for (std::size_t i=0; i<shellface_boundaries.size(); i++)
      unique_side_boundaries.push_back(shellface_boundaries[i]);
  }
  mesh.get_boundary_info().build_node_boundary_ids(unique_node_boundaries);

  num_side_sets = cast_int<int>(unique_side_boundaries.size());
  num_node_sets = cast_int<int>(unique_node_boundaries.size());

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
  num_elem_blk = cast_int<int>(subdomain_map.size());

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



void ExodusII_IO_Helper::write_nodal_coordinates(const MeshBase & mesh, bool use_discontinuous)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Clear existing data from any previous calls.
  x.clear();
  y.clear();
  z.clear();
  node_num_map.clear();

  // Reserve space in the nodal coordinate vectors.  num_nodes is
  // exact, this just allows us to do away with one potentially
  // error-inducing loop index.
  x.reserve(num_nodes);
  y.reserve(num_nodes);
  z.reserve(num_nodes);

  // And in the node_num_map - since the nodes aren't organized in
  // blocks, libmesh will always write out the identity map
  // here... unless there has been some refinement and coarsening, or
  // node deletion without a corresponding call to contract(). You
  // need to write this any time there could be 'holes' in the node
  // numbering, so we write it every time.
  node_num_map.reserve(num_nodes);

  // Clear out any previously-mapped node IDs.
  libmesh_node_num_to_exodus.clear();

  if (!use_discontinuous)
    {
      MeshBase::const_node_iterator it = mesh.nodes_begin();
      const MeshBase::const_node_iterator end = mesh.nodes_end();
      for (; it != end; ++it)
        {
          const Node & node = **it;

          x.push_back(node(0) + _coordinate_offset(0));

#if LIBMESH_DIM > 1
          y.push_back(node(1) + _coordinate_offset(1));
#else
          y.push_back(0.);
#endif
#if LIBMESH_DIM > 2
          z.push_back(node(2) + _coordinate_offset(2));
#else
          z.push_back(0.);
#endif

          // Fill in node_num_map entry with the proper (1-based) node id
          node_num_map.push_back(node.id() + 1);

          // Also map the zero-based libmesh node id to the 1-based
          // Exodus ID it will be assigned (this is equivalent to the
          // current size of the x vector).
          libmesh_node_num_to_exodus[ cast_int<int>(node.id()) ] = cast_int<int>(x.size());
        }
    }
  else
    {
      MeshBase::const_element_iterator it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end();

      for (; it!=end; ++it)
        {
          const Elem * elem = *it;

          for (unsigned int n=0; n<elem->n_nodes(); n++)
            {
              x.push_back(elem->point(n)(0));
#if LIBMESH_DIM > 1
              y.push_back(elem->point(n)(1));
#else
              y.push_back(0.);
#endif
#if LIBMESH_DIM > 2
              z.push_back(elem->point(n)(2));
#else
              z.push_back(0.);
#endif

              // Let's skip the node_num_map in the discontinuous
              // case, since we're effectively duplicating nodes for
              // the sake of discontinuous visualization, so it isn't
              // clear how to deal with node_num_map here. This means
              // that writing discontinuous meshes won't work with
              // element numberings that have "holes".
            }
        }
    }

  if (_single_precision)
    {
      std::vector<float>
        x_single(x.begin(), x.end()),
        y_single(y.begin(), y.end()),
        z_single(z.begin(), z.end());

      ex_err = exII::ex_put_coord(ex_id,
                                  x_single.empty() ? libmesh_nullptr : &x_single[0],
                                  y_single.empty() ? libmesh_nullptr : &y_single[0],
                                  z_single.empty() ? libmesh_nullptr : &z_single[0]);
    }
  else
    {
      ex_err = exII::ex_put_coord(ex_id,
                                  x.empty() ? libmesh_nullptr : &x[0],
                                  y.empty() ? libmesh_nullptr : &y[0],
                                  z.empty() ? libmesh_nullptr : &z[0]);
    }


  EX_CHECK_ERR(ex_err, "Error writing coordinates to Exodus file.");

  if (!use_discontinuous)
    {
      // Also write the (1-based) node_num_map to the file.
      ex_err = exII::ex_put_node_num_map(ex_id, &node_num_map[0]);
      EX_CHECK_ERR(ex_err, "Error writing node_num_map");
    }
}



void ExodusII_IO_Helper::write_elements(const MeshBase & mesh, bool use_discontinuous)
{
  // n_active_elem() is a parallel_only function
  unsigned int n_active_elem = mesh.n_active_elem();

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Map from block ID to a vector of element IDs in that block.  Element
  // IDs are now of type dof_id_type, subdomain IDs are of type subdomain_id_type.
  typedef std::map<subdomain_id_type, std::vector<dof_id_type> > subdomain_map_type;
  subdomain_map_type subdomain_map;

  MeshBase::const_element_iterator mesh_it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  //loop through element and map between block and element vector
  for (; mesh_it!=end; ++mesh_it)
    {
      const Elem * elem = *mesh_it;
      subdomain_map[ elem->subdomain_id() ].push_back(elem->id());
    }

  // element map vector
  num_elem_blk = cast_int<int>(subdomain_map.size());
  block_ids.resize(num_elem_blk);
  elem_num_map.resize(n_active_elem);
  std::vector<int>::iterator curr_elem_map_end = elem_num_map.begin();

  // Note: It appears that there is a bug in exodusII::ex_put_name where
  // the index returned from the ex_id_lkup is erronously used.  For now
  // the work around is to use the alternative function ex_put_names, but
  // this function requires a char ** datastructure.
  NamesData names_table(num_elem_blk, MAX_STR_LENGTH);

  // This counter is used to fill up the libmesh_elem_num_to_exodus map in the loop below.
  unsigned libmesh_elem_num_to_exodus_counter = 0;

  // counter indexes into the block_ids vector
  unsigned int counter = 0;

  // node counter for discontinuous plotting
  unsigned int node_counter = 0;

  for (subdomain_map_type::iterator it=subdomain_map.begin(); it!=subdomain_map.end(); ++it)
    {
      block_ids[counter] = (*it).first;
      names_table.push_back_entry(mesh.subdomain_name((*it).first));

      // Get a reference to a vector of element IDs for this subdomain.
      subdomain_map_type::mapped_type & tmp_vec = (*it).second;

      ExodusII_IO_Helper::ElementMaps em;
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
      // Skip infinite element-blocks; they can not be viewed in most visualization software
      // as paraview.
      if (mesh.elem_ref(tmp_vec[0]).infinite())
        continue;
#endif

      //Use the first element in this block to get representative information.
      //Note that Exodus assumes all elements in a block are of the same type!
      //We are using that same assumption here!
      const ExodusII_IO_Helper::Conversion conv =
        em.assign_conversion(mesh.elem_ref(tmp_vec[0]).type());
      num_nodes_per_elem = mesh.elem_ref(tmp_vec[0]).n_nodes();

      ex_err = exII::ex_put_elem_block(ex_id,
                                       (*it).first,
                                       conv.exodus_elem_type().c_str(),
                                       tmp_vec.size(),
                                       num_nodes_per_elem,
                                       /*num_attr=*/0);

      EX_CHECK_ERR(ex_err, "Error writing element block.");

      connect.resize(tmp_vec.size()*num_nodes_per_elem);

      for (std::size_t i=0; i<tmp_vec.size(); i++)
        {
          unsigned int elem_id = tmp_vec[i];
          libmesh_elem_num_to_exodus[elem_id] = ++libmesh_elem_num_to_exodus_counter; // 1-based indexing for Exodus

          const Elem & elem = mesh.elem_ref(elem_id);

          // We *might* be able to get away with writing mixed element
          // types which happen to have the same number of nodes, but
          // do we actually *want* to get away with that?
          // .) No visualization software would be able to handle it.
          // .) There'd be no way for us to read it back in reliably.
          // .) Even elements with the same number of nodes may have different connectivities (?)

          // This needs to be more than an assert so we don't fail
          // with a mysterious segfault while trying to write mixed
          // element meshes in optimized mode.
          if (elem.type() != conv.get_canonical_type())
            libmesh_error_msg("Error: Exodus requires all elements with a given subdomain ID to be the same type.\n" \
                              << "Can't write both "                  \
                              << Utility::enum_to_string(elem.type()) \
                              << " and "                              \
                              << Utility::enum_to_string(conv.get_canonical_type()) \
                              << " in the same block!");


          for (unsigned int j=0; j<static_cast<unsigned int>(num_nodes_per_elem); ++j)
            {
              unsigned connect_index   = (i*num_nodes_per_elem)+j;
              unsigned elem_node_index = conv.get_inverse_node_map(j); // inverse node map is for writing.
              if (verbose)
                {
                  libMesh::out << "Exodus node index " << j
                               << " = LibMesh node index " << elem_node_index << std::endl;
                }

              if (!use_discontinuous)
                {
                  // The global id for the current node in libmesh.
                  dof_id_type libmesh_node_id = elem.node_id(elem_node_index);

                  // Find the zero-based libmesh id in the map, this
                  // should be faster than doing linear searches on
                  // the node_num_map.
                  std::map<int, int>::iterator pos =
                    libmesh_node_num_to_exodus.find(cast_int<int>(libmesh_node_id));

                  // Make sure it was found.
                  if (pos == libmesh_node_num_to_exodus.end())
                    libmesh_error_msg("libmesh node id " << libmesh_node_id << " not found in node_num_map.");

                  // Write the Exodus global node id associated with
                  // this libmesh node number to the connectivity
                  // array.
                  connect[connect_index] = pos->second;
                }
              else
                {
                  // FIXME: We are hard-coding the 1-based node
                  // numbering assumption here, so writing
                  // "discontinuous" Exodus files won't work with node
                  // numberings that have "holes".
                  connect[connect_index] = node_counter + elem_node_index + 1;
                }
            }

          node_counter += num_nodes_per_elem;
        }

      ex_err = exII::ex_put_elem_conn(ex_id, (*it).first, &connect[0]);
      EX_CHECK_ERR(ex_err, "Error writing element connectivities");

      // This transform command stores its result in a range that begins at the third argument,
      // so this command is adding values to the elem_num_map vector starting from curr_elem_map_end.
      curr_elem_map_end = std::transform(tmp_vec.begin(),
                                         tmp_vec.end(),
                                         curr_elem_map_end,
                                         std::bind2nd(std::plus<subdomain_map_type::mapped_type::value_type>(), 1));  // Adds one to each id to make a 1-based exodus file!

      // But if we don't want to add one, we just want to put the values
      // of tmp_vec into elem_map in the right location, we can use
      // std::copy().
      // curr_elem_map_end = std::copy(tmp_vec.begin(), tmp_vec.end(), curr_elem_map_end);

      counter++;
    }

  // write out the element number map that we created
  ex_err = exII::ex_put_elem_num_map(ex_id, &elem_num_map[0]);
  EX_CHECK_ERR(ex_err, "Error writing element map");

  // Write out the block names
  if (num_elem_blk > 0)
    {
      ex_err = exII::ex_put_names(ex_id, exII::EX_ELEM_BLOCK, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing element names");
    }
}




void ExodusII_IO_Helper::write_sidesets(const MeshBase & mesh)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  ExodusII_IO_Helper::ElementMaps em;

  // Maps from sideset id to the element and sides
  std::map<int, std::vector<int> > elem;
  std::map<int, std::vector<int> > side;
  std::vector<boundary_id_type> side_boundary_ids;

  {
    std::vector< dof_id_type > el;
    std::vector< unsigned short int > sl;
    std::vector< boundary_id_type > il;

    mesh.get_boundary_info().build_side_list(el, sl, il);

    // Accumulate the vectors to pass into ex_put_side_set
    for (std::size_t i=0; i<el.size(); i++)
      {
        std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
        /**
         * We need to build up active elements if AMR is enabled and add
         * them to the exodus sidesets instead of the potentially inactive "parent" elements
         */
        mesh.elem_ref(el[i]).active_family_tree_by_side(family, sl[i], false);
#else
        family.push_back(mesh.elem_ptr(el[i]));
#endif

        for (std::size_t j=0; j<family.size(); ++j)
          {
            const ExodusII_IO_Helper::Conversion conv =
              em.assign_conversion(mesh.elem_ptr(family[j]->id())->type());

            // Use the libmesh to exodus datastructure map to get the proper sideset IDs
            // The datastructure contains the "collapsed" contiguous ids
            elem[il[i]].push_back(libmesh_elem_num_to_exodus[family[j]->id()]);
            side[il[i]].push_back(conv.get_inverse_side_map(sl[i]));
          }
      }

    mesh.get_boundary_info().build_side_boundary_ids(side_boundary_ids);
  }

  {
    // add data for shell faces, if needed

    std::vector< dof_id_type > el;
    std::vector< unsigned short int > sl;
    std::vector< boundary_id_type > il;

    mesh.get_boundary_info().build_shellface_list(el, sl, il);

    // Accumulate the vectors to pass into ex_put_side_set
    for (std::size_t i=0; i<el.size(); i++)
      {
        std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
        /**
         * We need to build up active elements if AMR is enabled and add
         * them to the exodus sidesets instead of the potentially inactive "parent" elements
         */
        mesh.elem_ref(el[i]).active_family_tree_by_side(family, sl[i], false);
#else
        family.push_back(mesh.elem_ptr(el[i]));
#endif

        for (std::size_t j=0; j<family.size(); ++j)
          {
            const ExodusII_IO_Helper::Conversion conv =
              em.assign_conversion(mesh.elem_ptr(family[j]->id())->type());

            // Use the libmesh to exodus datastructure map to get the proper sideset IDs
            // The datastructure contains the "collapsed" contiguous ids
            elem[il[i]].push_back(libmesh_elem_num_to_exodus[family[j]->id()]);
            side[il[i]].push_back(conv.get_inverse_shellface_map(sl[i]));
          }
      }

    std::vector<boundary_id_type> shellface_boundary_ids;
    mesh.get_boundary_info().build_shellface_boundary_ids(shellface_boundary_ids);
    for (std::size_t i=0; i<shellface_boundary_ids.size(); i++)
      side_boundary_ids.push_back(shellface_boundary_ids[i]);
  }

  // Write out the sideset names, but only if there is something to write
  if (side_boundary_ids.size() > 0)
    {
      NamesData names_table(side_boundary_ids.size(), MAX_STR_LENGTH);

      for (std::size_t i=0; i<side_boundary_ids.size(); i++)
        {
          boundary_id_type ss_id = side_boundary_ids[i];
          int actual_id = ss_id;

          names_table.push_back_entry(mesh.get_boundary_info().get_sideset_name(ss_id));

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

  mesh.get_boundary_info().build_node_list(nl, il);

  // Maps from nodeset id to the nodes
  std::map<boundary_id_type, std::vector<int> > node;

  // Accumulate the vectors to pass into ex_put_node_set
  for (std::size_t i=0; i<nl.size(); i++)
    node[il[i]].push_back(nl[i]+1);

  std::vector<boundary_id_type> node_boundary_ids;
  mesh.get_boundary_info().build_node_boundary_ids(node_boundary_ids);

  // Write out the nodeset names, but only if there is something to write
  if (node_boundary_ids.size() > 0)
    {
      NamesData names_table(node_boundary_ids.size(), MAX_STR_LENGTH);

      for (std::size_t i=0; i<node_boundary_ids.size(); i++)
        {
          boundary_id_type nodeset_id = node_boundary_ids[i];

          int actual_id = nodeset_id;

          names_table.push_back_entry(mesh.get_boundary_info().get_nodeset_name(nodeset_id));

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



void ExodusII_IO_Helper::initialize_element_variables(std::vector<std::string> names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Quick return if there are no element variables to write
  if (names.size() == 0)
    return;

  // Quick return if we have already called this function
  if (_elem_vars_initialized)
    return;

  // Be sure that variables in the file match what we are asking for
  if (num_elem_vars > 0)
    {
      this->check_existing_vars(ELEMENTAL, names, this->elem_var_names);
      return;
    }

  // Set the flag so we can skip this stuff on subsequent calls to
  // initialize_element_variables()
  _elem_vars_initialized = true;

  this->write_var_names(ELEMENTAL, names);

  // Form the element variable truth table and send to Exodus.
  // This tells which variables are written to which blocks,
  // and can dramatically speed up writing element variables
  //
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
}



void ExodusII_IO_Helper::initialize_nodal_variables(std::vector<std::string> names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Quick return if there are no nodal variables to write
  if (names.size() == 0)
    return;

  // Quick return if we have already called this function
  if (_nodal_vars_initialized)
    return;

  // Be sure that variables in the file match what we are asking for
  if (num_nodal_vars > 0)
    {
      this->check_existing_vars(NODAL, names, this->nodal_var_names);
      return;
    }

  // Set the flag so we can skip the rest of this function on subsequent calls.
  _nodal_vars_initialized = true;

  this->write_var_names(NODAL, names);
}



void ExodusII_IO_Helper::initialize_global_variables(std::vector<std::string> names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Quick return if there are no global variables to write
  if (names.size() == 0)
    return;

  if (_global_vars_initialized)
    return;

  // Be sure that variables in the file match what we are asking for
  if (num_global_vars > 0)
    {
      this->check_existing_vars(GLOBAL, names, this->global_var_names);
      return;
    }

  _global_vars_initialized = true;

  this->write_var_names(GLOBAL, names);
}



void ExodusII_IO_Helper::check_existing_vars(ExodusVarType type,
                                             std::vector<std::string> & names,
                                             std::vector<std::string> & names_from_file)
{
  // There may already be global variables in the file (for example,
  // if we're appending) and in that case, we
  // 1.) Cannot initialize them again.
  // 2.) Should check to be sure that the global variable names are the same.

  // Fills up names_from_file for us
  this->read_var_names(type);

  // Both the names of the global variables and their order must match
  if (names_from_file != names)
    {
      libMesh::err << "Error! The Exodus file already contains the variables:" << std::endl;
      for (std::size_t i=0; i<names_from_file.size(); ++i)
        libMesh::out << names_from_file[i] << std::endl;

      libMesh::err << "And you asked to write:" << std::endl;
      for (std::size_t i=0; i<names.size(); ++i)
        libMesh::out << names[i] << std::endl;

      libmesh_error_msg("Cannot overwrite existing variables in Exodus II file.");
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



void ExodusII_IO_Helper::write_element_values(const MeshBase & mesh, const std::vector<Real> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Loop over the element blocks and write the data one block at a time
  std::map<unsigned int, std::vector<unsigned int> > subdomain_map;

  // Ask the file how many element vars it has, store it in the num_elem_vars variable.
  ex_err = exII::ex_get_var_param(ex_id, "e", &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  MeshBase::const_element_iterator mesh_it = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();

  // loop through element and map between block and element vector
  for (; mesh_it!=end; ++mesh_it)
    {
      const Elem * elem = *mesh_it;
      subdomain_map[elem->subdomain_id()].push_back(elem->id());
    }

  // Use mesh.n_elem() to access into the values vector rather than
  // the number of elements the Exodus writer thinks the mesh has,
  // which may not include inactive elements.
  dof_id_type n_elem = mesh.n_elem();

  // For each variable, create a 'data' array which holds all the elemental variable
  // values *for a given block* on this processor, then write that data vector to file
  // before moving onto the next block.
  for (unsigned int i=0; i<static_cast<unsigned>(num_elem_vars); ++i)
    {
      // The size of the subdomain map is the number of blocks.
      std::map<unsigned int, std::vector<unsigned int> >::iterator it = subdomain_map.begin();

      for (unsigned int j=0; it!=subdomain_map.end(); ++it, ++j)
        {
          const std::vector<unsigned int> & elem_nums = (*it).second;
          const unsigned int num_elems_this_block =
            cast_int<unsigned int>(elem_nums.size());
          std::vector<Real> data(num_elems_this_block);

          for (unsigned int k=0; k<num_elems_this_block; ++k)
            data[k] = values[i*n_elem + elem_nums[k]];

          if (_single_precision)
            {
              std::vector<float> cast_data(data.begin(), data.end());

              ex_err = exII::ex_put_elem_var(ex_id,
                                             timestep,
                                             i+1,
                                             this->get_block_id(j),
                                             num_elems_this_block,
                                             &cast_data[0]);
            }
          else
            {
              ex_err = exII::ex_put_elem_var(ex_id,
                                             timestep,
                                             i+1,
                                             this->get_block_id(j),
                                             num_elems_this_block,
                                             &data[0]);
            }
          EX_CHECK_ERR(ex_err, "Error writing element values.");
        }
    }

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_nodal_values(int var_id, const std::vector<Real> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_single_precision)
    {
      std::vector<float> cast_values(values.begin(), values.end());
      ex_err = exII::ex_put_nodal_var(ex_id, timestep, var_id, num_nodes, &cast_values[0]);
    }
  else
    {
      ex_err = exII::ex_put_nodal_var(ex_id, timestep, var_id, num_nodes, &values[0]);
    }
  EX_CHECK_ERR(ex_err, "Error writing nodal values.");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::write_information_records(const std::vector<std::string> & records)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // There may already be information records in the file (for
  // example, if we're appending) and in that case, according to the
  // Exodus documentation, writing more information records is not
  // supported.
  int num_info = inquire(exII::EX_INQ_INFO, "Error retrieving the number of information records from file!");
  if (num_info > 0)
    {
      libMesh::err << "Warning! The Exodus file already contains information records.\n"
                   << "Exodus does not support writing additional records in this situation."
                   << std::endl;
      return;
    }

  int num_records = cast_int<int>(records.size());

  if (num_records > 0)
    {
      NamesData info(num_records, MAX_LINE_LENGTH);

      // If an entry is longer than MAX_LINE_LENGTH characters it's not an error, we just
      // write the first MAX_LINE_LENGTH characters to the file.
      for (std::size_t i=0; i<records.size(); ++i)
        info.push_back_entry(records[i]);

      ex_err = exII::ex_put_info(ex_id, num_records, info.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing global values.");

      ex_err = exII::ex_update(ex_id);
      EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
    }
}



void ExodusII_IO_Helper::write_global_values(const std::vector<Real> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_single_precision)
    {
      std::vector<float> cast_values(values.begin(), values.end());
      ex_err = exII::ex_put_glob_vars(ex_id, timestep, num_global_vars, &cast_values[0]);
    }
  else
    {
      ex_err = exII::ex_put_glob_vars(ex_id, timestep, num_global_vars, &values[0]);
    }
  EX_CHECK_ERR(ex_err, "Error writing global values.");

  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::use_mesh_dimension_instead_of_spatial_dimension(bool val)
{
  _use_mesh_dimension_instead_of_spatial_dimension = val;
}



void ExodusII_IO_Helper::write_as_dimension(unsigned dim)
{
  _write_as_dimension = dim;
}



void ExodusII_IO_Helper::set_coordinate_offset(Point p)
{
  _coordinate_offset = p;
}


std::vector<std::string> ExodusII_IO_Helper::get_complex_names(const std::vector<std::string> & names) const
{
  std::vector<std::string>::const_iterator names_it = names.begin();
  std::vector<std::string>::const_iterator names_end = names.end();

  std::vector<std::string> complex_names;

  // This will loop over all names and create new "complex" names
  // (i.e. names that start with r_, i_ or a_
  for (; names_it != names_end; ++names_it)
    {
      std::stringstream name_real, name_imag, name_abs;
      name_real << "r_" << *names_it;
      name_imag << "i_" << *names_it;
      name_abs << "a_" << *names_it;

      complex_names.push_back(name_real.str());
      complex_names.push_back(name_imag.str());
      complex_names.push_back(name_abs.str());
    }

  return complex_names;
}



// ------------------------------------------------------------
// ExodusII_IO_Helper::Conversion class members
ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(std::string type_str)
{
  init_element_equivalence_map();

  // Do only upper-case comparisons
  std::transform(type_str.begin(), type_str.end(), type_str.begin(), ::toupper);

  std::map<std::string, ElemType>::iterator it =
    element_equivalence_map.find(type_str);

  if (it != element_equivalence_map.end())
    return assign_conversion( it->second );
  else
    libmesh_error_msg("ERROR! Unrecognized element type_str: " << type_str);

  libmesh_error_msg("We'll never get here!");
  return assign_conversion (EDGE2);
}



ExodusII_IO_Helper::Conversion ExodusII_IO_Helper::ElementMaps::assign_conversion(const ElemType type)
{
  switch (type)
    {
    case NODEELEM:
      {
        const Conversion conv(nodeelem_node_map,
                              ARRAY_LENGTH(nodeelem_node_map),
                              nodeelem_node_map, // inverse node map same as forward node map
                              ARRAY_LENGTH(nodeelem_node_map),
                              libmesh_nullptr, // NODELEM doesn't have any edges
                              0,
                              libmesh_nullptr,
                              0,
                              NODEELEM, "SPHERE");
        return conv;
      }

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

    case QUADSHELL4:
      {
        return Conversion(quad4_node_map,
                          ARRAY_LENGTH(quad4_node_map), // node mapping is the same as for quad4
                          quad4_node_map,
                          ARRAY_LENGTH(quad4_node_map),
                          quadshell4_edge_map,
                          ARRAY_LENGTH(quadshell4_edge_map),
                          quadshell4_inverse_edge_map,
                          ARRAY_LENGTH(quadshell4_inverse_edge_map),
                          quadshell4_shellface_map,
                          ARRAY_LENGTH(quadshell4_shellface_map),
                          quadshell4_inverse_shellface_map,
                          ARRAY_LENGTH(quadshell4_inverse_shellface_map),
                          2, // the side index offset for QUADSHELL4 is 2
                          QUADSHELL4,
                          "SHELL4");
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
        return Conversion(tri3_node_map,
                          ARRAY_LENGTH(tri3_node_map),
                          tri3_node_map, // inverse node map same as forward node map
                          ARRAY_LENGTH(tri3_node_map),
                          tri_edge_map,
                          ARRAY_LENGTH(tri_edge_map),
                          tri_inverse_edge_map,
                          ARRAY_LENGTH(tri_inverse_edge_map),
                          TRI3,
                          "TRI3");
      }

    case TRISHELL3:
      {
        return Conversion(tri3_node_map,
                          ARRAY_LENGTH(tri3_node_map), // node mapping is the same as for tri3
                          tri3_node_map,
                          ARRAY_LENGTH(tri3_node_map),
                          trishell3_edge_map,
                          ARRAY_LENGTH(trishell3_edge_map),
                          trishell3_inverse_edge_map,
                          ARRAY_LENGTH(trishell3_inverse_edge_map),
                          trishell3_shellface_map,
                          ARRAY_LENGTH(trishell3_shellface_map),
                          trishell3_inverse_shellface_map,
                          ARRAY_LENGTH(trishell3_inverse_shellface_map),
                          2, // the side index offset for TRISHELL4 is 2
                          TRISHELL3,
                          "TRISHELL3");
      }

    case TRI3SUBDIVISION:
      {
        const Conversion conv(tri3_node_map,
                              ARRAY_LENGTH(tri3_node_map),
                              tri3_node_map, // inverse node map same as forward node map
                              ARRAY_LENGTH(tri3_node_map),
                              tri_edge_map,
                              ARRAY_LENGTH(tri_edge_map),
                              tri_inverse_edge_map,
                              ARRAY_LENGTH(tri_inverse_edge_map),
                              TRI3SUBDIVISION,
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

    case PYRAMID13:
      {
        const Conversion conv(pyramid13_node_map,
                              ARRAY_LENGTH(pyramid13_node_map),
                              pyramid13_node_map, // inverse node map same as forward node map
                              ARRAY_LENGTH(pyramid13_node_map),
                              pyramid_face_map,
                              ARRAY_LENGTH(pyramid_face_map),
                              pyramid_inverse_face_map,
                              ARRAY_LENGTH(pyramid_inverse_face_map),
                              PYRAMID13,
                              "PYRAMID13");
        return conv;
      }

    case PYRAMID14:
      {
        const Conversion conv(pyramid14_node_map,
                              ARRAY_LENGTH(pyramid14_node_map),
                              pyramid14_node_map, // inverse node map same as forward node map
                              ARRAY_LENGTH(pyramid14_node_map),
                              pyramid_face_map,
                              ARRAY_LENGTH(pyramid_face_map),
                              pyramid_inverse_face_map,
                              ARRAY_LENGTH(pyramid_inverse_face_map),
                              PYRAMID14,
                              "PYRAMID14");
        return conv;
      }

    default:
      libmesh_error_msg("Unsupported element type: " << type);
    }

  libmesh_error_msg("We'll never get here!");
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



int ExodusII_IO_Helper::Conversion::get_side_map(int i) const
{
  // If we asked for a side that doesn't exist, return an invalid_id
  // and allow higher-level code to handle it.
  if (static_cast<size_t>(i) >= side_map_size)
    return invalid_id;

  return side_map[i];
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



char ** ExodusII_IO_Helper::NamesData::get_char_star_star()
{
  return &data_table_pointers[0];
}



char * ExodusII_IO_Helper::NamesData::get_char_star(int i)
{
  if (static_cast<unsigned>(i) >= table_size)
    libmesh_error_msg("Requested char * " << i << " but only have " << table_size << "!");

  else
    return &(data_table[i][0]);
}


} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_EXODUS_API
