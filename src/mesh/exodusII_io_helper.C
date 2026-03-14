// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh includes
#include "libmesh/boundary_info.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fpe_disabler.h"
#include "libmesh/remote_elem.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"
#include "libmesh/libmesh_logging.h"

#ifdef DEBUG
#include "libmesh/mesh_tools.h"  // for elem_types warning
#endif

#include <libmesh/ignore_warnings.h>
namespace exII {
extern "C" {
#include "exodusII.h" // defines MAX_LINE_LENGTH, MAX_STR_LENGTH used later
}
}
#include <libmesh/restore_warnings.h>

// C++ includes
#include <algorithm>
#include <cfenv> // workaround for HDF5 bug
#include <cstdlib> // std::strtol
#include <sstream>
#include <unordered_map>

// Anonymous namespace for file local data and helper functions
namespace
{

// ExodusII defaults to 32 bytes names, but we've had user complaints
// about truncation with those.
// It looks like the maximum they'll support is 80 byte names.
static constexpr int libmesh_max_str_length = MAX_LINE_LENGTH;

using namespace libMesh;

// File scope constant node/edge/face mapping arrays.
// 2D inverse face map definitions.
// These take a libMesh ID and turn it into an Exodus ID
const std::vector<int> trishell3_inverse_edge_map = {3, 4, 5};
const std::vector<int> quadshell4_inverse_edge_map = {3, 4, 5, 6};

// 3D node map definitions
// The hex27, prism20-21, and tet14 appear to be the only elements
// with a non-identity mapping between Exodus' node numbering and
// libmesh's.  Exodus doesn't even number prisms hierarchically!
const std::vector<int> hex27_node_map = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  // Mid-face nodes and center node
  21, 25, 24, 26, 23, 22, 20};
//20  21  22  23  24  25  26 // LibMesh indices

const std::vector<int> hex27_inverse_node_map = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  // Mid-face nodes and center node
  26, 20, 25, 24, 22, 21, 23};
//20  21  22  23  24  25  26

const std::vector<int> prism20_node_map = {
  // Vertices
  0, 1, 2, 3, 4, 5,
  // Matching mid-edge nodes
  6, 7, 8, 9, 10, 11, 12, 13, 14,
  // Non-matching nodes
  19, 17, 18, 15, 16};
//15  16  17  18  19 // LibMesh indices

const std::vector<int> prism20_inverse_node_map = {
  // Vertices
  0, 1, 2, 3, 4, 5,
  // Matching mid-edge nodes
  6, 7, 8, 9, 10, 11, 12, 13, 14,
  // Non-matching nodes
  18, 19, 16, 17, 15};
//15  16  17  18  19

const std::vector<int> prism21_node_map = {
  // Vertices
  0, 1, 2, 3, 4, 5,
  // Matching mid-edge nodes
  6, 7, 8, 9, 10, 11, 12, 13, 14,
  // Non-matching nodes
  20, 18, 19, 16, 17, 15};
//15  16  17  18  19  20 // LibMesh indices

const std::vector<int> prism21_inverse_node_map = {
  // Vertices
  0, 1, 2, 3, 4, 5,
  // Matching mid-edge nodes
  6, 7, 8, 9, 10, 11, 12, 13, 14,
  // Non-matching nodes
  20, 18, 19, 16, 17, 15};
//15  16  17  18  19  20

const std::vector<int> tet14_node_map = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
  // Mid-face nodes
  10, 13, 11, 12};
//10  11  12  13 // LibMesh indices

const std::vector<int> tet14_inverse_node_map = {
  // Vertex and mid-edge nodes
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
  // Mid-face nodes
  10, 12, 13, 11};
//10  11  12  13


// 3D face map definitions
const std::vector<int> tet_face_map = {1, 2, 3, 0};
const std::vector<int> hex_face_map = {1, 2, 3, 4, 0, 5};
const std::vector<int> prism_face_map = {1, 2, 3, 0, 4};

// These take a libMesh ID and turn it into an Exodus ID
const std::vector<int> tet_inverse_face_map = {4, 1, 2, 3};
const std::vector<int> hex_inverse_face_map = {5, 1, 2, 3, 4, 6};
const std::vector<int> prism_inverse_face_map = {4, 1, 2, 3, 5};

// 3D element edge maps. Map 0-based Exodus id -> libMesh id.
// Commented out until we have code that needs it, to keep compiler
// warnings happy.
// const std::vector<int> hex_edge_map =
  // {0,1,2,3,8,9,10,11,4,5,7,6};

// 3D inverse element edge maps. Map libmesh edge ids to 1-based Exodus edge ids.
// Commented out until we have code that needs it, to keep compiler
// warnings happy.
// const std::vector<int> hex_inverse_edge_map =
  // {1,2,3,4,9,10,12,11,5,6,7,8};

  /**
   * \returns The value obtained from a generic exII::ex_inquire() call.
   */
  int inquire(libMesh::ExodusII_IO_Helper & e2h, exII::ex_inquiry req_info_in, std::string error_msg="")
  {
    int ret_int = 0;
    char ret_char = 0;
    float ret_float = 0.;

    e2h.ex_err = exII::ex_inquire(e2h.ex_id,
                                  req_info_in,
                                  &ret_int,
                                  &ret_float,
                                  &ret_char);

    EX_CHECK_ERR(e2h.ex_err, error_msg);

    return ret_int;
  }

  // Bezier Extraction test: if we see BEx data we had better be in a
  // Bezier element block
  inline bool is_bezier_elem(const char * elem_type_str)
  {
    // Reading Bezier Extraction from Exodus files requires ExodusII v8
#if EX_API_VERS_NODOT < 800
    libMesh::libmesh_ignore(elem_type_str);
    return false;
#else
    if (strlen(elem_type_str) <= 4)
      return false;
    return (std::string(elem_type_str, elem_type_str+4) == "BEX_");
#endif
  }

} // end anonymous namespace



namespace libMesh
{

// ExodusII_IO_Helper::Conversion static data
const int ExodusII_IO_Helper::Conversion::invalid_id = std::numeric_limits<int>::max();

ExodusII_IO_Helper::ExodusII_IO_Helper(const ParallelObject & parent,
                                       bool v,
                                       bool run_only_on_proc0,
                                       bool single_precision) :
  ParallelObject(parent),
  ex_id(0),
  ex_err(0),
  header_info(), // zero-initialize
  title(header_info.title),
  num_dim(header_info.num_dim),
  num_nodes(header_info.num_nodes),
  num_elem(header_info.num_elem),
  num_elem_blk(header_info.num_elem_blk),
  num_edge(header_info.num_edge),
  num_edge_blk(header_info.num_edge_blk),
  num_node_sets(header_info.num_node_sets),
  num_side_sets(header_info.num_side_sets),
  num_elem_sets(header_info.num_elem_sets),
  num_global_vars(0),
  num_sideset_vars(0),
  num_nodeset_vars(0),
  num_elemset_vars(0),
  num_elem_this_blk(0),
  num_nodes_per_elem(0),
  num_attr(0),
  num_elem_all_sidesets(0),
  num_elem_all_elemsets(0),
  bex_num_elem_cvs(0),
  num_time_steps(0),
  num_nodal_vars(0),
  num_elem_vars(0),
  verbose(v),
  set_unique_ids_from_maps(false),
  opened_for_writing(false),
  opened_for_reading(false),
  _run_only_on_proc0(run_only_on_proc0),
  _opened_by_create(false),
  _elem_vars_initialized(false),
  _global_vars_initialized(false),
  _nodal_vars_initialized(false),
  _use_mesh_dimension_instead_of_spatial_dimension(false),
  _write_hdf5(true),
  _max_name_length(32),
  _end_elem_id(0),
  _write_as_dimension(0),
  _single_precision(single_precision)
{
  title.resize(MAX_LINE_LENGTH+1);
  elem_type.resize(libmesh_max_str_length);
  init_element_equivalence_map();
  init_conversion_map();
}



ExodusII_IO_Helper::~ExodusII_IO_Helper() = default;



int ExodusII_IO_Helper::get_exodus_version()
{
  return EX_API_VERS_NODOT;
}



// Initialization function for conversion_map object
void ExodusII_IO_Helper::init_conversion_map()
{
  auto convert_type = [this](ElemType type,
                             std::string_view exodus_type,
                             const std::vector<int> * node_map = nullptr,
                             const std::vector<int> * inverse_node_map = nullptr,
                             const std::vector<int> * side_map = nullptr,
                             const std::vector<int> * inverse_side_map = nullptr,
                             const std::vector<int> * shellface_map = nullptr,
                             const std::vector<int> * inverse_shellface_map = nullptr,
                             size_t shellface_index_offset = 0)
  {
    std::unique_ptr<Elem> elem = Elem::build(type);
    auto & conv = conversion_map[elem->dim()][type];
    conv.libmesh_type = type;
    conv.exodus_type = exodus_type;
    conv.node_map = node_map;
    conv.inverse_node_map = inverse_node_map;
    conv.side_map = side_map;
    conv.inverse_side_map = inverse_side_map;
    conv.shellface_map = shellface_map;
    conv.inverse_shellface_map = inverse_shellface_map;
    conv.shellface_index_offset = shellface_index_offset;
    conv.n_nodes = elem->n_nodes();
    for (int d = elem->dim()+1; d <= 3; ++d)
      conversion_map[d][type] = conv;
  };

  convert_type(NODEELEM, "SPHERE");
  convert_type(EDGE2, "EDGE2");
  convert_type(EDGE3, "EDGE3");
  convert_type(EDGE4, "EDGE4");
  convert_type(QUAD4, "QUAD4");
  convert_type(QUAD8, "QUAD8");
  convert_type(QUAD9, "QUAD9");
  convert_type(QUADSHELL4, "SHELL4", nullptr, nullptr, nullptr,
               /* inverse_side_map = */ &quadshell4_inverse_edge_map,
               nullptr, nullptr, /* shellface_index_offset = */ 2);
  convert_type(QUADSHELL8, "SHELL8", nullptr, nullptr, nullptr,
               /* inverse_side_map = */ &quadshell4_inverse_edge_map,
               nullptr, nullptr, /* shellface_index_offset = */ 2);
  convert_type(QUADSHELL9, "SHELL9", nullptr, nullptr, nullptr,
               /* inverse_side_map = */ &quadshell4_inverse_edge_map,
               nullptr, nullptr, /* shellface_index_offset = */ 2);

  convert_type(TRI3, "TRI3");
  convert_type(TRI6, "TRI6");
  convert_type(TRI7, "TRI7");
  // Exodus does weird things to triangle side mapping in 3D.  See
  // https://sandialabs.github.io/seacas-docs/html/element_types.html#tri
  conversion_map[3][TRI3].inverse_side_map = &trishell3_inverse_edge_map;
  conversion_map[3][TRI3].shellface_index_offset = 2;
  conversion_map[3][TRI6].inverse_side_map = &trishell3_inverse_edge_map;
  conversion_map[3][TRI6].shellface_index_offset = 2;
  conversion_map[3][TRI7].inverse_side_map = &trishell3_inverse_edge_map;
  conversion_map[3][TRI7].shellface_index_offset = 2;

  convert_type(TRISHELL3, "TRISHELL3", nullptr, nullptr, nullptr,
               /* inverse_side_map = */ &trishell3_inverse_edge_map,
               nullptr, nullptr, /* shellface_index_offset = */ 2);
  convert_type(TRI3SUBDIVISION, "TRI3");
  convert_type(HEX8, "HEX8", nullptr, nullptr,
               &hex_face_map, &hex_inverse_face_map);
  convert_type(HEX20, "HEX20", nullptr, nullptr,
               &hex_face_map, &hex_inverse_face_map);
  convert_type(HEX27, "HEX27", &hex27_node_map,
               &hex27_inverse_node_map,
               &hex_face_map, &hex_inverse_face_map);
  convert_type(TET4, "TETRA4", nullptr, nullptr,
               &tet_face_map, &tet_inverse_face_map);
  convert_type(TET10, "TETRA10", nullptr, nullptr,
               &tet_face_map, &tet_inverse_face_map);
  convert_type(TET14, "TETRA14", &tet14_node_map,
               &tet14_inverse_node_map,
               &tet_face_map, &tet_inverse_face_map);
  convert_type(PRISM6, "WEDGE", nullptr, nullptr,
               &prism_face_map, &prism_inverse_face_map);
  convert_type(PRISM15, "WEDGE15", nullptr, nullptr,
               &prism_face_map, &prism_inverse_face_map);
  convert_type(PRISM18, "WEDGE18", nullptr, nullptr,
               &prism_face_map, &prism_inverse_face_map);
  convert_type(PRISM20, "WEDGE20", &prism20_node_map,
               &prism20_inverse_node_map,
               &prism_face_map, &prism_inverse_face_map);
  convert_type(PRISM21, "WEDGE21", &prism21_node_map,
               &prism21_inverse_node_map,
               &prism_face_map, &prism_inverse_face_map);
  convert_type(PYRAMID5, "PYRAMID5");
  convert_type(PYRAMID13, "PYRAMID13");
  convert_type(PYRAMID14, "PYRAMID14");
  convert_type(PYRAMID18, "PYRAMID18");
}



// This function initializes the element_equivalence_map the first time it
// is called, and returns early all other times.
void ExodusII_IO_Helper::init_element_equivalence_map()
{
  // We use an ExodusII SPHERE element to represent a NodeElem
  element_equivalence_map["SPHERE"] = NODEELEM;

  // EDGE2 equivalences
  element_equivalence_map["EDGE"]   = EDGE2;
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

  // EDGE4 equivalences
  element_equivalence_map["EDGE4"]  = EDGE4;
  element_equivalence_map["TRUSS4"] = EDGE4;
  element_equivalence_map["BEAM4"]  = EDGE4;
  element_equivalence_map["BAR4"]   = EDGE4;

  // This whole design is going to need to be refactored whenever we
  // support higher-order IGA, with one element type having variable
  // polynomiaal degree...
  element_equivalence_map["BEX_CURVE"] = EDGE3;

  // QUAD4 equivalences
  element_equivalence_map["QUAD"]   = QUAD4;
  element_equivalence_map["QUAD4"]  = QUAD4;

  // QUADSHELL4 equivalences
  element_equivalence_map["SHELL"]  = QUADSHELL4;
  element_equivalence_map["SHELL4"] = QUADSHELL4;

  // QUAD8 equivalences
  element_equivalence_map["QUAD8"]  = QUAD8;

  // QUADSHELL8 equivalences
  element_equivalence_map["SHELL8"] = QUADSHELL8;

  // QUAD9 equivalences
  element_equivalence_map["QUAD9"]  = QUAD9;
  // This only supports p==2 IGA:
  element_equivalence_map["BEX_QUAD"]  = QUAD9;

  // QUADSHELL9 equivalences
  element_equivalence_map["SHELL9"] = QUADSHELL9;

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
  // This only supports p==2 IGA:
  element_equivalence_map["BEX_TRIANGLE"] = TRI6;

  // TRI7 equivalences
  element_equivalence_map["TRI7"]      = TRI7;

  // HEX8 equivalences
  element_equivalence_map["HEX"]  = HEX8;
  element_equivalence_map["HEX8"] = HEX8;

  // HEX20 equivalences
  element_equivalence_map["HEX20"] = HEX20;

  // HEX27 equivalences
  element_equivalence_map["HEX27"] = HEX27;
  // This only supports p==2 IGA:
  element_equivalence_map["BEX_HEX"] = HEX27;

  // TET4 equivalences
  element_equivalence_map["TETRA"]  = TET4;
  element_equivalence_map["TETRA4"] = TET4;

  // TET10 equivalences
  element_equivalence_map["TETRA10"] = TET10;
  // This only supports p==2 IGA:
  element_equivalence_map["BEX_TETRA"] = TET10;

  // TET14 (in Exodus 8) equivalence
  element_equivalence_map["TETRA14"] = TET14;

  // PRISM6 equivalences
  element_equivalence_map["WEDGE"] = PRISM6;
  element_equivalence_map["WEDGE6"] = PRISM6;

  // PRISM15 equivalences
  element_equivalence_map["WEDGE15"] = PRISM15;

  // PRISM18 equivalences
  element_equivalence_map["WEDGE18"] = PRISM18;
  // This only supports p==2 IGA:
  element_equivalence_map["BEX_WEDGE"] = PRISM18;

  // PRISM20 equivalences
  element_equivalence_map["WEDGE20"] = PRISM20;

  // PRISM21 equivalences
  element_equivalence_map["WEDGE21"] = PRISM21;

  // PYRAMID equivalences
  element_equivalence_map["PYRAMID"]  = PYRAMID5;
  element_equivalence_map["PYRAMID5"] = PYRAMID5;
  element_equivalence_map["PYRAMID13"] = PYRAMID13;
  element_equivalence_map["PYRAMID14"] = PYRAMID14;
  element_equivalence_map["PYRAMID18"] = PYRAMID18;
}

const ExodusII_IO_Helper::Conversion &
ExodusII_IO_Helper::get_conversion(const ElemType type) const
{
  auto & maps_for_dim = libmesh_map_find(conversion_map, this->num_dim);
  return libmesh_map_find(maps_for_dim, type);
}

const ExodusII_IO_Helper::Conversion &
ExodusII_IO_Helper::get_conversion(std::string type_str) const
{
  // Do only upper-case comparisons
  std::transform(type_str.begin(), type_str.end(), type_str.begin(), ::toupper);
  return get_conversion (libmesh_map_find(element_equivalence_map, type_str));
}

const char * ExodusII_IO_Helper::get_elem_type() const
{
  return elem_type.data();
}



void ExodusII_IO_Helper::message(std::string_view msg)
{
  if (verbose) libMesh::out << msg << std::endl;
}



void ExodusII_IO_Helper::message(std::string_view msg, int i)
{
  if (verbose) libMesh::out << msg << i << "." << std::endl;
}


ExodusII_IO_Helper::MappedOutputVector::
MappedOutputVector(const std::vector<Real> & our_data_in,
                   bool single_precision_in)
  : our_data(our_data_in),
    single_precision(single_precision_in)
{
  if (single_precision)
    {
      if (sizeof(Real) != sizeof(float))
        {
          float_vec.resize(our_data.size());
          // boost float128 demands explicit downconversions
          for (std::size_t i : index_range(our_data))
            float_vec[i] = float(our_data[i]);
        }
    }

  else if (sizeof(Real) != sizeof(double))
    {
      double_vec.resize(our_data.size());
      // boost float128 demands explicit downconversions
      for (std::size_t i : index_range(our_data))
        double_vec[i] = double(our_data[i]);
    }
}

void *
ExodusII_IO_Helper::MappedOutputVector::data()
{
  if (single_precision)
    {
      if (sizeof(Real) != sizeof(float))
        return static_cast<void*>(float_vec.data());
    }

  else if (sizeof(Real) != sizeof(double))
    return static_cast<void*>(double_vec.data());

  // Otherwise return a (suitably casted) pointer to the original underlying data.
  return const_cast<void *>(static_cast<const void *>(our_data.data()));
}

ExodusII_IO_Helper::MappedInputVector::
MappedInputVector(std::vector<Real> & our_data_in,
                  bool single_precision_in)
  : our_data(our_data_in),
    single_precision(single_precision_in)
{
  // Allocate temporary space to store enough floats/doubles, if required.
  if (single_precision)
    {
      if (sizeof(Real) != sizeof(float))
        float_vec.resize(our_data.size());
    }
  else if (sizeof(Real) != sizeof(double))
    double_vec.resize(our_data.size());
}

ExodusII_IO_Helper::MappedInputVector::
~MappedInputVector()
{
  if (single_precision)
    {
      if (sizeof(Real) != sizeof(float))
        our_data.assign(float_vec.begin(), float_vec.end());
    }
  else if (sizeof(Real) != sizeof(double))
    our_data.assign(double_vec.begin(), double_vec.end());
}

void *
ExodusII_IO_Helper::MappedInputVector::data()
{
  if (single_precision)
    {
      if (sizeof(Real) != sizeof(float))
        return static_cast<void*>(float_vec.data());
    }

  else if (sizeof(Real) != sizeof(double))
    return static_cast<void*>(double_vec.data());

  // Otherwise return a (suitably casted) pointer to the original underlying data.
  return static_cast<void *>(our_data.data());
}

void ExodusII_IO_Helper::open(const char * filename, bool read_only)
{
  // Version of Exodus you are using
  float ex_version = 0.;

  int comp_ws = 0;

  if (_single_precision)
    comp_ws = cast_int<int>(sizeof(float));

  // Fall back on double precision when necessary since ExodusII
  // doesn't seem to support long double
  else
    comp_ws = cast_int<int>(std::min(sizeof(Real), sizeof(double)));

  // Word size in bytes of the floating point data as they are stored
  // in the ExodusII file.  "If this argument is 0, the word size of the
  // floating point data already stored in the file is returned"
  int io_ws = 0;

  {
    FPEDisabler disable_fpes;
    ex_id = exII::ex_open(filename,
                          read_only ? EX_READ : EX_WRITE,
                          &comp_ws,
                          &io_ws,
                          &ex_version);
  }

  std::string err_msg = std::string("Error opening ExodusII mesh file: ") + std::string(filename);
  EX_CHECK_ERR(ex_id, err_msg);
  if (verbose) libMesh::out << "File opened successfully." << std::endl;

  // If we're writing then we'll want to use the specified length;
  // if we're reading then we'll override this by what's in the file.
  int max_name_length_to_set = _max_name_length;

  if (read_only)
  {
    opened_for_reading = true;

    // ExodusII reads truncate to 32-char strings by default; we'd
    // like to support whatever's in the file, so as early as possible
    // let's find out what that is.
    int max_name_length = exII::ex_inquire_int(ex_id, exII::EX_INQ_DB_MAX_USED_NAME_LENGTH);

    libmesh_error_msg_if(max_name_length > MAX_LINE_LENGTH,
                         "Unexpected maximum name length of " <<
                         max_name_length << " in file " << filename <<
                         " exceeds expected " << MAX_LINE_LENGTH);

    // I don't think the 32 here should be necessary, but let's make
    // sure we don't accidentally make things *worse* for anyone.
    max_name_length_to_set = std::max(max_name_length, 32);
  }
  else
    opened_for_writing = true;

  ex_err = exII::ex_set_max_name_length(ex_id, max_name_length_to_set);
  EX_CHECK_ERR(ex_err, "Error setting max ExodusII name length.");

  current_filename = std::string(filename);
}



ExodusHeaderInfo
ExodusII_IO_Helper::read_header() const
{
  // Read init params using newer API that reads into a struct.  For
  // backwards compatibility, assign local member values from struct
  // afterwards. Note: using the new API allows us to automatically
  // read edge and face block/set information if it's present in the
  // file.
  exII::ex_init_params params = {};
  int err_flag = exII::ex_get_init_ext(ex_id, &params);
  EX_CHECK_ERR(err_flag, "Error retrieving header info.");

  // Extract required data into our struct
  ExodusHeaderInfo h;
  h.title.assign(params.title, params.title + MAX_LINE_LENGTH);
  h.num_dim = params.num_dim;
  h.num_nodes = params.num_nodes;
  h.num_elem = params.num_elem;
  h.num_elem_blk = params.num_elem_blk;
  h.num_node_sets = params.num_node_sets;
  h.num_side_sets = params.num_side_sets;
  h.num_elem_sets = params.num_elem_sets;
  h.num_edge_blk = params.num_edge_blk;
  h.num_edge = params.num_edge;

  // And return it
  return h;
}



void ExodusII_IO_Helper::read_and_store_header_info()
{
  // Read header params from file, storing them in this class's
  // ExodusHeaderInfo struct.  This automatically updates the local
  // num_dim, num_elem, etc. references.
  this->header_info = this->read_header();

  // Read the number of timesteps which are present in the file
  this->read_num_time_steps();

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_NODAL, &num_nodal_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of nodal variables.");

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_ELEM_BLOCK, &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_GLOBAL, &num_global_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of global variables.");

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_SIDE_SET, &num_sideset_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of sideset variables.");

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_NODE_SET, &num_nodeset_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of nodeset variables.");

  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_ELEM_SET, &num_elemset_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemset variables.");

  message("Exodus header info retrieved successfully.");
}




void ExodusII_IO_Helper::read_qa_records()
{
  // The QA records are four MAX_STR_LENGTH-byte character strings.
  int num_qa_rec =
    inquire(*this, exII::EX_INQ_QA, "Error retrieving number of QA records");

  if (verbose)
    libMesh::out << "Found "
                 << num_qa_rec
                 << " QA record(s) in the Exodus file."
                 << std::endl;

  if (num_qa_rec > 0)
    {
      // Actual (num_qa_rec x 4) storage for strings. The object we
      // pass to the Exodus API will just contain pointers into the
      // qa_storage object, which will have all automatic memory
      // management.
      std::vector<std::vector<std::vector<char>>> qa_storage(num_qa_rec);
      for (auto i : make_range(num_qa_rec))
        {
          qa_storage[i].resize(4);
          for (auto j : make_range(4))
            qa_storage[i][j].resize(libmesh_max_str_length+1);
        }

      // inner_array_t is a fixed-size array of 4 strings
      typedef char * inner_array_t[4];

      // There is at least one compiler (Clang 12.0.1) that complains about
      // "a non-scalar type used in a pseudo-destructor expression" when
      // we try to instantiate a std::vector of inner_array_t objects as in:
      // std::vector<inner_array_t> qa_record(num_qa_rec);
      // So, we instead attempt to achieve the same effect with a std::unique_ptr.
      auto qa_record = std::make_unique<inner_array_t[]>(num_qa_rec);

      // Create data structure to be passed to Exodus API by setting
      // pointers to the actual strings which are in qa_storage.
      for (auto i : make_range(num_qa_rec))
        for (auto j : make_range(4))
          qa_record[i][j] = qa_storage[i][j].data();

      ex_err = exII::ex_get_qa (ex_id, qa_record.get());
      EX_CHECK_ERR(ex_err, "Error reading the QA records.");

      // Print the QA records
      if (verbose)
        {
          for (auto i : make_range(num_qa_rec))
            {
              libMesh::out << "QA Record: " << i << std::endl;
              for (auto j : make_range(4))
                libMesh::out << qa_record[i][j] << std::endl;
            }
        }
    }
}




void ExodusII_IO_Helper::print_header()
{
  if (verbose)
    libMesh::out << "Title: \t" << title.data() << std::endl
                 << "Mesh Dimension: \t"   << num_dim << std::endl
                 << "Number of Nodes: \t" << num_nodes << std::endl
                 << "Number of elements: \t" << num_elem << std::endl
                 << "Number of elt blocks: \t" << num_elem_blk << std::endl
                 << "Number of node sets: \t" << num_node_sets << std::endl
                 << "Number of side sets: \t" << num_side_sets << std::endl
                 << "Number of elem sets: \t" << num_elem_sets << std::endl;
}



void ExodusII_IO_Helper::read_nodes()
{
  LOG_SCOPE("read_nodes()", "ExodusII_IO_Helper");

  x.resize(num_nodes);
  y.resize(num_nodes);
  z.resize(num_nodes);

  if (num_nodes)
    {
      ex_err = exII::ex_get_coord
        (ex_id,
         MappedInputVector(x, _single_precision).data(),
         MappedInputVector(y, _single_precision).data(),
         MappedInputVector(z, _single_precision).data());

      EX_CHECK_ERR(ex_err, "Error retrieving nodal data.");
      message("Nodal data retrieved successfully.");
    }

  // If a nodal attribute bex_weight exists, we get spline weights
  // from it
  int n_nodal_attr = 0;
  ex_err = exII::ex_get_attr_param(ex_id, exII::EX_NODAL, 0, & n_nodal_attr);
  EX_CHECK_ERR(ex_err, "Error getting number of nodal attributes.");

  if (n_nodal_attr > 0)
    {
      std::vector<std::vector<char>> attr_name_data
        (n_nodal_attr, std::vector<char>(libmesh_max_str_length + 1));
      std::vector<char *> attr_names(n_nodal_attr);
      for (auto i : index_range(attr_names))
        attr_names[i] = attr_name_data[i].data();

      ex_err = exII::ex_get_attr_names(ex_id, exII::EX_NODAL, 0, attr_names.data());
      EX_CHECK_ERR(ex_err, "Error getting nodal attribute names.");

      for (auto i : index_range(attr_names))
        if (std::string("bex_weight") == attr_names[i])
          {
            w.resize(num_nodes);
            ex_err =
              exII::ex_get_one_attr (ex_id, exII::EX_NODAL, 0, i+1,
                                     MappedInputVector(w, _single_precision).data());
            EX_CHECK_ERR(ex_err, "Error getting Bezier Extraction nodal weights");
          }
    }
}



void ExodusII_IO_Helper::read_node_num_map ()
{
  node_num_map.resize(num_nodes);

  // Note: we cannot use the exII::ex_get_num_map() here because it
  // (apparently) does not behave like ex_get_node_num_map() when
  // there is no node number map in the file: it throws an error
  // instead of returning a default identity array (1,2,3,...).
  ex_err = exII::ex_get_node_num_map
    (ex_id, node_num_map.empty() ? nullptr : node_num_map.data());

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


void ExodusII_IO_Helper::read_bex_cv_blocks()
{
  // If a bex blob exists, we look for Bezier Extraction coefficient
  // data there.

  // These APIs require newer Exodus than 5.22
#if EX_API_VERS_NODOT >= 800
  int n_blobs = exII::ex_inquire_int(ex_id, exII::EX_INQ_BLOB);

  if (n_blobs > 0)
    {
      std::vector<exII::ex_blob> blobs(n_blobs);
      std::vector<std::vector<char>> blob_names(n_blobs);
      for (auto i : make_range(n_blobs))
        {
          blob_names[i].resize(libmesh_max_str_length+1);
          blobs[i].name = blob_names[i].data();
        }

      ex_err = exII::ex_get_blobs(ex_id, blobs.data());
      EX_CHECK_ERR(ex_err, "Error getting blobs.");

      bool found_blob = false;
      const exII::ex_blob * my_blob = &blobs[0];
      for (const auto & blob : blobs)
        {
          if (std::string("bex_cv_blob") == blob.name)
            {
              found_blob = true;
              my_blob = &blob;
            }
        }

      if (!found_blob)
        libmesh_error_msg("Found no bex_cv_blob for bezier elements");

      const int n_blob_attr =
        exII::ex_get_attribute_count(ex_id, exII::EX_BLOB,
                                     my_blob->id);

      std::vector<exII::ex_attribute> attributes(n_blob_attr);
      ex_err = exII::ex_get_attribute_param(ex_id, exII::EX_BLOB,
                                            my_blob->id,
                                            attributes.data());
      EX_CHECK_ERR(ex_err, "Error getting bex blob attribute parameters.");

      int bex_num_dense_cv_blocks = 0;
      std::vector<int> bex_dense_cv_info;
      for (auto & attr : attributes)
        {
          if (std::string("bex_dense_cv_info") == attr.name)
            {
              const std::size_t value_count = attr.value_count;
              if (value_count % 2)
                libmesh_error_msg("Found odd number of bex_dense_cv_info");

              bex_dense_cv_info.resize(value_count);
              attr.values = bex_dense_cv_info.data();
              exII::ex_get_attribute(ex_id, &attr);

              bex_num_dense_cv_blocks = value_count / 2;

              libmesh_error_msg_if(bex_num_dense_cv_blocks > 1,
                                   "Found more than 1 dense bex CV block; unsure how to handle that");
            }
        }

      if (bex_dense_cv_info.empty())
        libmesh_error_msg("No bex_dense_cv_info found");

      int n_blob_vars;
      exII::ex_get_variable_param(ex_id, exII::EX_BLOB, &n_blob_vars);
      std::vector<char> var_name (libmesh_max_str_length + 1);
      for (auto v_id : make_range(1,n_blob_vars+1))
        {
          ex_err = exII::ex_get_variable_name(ex_id, exII::EX_BLOB, v_id, var_name.data());
          EX_CHECK_ERR(ex_err, "Error reading bex blob var name.");

          if (std::string("bex_dense_cv_blocks") == var_name.data())
            {
              std::vector<double> bex_dense_cv_blocks(my_blob->num_entry);

              ex_err = exII::ex_get_var(ex_id, 1, exII::EX_BLOB, v_id,
                                        my_blob->id, my_blob->num_entry,
                                        bex_dense_cv_blocks.data());
              EX_CHECK_ERR(ex_err, "Error reading bex_dense_cv_blocks.");

              bex_dense_constraint_vecs.clear();
              bex_dense_constraint_vecs.resize(bex_num_dense_cv_blocks);

              std::size_t offset = 0;
              for (auto i : IntRange<std::size_t>(0, bex_num_dense_cv_blocks))
                {
                  bex_dense_constraint_vecs[i].resize(bex_dense_cv_info[2*i]);
                  const int vecsize = bex_dense_cv_info[2*i+1];
                  for (auto & vec : bex_dense_constraint_vecs[i])
                    {
                      vec.resize(vecsize);
                      std::copy(std::next(bex_dense_cv_blocks.begin(), offset),
                                std::next(bex_dense_cv_blocks.begin(), offset + vecsize),
                                vec.begin());
                      offset += vecsize;
                    }
                }
              libmesh_assert(offset == bex_dense_cv_blocks.size());
            }
        }
    }
#endif // EX_API_VERS_NODOT >= 800
}


void ExodusII_IO_Helper::print_nodes(std::ostream & out_stream)
{
  for (int i=0; i<num_nodes; i++)
    out_stream << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
}



void ExodusII_IO_Helper::read_block_info()
{
  if (num_elem_blk)
    {
      // Read all element block IDs.
      block_ids.resize(num_elem_blk);
      ex_err = exII::ex_get_ids(ex_id,
                                exII::EX_ELEM_BLOCK,
                                block_ids.data());

      EX_CHECK_ERR(ex_err, "Error getting block IDs.");
      message("All block IDs retrieved successfully.");

      char name_buffer[libmesh_max_str_length+1];
      for (int i=0; i<num_elem_blk; ++i)
        {
          ex_err = exII::ex_get_name(ex_id, exII::EX_ELEM_BLOCK,
                                     block_ids[i], name_buffer);
          EX_CHECK_ERR(ex_err, "Error getting block name.");
          id_to_block_names[block_ids[i]] = name_buffer;
        }
      message("All block names retrieved successfully.");
    }

  if (num_edge_blk)
    {
      // Read all edge block IDs.
      edge_block_ids.resize(num_edge_blk);
      ex_err = exII::ex_get_ids(ex_id,
                                exII::EX_EDGE_BLOCK,
                                edge_block_ids.data());

      EX_CHECK_ERR(ex_err, "Error getting edge block IDs.");
      message("All edge block IDs retrieved successfully.");

      // Read in edge block names
      char name_buffer[libmesh_max_str_length+1];
      for (int i=0; i<num_edge_blk; ++i)
        {
          ex_err = exII::ex_get_name(ex_id, exII::EX_EDGE_BLOCK,
                                     edge_block_ids[i], name_buffer);
          EX_CHECK_ERR(ex_err, "Error getting block name.");
          id_to_edge_block_names[edge_block_ids[i]] = name_buffer;
        }
      message("All edge block names retrieved successfully.");
    }
}



int ExodusII_IO_Helper::get_block_id(int index)
{
  libmesh_assert_less (index, block_ids.size());

  return block_ids[index];
}



std::string ExodusII_IO_Helper::get_block_name(int index)
{
  libmesh_assert_less (index, block_ids.size());

  return id_to_block_names[block_ids[index]];
}



int ExodusII_IO_Helper::get_side_set_id(int index)
{
  libmesh_assert_less (index, ss_ids.size());

  return ss_ids[index];
}



std::string ExodusII_IO_Helper::get_side_set_name(int index)
{
  libmesh_assert_less (index, ss_ids.size());

  return id_to_ss_names[ss_ids[index]];
}



int ExodusII_IO_Helper::get_node_set_id(int index)
{
  libmesh_assert_less (index, nodeset_ids.size());

  return nodeset_ids[index];
}



std::string ExodusII_IO_Helper::get_node_set_name(int index)
{
  libmesh_assert_less (index, nodeset_ids.size());

  return id_to_ns_names[nodeset_ids[index]];
}




void ExodusII_IO_Helper::read_elem_in_block(int block)
{
  LOG_SCOPE("read_elem_in_block()", "ExodusII_IO_Helper");

  libmesh_assert_less (block, block_ids.size());

  // Unlike the other "extended" APIs, this one does not use a parameter struct.
  int num_edges_per_elem = 0;
  int num_faces_per_elem = 0;
  int num_node_data_per_elem = 0;
  ex_err = exII::ex_get_block(ex_id,
                              exII::EX_ELEM_BLOCK,
                              block_ids[block],
                              elem_type.data(),
                              &num_elem_this_blk,
                              &num_node_data_per_elem,
                              &num_edges_per_elem, // 0 or -1 if no "extended" block info
                              &num_faces_per_elem, // 0 or -1 if no "extended" block info
                              &num_attr);

  EX_CHECK_ERR(ex_err, "Error getting block info.");
  message("Info retrieved successfully for block: ", block);

  // Warn that we don't currently support reading blocks with extended info.
  // Note: the docs say -1 will be returned for this but I found that it was
  // actually 0, so not sure which it will be in general.
  if (!(num_edges_per_elem == 0) && !(num_edges_per_elem == -1))
    libmesh_warning("Exodus files with extended edge connectivity not currently supported.");
  if (!(num_faces_per_elem == 0) && !(num_faces_per_elem == -1))
    libmesh_warning("Exodus files with extended face connectivity not currently supported.");

  // If we have a Bezier element here, then we've packed constraint
  // vector connectivity at the end of the nodal connectivity, and
  // num_nodes_per_elem reflected both.
  const bool is_bezier = is_bezier_elem(elem_type.data());
  if (is_bezier)
    {
      const auto & conv = get_conversion(std::string(elem_type.data()));
      num_nodes_per_elem = conv.n_nodes;
    }
  else
    num_nodes_per_elem = num_node_data_per_elem;

  if (verbose)
    libMesh::out << "Read a block of " << num_elem_this_blk
                 << " " << elem_type.data() << "(s)"
                 << " having " << num_nodes_per_elem
                 << " nodes per element." << std::endl;

  // Read in the connectivity of the elements of this block,
  // watching out for the case where we actually have no
  // elements in this block (possible with parallel files)
  connect.resize(num_node_data_per_elem*num_elem_this_blk);

  if (!connect.empty())
    {
      ex_err = exII::ex_get_conn(ex_id,
                                 exII::EX_ELEM_BLOCK,
                                 block_ids[block],
                                 connect.data(), // node_conn
                                 nullptr,        // elem_edge_conn (unused)
                                 nullptr);       // elem_face_conn (unused)

      EX_CHECK_ERR(ex_err, "Error reading block connectivity.");
      message("Connectivity retrieved successfully for block: ", block);
    }

  // If we had any attributes for this block, check to see if some of
  // them were Bezier-extension attributes.

  // num_attr above is zero, not actually the number of block attributes?
  // ex_get_attr_param *also* gives me zero?  Really, Exodus?
#if EX_API_VERS_NODOT >= 800
  int real_n_attr = exII::ex_get_attribute_count(ex_id, exII::EX_ELEM_BLOCK, block_ids[block]);
  EX_CHECK_ERR(real_n_attr, "Error getting number of element block attributes.");

  if (real_n_attr > 0)
    {
      std::vector<exII::ex_attribute> attributes(real_n_attr);

      ex_err = exII::ex_get_attribute_param(ex_id, exII::EX_ELEM_BLOCK, block_ids[block], attributes.data());
      EX_CHECK_ERR(ex_err, "Error getting element block attribute parameters.");

      ex_err = exII::ex_get_attributes(ex_id, real_n_attr, attributes.data());
      EX_CHECK_ERR(ex_err, "Error getting element block attribute values.");

      for (auto attr : attributes)
        {
          if (std::string("bex_elem_degrees") == attr.name)
            {
              if (attr.type != exII::EX_INTEGER)
                libmesh_error_msg("Found non-integer bex_elem_degrees");

              if (attr.value_count > 3)
                libmesh_error_msg("Looking for at most 3 bex_elem_degrees; found " << attr.value_count);

              libmesh_assert(is_bezier);

              std::vector<int> bex_elem_degrees(3); // max dim

              const int * as_int = static_cast<int *>(attr.values);
              std::copy(as_int, as_int+attr.value_count, bex_elem_degrees.begin());


              // Right now Bezier extraction elements aren't possible
              // for p>2 and aren't useful for p<2, and we don't
              // support anisotropic p...
#ifndef NDEBUG
              const auto & conv = get_conversion(std::string(elem_type.data()));

              for (auto d : IntRange<int>(0, conv.dim))
                libmesh_assert_equal_to(bex_elem_degrees[d], 2);
#endif
            }
            // ex_get_attributes did a values=calloc(); free() is our job.
            if (attr.values)
              free(attr.values);
        }
    }

  if (is_bezier)
    {
      // We'd better have the number of cvs we expect
      if( num_node_data_per_elem > num_nodes_per_elem )
        bex_num_elem_cvs = num_node_data_per_elem / 2;
      else
        bex_num_elem_cvs = num_nodes_per_elem;
      libmesh_assert_greater_equal(bex_num_elem_cvs, 0);

      // The old connect vector is currently a mix of the expected
      // connectivity and any Bezier extraction connectivity;
      // disentangle that, if necessary.
      bex_cv_conn.resize(num_elem_this_blk);
      if (num_node_data_per_elem > num_nodes_per_elem)
        {
          std::vector<int> old_connect(bex_num_elem_cvs * num_elem_this_blk);
          old_connect.swap(connect);
          auto src = old_connect.data();
          auto dst = connect.data();
          for (auto e : IntRange<std::size_t>(0, num_elem_this_blk))
            {
              std::copy(src, src + bex_num_elem_cvs, dst);
              src += bex_num_elem_cvs;
              dst += bex_num_elem_cvs;

              bex_cv_conn[e].resize(bex_num_elem_cvs);
              std::copy(src, src + bex_num_elem_cvs,
                        bex_cv_conn[e].begin());
              src += bex_num_elem_cvs;
            }
        }
    }

#endif // EX_API_VERS_NODOT >= 800
}



void ExodusII_IO_Helper::read_edge_blocks(MeshBase & mesh)
{
  LOG_SCOPE("read_edge_blocks()", "ExodusII_IO_Helper");

  // Check for quick return if there are no edge blocks.
  if (num_edge_blk == 0)
    return;

  // Build data structure that we can quickly search for edges
  // and then add required BoundaryInfo information. This is a
  // map from edge->key() to a list of (elem_id, edge_id) pairs
  // for the Edge in question. Since edge->key() is edge orientation
  // invariant, this map does not distinguish different orientations
  // of the same Edge. Since edge->key() is also not guaranteed to be
  // unique (though it is very unlikely for two distinct edges to have
  // the same key()), when we later look up an (elem_id, edge_id) pair
  // in the edge_map, we need to verify that the edge indeed matches
  // the searched edge by doing some further checks.
  typedef std::pair<dof_id_type, unsigned int> ElemEdgePair;
  std::unordered_map<dof_id_type, std::vector<ElemEdgePair>> edge_map;
  std::unique_ptr<Elem> edge_ptr;
  for (const auto & elem : mesh.element_ptr_range())
    for (auto e : elem->edge_index_range())
      {
        elem->build_edge_ptr(edge_ptr, e);
        dof_id_type edge_key = edge_ptr->key();

        // Creates vector if not already there
        auto & vec = edge_map[edge_key];
        vec.emplace_back(elem->id(), e);

        // If edge_ptr is a higher-order Elem (EDGE3 or higher) then also add
        // a map entry for the lower-order (EDGE2) element which has matching
        // vertices. This allows us to match lower-order edge blocks to edges
        // of higher-order 3D elems (e.g. HEX20, TET10) and simplifies the
        // definition of edge blocks.
        if (edge_ptr->default_order() != FIRST)
          {
            // Construct a temporary low-order edge so that we can compute its key()
            auto low_order_edge =
              Elem::build(Elem::first_order_equivalent_type(edge_ptr->type()));

            // Assign node pointers to low-order edge
            for (unsigned int v=0; v<edge_ptr->n_vertices(); ++v)
              low_order_edge->set_node(v, edge_ptr->node_ptr(v));

            // Compute the key for the temporary low-order edge we just built
            dof_id_type low_order_edge_key = low_order_edge->key();

            // Add this key to the map associated with the same (elem,
            // edge) pair as the higher-order edge
            auto & low_order_vec = edge_map[low_order_edge_key];
            low_order_vec.emplace_back(elem->id(), e);
          }
      }

  // Get reference to the mesh's BoundaryInfo object, as we will be
  // adding edges to this below.
  BoundaryInfo & bi = mesh.get_boundary_info();

  for (const auto & edge_block_id : edge_block_ids)
    {
      // exII::ex_get_block() output parameters.  Unlike the other
      // "extended" APIs, exII::ex_get_block() does not use a
      // parameter struct.
      int num_edge_this_blk = 0;
      int num_nodes_per_edge = 0;
      int num_edges_per_edge = 0;
      int num_faces_per_edge = 0;
      int num_attr_per_edge = 0;
      ex_err = exII::ex_get_block(ex_id,
                                  exII::EX_EDGE_BLOCK,
                                  edge_block_id,
                                  elem_type.data(),
                                  &num_edge_this_blk,
                                  &num_nodes_per_edge,
                                  &num_edges_per_edge, // 0 or -1 for edge blocks
                                  &num_faces_per_edge, // 0 or -1 for edge blocks
                                  &num_attr_per_edge);

      EX_CHECK_ERR(ex_err, "Error getting edge block info.");
      message("Info retrieved successfully for block: ", edge_block_id);

      // Read in the connectivity of the edges of this block,
      // watching out for the case where we actually have no
      // elements in this block (possible with parallel files)
      connect.resize(num_nodes_per_edge * num_edge_this_blk);

      if (!connect.empty())
        {
          ex_err = exII::ex_get_conn(ex_id,
                                     exII::EX_EDGE_BLOCK,
                                     edge_block_id,
                                     connect.data(), // node_conn
                                     nullptr,        // elem_edge_conn (unused)
                                     nullptr);       // elem_face_conn (unused)

          EX_CHECK_ERR(ex_err, "Error reading block connectivity.");
          message("Connectivity retrieved successfully for block: ", edge_block_id);

          // All edge types have an identity mapping from the corresponding
          // Exodus type, so we don't need to bother with mapping ids, but
          // we do need to know what kind of elements to build.
          const auto & conv = get_conversion(std::string(elem_type.data()));

          // Loop over indices in connectivity array, build edge elements,
          // look them up in the edge_map.
          for (auto [i, sz] = std::make_tuple(0u, connect.size()); i<sz; i+=num_nodes_per_edge)
            {
              auto edge = Elem::build(conv.libmesh_elem_type());
              for (int n=0; n<num_nodes_per_edge; ++n)
                {
                  auto exodus_node_id = this->connect[i+n];
                  dof_id_type libmesh_node_id = this->get_libmesh_node_id(exodus_node_id);
                  edge->set_node(n, mesh.node_ptr(libmesh_node_id));
                }

              // Compute key for the edge Elem we just built.
              dof_id_type edge_key = edge->key();

              // If this key is not found in the edge_map, which is
              // supposed to include every edge in the Mesh, then we
              // will throw an error now.
              auto & elem_edge_pair_vec =
                libmesh_map_find(edge_map, edge_key);

              for (const auto & elem_edge_pair : elem_edge_pair_vec)
                {
                  // We only want to match edges which have the same
                  // nodes (possibly with different orientation) to the one in the
                  // Exodus file, otherwise we ignore this elem_edge_pair.
                  //
                  // Note: this also handles the situation where two
                  // edges have the same key (hash collision) as then
                  // this check avoids a false positive.

                  // Build edge indicated by elem_edge_pair
                  mesh.elem_ptr(elem_edge_pair.first)->
                    build_edge_ptr(edge_ptr, elem_edge_pair.second);

                  // Determine whether this candidate edge is a "real" match,
                  // i.e. has the same nodes with a possibly different
                  // orientation. Note that here we only check that
                  // the vertices match regardless of how many nodes
                  // the edge has, which allows us to match a
                  // lower-order edge to a higher-order Elem.
                  bool is_match =
                    ((edge_ptr->node_id(0) == edge->node_id(0)) && (edge_ptr->node_id(1) == edge->node_id(1))) ||
                    ((edge_ptr->node_id(0) == edge->node_id(1)) && (edge_ptr->node_id(1) == edge->node_id(0)));

                  if (is_match)
                    {
                      // Add this (elem, edge, id) combo to the BoundaryInfo object.
                      bi.add_edge(elem_edge_pair.first,
                                  elem_edge_pair.second,
                                  edge_block_id);
                    }
                } // end loop over elem_edge_pairs
            } // end loop over connectivity array

          // Set edgeset name in the BoundaryInfo object.
          bi.edgeset_name(edge_block_id) = id_to_edge_block_names[edge_block_id];
        } // end if !connect.empty()
    } // end for edge_block_id : edge_block_ids
}



void ExodusII_IO_Helper::read_elem_num_map ()
{
  elem_num_map.resize(num_elem);

  // Note: we cannot use the exII::ex_get_num_map() here because it
  // (apparently) does not behave like ex_get_elem_num_map() when
  // there is no elem number map in the file: it throws an error
  // instead of returning a default identity array (1,2,3,...).
  ex_err = exII::ex_get_elem_num_map
    (ex_id, elem_num_map.empty() ? nullptr : elem_num_map.data());

  EX_CHECK_ERR(ex_err, "Error retrieving element number map.");
  message("Element numbering map retrieved successfully.");

  if (num_elem)
    {
      // The elem_num_map may contain ids larger than num_elem.  In
      // other words, the elem_num_map is not necessarily just a
      // permutation of the "trivial" 1,2,3,...  mapping, it can
      // contain effectively "any" numbers. Therefore, to get
      // "_end_elem_id", we need to check what the max entry in the
      // elem_num_map is.
      auto it = std::max_element(elem_num_map.begin(), elem_num_map.end());
      _end_elem_id = *it;
    }
  else
    _end_elem_id = 0;

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
      ex_err = exII::ex_get_ids(ex_id,
                                exII::EX_SIDE_SET,
                                ss_ids.data());
      EX_CHECK_ERR(ex_err, "Error retrieving sideset information.");
      message("All sideset information retrieved successfully.");

      // Resize appropriate data structures -- only do this once outside the loop
      num_sides_per_set.resize(num_side_sets);
      num_df_per_set.resize(num_side_sets);

      // Inquire about the length of the concatenated side sets element list
      num_elem_all_sidesets = inquire(*this, exII::EX_INQ_SS_ELEM_LEN, "Error retrieving length of the concatenated side sets element list!");

      elem_list.resize (num_elem_all_sidesets);
      side_list.resize (num_elem_all_sidesets);
      id_list.resize   (num_elem_all_sidesets);
    }

  char name_buffer[libmesh_max_str_length+1];
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
      ex_err = exII::ex_get_ids(ex_id,
                                exII::EX_NODE_SET,
                                nodeset_ids.data());
      EX_CHECK_ERR(ex_err, "Error retrieving nodeset information.");
      message("All nodeset information retrieved successfully.");

      // Resize appropriate data structures -- only do this once outside the loop
      num_nodes_per_set.resize(num_node_sets);
      num_node_df_per_set.resize(num_node_sets);
    }

  char name_buffer[libmesh_max_str_length+1];
  for (int i=0; i<num_node_sets; ++i)
    {
      ex_err = exII::ex_get_name(ex_id, exII::EX_NODE_SET,
                                 nodeset_ids[i], name_buffer);
      EX_CHECK_ERR(ex_err, "Error getting node set name.");
      id_to_ns_names[nodeset_ids[i]] = name_buffer;
    }
  message("All node set names retrieved successfully.");
}



void ExodusII_IO_Helper::read_elemset_info()
{
  elemset_ids.resize(num_elem_sets);
  if (num_elem_sets > 0)
    {
      ex_err = exII::ex_get_ids(ex_id,
                                exII::EX_ELEM_SET,
                                elemset_ids.data());
      EX_CHECK_ERR(ex_err, "Error retrieving elemset information.");
      message("All elemset information retrieved successfully.");

      // Resize appropriate data structures -- only do this once outside the loop
      num_elems_per_set.resize(num_elem_sets);
      num_elem_df_per_set.resize(num_elem_sets);

      // Inquire about the length of the concatenated elemset list
      num_elem_all_elemsets =
        inquire(*this, exII::EX_INQ_ELS_LEN,
                "Error retrieving length of the concatenated elem sets element list!");

      elemset_list.resize(num_elem_all_elemsets);
      elemset_id_list.resize(num_elem_all_elemsets);

      // Debugging
      // libMesh::out << "num_elem_all_elemsets = " << num_elem_all_elemsets << std::endl;
    }

  char name_buffer[libmesh_max_str_length+1];
  for (int i=0; i<num_elem_sets; ++i)
    {
      ex_err = exII::ex_get_name(ex_id, exII::EX_ELEM_SET,
                                 elemset_ids[i], name_buffer);
      EX_CHECK_ERR(ex_err, "Error getting node set name.");
      id_to_elemset_names[elemset_ids[i]] = name_buffer;
    }
  message("All elem set names retrieved successfully.");
}



void ExodusII_IO_Helper::read_sideset(int id, int offset)
{
  LOG_SCOPE("read_sideset()", "ExodusII_IO_Helper");

  libmesh_assert_less (id, ss_ids.size());
  libmesh_assert_less (id, num_sides_per_set.size());
  libmesh_assert_less (id, num_df_per_set.size());
  libmesh_assert_less_equal (offset, elem_list.size());
  libmesh_assert_less_equal (offset, side_list.size());

  ex_err = exII::ex_get_set_param(ex_id,
                                  exII::EX_SIDE_SET,
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


  // Don't call ex_get_set unless there are actually sides there to get.
  // Exodus prints an annoying warning in DEBUG mode otherwise...
  if (num_sides_per_set[id] > 0)
    {
      ex_err = exII::ex_get_set(ex_id,
                                exII::EX_SIDE_SET,
                                ss_ids[id],
                                &elem_list[offset],
                                &side_list[offset]);
      EX_CHECK_ERR(ex_err, "Error retrieving sideset data.");
      message("Data retrieved successfully for sideset: ", id);

      for (int i=0; i<num_sides_per_set[id]; i++)
        id_list[i+offset] = ss_ids[id];
    }
}



void ExodusII_IO_Helper::read_elemset(int id, int offset)
{
  LOG_SCOPE("read_elemset()", "ExodusII_IO_Helper");

  libmesh_assert_less (id, elemset_ids.size());
  libmesh_assert_less (id, num_elems_per_set.size());
  libmesh_assert_less (id, num_elem_df_per_set.size());
  libmesh_assert_less_equal (offset, elemset_list.size());

  ex_err = exII::ex_get_set_param(ex_id,
                                  exII::EX_ELEM_SET,
                                  elemset_ids[id],
                                  &num_elems_per_set[id],
                                  &num_elem_df_per_set[id]);
  EX_CHECK_ERR(ex_err, "Error retrieving elemset parameters.");
  message("Parameters retrieved successfully for elemset: ", id);


  // It's OK for offset==elemset_list.size() as long as num_elems_per_set[id]==0
  // because in that case we don't actually read anything...
  #ifdef DEBUG
  if (static_cast<unsigned int>(offset) == elemset_list.size())
    libmesh_assert_equal_to (num_elems_per_set[id], 0);
  #endif

  // Don't call ex_get_set() unless there are actually elems there to get.
  // Exodus prints an annoying warning in DEBUG mode otherwise...
  if (num_elems_per_set[id] > 0)
    {
      ex_err = exII::ex_get_set(ex_id,
                                exII::EX_ELEM_SET,
                                elemset_ids[id],
                                &elemset_list[offset],
                                /*set_extra_list=*/nullptr);
      EX_CHECK_ERR(ex_err, "Error retrieving elemset data.");
      message("Data retrieved successfully for elemset: ", id);

      // Create vector containing elemset ids for each element in the set
      for (int i=0; i<num_elems_per_set[id]; i++)
        elemset_id_list[i+offset] = elemset_ids[id];
    }
}



void ExodusII_IO_Helper::read_all_nodesets()
{
  LOG_SCOPE("read_all_nodesets()", "ExodusII_IO_Helper");

  // Figure out how many nodesets there are in the file so we can
  // properly resize storage as necessary.
  num_node_sets =
    inquire
    (*this, exII::EX_INQ_NODE_SETS,
     "Error retrieving number of node sets");

  // Figure out how many nodes there are in all the nodesets.
  int total_nodes_in_all_sets =
    inquire
    (*this, exII::EX_INQ_NS_NODE_LEN,
     "Error retrieving number of nodes in all node sets.");

  // Figure out how many distribution factors there are in all the nodesets.
  int total_df_in_all_sets =
    inquire
    (*this, exII::EX_INQ_NS_DF_LEN,
     "Error retrieving number of distribution factors in all node sets.");

  // If there are no nodesets, there's nothing to read in.
  if (num_node_sets == 0)
    return;

  // Allocate space to read all the nodeset data.
  // Use existing class members where possible to avoid shadowing
  nodeset_ids.clear();          nodeset_ids.resize(num_node_sets);
  num_nodes_per_set.clear();    num_nodes_per_set.resize(num_node_sets);
  num_node_df_per_set.clear();  num_node_df_per_set.resize(num_node_sets);
  node_sets_node_index.clear(); node_sets_node_index.resize(num_node_sets);
  node_sets_dist_index.clear(); node_sets_dist_index.resize(num_node_sets);
  node_sets_node_list.clear();  node_sets_node_list.resize(total_nodes_in_all_sets);
  node_sets_dist_fact.clear();  node_sets_dist_fact.resize(total_df_in_all_sets);

  // Handle single-precision files
  MappedInputVector mapped_node_sets_dist_fact(node_sets_dist_fact, _single_precision);

  // Build exII::ex_set_spec struct
  exII::ex_set_specs set_specs = {};
  set_specs.sets_ids            = nodeset_ids.data();
  set_specs.num_entries_per_set = num_nodes_per_set.data();
  set_specs.num_dist_per_set    = num_node_df_per_set.data();
  set_specs.sets_entry_index    = node_sets_node_index.data();
  set_specs.sets_dist_index     = node_sets_dist_index.data();
  set_specs.sets_entry_list     = node_sets_node_list.data();
  set_specs.sets_extra_list     = nullptr;
  set_specs.sets_dist_fact      = total_df_in_all_sets ? mapped_node_sets_dist_fact.data() : nullptr;

  ex_err = exII::ex_get_concat_sets(ex_id, exII::EX_NODE_SET, &set_specs);
  EX_CHECK_ERR(ex_err, "Error reading concatenated nodesets");

  // Read the nodeset names from file!
  char name_buffer[libmesh_max_str_length+1];
  for (int i=0; i<num_node_sets; ++i)
    {
      ex_err = exII::ex_get_name
        (ex_id,
         exII::EX_NODE_SET,
         nodeset_ids[i],
         name_buffer);
      EX_CHECK_ERR(ex_err, "Error getting node set name.");
      id_to_ns_names[nodeset_ids[i]] = name_buffer;
    }
}



void ExodusII_IO_Helper::close() noexcept
{
  // Call ex_close on every processor that did ex_open or ex_create;
  // newer Exodus versions error if we try to reopen a file that
  // hasn't been officially closed.  Don't close the file if we didn't
  // open it; this also raises an Exodus error.

  // We currently do read-only ex_open on every proc (to do read
  // operations on every proc), but we do ex_open and ex_create for
  // writes on every proc only with Nemesis files.
  if (!(_opened_by_create || opened_for_writing) ||
      (this->processor_id() == 0) ||
      (!_run_only_on_proc0))
    {
      if (opened_for_writing || opened_for_reading)
        {
          ex_err = exII::ex_close(ex_id);
          // close() is called from the destructor, so it may be called e.g.
          // during stack unwinding while processing an exception. In that case
          // we don't want to throw another exception or immediately terminate
          // the code, since that would prevent any possible recovery from the
          // exception in question. So we just log the error closing the file
          // and continue.
          if (ex_err < 0)
            message("Error closing Exodus file.");
          else
            message("Exodus file closed successfully.");
        }
    }

  // Now that the file is closed, it's no longer opened for
  // reading or writing.
  opened_for_writing = false;
  opened_for_reading = false;
  _opened_by_create = false;
}



void ExodusII_IO_Helper::read_time_steps()
{
  // Make sure we have an up-to-date count of the number of time steps in the file.
  this->read_num_time_steps();

  if (num_time_steps > 0)
    {
      time_steps.resize(num_time_steps);
      ex_err = exII::ex_get_all_times
        (ex_id,
         MappedInputVector(time_steps, _single_precision).data());
      EX_CHECK_ERR(ex_err, "Error reading timesteps!");
    }
}



void ExodusII_IO_Helper::read_num_time_steps()
{
  num_time_steps =
    inquire(*this, exII::EX_INQ_TIME, "Error retrieving number of time steps");
}



void ExodusII_IO_Helper::read_nodal_var_values(std::string nodal_var_name, int time_step)
{
  LOG_SCOPE("read_nodal_var_values()", "ExodusII_IO_Helper");

  // Read the nodal variable names from file, so we can see if we have the one we're looking for
  this->read_var_names(NODAL);

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
      libMesh::err << "Available variables: " << std::endl;
      for (const auto & var_name : nodal_var_names)
        libMesh::err << var_name << std::endl;

      libmesh_error_msg("Unable to locate variable named: " << nodal_var_name);
    }

  // Clear out any previously read nodal variable values
  this->nodal_var_values.clear();

  std::vector<Real> unmapped_nodal_var_values(num_nodes);

  // Call the Exodus API to read the nodal variable values
  ex_err = exII::ex_get_var
    (ex_id,
     time_step,
     exII::EX_NODAL,
     var_index+1,
     1, // exII::ex_entity_id, not sure exactly what this is but in the ex_get_nodal_var.c shim, they pass 1
     num_nodes,
     MappedInputVector(unmapped_nodal_var_values, _single_precision).data());
  EX_CHECK_ERR(ex_err, "Error reading nodal variable values!");

  for (auto i : make_range(num_nodes))
    {
      // Determine the libmesh node id implied by "i". The
      // get_libmesh_node_id() helper function expects a 1-based
      // Exodus node id, so we construct the "implied" Exodus node id
      // from "i" by adding 1.
      //
      // If the user has set the "set_unique_ids_from_maps" flag to
      // true, then calling get_libmesh_node_id(i+1) will just return
      // i, otherwise it will determine the value (with error
      // checking) using this->node_num_map.
      auto libmesh_node_id = this->get_libmesh_node_id(/*exodus_node_id=*/i+1);

      // Store the nodal value in the map.
      this->nodal_var_values[libmesh_node_id] = unmapped_nodal_var_values[i];
    }
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
    case SIDESET:
      this->read_var_names_impl("s", num_sideset_vars, sideset_var_names);
      break;
    case NODESET:
      this->read_var_names_impl("m", num_nodeset_vars, nodeset_var_names);
      break;
    case ELEMSET:
      this->read_var_names_impl("t", num_elemset_vars, elemset_var_names);
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
  NamesData names_table(count, libmesh_max_str_length);

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




void
ExodusII_IO_Helper::write_var_names(ExodusVarType type,
                                    const std::vector<std::string> & names)
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
    case SIDESET:
      {
        // Note: calling this function *sets* num_sideset_vars to the
        // number of entries in the 'names' vector, num_sideset_vars
        // does not already need to be set before calling this.
        this->write_var_names_impl("s", num_sideset_vars, names);
        break;
      }
    case NODESET:
      {
        this->write_var_names_impl("m", num_nodeset_vars, names);
        break;
      }
    case ELEMSET:
      {
        this->write_var_names_impl("t", num_elemset_vars, names);
        break;
      }
    default:
      libmesh_error_msg("Unrecognized ExodusVarType " << type);
    }
}



void
ExodusII_IO_Helper::write_var_names_impl(const char * var_type,
                                         int & count,
                                         const std::vector<std::string> & names)
{
  // Update the count variable so that it's available to other parts of the class.
  count = cast_int<int>(names.size());

  // Write that number of variables to the file.
  ex_err = exII::ex_put_var_param(ex_id, var_type, count);
  EX_CHECK_ERR(ex_err, "Error setting number of vars.");

  // Nemesis doesn't like trying to write nodal variable names in
  // files with no nodes.
  if (!this->num_nodes)
    return;

  if (count > 0)
    {
      NamesData names_table(count, _max_name_length);

      // Store the input names in the format required by Exodus.
      for (int i=0; i != count; ++i)
        {
          if(names[i].length() > _max_name_length)
            libmesh_warning(
              "*** Warning, Exodus variable name \"" <<
              names[i] << "\" too long (current max " <<
              _max_name_length << "/" << libmesh_max_str_length <<
              " characters). Name will be truncated. ");
          names_table.push_back_entry(names[i]);
        }

      if (verbose)
        {
          libMesh::out << "Writing variable name(s) to file: " << std::endl;
          for (int i=0; i != count; ++i)
            libMesh::out << names_table.get_char_star(i) << std::endl;
        }

      ex_err = exII::ex_put_var_names(ex_id,
                                      var_type,
                                      count,
                                      names_table.get_char_star_star()
                                      );

      EX_CHECK_ERR(ex_err, "Error writing variable names.");
    }
}



void ExodusII_IO_Helper::read_elemental_var_values(std::string elemental_var_name,
                                                   int time_step,
                                                   std::map<dof_id_type, Real> & elem_var_value_map)
{
  LOG_SCOPE("read_elemental_var_values()", "ExodusII_IO_Helper");

  this->read_var_names(ELEMENTAL);

  // See if we can find the variable we are looking for
  unsigned int var_index = 0;
  bool found = false;

  // Do a linear search for elem_var_name in elemental_var_names
  for (; var_index != elem_var_names.size(); ++var_index)
    if (elem_var_names[var_index] == elemental_var_name)
      {
        found = true;
        break;
      }

  if (!found)
    {
      libMesh::err << "Available variables: " << std::endl;
      for (const auto & var_name : elem_var_names)
        libMesh::err << var_name << std::endl;

      libmesh_error_msg("Unable to locate variable named: " << elemental_var_name);
    }

  // Sequential index which we can use to look up the element ID in the elem_num_map.
  unsigned ex_el_num = 0;

  // Element variable truth table
  std::vector<int> var_table(block_ids.size() * elem_var_names.size());
  exII::ex_get_truth_table(ex_id, exII::EX_ELEM_BLOCK, block_ids.size(), elem_var_names.size(), var_table.data());

  for (unsigned i=0; i<static_cast<unsigned>(num_elem_blk); i++)
    {
      ex_err = exII::ex_get_block(ex_id,
                                  exII::EX_ELEM_BLOCK,
                                  block_ids[i],
                                  /*elem_type=*/nullptr,
                                  &num_elem_this_blk,
                                  /*num_nodes_per_entry=*/nullptr,
                                  /*num_edges_per_entry=*/nullptr,
                                  /*num_faces_per_entry=*/nullptr,
                                  /*num_attr=*/nullptr);
      EX_CHECK_ERR(ex_err, "Error getting number of elements in block.");

      // If the current variable isn't active on this subdomain, advance
      // the index by the number of elements on this block and go to the
      // next loop iteration.
      if (!var_table[elem_var_names.size()*i + var_index])
        {
          ex_el_num += num_elem_this_blk;
          continue;
        }

      std::vector<Real> block_elem_var_values(num_elem_this_blk);

      ex_err = exII::ex_get_var
        (ex_id,
         time_step,
         exII::EX_ELEM_BLOCK,
         var_index+1,
         block_ids[i],
         num_elem_this_blk,
         MappedInputVector(block_elem_var_values, _single_precision).data());
      EX_CHECK_ERR(ex_err, "Error getting elemental values.");

      for (unsigned j=0; j<static_cast<unsigned>(num_elem_this_blk); j++)
        {
          // Determine the libmesh id of the element with zero-based
          // index "ex_el_num".  This function expects a one-based
          // index, so we add 1 to ex_el_num when we pass it in.
          auto libmesh_elem_id =
            this->get_libmesh_elem_id(ex_el_num + 1);

          // Store the elemental value in the map.
          elem_var_value_map[libmesh_elem_id] = block_elem_var_values[j];

          // Go to the next sequential element ID.
          ex_el_num++;
        }
    }
}



std::vector<std::map<dof_id_type, dof_id_type>>
ExodusII_IO_Helper::read_extra_integers
  (const std::vector<std::string> & extra_integer_var_names)
{
  std::vector<std::map<dof_id_type, dof_id_type>>
    extra_integer_values(extra_integer_var_names.size());

  if (extra_integer_var_names.empty())
    return extra_integer_values;

  // Make sure we have an up-to-date count of the number of time steps in the file.
  this->read_num_time_steps();

  // Prepare to check if each real number is outside of the
  // range we can convert exactly

  const int exodus_digits = _single_precision ?
                            std::numeric_limits<float>::digits :
                            std::numeric_limits<double>::digits;

  const int shift = std::min(std::numeric_limits<Real>::digits,
                             exodus_digits);

  const long long max_representation = 1LL << shift;

  for (auto i : index_range(extra_integer_var_names))
    {
      std::map<dof_id_type, Real> raw_vals;

      // Read element extra "integers" as doubles from the last time step
      this->read_elemental_var_values(extra_integer_var_names[i], this->num_time_steps, raw_vals);

      // Convert doubles to actual integers, at least within the
      // largest convex subset of doubles where this is a bijection.
      auto & values = extra_integer_values[i];

      for (auto [elem_id, extra_val] : raw_vals)
        {
          if (extra_val == Real(-1))
            {
              values[elem_id] = DofObject::invalid_id;
              continue;
            }

          // Ignore FE_INVALID here even if we've enabled FPEs; a
          // thrown exception is preferred over an FPE signal.
          FPEDisabler disable_fpes;
          const long long int_val = std::llround(extra_val);

          libmesh_error_msg_if(int_val > max_representation,
                               "Error! An element integer value higher than "
                               << max_representation
                               << " was found! Exodus uses real numbers for storing element "
                               " integers, which can only represent integers from 0 to "
                               << max_representation << ".");

          libmesh_error_msg_if(int_val < 0,
                               "Error! An element integer value less than -1"
                               << " was found! Exodus uses real numbers for storing element "
                               " integers, which can only represent integers from 0 to "
                               << max_representation << ".");

          values[elem_id] = cast_int<dof_id_type>(int_val);
        }
    }

  return extra_integer_values;
}



dof_id_type ExodusII_IO_Helper::get_libmesh_node_id(int exodus_node_id)
{
  return this->get_libmesh_id(exodus_node_id, this->node_num_map);
}

dof_id_type ExodusII_IO_Helper::get_libmesh_elem_id(int exodus_elem_id)
{
  return this->get_libmesh_id(exodus_elem_id, this->elem_num_map);
}

dof_id_type
ExodusII_IO_Helper::get_libmesh_id(int exodus_id,
                                   const std::vector<int> & num_map)
{
  // The input exodus_id is assumed to be a (1-based) index into
  // the {node,elem}_num_map, so in order to use exodus_id as an index
  // in C++, we need to first make it zero-based.
  auto exodus_id_zero_based =
    cast_int<dof_id_type>(exodus_id - 1);

  // Throw an informative error message rather than accessing past the
  // end of the provided num_map. If we are setting Elem unique_ids
  // based on the num_map, we don't need to do this check.
  if (!this->set_unique_ids_from_maps)
    libmesh_error_msg_if(exodus_id_zero_based >= num_map.size(),
                         "Cannot get LibMesh id for Exodus id: " << exodus_id);

  // If the user set the flag which stores Exodus node/elem ids as
  // unique_ids instead of regular ids, then the libmesh id we are
  // looking for is actually just "exodus_id_zero_based". Otherwise,
  // we need to look up the Node/Elem's id in the provided num_map,
  // *and* then subtract 1 from that because the entries in the
  // num_map are also 1-based.
  dof_id_type libmesh_id =
    this->set_unique_ids_from_maps ?
    cast_int<dof_id_type>(exodus_id_zero_based) :
    cast_int<dof_id_type>(num_map[exodus_id_zero_based] - 1);

  return libmesh_id;
}



void
ExodusII_IO_Helper::
conditionally_set_node_unique_id(MeshBase & mesh, Node * node, int zero_based_node_num_map_index)
{
  this->set_dof_object_unique_id(mesh, node, libmesh_vector_at(this->node_num_map, zero_based_node_num_map_index));
}

void
ExodusII_IO_Helper::
conditionally_set_elem_unique_id(MeshBase & mesh, Elem * elem, int zero_based_elem_num_map_index)
{
  this->set_dof_object_unique_id(mesh, elem, libmesh_vector_at(this->elem_num_map, zero_based_elem_num_map_index));
}

void
ExodusII_IO_Helper::set_dof_object_unique_id(
  MeshBase & mesh,
  DofObject * dof_object,
  int exodus_mapped_id)
{
  if (this->set_unique_ids_from_maps)
  {
    // Exodus ids are always 1-based while libmesh ids are always
    // 0-based, so to make a libmesh unique_id here, we subtract 1
    // from the exodus_mapped_id to make it 0-based.
    auto exodus_mapped_id_zero_based =
      cast_int<dof_id_type>(exodus_mapped_id - 1);

    // Set added_node's unique_id to "exodus_mapped_id_zero_based".
    dof_object->set_unique_id(cast_int<unique_id_type>(exodus_mapped_id_zero_based));

    // Normally the Mesh is responsible for setting the unique_ids
    // of Nodes/Elems in a consistent manner, so when we set the unique_id
    // of a Node/Elem manually based on the {node,elem}_num_map, we need to
    // make sure that the "next" unique id assigned by the Mesh
    // will still be valid. We do this by making sure that the
    // next_unique_id is greater than the one we set manually. The
    // APIs for doing this are only defined when unique ids are
    // enabled.
#ifdef LIBMESH_ENABLE_UNIQUE_ID
    unique_id_type next_unique_id = mesh.next_unique_id();
    mesh.set_next_unique_id(std::max(next_unique_id, static_cast<unique_id_type>(exodus_mapped_id_zero_based + 1)));
#else
    // Avoid compiler warnings about the unused variable
    libmesh_ignore(mesh);
#endif
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

      if (_single_precision)
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

      // By default we just open the Exodus file in "EX_CLOBBER" mode,
      // which, according to "ncdump -k", writes the file in "64-bit
      // offset" mode, which is a NETCDF3 file format.
      int mode = EX_CLOBBER;

      // If HDF5 is available, by default we will write Exodus files
      // in a more modern NETCDF4-compatible format. For this file
      // type, "ncdump -k" will report "netCDF-4".
#ifdef LIBMESH_HAVE_HDF5
      if (this->_write_hdf5)
        {
          mode |= EX_NETCDF4;
          mode |= EX_NOCLASSIC;
        }
#endif

      {
        FPEDisabler disable_fpes;
        ex_id = exII::ex_create(filename.c_str(), mode, &comp_ws, &io_ws);
      }

      EX_CHECK_ERR(ex_id, "Error creating ExodusII/Nemesis mesh file.");

      // We don't have access to the names we might be writing until we
      // write them, so we can't set a guaranteed max name length here.
      // But it looks like the most ExodusII can support is 80, so we'll
      // just waste 48 bytes here and there.
      ex_err = exII::ex_set_max_name_length(ex_id, _max_name_length);
      EX_CHECK_ERR(ex_err, "Error setting max ExodusII name length.");

      if (verbose)
        libMesh::out << "File created successfully." << std::endl;
    }

  opened_for_writing = true;
  _opened_by_create = true;
  current_filename = filename;
}



void ExodusII_IO_Helper::initialize(std::string str_title, const MeshBase & mesh, bool use_discontinuous)
{
  // The majority of this function only executes on processor 0, so any functions
  // which are collective, like n_active_elem() or n_edge_conds() must be called
  // before the processors' execution paths diverge.
  libmesh_parallel_only(mesh.comm());

  unsigned int n_active_elem = mesh.n_active_elem();
  const BoundaryInfo & bi = mesh.get_boundary_info();
  num_edge = bi.n_edge_conds();

  // We need to know about all processors' subdomains
  build_subdomain_map(mesh, false);

  num_elem = n_active_elem;
  num_nodes = 0;

  // If we're adding face elements they'll need copies of their nodes,
  // and we'll need to manage the extra copies.
  this->calculate_added_side_node_offsets(mesh);

  // If _write_as_dimension is nonzero, use it to set num_dim in the Exodus file.
  if (_write_as_dimension)
    num_dim = _write_as_dimension;
  else if (_use_mesh_dimension_instead_of_spatial_dimension)
    num_dim = mesh.mesh_dimension();
  else
    num_dim = mesh.spatial_dimension();

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (!use_discontinuous)
    {
      // Don't rely on mesh.n_nodes() here.  If ReplicatedMesh nodes
      // have been deleted without renumbering after, it will be
      // incorrect.
      num_nodes += cast_int<int>(std::distance(mesh.nodes_begin(),
                                               mesh.nodes_end()));
    }
  else
    {
      for (const auto & elem : mesh.active_element_ptr_range())
        num_nodes += elem->n_nodes();
    }

  std::set<boundary_id_type> unique_side_boundaries;
  std::vector<boundary_id_type> unique_node_boundaries;

  // Build set of unique sideset (+shellface) ids
  {
    // Start with "side" boundaries (i.e. of 3D elements)
    std::vector<boundary_id_type> side_boundaries;
    bi.build_side_boundary_ids(side_boundaries);
    unique_side_boundaries.insert(side_boundaries.begin(), side_boundaries.end());

    // Add shell face boundaries to the list of side boundaries, since ExodusII
    // treats these the same way.
    std::vector<boundary_id_type> shellface_boundaries;
    bi.build_shellface_boundary_ids(shellface_boundaries);
    unique_side_boundaries.insert(shellface_boundaries.begin(), shellface_boundaries.end());

    // Add any empty-but-named side boundary ids
    for (const auto & pr : bi.get_sideset_name_map())
      unique_side_boundaries.insert(pr.first);
  }

  // Build set of unique nodeset ids
  bi.build_node_boundary_ids(unique_node_boundaries);
  for (const auto & pair : bi.get_nodeset_name_map())
    {
      const boundary_id_type id = pair.first;

      if (std::find(unique_node_boundaries.begin(),
                    unique_node_boundaries.end(), id)
            == unique_node_boundaries.end())
        unique_node_boundaries.push_back(id);
    }

  num_side_sets = cast_int<int>(unique_side_boundaries.size());
  num_node_sets = cast_int<int>(unique_node_boundaries.size());

  num_elem_blk = cast_int<int>(this->_subdomain_map.size());

  if (str_title.size() > MAX_LINE_LENGTH)
    {
      libMesh::err << "Warning, Exodus files cannot have titles longer than "
                   << MAX_LINE_LENGTH
                   << " characters.  Your title will be truncated."
                   << std::endl;
      str_title.resize(MAX_LINE_LENGTH);
    }

  // Edge BCs are handled a bit differently than sidesets and nodesets.
  // They are written as separate "edge blocks", and then edge variables
  // can be defined on those blocks. That is, they are not written as
  // edge sets, since edge sets must refer to edges stored elsewhere.
  // We write a separate edge block for each unique boundary id that
  // we have.
  num_edge_blk = bi.get_edge_boundary_ids().size();

  // Check whether the Mesh Elems have an extra_integer called "elemset_code".
  // If so, this means that the mesh defines elemsets via the
  // extra_integers capability of Elems.
  if (mesh.has_elem_integer("elemset_code"))
    {
      // unsigned int elemset_index =
      //   mesh.get_elem_integer_index("elemset_code");

      // Debugging
      // libMesh::out << "Mesh defines an elemset_code at index " << elemset_index << std::endl;

      // Store the number of elemsets in the exo file header.
      num_elem_sets = mesh.n_elemsets();
    }

  // Build an ex_init_params() structure that is to be passed to the
  // newer ex_put_init_ext() API. The new API will eventually allow us
  // to store edge and face data in the Exodus file.
  //
  // Notes:
  // * We use C++11 zero initialization syntax to make sure that all
  //   members of the struct (including ones we aren't using) are
  //   given sensible values.
  // * For the "title" field, we manually do a null-terminated string
  //   copy since std::string does not null-terminate but it does
  //   return the number of characters successfully copied.
  exII::ex_init_params params = {};
  params.title[str_title.copy(params.title, MAX_LINE_LENGTH)] = '\0';
  params.num_dim = num_dim;
  params.num_nodes = num_nodes;
  params.num_elem = num_elem;
  params.num_elem_blk = num_elem_blk;
  params.num_node_sets = num_node_sets;
  params.num_side_sets = num_side_sets;
  params.num_elem_sets = num_elem_sets;
  params.num_edge_blk = num_edge_blk;
  params.num_edge = num_edge;

  ex_err = exII::ex_put_init_ext(ex_id, &params);
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

  auto push_node = [this](const Point & p) {
    x.push_back(p(0) + _coordinate_offset(0));

#if LIBMESH_DIM > 1
    y.push_back(p(1) + _coordinate_offset(1));
#else
    y.push_back(0.);
#endif
#if LIBMESH_DIM > 2
    z.push_back(p(2) + _coordinate_offset(2));
#else
    z.push_back(0.);
#endif
  };

  // And in the node_num_map. If the user has set the
  // _set_unique_ids_from_maps flag, then we will write the Node
  // unique_ids to the node_num_map, otherwise we will just write a
  // trivial node_num_map, since in that we don't write the unique_id
  // information to the Exodus file. In other words, set the
  // _set_unique_ids_from_maps flag to true on both the reading and
  // writing ExodusII_IO objects if you want to preserve the
  // node_num_map information without actually renumbering the Nodes
  // in libmesh according to the node_num_map.
  //
  // One reason why you might not want to actually renumber the Nodes
  // in libmesh according to the node_num_map is that it can introduce
  // undesirable large "gaps" in the numbering, e.g. Nodes numbered
  // [0, 1, 1000, 10001] which is not ideal for the ReplicatedMesh
  // _nodes data structure, which stores the Nodes in a contiguous
  // array based on Node id.

  // Let's skip the node_num_map in the discontinuous and add_sides
  // cases, since we're effectively duplicating nodes for the sake of
  // discontinuous visualization, so it isn't clear how to deal with
  // node_num_map here. This means that writing meshes in such a way
  // won't work with element numberings that have id "holes".

  if (!use_discontinuous && !_add_sides)
    node_num_map.reserve(num_nodes);

  // Clear out any previously-mapped node IDs.
  libmesh_node_num_to_exodus.clear();

  if (!use_discontinuous)
    {
      for (const auto & node_ptr : mesh.node_ptr_range())
        {
          const Node & node = *node_ptr;

          push_node(node);

          // Fill in node_num_map entry with the proper (1-based) node
          // id, unless we're not going to be able to keep the map up
          // later. If the user has chosen to _set_unique_ids_from_maps,
          // then we fill up the node_num_map with (1-based) unique
          // ids rather than node ids.
          if (!_add_sides)
            {
              if (this->set_unique_ids_from_maps)
                node_num_map.push_back(node.unique_id() + 1);
              else
                node_num_map.push_back(node.id() + 1);
            }

          // Also map the zero-based libmesh node id to the (1-based)
          // index in the node_num_map it corresponds to
          // (this is equivalent to the current size of the "x" vector,
          // so we just use x.size()). This map is used to look up
          // an Exodus Node id given a libMesh Node id, so it does
          // involve unique_ids.
          libmesh_node_num_to_exodus[ cast_int<int>(node.id()) ] = cast_int<int>(x.size());
        } // end for (node_ptr)
    }
  else // use_discontinuous
    {
      for (const auto & elem : mesh.active_element_ptr_range())
        for (const Node & node : elem->node_ref_range())
          {
            push_node(node);

            // Let's skip the node_num_map in the discontinuous
            // case, since we're effectively duplicating nodes for
            // the sake of discontinuous visualization, so it isn't
            // clear how to deal with node_num_map here. This means
            // that writing discontinuous meshes won't work with
            // element numberings that have "holes".
          }
    }

  if (_add_sides)
    {
      // To match the numbering of parallel-generated nodal solutions
      // on fake side nodes, we need to loop through elements from
      // earlier ranks first.
      std::vector<std::vector<const Elem *>>
        elems_by_pid(mesh.n_processors());

      for (const auto & elem : mesh.active_element_ptr_range())
        elems_by_pid[elem->processor_id()].push_back(elem);

      for (auto p : index_range(elems_by_pid))
        for (const Elem * elem : elems_by_pid[p])
          for (auto s : elem->side_index_range())
            {
              if (EquationSystems::redundant_added_side(*elem,s))
                continue;

              const std::vector<unsigned int> side_nodes =
                elem->nodes_on_side(s);

              for (auto n : side_nodes)
                push_node(elem->point(n));
            }

      // Node num maps just don't make sense if we're adding a bunch
      // of visualization nodes that are independent copies of the
      // same libMesh node.
      node_num_map.clear();
    }

  ex_err = exII::ex_put_coord
    (ex_id,
     x.empty() ? nullptr : MappedOutputVector(x, _single_precision).data(),
     y.empty() ? nullptr : MappedOutputVector(y, _single_precision).data(),
     z.empty() ? nullptr : MappedOutputVector(z, _single_precision).data());

  EX_CHECK_ERR(ex_err, "Error writing coordinates to Exodus file.");

  if (!use_discontinuous && !_add_sides)
    {
      // Also write the (1-based) node_num_map to the file.
      ex_err = exII::ex_put_node_num_map(ex_id, node_num_map.data());
      EX_CHECK_ERR(ex_err, "Error writing node_num_map");
    }
}



void ExodusII_IO_Helper::write_elements(const MeshBase & mesh, bool use_discontinuous)
{
  LOG_SCOPE("write_elements()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // element map vector
  num_elem_blk = cast_int<int>(this->_subdomain_map.size());
  block_ids.resize(num_elem_blk);

  std::vector<int> elem_blk_id;
  std::vector<int> num_elem_this_blk_vec;
  std::vector<int> num_nodes_per_elem_vec;
  std::vector<int> num_edges_per_elem_vec;
  std::vector<int> num_faces_per_elem_vec;
  std::vector<int> num_attr_vec;
  NamesData elem_type_table(num_elem_blk, _max_name_length);

  // Note: It appears that there is a bug in exodusII::ex_put_name where
  // the index returned from the ex_id_lkup is erroneously used.  For now
  // the work around is to use the alternative function ex_put_names, but
  // this function requires a char ** data structure.
  NamesData names_table(num_elem_blk, _max_name_length);

  num_elem = 0;

  // counter indexes into the block_ids vector
  unsigned int counter = 0;
  for (auto & [subdomain_id, element_id_vec] : this->_subdomain_map)
    {
      block_ids[counter] = subdomain_id;

      const ElemType elem_t = (subdomain_id >= this->_subdomain_id_end) ?
        ElemType(subdomain_id - this->_subdomain_id_end) :
        mesh.elem_ref(element_id_vec[0]).type();

      if (subdomain_id >= this->_subdomain_id_end)
        {
          libmesh_assert(_add_sides);
          libmesh_assert(element_id_vec.size() == 1);
          num_elem_this_blk_vec.push_back
            (cast_int<int>(element_id_vec[0]));
          names_table.push_back_entry
            (Utility::enum_to_string<ElemType>(elem_t));
        }
      else
        {
          libmesh_assert(!element_id_vec.empty());
          num_elem_this_blk_vec.push_back
            (cast_int<int>(element_id_vec.size()));
          names_table.push_back_entry
            (mesh.subdomain_name(subdomain_id));
        }

      num_elem += num_elem_this_blk_vec.back();

      // Use the first element in this block to get representative information.
      // Note that Exodus assumes all elements in a block are of the same type!
      // We are using that same assumption here!
      const auto & conv = get_conversion(elem_t);
      num_nodes_per_elem = Elem::type_to_n_nodes_map[elem_t];
      if (Elem::type_to_n_nodes_map[elem_t] == invalid_uint)
        libmesh_not_implemented_msg("Support for Polygons/Polyhedra not yet implemented");

      elem_blk_id.push_back(subdomain_id);
      elem_type_table.push_back_entry(conv.exodus_elem_type().c_str());
      num_nodes_per_elem_vec.push_back(num_nodes_per_elem);
      num_attr_vec.push_back(0); // we don't currently use elem block attributes.
      num_edges_per_elem_vec.push_back(0); // We don't currently store any edge blocks
      num_faces_per_elem_vec.push_back(0); // We don't currently store any face blocks
      ++counter;
    }

  // Here we reserve() space so that we can push_back() onto the
  // elem_num_map in the loops below.
  this->elem_num_map.reserve(num_elem);

  // In the case of discontinuous plotting we initialize a map from
  // (element, node) pairs to the corresponding discontinuous node index.
  // This ordering must match the ordering used in write_nodal_coordinates.
  //
  // Note: This map takes the place of the libmesh_node_num_to_exodus map in
  // the discontinuous case.
  std::map<std::pair<dof_id_type, unsigned int>, dof_id_type> discontinuous_node_indices;
  dof_id_type node_counter = 1; // Exodus numbering is 1-based
  if (use_discontinuous)
  {
    for (const auto & elem : mesh.active_element_ptr_range())
      for (auto n : elem->node_index_range())
        discontinuous_node_indices[std::make_pair(elem->id(),n)] =
          node_counter++;
  }
  else
    node_counter = mesh.max_node_id() + 1; // Exodus numbering is 1-based

  if (_add_sides)
    {
      for (const Elem * elem : mesh.active_element_ptr_range())
        {
          // We'll use "past-the-end" indices to indicate side node
          // copies
          unsigned int local_node_index = elem->n_nodes();

          for (auto s : elem->side_index_range())
            {
              if (EquationSystems::redundant_added_side(*elem,s))
                continue;

              const std::vector<unsigned int> side_nodes =
                elem->nodes_on_side(s);

              for (auto n : index_range(side_nodes))
                {
                  libmesh_ignore(n);
                  discontinuous_node_indices
                    [std::make_pair(elem->id(),local_node_index++)] =
                    node_counter++;
                }
            }
        }
    }

  // Reference to the BoundaryInfo object for convenience.
  const BoundaryInfo & bi = mesh.get_boundary_info();

  // Build list of (elem, edge, id) triples
  std::vector<BoundaryInfo::BCTuple> edge_tuples = bi.build_edge_list();

  // Build the connectivity array for each edge block. The connectivity array
  // is a vector<int> with "num_edges * num_nodes_per_edge" entries. We write
  // the Exodus node numbers to the connectivity arrays so that they can
  // be used directly in the calls to exII::ex_put_conn() below. We also keep
  // track of the ElemType and the number of nodes for each boundary_id. All
  // edges with a given boundary_id must be of the same type.
  std::map<boundary_id_type, std::vector<int>> edge_id_to_conn;
  std::map<boundary_id_type, std::pair<ElemType, unsigned int>> edge_id_to_elem_type;

  std::unique_ptr<const Elem> edge;
  for (const auto & t : edge_tuples)
    {
      dof_id_type elem_id = std::get<0>(t);
      unsigned int edge_id = std::get<1>(t);
      boundary_id_type b_id = std::get<2>(t);

      // Build the edge in question
      mesh.elem_ptr(elem_id)->build_edge_ptr(edge, edge_id);

      // Error checking: make sure that all edges in this block are
      // the same geometric type.
      if (const auto check_it = edge_id_to_elem_type.find(b_id);
          check_it == edge_id_to_elem_type.end())
        {
          // Keep track of the ElemType and number of nodes in this boundary id.
          edge_id_to_elem_type[b_id] = std::make_pair(edge->type(), edge->n_nodes());
        }
      else
        {
          // Make sure the existing data is consistent
          const auto & val_pair = check_it->second;
          libmesh_error_msg_if(val_pair.first != edge->type() || val_pair.second != edge->n_nodes(),
                               "All edges in a block must have same geometric type.");
        }

      // Get reference to the connectivity array for this block
      auto & conn = edge_id_to_conn[b_id];

      // For each node on the edge, look up the exodus node id and
      // store it in the conn array. Note: all edge types have
      // identity node mappings so we don't bother with Conversion
      // objects here.
      for (auto n : edge->node_index_range())
        {
          // We look up Exodus node numbers differently if we are
          // writing a discontinuous Exodus file.
          int exodus_node_id = -1;

          if (!use_discontinuous)
            {
              dof_id_type libmesh_node_id = edge->node_ptr(n)->id();
              exodus_node_id = libmesh_map_find
                (libmesh_node_num_to_exodus, cast_int<int>(libmesh_node_id));
            }
          else
            {
              // Get the node on the element containing this edge
              // which corresponds to edge node n. Then use that id to look up
              // the exodus_node_id in the discontinuous_node_indices map.
              unsigned int pn = mesh.elem_ptr(elem_id)->local_edge_node(edge_id, n);
              exodus_node_id = libmesh_map_find
                (discontinuous_node_indices, std::make_pair(elem_id, pn));
            }

          conn.push_back(exodus_node_id);
        }
    }

  // Make sure we have the same number of edge ids that we thought we would.
  libmesh_assert(static_cast<int>(edge_id_to_conn.size()) == num_edge_blk);

  // Build data structures describing edge blocks. This information must be
  // be passed to exII::ex_put_concat_all_blocks() at the same time as the
  // information about elem blocks.
  std::vector<int> edge_blk_id;
  NamesData edge_type_table(num_edge_blk, _max_name_length);
  std::vector<int> num_edge_this_blk_vec;
  std::vector<int> num_nodes_per_edge_vec;
  std::vector<int> num_attr_edge_vec;

  // We also build a data structure of edge block names which can
  // later be passed to exII::ex_put_names().
  NamesData edge_block_names_table(num_edge_blk, _max_name_length);

  // Note: We are going to use the edge **boundary** ids as **block** ids.
  for (const auto & pr : edge_id_to_conn)
    {
      // Store the edge block id in the array to be passed to Exodus.
      boundary_id_type id = pr.first;
      edge_blk_id.push_back(id);

      // Set Exodus element type and number of nodes for this edge block.
      const auto & elem_type_node_count = edge_id_to_elem_type[id];
      const auto & conv = get_conversion(elem_type_node_count.first);
      edge_type_table.push_back_entry(conv.exodus_type.c_str());
      num_nodes_per_edge_vec.push_back(elem_type_node_count.second);

      // The number of edges is the number of entries in the connectivity
      // array divided by the number of nodes per edge.
      num_edge_this_blk_vec.push_back(pr.second.size() / elem_type_node_count.second);

      // We don't store any attributes currently
      num_attr_edge_vec.push_back(0);

      // Store the name of this edge block
      edge_block_names_table.push_back_entry(bi.get_edgeset_name(id));
    }

  // Zero-initialize and then fill in an exII::ex_block_params struct
  // with the data we have collected. This new API replaces the old
  // exII::ex_put_concat_elem_block() API, and will eventually allow
  // us to also allocate space for edge/face blocks if desired.
  //
  // TODO: It seems like we should be able to take advantage of the
  // optimization where you set define_maps==1, but when I tried this
  // I got the error: "failed to find node map size". I think the
  // problem is that we need to first specify a nonzero number of
  // node/elem maps during the call to ex_put_init_ext() in order for
  // this to work correctly.
  exII::ex_block_params params = {};

  // Set pointers for information about elem blocks.
  params.elem_blk_id = elem_blk_id.data();
  params.elem_type = elem_type_table.get_char_star_star();
  params.num_elem_this_blk = num_elem_this_blk_vec.data();
  params.num_nodes_per_elem = num_nodes_per_elem_vec.data();
  params.num_edges_per_elem = num_edges_per_elem_vec.data();
  params.num_faces_per_elem = num_faces_per_elem_vec.data();
  params.num_attr_elem = num_attr_vec.data();
  params.define_maps = 0;

  // Set pointers to edge block information only if we actually have some.
  if (num_edge_blk)
    {
      params.edge_blk_id = edge_blk_id.data();
      params.edge_type = edge_type_table.get_char_star_star();
      params.num_edge_this_blk = num_edge_this_blk_vec.data();
      params.num_nodes_per_edge = num_nodes_per_edge_vec.data();
      params.num_attr_edge = num_attr_edge_vec.data();
    }

  ex_err = exII::ex_put_concat_all_blocks(ex_id, &params);
  EX_CHECK_ERR(ex_err, "Error writing element blocks.");

  // This counter is used to fill up the libmesh_elem_num_to_exodus map in the loop below.
  unsigned libmesh_elem_num_to_exodus_counter = 0;

  // We need these later if we're adding fake sides, but we don't need
  // to recalculate it.
  auto num_elem_this_blk_it = num_elem_this_blk_vec.begin();

  // We write "fake" ids to the elem_num_map when adding fake sides.
  // I don't think it's too important exactly what fake ids are used,
  // as long as they don't conflict with any other ids that are
  // already in the elem_num_map.
  auto next_fake_id = mesh.max_elem_id() + 1; // 1-based numbering in Exodus
#ifdef LIBMESH_ENABLE_UNIQUE_ID
  if (this->set_unique_ids_from_maps)
    next_fake_id = mesh.next_unique_id();
#endif

  for (auto & [subdomain_id, element_id_vec] : this->_subdomain_map)
    {
      // Use the first element in the block to get representative
      // information for a "real" block.  Note that Exodus assumes all
      // elements in a block are of the same type!  We are using that
      // same assumption here!
      const ElemType elem_t = (subdomain_id >= this->_subdomain_id_end) ?
        ElemType(subdomain_id - this->_subdomain_id_end) :
        mesh.elem_ref(element_id_vec[0]).type();

      const auto & conv = get_conversion(elem_t);
      num_nodes_per_elem = Elem::type_to_n_nodes_map[elem_t];
      if (Elem::type_to_n_nodes_map[elem_t] == invalid_uint)
        libmesh_not_implemented_msg("Support for Polygons/Polyhedra not yet implemented");

      // If this is a *real* block, we just loop over vectors of
      // element ids to add.
      if (subdomain_id < this->_subdomain_id_end)
      {
        connect.resize(element_id_vec.size()*num_nodes_per_elem);

        for (auto i : index_range(element_id_vec))
        {
          unsigned int elem_id = element_id_vec[i];
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
          libmesh_error_msg_if(elem.type() != conv.libmesh_elem_type(),
                               "Error: Exodus requires all elements with a given subdomain ID to be the same type.\n"
                               << "Can't write both "
                               << Utility::enum_to_string(elem.type())
                               << " and "
                               << Utility::enum_to_string(conv.libmesh_elem_type())
                               << " in the same block!");

          for (unsigned int j=0; j<static_cast<unsigned int>(num_nodes_per_elem); ++j)
            {
              unsigned int connect_index   = cast_int<unsigned int>((i*num_nodes_per_elem)+j);
              unsigned elem_node_index = conv.get_inverse_node_map(j); // inverse node map is for writing.
              if (!use_discontinuous)
                {
                  // The global id for the current node in libmesh.
                  dof_id_type libmesh_node_id = elem.node_id(elem_node_index);

                  // Write the Exodus global node id associated with
                  // this libmesh node number to the connectivity
                  // array, or throw an error if it's not found.
                  connect[connect_index] =
                    libmesh_map_find(libmesh_node_num_to_exodus,
                                     cast_int<int>(libmesh_node_id));
                }
              else
                {
                  // Look up the (elem_id, elem_node_index) pair in the map.
                  connect[connect_index] =
                    libmesh_map_find(discontinuous_node_indices,
                                     std::make_pair(elem_id, elem_node_index));
                }
            } // end for(j)

          // push_back() either elem_id+1 or the current Elem's
          // unique_id+1 into the elem_num_map, depending on the value
          // of the set_unique_ids_from_maps flag.
          if (this->set_unique_ids_from_maps)
            this->elem_num_map.push_back(elem.unique_id() + 1);
          else
            this->elem_num_map.push_back(elem_id + 1);

        } // end for(i)
      }
      else // subdomain_id >= subdomain_id_end
      {
        // If this is a "fake" block of added sides, we build those as
        // we go.
        libmesh_assert(_add_sides);

        libmesh_assert(num_elem_this_blk_it != num_elem_this_blk_vec.end());
        num_elem_this_blk = *num_elem_this_blk_it;

        connect.resize(num_elem_this_blk*num_nodes_per_elem);

        std::size_t connect_index = 0;
        for (const auto & elem : mesh.active_element_ptr_range())
          {
            unsigned int local_node_index = elem->n_nodes();

            for (auto s : elem->side_index_range())
              {
                if (EquationSystems::redundant_added_side(*elem,s))
                  continue;

                if (elem->side_type(s) != elem_t)
                  continue;

                const std::vector<unsigned int> side_nodes =
                  elem->nodes_on_side(s);

                for (auto n : index_range(side_nodes))
                  {
                    libmesh_ignore(n);
                    const int exodus_node_id = libmesh_map_find
                      (discontinuous_node_indices,
                       std::make_pair(elem->id(), local_node_index++));
                    libmesh_assert_less(connect_index, connect.size());
                    connect[connect_index++] = exodus_node_id;
                  }
              }
          }

        // Store num_elem_this_blk "fake" ids into the
        // elem_num_map. Use a traditional for-loop to avoid unused
        // variable warnings about the loop counter.
        for (int i=0; i<num_elem_this_blk; ++i)
          this->elem_num_map.push_back(next_fake_id++);
      }

      ++num_elem_this_blk_it;

      ex_err = exII::ex_put_conn
        (ex_id,
         exII::EX_ELEM_BLOCK,
         subdomain_id,
         connect.data(), // node_conn
         nullptr,        // elem_edge_conn (unused)
         nullptr);       // elem_face_conn (unused)
      EX_CHECK_ERR(ex_err, "Error writing element connectivities");
    } // end for (auto & [subdomain_id, element_id_vec] : this->_subdomain_map)

  // write out the element number map that we created
  ex_err = exII::ex_put_elem_num_map(ex_id, elem_num_map.data());
  EX_CHECK_ERR(ex_err, "Error writing element map");

  // Write out the block names
  if (num_elem_blk > 0)
    {
      ex_err = exII::ex_put_names(ex_id, exII::EX_ELEM_BLOCK, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing element block names");
    }

  // Write out edge blocks if we have any
  for (const auto & pr : edge_id_to_conn)
    {
      ex_err = exII::ex_put_conn
        (ex_id,
         exII::EX_EDGE_BLOCK,
         pr.first,
         pr.second.data(), // node_conn
         nullptr,          // elem_edge_conn (unused)
         nullptr);         // elem_face_conn (unused)
      EX_CHECK_ERR(ex_err, "Error writing element connectivities");
    }

  // Write out the edge block names, if any.
  if (num_edge_blk > 0)
    {
      ex_err = exII::ex_put_names
        (ex_id,
         exII::EX_EDGE_BLOCK,
         edge_block_names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing edge block names");
    }
}




void ExodusII_IO_Helper::write_sidesets(const MeshBase & mesh)
{
  LOG_SCOPE("write_sidesets()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Maps from sideset id to lists of corresponding element ids and side ids
  std::map<int, std::vector<int>> elem_lists;
  std::map<int, std::vector<int>> side_lists;
  std::set<boundary_id_type> side_boundary_ids;

  {
    // Accumulate the vectors to pass into ex_put_side_set
    // build_side_lists() returns a vector of (elem, side, bc) tuples.
    for (const auto & t : mesh.get_boundary_info().build_side_list())
      {
        std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
        /**
         * We need to build up active elements if AMR is enabled and add
         * them to the exodus sidesets instead of the potentially inactive "parent" elements
         */
        mesh.elem_ref(std::get<0>(t)).active_family_tree_by_side(family, std::get<1>(t), false);
#else
        family.push_back(mesh.elem_ptr(std::get<0>(t)));
#endif

        for (const auto & f : family)
          {
            const auto & conv = get_conversion(mesh.elem_ptr(f->id())->type());

            // Use the libmesh to exodus data structure map to get the proper sideset IDs
            // The data structure contains the "collapsed" contiguous ids
            elem_lists[std::get<2>(t)].push_back(libmesh_elem_num_to_exodus[f->id()]);
            side_lists[std::get<2>(t)].push_back(conv.get_inverse_side_map(std::get<1>(t)));
          }
      }

    std::vector<boundary_id_type> tmp;
    mesh.get_boundary_info().build_side_boundary_ids(tmp);
    side_boundary_ids.insert(tmp.begin(), tmp.end());
  }

  {
    // add data for shell faces, if needed

    // Accumulate the vectors to pass into ex_put_side_set
    for (const auto & t : mesh.get_boundary_info().build_shellface_list())
      {
        std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
        /**
         * We need to build up active elements if AMR is enabled and add
         * them to the exodus sidesets instead of the potentially inactive "parent" elements
         */
        mesh.elem_ref(std::get<0>(t)).active_family_tree_by_side(family, std::get<1>(t), false);
#else
        family.push_back(mesh.elem_ptr(std::get<0>(t)));
#endif

        for (const auto & f : family)
          {
            const auto & conv = get_conversion(mesh.elem_ptr(f->id())->type());

            // Use the libmesh to exodus data structure map to get the proper sideset IDs
            // The data structure contains the "collapsed" contiguous ids
            elem_lists[std::get<2>(t)].push_back(libmesh_elem_num_to_exodus[f->id()]);
            side_lists[std::get<2>(t)].push_back(conv.get_inverse_shellface_map(std::get<1>(t)));
          }
      }

    std::vector<boundary_id_type> tmp;
    mesh.get_boundary_info().build_shellface_boundary_ids(tmp);
    side_boundary_ids.insert(tmp.begin(), tmp.end());
  }

  // Add any empty-but-named side boundary ids
  for (const auto & pr : mesh.get_boundary_info().get_sideset_name_map())
    side_boundary_ids.insert(pr.first);

  // Write out the sideset names, but only if there is something to write
  if (side_boundary_ids.size() > 0)
    {
      NamesData names_table(side_boundary_ids.size(), _max_name_length);

      std::vector<exII::ex_set> sets(side_boundary_ids.size());

      // Loop over "side_boundary_ids" and "sets" simultaneously
      for (auto [i, it] = std::tuple{0u, side_boundary_ids.begin()}; i<sets.size(); ++i, ++it)
        {
          boundary_id_type ss_id = *it;
          names_table.push_back_entry(mesh.get_boundary_info().get_sideset_name(ss_id));

          sets[i].id = ss_id;
          sets[i].type = exII::EX_SIDE_SET;
          sets[i].num_distribution_factor = 0;
          sets[i].distribution_factor_list = nullptr;

          if (const auto elem_it = elem_lists.find(ss_id);
              elem_it == elem_lists.end())
            {
              sets[i].num_entry = 0;
              sets[i].entry_list = nullptr;
              sets[i].extra_list = nullptr;
            }
          else
            {
              sets[i].num_entry = elem_it->second.size();
              sets[i].entry_list = elem_it->second.data();
              sets[i].extra_list = libmesh_map_find(side_lists, ss_id).data();
            }
        }

      ex_err = exII::ex_put_sets(ex_id, side_boundary_ids.size(), sets.data());
      EX_CHECK_ERR(ex_err, "Error writing sidesets");

      ex_err = exII::ex_put_names(ex_id, exII::EX_SIDE_SET, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing sideset names");
    }
}



void ExodusII_IO_Helper::write_nodesets(const MeshBase & mesh)
{
  LOG_SCOPE("write_nodesets()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // build_node_list() builds a sorted list of (node-id, bc-id) tuples
  // that is sorted by node-id, but we actually want it to be sorted
  // by bc-id, i.e. the second argument of the tuple.
  auto bc_tuples =
    mesh.get_boundary_info().build_node_list();

  // We use std::stable_sort() here so that the entries within a
  // single nodeset remain sorted in node-id order, but now the
  // smallest boundary id's nodes appear first in the list, followed
  // by the second smallest, etc. That is, we are purposely doing two
  // different sorts here, with the first one being within the
  // build_node_list() call itself.
  std::stable_sort(bc_tuples.begin(), bc_tuples.end(),
                   [](const BoundaryInfo::NodeBCTuple & t1,
                      const BoundaryInfo::NodeBCTuple & t2)
                   { return std::get<1>(t1) < std::get<1>(t2); });

  std::vector<boundary_id_type> node_boundary_ids;
  mesh.get_boundary_info().build_node_boundary_ids(node_boundary_ids);

  // Add any empty-but-named node boundary ids
  for (const auto & pair : mesh.get_boundary_info().get_nodeset_name_map())
    {
      const boundary_id_type id = pair.first;

      if (std::find(node_boundary_ids.begin(),
                    node_boundary_ids.end(), id)
            == node_boundary_ids.end())
        node_boundary_ids.push_back(id);
    }

  // Write out the nodeset names, but only if there is something to write
  if (node_boundary_ids.size() > 0)
    {
      NamesData names_table(node_boundary_ids.size(), _max_name_length);

      // Vectors to be filled and passed to exII::ex_put_concat_sets()
      // Use existing class members and avoid variable shadowing.
      nodeset_ids.clear();
      num_nodes_per_set.clear();
      num_node_df_per_set.clear();
      node_sets_node_index.clear();
      node_sets_node_list.clear();

      // Pre-allocate space
      nodeset_ids.reserve(node_boundary_ids.size());
      num_nodes_per_set.reserve(node_boundary_ids.size());
      num_node_df_per_set.resize(node_boundary_ids.size()); // all zeros
      node_sets_node_index.reserve(node_boundary_ids.size());
      node_sets_node_list.reserve(bc_tuples.size());

      // Assign entries to node_sets_node_list, keeping track of counts as we go.
      std::map<boundary_id_type, unsigned int> nodeset_counts;
      for (auto id : node_boundary_ids)
        nodeset_counts[id] = 0;

      for (const auto & t : bc_tuples)
        {
          const dof_id_type & node_id = std::get<0>(t) + 1; // Note: we use 1-based node ids in Exodus!
          const boundary_id_type & nodeset_id = std::get<1>(t);
          node_sets_node_list.push_back(node_id);
          nodeset_counts[nodeset_id] += 1;
        }

      // Fill in other indexing vectors needed by Exodus
      unsigned int running_sum = 0;
      for (const auto & pr : nodeset_counts)
        {
          nodeset_ids.push_back(pr.first);
          num_nodes_per_set.push_back(pr.second);
          node_sets_node_index.push_back(running_sum);
          names_table.push_back_entry(mesh.get_boundary_info().get_nodeset_name(pr.first));
          running_sum += pr.second;
        }

      // Fill in an exII::ex_set_specs object which can then be passed to
      // the ex_put_concat_sets() function.
      exII::ex_set_specs set_data = {};
      set_data.sets_ids = nodeset_ids.data();
      set_data.num_entries_per_set = num_nodes_per_set.data();
      set_data.num_dist_per_set = num_node_df_per_set.data(); // zeros
      set_data.sets_entry_index = node_sets_node_index.data();
      set_data.sets_dist_index = node_sets_node_index.data(); // dummy value
      set_data.sets_entry_list = node_sets_node_list.data();

      // Write all nodesets together.
      ex_err = exII::ex_put_concat_sets(ex_id, exII::EX_NODE_SET, &set_data);
      EX_CHECK_ERR(ex_err, "Error writing concatenated nodesets");

      // Write out the nodeset names
      ex_err = exII::ex_put_names(ex_id, exII::EX_NODE_SET, names_table.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing nodeset names");
    }
}



void ExodusII_IO_Helper::initialize_element_variables(std::vector<std::string> names,
                                                      const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Quick return if there are no element variables to write
  if (names.size() == 0)
    return;

  // Be sure that variables in the file match what we are asking for
  if (num_elem_vars > 0)
    {
      this->check_existing_vars(ELEMENTAL, names, this->elem_var_names);
      return;
    }

  // Quick return if we have already called this function
  if (_elem_vars_initialized)
    return;

  // Set the flag so we can skip this stuff on subsequent calls to
  // initialize_element_variables()
  _elem_vars_initialized = true;

  this->write_var_names(ELEMENTAL, names);

  // Use the truth table to indicate which subdomain/variable pairs are
  // active according to vars_active_subdomains.
  std::vector<int> truth_tab(num_elem_blk*num_elem_vars, 0);
  for (auto var_num : index_range(vars_active_subdomains))
    {
      // If the list of active subdomains is empty, it is interpreted as being
      // active on *all* subdomains.
      std::set<subdomain_id_type> current_set;
      if (vars_active_subdomains[var_num].empty())
        for (auto block_id : block_ids)
          current_set.insert(cast_int<subdomain_id_type>(block_id));
      else
        current_set = vars_active_subdomains[var_num];

      // Find index into the truth table for each id in current_set.
      for (auto block_id : current_set)
        {
          auto it = std::find(block_ids.begin(), block_ids.end(), block_id);
          libmesh_error_msg_if(it == block_ids.end(),
                               "ExodusII_IO_Helper: block id " << block_id << " not found in block_ids.");

          std::size_t block_index =
            std::distance(block_ids.begin(), it);

          std::size_t truth_tab_index = block_index*num_elem_vars + var_num;
          truth_tab[truth_tab_index] = 1;
        }
    }

  ex_err = exII::ex_put_truth_table
    (ex_id,
     exII::EX_ELEM_BLOCK,
     num_elem_blk,
     num_elem_vars,
     truth_tab.data());
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

  // Both the number of variables and their names (up to the first
  // MAX_STR_LENGTH characters) must match for the names we are
  // planning to write and the names already in the file.
  bool match =
    std::equal(names.begin(), names.end(),
               names_from_file.begin(),
               [this](const std::string & a,
                      const std::string & b) -> bool
               {
                 return a.compare(/*pos=*/0, /*len=*/_max_name_length, b) == 0;
               });

  if (!match)
    {
      libMesh::err << "Error! The Exodus file already contains the variables:" << std::endl;
      for (const auto & name : names_from_file)
        libMesh::err << name << std::endl;

      libMesh::err << "And you asked to write:" << std::endl;
      for (const auto & name : names)
        libMesh::err << name << std::endl;

      libmesh_error_msg("Cannot overwrite existing variables in Exodus II file.");
    }
}



void ExodusII_IO_Helper::write_timestep(int timestep, Real time)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (_single_precision)
    {
      float cast_time = float(time);
      ex_err = exII::ex_put_time(ex_id, timestep, &cast_time);
    }
  else
    {
      double cast_time = double(time);
      ex_err = exII::ex_put_time(ex_id, timestep, &cast_time);
    }
  EX_CHECK_ERR(ex_err, "Error writing timestep.");

  this->update();
}



void
ExodusII_IO_Helper::write_elemsets(const MeshBase & mesh)
{
  LOG_SCOPE("write_elemsets()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // TODO: Add support for named elemsets
  // NamesData names_table(elemsets.size(), _max_name_length);

  // We only need to write elemsets if the Mesh has an extra elem
  // integer called "elemset_code" defined on it.
  if (mesh.has_elem_integer("elemset_code"))
    {
      std::map<elemset_id_type, std::vector<int>> exodus_elemsets;

      unsigned int elemset_index =
        mesh.get_elem_integer_index("elemset_code");

      // Catch ids returned from MeshBase::get_elemsets() calls
      MeshBase::elemset_type set_ids;
      for (const auto & elem : mesh.element_ptr_range())
        {
          dof_id_type elemset_code =
            elem->get_extra_integer(elemset_index);

          // Look up which element set ids (if any) this elemset_code corresponds to.
          mesh.get_elemsets(elemset_code, set_ids);

          // Debugging
          // libMesh::out << "elemset_code = " << elemset_code << std::endl;
          // for (const auto & set_id : set_ids)
          //   libMesh::out << set_id << " ";
          // libMesh::out << std::endl;

          // Store this Elem id in every set to which it belongs.
          for (const auto & set_id : set_ids)
            exodus_elemsets[set_id].push_back(libmesh_elem_num_to_exodus[elem->id()]);
        }

      // Debugging: print contents of exodus_elemsets map
      // for (const auto & [set_id, elem_ids] : exodus_elemsets)
      //   {
      //     libMesh::out << "elemset " << set_id << ": ";
      //     for (const auto & elem_id : elem_ids)
      //       libMesh::out << elem_id << " ";
      //     libMesh::out << std::endl;
      //   }

      // Only continue if we actually had some elements in sets
      if (!exodus_elemsets.empty())
        {
          // Reserve space, loop over newly-created map, construct
          // exII::ex_set objects to be passed to exII::ex_put_sets(). Note:
          // we do non-const iteration since Exodus requires non-const pointers
          // to be passed to its APIs.
          std::vector<exII::ex_set> sets;
          sets.reserve(exodus_elemsets.size());

          for (auto & [elem_set_id, ids_vec] : exodus_elemsets)
            {
              // TODO: Add support for named elemsets
              // names_table.push_back_entry(mesh.get_elemset_name(elem_set_id));

              exII::ex_set & current_set = sets.emplace_back();
              current_set.id = elem_set_id;
              current_set.type = exII::EX_ELEM_SET;
              current_set.num_entry = ids_vec.size();
              current_set.num_distribution_factor = 0;
              current_set.entry_list = ids_vec.data();
              current_set.extra_list = nullptr; // extra_list is used for sidesets, not needed for elemsets
              current_set.distribution_factor_list = nullptr; // not used for elemsets
            }

          // Sanity check: make sure the number of elemsets we already wrote to the header
          // matches the number of elemsets we just constructed by looping over the Mesh.
          libmesh_assert_msg(num_elem_sets == cast_int<int>(exodus_elemsets.size()),
                             "Mesh has " << exodus_elemsets.size()
                             << " elemsets, but header was written with num_elem_sets == " << num_elem_sets);
          libmesh_assert_msg(num_elem_sets == cast_int<int>(mesh.n_elemsets()),
                             "mesh.n_elemsets() == " << mesh.n_elemsets()
                             << ", but header was written with num_elem_sets == " << num_elem_sets);

          ex_err = exII::ex_put_sets(ex_id, exodus_elemsets.size(), sets.data());
          EX_CHECK_ERR(ex_err, "Error writing elemsets");

          // TODO: Add support for named elemsets
          // ex_err = exII::ex_put_names(ex_id, exII::EX_ELEM_SET, names_table.get_char_star_star());
          // EX_CHECK_ERR(ex_err, "Error writing elemset names");
        } // end if (!exodus_elemsets.empty())
    } // end if (mesh.has_elem_integer("elemset_code"))
}



void
ExodusII_IO_Helper::
write_sideset_data(const MeshBase & mesh,
                   int timestep,
                   const std::vector<std::string> & var_names,
                   const std::vector<std::set<boundary_id_type>> & side_ids,
                   const std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals)
{
  LOG_SCOPE("write_sideset_data()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Write the sideset variable names to file. This function should
  // only be called once for SIDESET variables, repeated calls to
  // write_var_names overwrites/changes the order of names that were
  // there previously, and will mess up any data that has already been
  // written.
  this->write_var_names(SIDESET, var_names);

  // I hope that we are allowed to call read_sideset_info() even
  // though we are in the middle of writing? It seems to work provided
  // that you have already written the mesh itself... read_sideset_info()
  // fills in the following data members:
  // .) num_side_sets
  // .) ss_ids
  this->read_sideset_info();

  // Write "truth" table for sideset variables.  The function
  // exII::ex_put_variable_param() must be called before
  // exII::ex_put_truth_table(). For us, this happens during the call
  // to ExodusII_IO_Helper::write_var_names(). sset_var_tab is a logically
  // (num_side_sets x num_sset_var) integer array of 0s and 1s
  // indicating which sidesets a given sideset variable is defined on.
  std::vector<int> sset_var_tab(num_side_sets * var_names.size());

  // We now call read_sideset() once per sideset and write any sideset
  // variable values which are defined there.
  int offset=0;
  for (int ss=0; ss<num_side_sets; ++ss)
    {
      // We don't know num_sides_per_set for each set until we call
      // read_sideset(). The values for each sideset are stored (using
      // the offsets) into the 'elem_list' and 'side_list' arrays of
      // this class.
      offset += (ss > 0 ? num_sides_per_set[ss-1] : 0);
      this->read_sideset(ss, offset);

      // For each variable in var_names, write the values for the
      // current sideset, if any.
      for (auto var : index_range(var_names))
        {
          // If this var has no values on this sideset, go to the next one.
          if (!side_ids[var].count(ss_ids[ss]))
            continue;

          // Otherwise, fill in this entry of the sideset truth table.
          sset_var_tab[ss*var_names.size() + var] = 1;

          // Data vector that will eventually be passed to exII::ex_put_var().
          std::vector<Real> sset_var_vals(num_sides_per_set[ss]);

          // Get reference to the BCTuple -> Real map for this variable.
          const auto & data_map = bc_vals[var];

          // Loop over elem_list, side_list entries in current sideset.
          for (int i=0; i<num_sides_per_set[ss]; ++i)
            {
              // Get elem_id and side_id from the respective lists that
              // are filled in by calling read_sideset().
              //
              // Note: these are Exodus-specific ids, so we have to convert them
              // to libmesh ids, as that is what will be in the bc_tuples.
              //
              // TODO: we should probably consult the exodus_elem_num_to_libmesh
              // mapping in order to figure out which libmesh element id 'elem_id'
              // actually corresponds to here, instead of just assuming it will be
              // off by one. Unfortunately that data structure does not seem to
              // be used at the moment. If we assume that write_sideset_data() is
              // always called following write(), then this should be a fairly safe
              // assumption...
              dof_id_type elem_id = elem_list[i + offset] - 1;
              unsigned int side_id = side_list[i + offset] - 1;

              // Sanity check: make sure that the "off by one"
              // assumption we used above to set 'elem_id' is valid.
              libmesh_error_msg_if
                (libmesh_map_find(libmesh_elem_num_to_exodus, cast_int<int>(elem_id)) !=
                 cast_int<dof_id_type>(elem_list[i + offset]),
                 "Error mapping Exodus elem id to libmesh elem id.");

              // Map from Exodus side ids to libmesh side ids.
              const auto & conv = get_conversion(mesh.elem_ptr(elem_id)->type());

              // Map from Exodus side ids to libmesh side ids.
              unsigned int converted_side_id = conv.get_side_map(side_id);

              // Construct a key so we can quickly see whether there is any
              // data for this variable in the map.
              BoundaryInfo::BCTuple key = std::make_tuple
                (elem_id,
                 converted_side_id,
                 ss_ids[ss]);

              // Find the data for this (elem,side,id) tuple. Throw an
              // error if not found. Then store value in vector which
              // will be passed to Exodus.
              sset_var_vals[i] = libmesh_map_find(data_map, key);
            } // end for (i)

          // As far as I can tell, there is no "concat" version of writing
          // sideset data, you have to call ex_put_sset_var() once per (variable,
          // sideset) pair.
          if (sset_var_vals.size() > 0)
            {
              ex_err = exII::ex_put_var
                (ex_id,
                 timestep,
                 exII::EX_SIDE_SET,
                 var + 1, // 1-based variable index of current variable
                 ss_ids[ss],
                 num_sides_per_set[ss],
                 MappedOutputVector(sset_var_vals, _single_precision).data());
              EX_CHECK_ERR(ex_err, "Error writing sideset vars.");
            }
        } // end for (var)
    } // end for (ss)

  // Finally, write the sideset truth table.
  ex_err =
    exII::ex_put_truth_table(ex_id,
                             exII::EX_SIDE_SET,
                             num_side_sets,
                             cast_int<int>(var_names.size()),
                             sset_var_tab.data());
  EX_CHECK_ERR(ex_err, "Error writing sideset var truth table.");
}



void
ExodusII_IO_Helper::
read_sideset_data(const MeshBase & mesh,
                  int timestep,
                  std::vector<std::string> & var_names,
                  std::vector<std::set<boundary_id_type>> & side_ids,
                  std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals)
{
  LOG_SCOPE("read_sideset_data()", "ExodusII_IO_Helper");

  // This reads the sideset variable names into the local
  // sideset_var_names data structure.
  this->read_var_names(SIDESET);

  if (num_sideset_vars)
    {
      // Read the sideset data truth table
      std::vector<int> sset_var_tab(num_side_sets * num_sideset_vars);
      ex_err = exII::ex_get_truth_table
        (ex_id,
         exII::EX_SIDE_SET,
         num_side_sets,
         num_sideset_vars,
         sset_var_tab.data());
      EX_CHECK_ERR(ex_err, "Error reading sideset variable truth table.");

      // Set up/allocate space in incoming data structures.
      var_names = sideset_var_names;
      side_ids.resize(num_sideset_vars);
      bc_vals.resize(num_sideset_vars);

      // Read the sideset data.
      //
      // Note: we assume that read_sideset() has already been called
      // for each sideset, so the required values in elem_list and
      // side_list are already present.
      //
      // TODO: As a future optimization, we could read only the values
      // requested by the user by looking at the input parameter
      // var_names and checking whether it already has entries in
      // it. We could do the same thing with the input side_ids
      // container and only read values for requested sidesets.
      int offset=0;
      for (int ss=0; ss<num_side_sets; ++ss)
        {
          offset += (ss > 0 ? num_sides_per_set[ss-1] : 0);
          for (int var=0; var<num_sideset_vars; ++var)
            {
              int is_present = sset_var_tab[num_sideset_vars*ss + var];

              if (is_present)
                {
                  // Record the fact that this variable is defined on this sideset.
                  side_ids[var].insert(ss_ids[ss]);

                  // Note: the assumption here is that a previous call
                  // to this->read_sideset_info() has already set the
                  // values of num_sides_per_set, so we just use those values here.
                  std::vector<Real> sset_var_vals(num_sides_per_set[ss]);
                  ex_err = exII::ex_get_var
                    (ex_id,
                     timestep,
                     exII::EX_SIDE_SET,
                     var + 1, // 1-based sideset variable index!
                     ss_ids[ss],
                     num_sides_per_set[ss],
                     MappedInputVector(sset_var_vals, _single_precision).data());
                  EX_CHECK_ERR(ex_err, "Error reading sideset variable.");

                  for (int i=0; i<num_sides_per_set[ss]; ++i)
                    {
                      dof_id_type exodus_elem_id = elem_list[i + offset];
                      unsigned int exodus_side_id = side_list[i + offset];

                      // FIXME: We should use exodus_elem_num_to_libmesh for this,
                      // but it apparently is never set up, so just
                      // subtract 1 from the Exodus elem id.
                      dof_id_type converted_elem_id = exodus_elem_id - 1;

                      // Map Exodus side id to libmesh side id.
                      // Map from Exodus side ids to libmesh side ids.
                      const auto & conv = get_conversion(mesh.elem_ptr(converted_elem_id)->type());

                      // Map from Exodus side id to libmesh side id.
                      // Note: the mapping is defined on 0-based indices, so subtract
                      // 1 before doing the mapping.
                      unsigned int converted_side_id = conv.get_side_map(exodus_side_id - 1);

                      // Make a BCTuple key from the converted information.
                      BoundaryInfo::BCTuple key = std::make_tuple
                        (converted_elem_id,
                         converted_side_id,
                         ss_ids[ss]);

                      // Store (elem, side, b_id) tuples in bc_vals[var]
                      bc_vals[var].emplace(key, sset_var_vals[i]);
                    } // end for (i)
                } // end if (present)
            } // end for (var)
        } // end for (ss)
    } // end if (num_sideset_vars)
}


void
ExodusII_IO_Helper::
get_sideset_data_indices (const MeshBase & mesh,
                          std::map<BoundaryInfo::BCTuple, unsigned int> & bc_array_indices)
{
  // Clear any existing data, we are going to build this data structure from scratch
  bc_array_indices.clear();

  // Store the sideset data array indices.
  //
  // Note: we assume that read_sideset() has already been called
  // for each sideset, so the required values in elem_list and
  // side_list are already present.
  int offset=0;
  for (int ss=0; ss<num_side_sets; ++ss)
    {
      offset += (ss > 0 ? num_sides_per_set[ss-1] : 0);
      for (int i=0; i<num_sides_per_set[ss]; ++i)
        {
          dof_id_type exodus_elem_id = elem_list[i + offset];
          unsigned int exodus_side_id = side_list[i + offset];

          // FIXME: We should use exodus_elem_num_to_libmesh for this,
          // but it apparently is never set up, so just
          // subtract 1 from the Exodus elem id.
          dof_id_type converted_elem_id = exodus_elem_id - 1;

          // Conversion operator for this Elem type
          const auto & conv = get_conversion(mesh.elem_ptr(converted_elem_id)->type());

          // Map from Exodus side id to libmesh side id.
          // Note: the mapping is defined on 0-based indices, so subtract
          // 1 before doing the mapping.
          unsigned int converted_side_id = conv.get_side_map(exodus_side_id - 1);

          // Make a BCTuple key from the converted information.
          BoundaryInfo::BCTuple key = std::make_tuple
            (converted_elem_id,
             converted_side_id,
             ss_ids[ss]);

          // Store (elem, side, b_id) tuple with corresponding array index
          bc_array_indices.emplace(key, cast_int<unsigned int>(i));
        } // end for (i)
    } // end for (ss)
}



void ExodusII_IO_Helper::
write_nodeset_data (int timestep,
                    const std::vector<std::string> & var_names,
                    const std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                    const std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals)
{
  LOG_SCOPE("write_nodeset_data()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Write the nodeset variable names to file. This function should
  // only be called once for NODESET variables, repeated calls to
  // write_var_names() overwrites/changes the order of names that were
  // there previously, and will mess up any data that has already been
  // written.
  this->write_var_names(NODESET, var_names);

  // For all nodesets, reads and fills in the arrays:
  // nodeset_ids
  // num_nodes_per_set
  // node_sets_node_index - starting index for each nodeset in the node_sets_node_list vector
  // node_sets_node_list
  // Note: we need these arrays so that we know what data to write
  this->read_all_nodesets();

  // The "truth" table for nodeset variables. nset_var_tab is a
  // logically (num_node_sets x num_nset_var) integer array of 0s and
  // 1s indicating which nodesets a given nodeset variable is defined
  // on.
  std::vector<int> nset_var_tab(num_node_sets * var_names.size());

  for (int ns=0; ns<num_node_sets; ++ns)
  {
    // The offset into the node_sets_node_list for the current nodeset
    int offset = node_sets_node_index[ns];

    // For each variable in var_names, write the values for the
    // current nodeset, if any.
    for (auto var : index_range(var_names))
      {
        // If this var has no values on this nodeset, go to the next one.
        if (!node_boundary_ids[var].count(nodeset_ids[ns]))
          continue;

        // Otherwise, fill in this entry of the nodeset truth table.
        nset_var_tab[ns*var_names.size() + var] = 1;

        // Data vector that will eventually be passed to exII::ex_put_var().
        std::vector<Real> nset_var_vals(num_nodes_per_set[ns]);

        // Get reference to the NodeBCTuple -> Real map for this variable.
        const auto & data_map = bc_vals[var];

        // Loop over entries in current nodeset.
        for (int i=0; i<num_nodes_per_set[ns]; ++i)
          {
            // Here we convert Exodus node ids to libMesh node ids by
            // subtracting 1.  We should probably use the
            // exodus_node_num_to_libmesh data structure for this, but
            // I don't think it is set up at the time when
            // write_nodeset_data() would normally be called.
            dof_id_type libmesh_node_id = node_sets_node_list[i + offset] - 1;

            // Construct a key to look up values in data_map.
            BoundaryInfo::NodeBCTuple key =
              std::make_tuple(libmesh_node_id, nodeset_ids[ns]);

            // We require that the user provided either no values for
            // this (var, nodeset) combination (in which case we don't
            // reach this point) or a value for _every_ node in this
            // nodeset for this var, so we use the libmesh_map_find()
            // macro to check for this.
            nset_var_vals[i] = libmesh_map_find(data_map, key);
          } // end for (node in nodeset[ns])

        // Write nodeset values to Exodus file
        if (nset_var_vals.size() > 0)
          {
            ex_err = exII::ex_put_var
              (ex_id,
               timestep,
               exII::EX_NODE_SET,
               var + 1, // 1-based variable index of current variable
               nodeset_ids[ns],
               num_nodes_per_set[ns],
               MappedOutputVector(nset_var_vals, _single_precision).data());
            EX_CHECK_ERR(ex_err, "Error writing nodeset vars.");
          }
      } // end for (var in var_names)
  } // end for (ns)

  // Finally, write the nodeset truth table.
  ex_err =
    exII::ex_put_truth_table(ex_id,
                             exII::EX_NODE_SET,
                             num_node_sets,
                             cast_int<int>(var_names.size()),
                             nset_var_tab.data());
  EX_CHECK_ERR(ex_err, "Error writing nodeset var truth table.");
}



void
ExodusII_IO_Helper::
write_elemset_data (int timestep,
                    const std::vector<std::string> & var_names,
                    const std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                    const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals)
{
  LOG_SCOPE("write_elemset_data()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Write the elemset variable names to file. This function should
  // only be called once for ELEMSET variables, repeated calls to
  // write_var_names() overwrites/changes the order of names that were
  // there previously, and will mess up any data that has already been
  // written.
  this->write_var_names(ELEMSET, var_names);

  // We now call the API to read the elemset info even though we are
  // in the middle of writing. This is a bit counter-intuitive, but it
  // seems to work provided that you have already written the mesh
  // itself... read_elemset_info() fills in the following data
  // members:
  // .) id_to_elemset_names
  // .) num_elems_per_set
  // .) num_elem_df_per_set
  // .) elemset_list
  // .) elemset_id_list
  // .) id_to_elemset_names
  this->read_elemset_info();

  // The "truth" table for elemset variables. elemset_var_tab is a
  // logically (num_elem_sets x num_elemset_vars) integer array of 0s and
  // 1s indicating which elemsets a given elemset variable is defined
  // on.
  std::vector<int> elemset_var_tab(num_elem_sets * var_names.size());

  int offset=0;
  for (int es=0; es<num_elem_sets; ++es)
    {
      // Debugging
      // libMesh::out << "Writing elemset variable values for elemset "
      //              << es << ", elemset_id = " << elemset_ids[es]
      //              << std::endl;

      // We know num_elems_per_set because we called read_elemset_info() above.
      offset += (es > 0 ? num_elems_per_set[es-1] : 0);
      this->read_elemset(es, offset);

      // For each variable in var_names, write the values for the
      // current elemset, if any.
      for (auto var : index_range(var_names))
        {
          // Debugging
          // libMesh::out << "Writing elemset variable values for var " << var << std::endl;

          // If this var has no values on this elemset, go to the next one.
          if (!elemset_ids_in[var].count(elemset_ids[es]))
            continue;

          // Otherwise, fill in this entry of the nodeset truth table.
          elemset_var_tab[es*var_names.size() + var] = 1;

          // Data vector that will eventually be passed to exII::ex_put_var().
          std::vector<Real> elemset_var_vals(num_elems_per_set[es]);

          // Get reference to the (elem_id, elemset_id) -> Real map for this variable.
          const auto & data_map = elemset_vals[var];

          // Loop over entries in current elemset.
          for (int i=0; i<num_elems_per_set[es]; ++i)
            {
              // Here we convert Exodus elem ids to libMesh node ids
              // simply by subtracting 1.  We should probably use the
              // exodus_elem_num_to_libmesh data structure for this,
              // but I don't think it is set up at the time when this
              // function is normally called.
              dof_id_type libmesh_elem_id = elemset_list[i + offset] - 1;

              // Construct a key to look up values in data_map.
              std::pair<dof_id_type, elemset_id_type> key =
                std::make_pair(libmesh_elem_id, elemset_ids[es]);

              // Debugging:
              // libMesh::out << "Searching for key = (" << key.first << ", " << key.second << ")" << std::endl;

              // We require that the user provided either no values for
              // this (var, elemset) combination (in which case we don't
              // reach this point) or a value for _every_ elem in this
              // elemset for this var, so we use the libmesh_map_find()
              // macro to check for this.
              elemset_var_vals[i] = libmesh_map_find(data_map, key);
            } // end for (node in nodeset[ns])

          // Write elemset values to Exodus file
          if (elemset_var_vals.size() > 0)
            {
              ex_err = exII::ex_put_var
                (ex_id,
                 timestep,
                 exII::EX_ELEM_SET,
                 var + 1, // 1-based variable index of current variable
                 elemset_ids[es],
                 num_elems_per_set[es],
                 MappedOutputVector(elemset_var_vals, _single_precision).data());
              EX_CHECK_ERR(ex_err, "Error writing elemset vars.");
            }
        } // end for (var in var_names)
    } // end for (ns)

  // Finally, write the elemset truth table to file.
  ex_err =
    exII::ex_put_truth_table(ex_id,
                             exII::EX_ELEM_SET, // exII::ex_entity_type
                             num_elem_sets,
                             cast_int<int>(var_names.size()),
                             elemset_var_tab.data());
  EX_CHECK_ERR(ex_err, "Error writing elemset var truth table.");
}



void
ExodusII_IO_Helper::
read_elemset_data (int timestep,
                   std::vector<std::string> & var_names,
                   std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                   std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals)
{
  LOG_SCOPE("read_elemset_data()", "ExodusII_IO_Helper");

  // This reads the elemset variable names into the local
  // elemset_var_names data structure.
  this->read_var_names(ELEMSET);

  // Debugging
  // libMesh::out << "elmeset variable names:" << std::endl;
  // for (const auto & name : elemset_var_names)
  //   libMesh::out << name << " ";
  // libMesh::out << std::endl;

  if (num_elemset_vars)
    {
      // Debugging
      // std::cout << "Reading " << num_elem_sets
      //           << " elemsets and " << num_elemset_vars
      //           << " elemset variables." << std::endl;

      // Read the elemset data truth table.
      std::vector<int> elemset_var_tab(num_elem_sets * num_elemset_vars);
      exII::ex_get_truth_table(ex_id,
                               exII::EX_ELEM_SET, // exII::ex_entity_type
                               num_elem_sets,
                               num_elemset_vars,
                               elemset_var_tab.data());
      EX_CHECK_ERR(ex_err, "Error reading elemset variable truth table.");

      // Debugging
      // libMesh::out << "Elemset variable truth table:" << std::endl;
      // for (const auto & val : elemset_var_tab)
      //   libMesh::out << val << " ";
      // libMesh::out << std::endl;

      // Debugging
      // for (auto i : make_range(num_elem_sets))
      //   {
      //     for (auto j : make_range(num_elemset_vars))
      //       libMesh::out << elemset_var_tab[num_elemset_vars*i + j] << " ";
      //     libMesh::out << std::endl;
      //   }

      // Set up/allocate space in incoming data structures. All vectors are
      // num_elemset_vars in length.
      var_names = elemset_var_names;
      elemset_ids_in.resize(num_elemset_vars);
      elemset_vals.resize(num_elemset_vars);

      // Read the elemset data
      int offset=0;
      for (int es=0; es<num_elem_sets; ++es)
        {
          offset += (es > 0 ? num_elems_per_set[es-1] : 0);
          for (int var=0; var<num_elemset_vars; ++var)
            {
              int is_present = elemset_var_tab[num_elemset_vars*es + var];

              if (is_present)
                {
                  // Debugging
                  // libMesh::out << "Variable " << var << " is present on elemset " << es << std::endl;

                  // Record the fact that this variable is defined on this elemset.
                  elemset_ids_in[var].insert(elemset_ids[es]);

                  // Note: the assumption here is that a previous call
                  // to this->read_elemset_info() has already set the
                  // values of num_elems_per_set, so we just use those values here.
                  std::vector<Real> elemset_var_vals(num_elems_per_set[es]);
                  ex_err = exII::ex_get_var
                    (ex_id,
                     timestep,
                     exII::EX_ELEM_SET, // exII::ex_entity_type
                     var + 1, // 1-based sideset variable index!
                     elemset_ids[es],
                     num_elems_per_set[es],
                     MappedInputVector(elemset_var_vals, _single_precision).data());
                  EX_CHECK_ERR(ex_err, "Error reading elemset variable.");

                  for (int i=0; i<num_elems_per_set[es]; ++i)
                    {
                      dof_id_type exodus_elem_id = elemset_list[i + offset];

                      // FIXME: We should use exodus_elem_num_to_libmesh for this,
                      // but it apparently is never set up, so just
                      // subtract 1 from the Exodus elem id.
                      dof_id_type converted_elem_id = exodus_elem_id - 1;

                      // Make key based on the elem and set ids
                      auto key = std::make_pair(converted_elem_id,
                                                static_cast<elemset_id_type>(elemset_ids[es]));

                      // Store value in the map
                      elemset_vals[var].emplace(key, elemset_var_vals[i]);
                    } // end for (i)
                } // end if (present)
            } // end for (var)
        } // end for (es)
    } // end if (num_elemset_vars)
}



void ExodusII_IO_Helper::
get_elemset_data_indices (std::map<std::pair<dof_id_type, elemset_id_type>, unsigned int> & elemset_array_indices)
{
  // Clear existing data, we are going to build these data structures from scratch
  elemset_array_indices.clear();

  // Read the elemset data.
  //
  // Note: we assume that the functions
  // 1.) this->read_elemset_info() and
  // 2.) this->read_elemset()
  // have already been called, so that we already know e.g. how
  // many elems are in each set, their ids, etc.
  int offset=0;
  for (int es=0; es<num_elem_sets; ++es)
    {
      offset += (es > 0 ? num_elems_per_set[es-1] : 0);

      // Note: we don't actually call exII::ex_get_var() here because
      // we don't need the values. We only need the indices into that vector
      // for each (elem_id, elemset_id) tuple.
      for (int i=0; i<num_elems_per_set[es]; ++i)
        {
          dof_id_type exodus_elem_id = elemset_list[i + offset];

          // FIXME: We should use exodus_elem_num_to_libmesh for this,
          // but it apparently is never set up, so just
          // subtract 1 from the Exodus elem id.
          dof_id_type converted_elem_id = exodus_elem_id - 1;

          // Make key based on the elem and set ids
          // Make a NodeBCTuple key from the converted information.
          auto key = std::make_pair(converted_elem_id,
                                    static_cast<elemset_id_type>(elemset_ids[es]));

          // Store the array index of this (node, b_id) tuple
          elemset_array_indices.emplace(key, cast_int<unsigned int>(i));
        } // end for (i)
    } // end for (es)
}



void ExodusII_IO_Helper::
read_nodeset_data (int timestep,
                   std::vector<std::string> & var_names,
                   std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                   std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals)
{
  LOG_SCOPE("read_nodeset_data()", "ExodusII_IO_Helper");

  // This reads the sideset variable names into the local
  // sideset_var_names data structure.
  this->read_var_names(NODESET);

  if (num_nodeset_vars)
    {
      // Read the nodeset data truth table
      std::vector<int> nset_var_tab(num_node_sets * num_nodeset_vars);
      ex_err = exII::ex_get_truth_table
        (ex_id,
         exII::EX_NODE_SET,
         num_node_sets,
         num_nodeset_vars,
         nset_var_tab.data());
      EX_CHECK_ERR(ex_err, "Error reading nodeset variable truth table.");

      // Set up/allocate space in incoming data structures.
      var_names = nodeset_var_names;
      node_boundary_ids.resize(num_nodeset_vars);
      bc_vals.resize(num_nodeset_vars);

      // Read the nodeset data.
      //
      // Note: we assume that the functions
      // 1.) this->read_nodeset_info() and
      // 2.) this->read_all_nodesets()
      // have already been called, so that we already know e.g. how
      // many nodes are in each set, their ids, etc.
      //
      // TODO: As a future optimization, we could read only the values
      // requested by the user by looking at the input parameter
      // var_names and checking whether it already has entries in
      // it.
      int offset=0;
      for (int ns=0; ns<num_node_sets; ++ns)
        {
          offset += (ns > 0 ? num_nodes_per_set[ns-1] : 0);
          for (int var=0; var<num_nodeset_vars; ++var)
            {
              int is_present = nset_var_tab[num_nodeset_vars*ns + var];

              if (is_present)
                {
                  // Record the fact that this variable is defined on this nodeset.
                  node_boundary_ids[var].insert(nodeset_ids[ns]);

                  // Note: the assumption here is that a previous call
                  // to this->read_nodeset_info() has already set the
                  // values of num_nodes_per_set, so we just use those values here.
                  std::vector<Real> nset_var_vals(num_nodes_per_set[ns]);
                  ex_err = exII::ex_get_var
                    (ex_id,
                     timestep,
                     exII::EX_NODE_SET,
                     var + 1, // 1-based nodeset variable index!
                     nodeset_ids[ns],
                     num_nodes_per_set[ns],
                     MappedInputVector(nset_var_vals, _single_precision).data());
                  EX_CHECK_ERR(ex_err, "Error reading nodeset variable.");

                  for (int i=0; i<num_nodes_per_set[ns]; ++i)
                    {
                      // The read_all_nodesets() function now reads all the node ids into the
                      // node_sets_node_list vector, which is of length "total_nodes_in_all_sets"
                      // The old read_nodset() function is no longer called as far as I can tell,
                      // and should probably be removed? The "offset" that we are using only
                      // depends on the current nodeset index and the num_nodes_per_set vector,
                      // which gets filled in by the call to read_all_nodesets().
                      dof_id_type exodus_node_id = node_sets_node_list[i + offset];

                      // FIXME: We should use exodus_node_num_to_libmesh for this,
                      // but it apparently is never set up, so just
                      // subtract 1 from the Exodus node id.
                      dof_id_type converted_node_id = exodus_node_id - 1;

                      // Make a NodeBCTuple key from the converted information.
                      BoundaryInfo::NodeBCTuple key = std::make_tuple
                        (converted_node_id, nodeset_ids[ns]);

                      // Store (node, b_id) tuples in bc_vals[var]
                      bc_vals[var].emplace(key, nset_var_vals[i]);
                    } // end for (i)
                } // end if (present)
            } // end for (var)
        } // end for (ns)
    } // end if (num_nodeset_vars)
}



void
ExodusII_IO_Helper::
get_nodeset_data_indices (std::map<BoundaryInfo::NodeBCTuple, unsigned int> & bc_array_indices)
{
  // Clear existing data, we are going to build these data structures from scratch
  bc_array_indices.clear();

  // Read the nodeset data.
  //
  // Note: we assume that the functions
  // 1.) this->read_nodeset_info() and
  // 2.) this->read_all_nodesets()
  // have already been called, so that we already know e.g. how
  // many nodes are in each set, their ids, etc.
  int offset=0;
  for (int ns=0; ns<num_node_sets; ++ns)
    {
      offset += (ns > 0 ? num_nodes_per_set[ns-1] : 0);
      // Note: we don't actually call exII::ex_get_var() here because
      // we don't need the values. We only need the indices into that vector
      // for each (node_id, boundary_id) tuple.
      for (int i=0; i<num_nodes_per_set[ns]; ++i)
        {
          // The read_all_nodesets() function now reads all the node ids into the
          // node_sets_node_list vector, which is of length "total_nodes_in_all_sets"
          // The old read_nodset() function is no longer called as far as I can tell,
          // and should probably be removed? The "offset" that we are using only
          // depends on the current nodeset index and the num_nodes_per_set vector,
          // which gets filled in by the call to read_all_nodesets().
          dof_id_type exodus_node_id = node_sets_node_list[i + offset];

          // FIXME: We should use exodus_node_num_to_libmesh for this,
          // but it apparently is never set up, so just
          // subtract 1 from the Exodus node id.
          dof_id_type converted_node_id = exodus_node_id - 1;

          // Make a NodeBCTuple key from the converted information.
          BoundaryInfo::NodeBCTuple key = std::make_tuple
            (converted_node_id, nodeset_ids[ns]);

          // Store the array index of this (node, b_id) tuple
          bc_array_indices.emplace(key, cast_int<unsigned int>(i));
        } // end for (i)
    } // end for (ns)
}

void ExodusII_IO_Helper::write_element_values
(const MeshBase & mesh,
 const std::vector<Real> & values,
 int timestep,
 const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains)
{
  LOG_SCOPE("write_element_values()", "ExodusII_IO_Helper");

  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Ask the file how many element vars it has, store it in the num_elem_vars variable.
  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_ELEM_BLOCK, &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  // We might be appending values to an existing file, in which case
  // we haven't done an initialize() and we need to build the
  // subdomain map here.
  if (this->_subdomain_map.empty())
    this->build_subdomain_map(mesh, false);

  // Use mesh.n_elem() to access into the values vector rather than
  // the number of elements the Exodus writer thinks the mesh has,
  // which may not include inactive elements.
  dof_id_type n_elem = mesh.n_elem();

  // Sanity check: we must have an entry in vars_active_subdomains for
  // each variable that we are potentially writing out.
  libmesh_assert_equal_to
    (vars_active_subdomains.size(),
     static_cast<unsigned>(num_elem_vars));

  // For each variable, create a 'data' array which holds all the elemental variable
  // values *for a given block* on this processor, then write that data vector to file
  // before moving onto the next block.
  for (unsigned int var_id=0; var_id<static_cast<unsigned>(num_elem_vars); ++var_id)
    {
      // The size of the subdomain map is the number of blocks.
      auto it = this->_subdomain_map.begin();

      // Reference to the set of active subdomains for the current variable.
      const auto & active_subdomains
        = vars_active_subdomains[var_id];

      for (unsigned int j=0; it!=this->_subdomain_map.end(); ++it, ++j)
        {
          // Skip any variable/subdomain pairs that are inactive.
          // Note that if active_subdomains is empty, it is interpreted
          // as being active on *all* subdomains.
          if (!(active_subdomains.empty() || active_subdomains.count(it->first)))
            continue;

          // Get reference to list of elem ids which are in the
          // current subdomain and count, allocate storage to hold
          // data that will be written to file.
          const auto & elem_nums = it->second;
          const unsigned int num_elems_this_block =
            cast_int<unsigned int>(elem_nums.size());
          std::vector<Real> data(num_elems_this_block);

          // variable-major ordering is:
          // (u1, u2, u3, ..., uN), (v1, v2, v3, ..., vN), ...
          // where N is the number of elements.
          for (unsigned int k=0; k<num_elems_this_block; ++k)
            data[k] = values[var_id*n_elem + elem_nums[k]];

          ex_err = exII::ex_put_var
            (ex_id,
             timestep,
             exII::EX_ELEM_BLOCK,
             var_id+1,
             this->get_block_id(j),
             num_elems_this_block,
             MappedOutputVector(data, _single_precision).data());

          EX_CHECK_ERR(ex_err, "Error writing element values.");
        }
    }

  this->update();
}



void ExodusII_IO_Helper::write_element_values_element_major
(const MeshBase & mesh,
 const std::vector<Real> & values,
 int timestep,
 const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains,
 const std::vector<std::string> & derived_var_names,
 const std::map<subdomain_id_type, std::vector<std::string>> & subdomain_to_var_names)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // Ask the file how many element vars it has, store it in the num_elem_vars variable.
  ex_err = exII::ex_get_variable_param(ex_id, exII::EX_ELEM_BLOCK, &num_elem_vars);
  EX_CHECK_ERR(ex_err, "Error reading number of elemental variables.");

  // We will eventually loop over the element blocks (subdomains) and
  // write the data one block (subdomain) at a time. Build a data
  // structure that keeps track of how many elements are in each
  // subdomain. This will allow us to reserve space in the data vector
  // we are going to write.
  std::map<subdomain_id_type, unsigned int> subdomain_to_n_elem;
  for (const auto & elem : mesh.active_element_ptr_range())
    subdomain_to_n_elem[elem->subdomain_id()] += 1;

  // Sanity check: we must have an entry in vars_active_subdomains for
  // each variable that we are potentially writing out.
  libmesh_assert_equal_to
    (vars_active_subdomains.size(),
     static_cast<unsigned>(num_elem_vars));

  // The size of the subdomain map is the number of blocks.
  auto subdomain_to_n_elem_iter = subdomain_to_n_elem.begin();

  // Store range of active Elem pointers. We are going to loop over
  // the elements n_vars * n_subdomains times, so let's make sure
  // the predicated iterators aren't slowing us down too much.
  ConstElemRange elem_range
    (mesh.active_elements_begin(),
     mesh.active_elements_end());

  for (unsigned int sbd_idx=0;
       subdomain_to_n_elem_iter != subdomain_to_n_elem.end();
       ++subdomain_to_n_elem_iter, ++sbd_idx)
    for (unsigned int var_id=0; var_id<static_cast<unsigned>(num_elem_vars); ++var_id)
      {
        // Reference to the set of active subdomains for the current variable.
        const auto & active_subdomains
          = vars_active_subdomains[var_id];

        // If the vars_active_subdomains container passed to this function
        // has an empty entry, it means the variable really is not active on
        // _any_ subdomains, not that it is active on _all_ subdomains. This
        // is just due to the way that we build the vars_active_subdomains
        // container.
        if (!active_subdomains.count(subdomain_to_n_elem_iter->first))
          continue;

        // Vector to hold values that will be written to Exodus file.
        std::vector<Real> data;
        data.reserve(subdomain_to_n_elem_iter->second);

        unsigned int values_offset = 0;
        for (auto & elem : elem_range)
          {
            // We'll use the Elem's subdomain id in several places below.
            subdomain_id_type sbd_id = elem->subdomain_id();

            // Get reference to the list of variable names defining
            // the indexing for the current Elem's subdomain.
            auto subdomain_to_var_names_iter =
              subdomain_to_var_names.find(sbd_id);

            // It's possible, but unusual, for there to be an Elem
            // from a subdomain that has no active variables from the
            // set of variables we are currently writing. If that
            // happens, we can just go to the next Elem because we
            // don't need to advance the offset into the values
            // vector, etc.
            if (subdomain_to_var_names_iter == subdomain_to_var_names.end())
              continue;

            const auto & var_names_this_sbd
              = subdomain_to_var_names_iter->second;

            // Only extract values if Elem is in the current subdomain.
            if (sbd_id == subdomain_to_n_elem_iter->first)
              {
                // Location of current var_id in the list of all variables on this
                // subdomain. FIXME: linear search but it's over a typically relatively
                // short vector of active variable names on this subdomain. We could do
                // a nested std::map<string,index> instead of a std::vector where the
                // location of the string is implicitly the index..
                auto pos =
                  std::find(var_names_this_sbd.begin(),
                            var_names_this_sbd.end(),
                            derived_var_names[var_id]);

                libmesh_error_msg_if(pos == var_names_this_sbd.end(),
                                     "Derived name " << derived_var_names[var_id] << " not found!");

                // Find the current variable's location in the list of all variable
                // names on the current Elem's subdomain.
                auto true_index =
                  std::distance(var_names_this_sbd.begin(), pos);

                data.push_back(values[values_offset + true_index]);
              }

            // The "true" offset is how much we have to advance the index for each Elem
            // in this subdomain.
            auto true_offset = var_names_this_sbd.size();

            // Increment to the next Elem's values
            values_offset += true_offset;
          } // for elem

        // Now write 'data' to Exodus file, in single precision if requested.
        if (!data.empty())
          {
            ex_err = exII::ex_put_var
              (ex_id,
               timestep,
               exII::EX_ELEM_BLOCK,
               var_id+1,
               this->get_block_id(sbd_idx),
               data.size(),
               MappedOutputVector(data, _single_precision).data());

            EX_CHECK_ERR(ex_err, "Error writing element values.");
          }
      } // for each var_id

  this->update();
}



void
ExodusII_IO_Helper::write_nodal_values(int var_id,
                                       const std::vector<Real> & values,
                                       int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (!values.empty())
    {
      libmesh_assert_equal_to(values.size(), std::size_t(num_nodes));

      ex_err = exII::ex_put_var
        (ex_id,
         timestep,
         exII::EX_NODAL,
         var_id,
         1, // exII::ex_entity_id, not sure exactly what this is but in the ex_put_nodal_var.c shim, they pass 1
         num_nodes,
         MappedOutputVector(values, _single_precision).data());

      EX_CHECK_ERR(ex_err, "Error writing nodal values.");

      this->update();
    }
}



void ExodusII_IO_Helper::write_information_records(const std::vector<std::string> & records)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  // There may already be information records in the file (for
  // example, if we're appending) and in that case, according to the
  // Exodus documentation, writing more information records is not
  // supported.
  int num_info = inquire(*this, exII::EX_INQ_INFO, "Error retrieving the number of information records from file!");
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
      for (const auto & record : records)
        info.push_back_entry(record);

      ex_err = exII::ex_put_info(ex_id, num_records, info.get_char_star_star());
      EX_CHECK_ERR(ex_err, "Error writing global values.");

      this->update();
    }
}



void ExodusII_IO_Helper::write_global_values(const std::vector<Real> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  if (!values.empty())
    {
      ex_err = exII::ex_put_var
        (ex_id,
         timestep,
         exII::EX_GLOBAL,
         1, // var index
         0, // obj_id (not used)
         num_global_vars,
         MappedOutputVector(values, _single_precision).data());

      EX_CHECK_ERR(ex_err, "Error writing global values.");

      this->update();
    }
}



void ExodusII_IO_Helper::update()
{
  ex_err = exII::ex_update(ex_id);
  EX_CHECK_ERR(ex_err, "Error flushing buffers to file.");
}



void ExodusII_IO_Helper::read_global_values(std::vector<Real> & values, int timestep)
{
  if ((_run_only_on_proc0) && (this->processor_id() != 0))
    return;

  values.clear();
  values.resize(num_global_vars);
  ex_err = exII::ex_get_var
    (ex_id,
     timestep,
     exII::EX_GLOBAL,
     1, // var_index
     1, // obj_id
     num_global_vars,
     MappedInputVector(values, _single_precision).data());

  EX_CHECK_ERR(ex_err, "Error reading global values.");
}



void ExodusII_IO_Helper::use_mesh_dimension_instead_of_spatial_dimension(bool val)
{
  _use_mesh_dimension_instead_of_spatial_dimension = val;
}


void ExodusII_IO_Helper::set_hdf5_writing(bool write_hdf5)
{
  _write_hdf5 = write_hdf5;
}


void ExodusII_IO_Helper::set_max_name_length(unsigned int max_length)
{
  // Opt mode error, because this may be exposed to users
  libmesh_error_msg_if (max_length > libmesh_max_str_length,
                        "Exodus maximum name length is limited to " <<
                        libmesh_max_str_length << " characters");

  // Devel+dbg mode assertion, because developers should do better
  libmesh_assert(!opened_for_writing);

  _max_name_length = max_length;
}


void ExodusII_IO_Helper::write_as_dimension(unsigned dim)
{
  _write_as_dimension = dim;
}



void ExodusII_IO_Helper::set_coordinate_offset(Point p)
{
  _coordinate_offset = p;
}


std::vector<std::string>
ExodusII_IO_Helper::get_complex_names(const std::vector<std::string> & names,
                                      bool write_complex_abs) const
{
  std::vector<std::string> complex_names;

  // This will loop over all names and create new "complex" names
  // (i.e. names that start with r_, i_ or a_)
  for (const auto & name : names)
    {
      complex_names.push_back("r_" + name);
      complex_names.push_back("i_" + name);
      if (write_complex_abs)
        complex_names.push_back("a_" + name);
    }

  return complex_names;
}



std::vector<std::set<subdomain_id_type>>
ExodusII_IO_Helper::
get_complex_vars_active_subdomains
(const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains,
 bool write_complex_abs) const
{
  std::vector<std::set<subdomain_id_type>> complex_vars_active_subdomains;

  for (auto & s : vars_active_subdomains)
    {
      // Push back the same data enough times for the real, imag, (and
      // possibly modulus) for the complex-valued solution.
      complex_vars_active_subdomains.push_back(s);
      complex_vars_active_subdomains.push_back(s);
      if (write_complex_abs)
        complex_vars_active_subdomains.push_back(s);
    }

  return complex_vars_active_subdomains;
}



std::map<subdomain_id_type, std::vector<std::string>>
ExodusII_IO_Helper::
get_complex_subdomain_to_var_names
(const std::map<subdomain_id_type, std::vector<std::string>> & subdomain_to_var_names,
 bool write_complex_abs) const
{
  // Eventual return value
  std::map<subdomain_id_type, std::vector<std::string>> ret;

  unsigned int num_complex_outputs = write_complex_abs ? 3 : 2;

  for (const auto & pr : subdomain_to_var_names)
    {
      // Initialize entry for current subdomain
      auto & vec = ret[pr.first];

      // Get list of non-complex variable names active on this subdomain.
      const auto & varnames = pr.second;

      // Allocate space for the complex-valued entries
      vec.reserve(num_complex_outputs * varnames.size());

      // For each varname in the input map, write three variable names
      // to the output formed by prepending "r_", "i_", and "a_",
      // respectively.
      for (const auto & varname : varnames)
        {
          vec.push_back("r_" + varname);
          vec.push_back("i_" + varname);
          if (write_complex_abs)
            vec.push_back("a_" + varname);
        }
    }
  return ret;
}



void ExodusII_IO_Helper::build_subdomain_map(const MeshBase & mesh, bool local)
{
  // Start from scratch
  this->_subdomain_map.clear();
  this->_subdomain_id_end = 0;

  // If we've been asked to add side elements, those will go in
  // their own blocks.
  if (this->_add_sides)
    {
      std::set<subdomain_id_type> sbd_ids;
      mesh.subdomain_ids(sbd_ids);
      if (!sbd_ids.empty())
        this->_subdomain_id_end = *sbd_ids.rbegin()+1;
    }

  // Loop through element and map between block and element vector.
  const auto range = local ? mesh.active_local_element_ptr_range() :
      mesh.active_element_ptr_range();
  for (const auto & elem : range)
    {
      // We skip writing infinite elements to the Exodus file, so
      // don't put them in the subdomain_map. That way the number of
      // blocks should be correct.
      if (elem->infinite())
        continue;

      this->_subdomain_map[ elem->subdomain_id() ].push_back(elem->id());

      // If we've been asked to add side elements, those will go in their own
      // blocks.  We don't have any ids to list for elements that don't
      // explicitly exist in the mesh, but we do an entry to keep
      // track of the number of elements we'll add in each new block.
      if (this->_add_sides)
        for (auto s : elem->side_index_range())
          {
            if (EquationSystems::redundant_added_side(*elem,s))
              continue;

            auto & marker =
              this->_subdomain_map[this->_subdomain_id_end + elem->side_type(s)];
            if (marker.empty())
              marker.push_back(1);
            else
              ++marker[0];
          }
    }

  if (!this->_add_sides && !this->_subdomain_map.empty())
    this->_subdomain_id_end = this->_subdomain_map.rbegin()->first + 1;
}



void ExodusII_IO_Helper::calculate_added_side_node_offsets(const MeshBase & mesh)
{

  // If we're adding face elements they'll need copies of their nodes.
  // We also have to count of how many nodes (and gaps between nodes!)
  // are on each processor, to calculate offsets for any nodal data
  // writing later.
  _added_side_node_offsets.clear();
  if (!_add_sides)
    return;

  dof_id_type num_side_elem = 0;
  dof_id_type num_local_side_nodes = 0;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      for (auto s : elem->side_index_range())
        {
          if (EquationSystems::redundant_added_side(*elem,s))
            continue;

          num_side_elem++;
          num_local_side_nodes += elem->nodes_on_side(s).size();
        }
    }

  mesh.comm().sum(num_side_elem);
  num_elem += num_side_elem;

  mesh.comm().allgather(num_local_side_nodes, _added_side_node_offsets);
  const processor_id_type n_proc = mesh.n_processors();
  libmesh_assert_equal_to(n_proc, _added_side_node_offsets.size());

  for (auto p : make_range(n_proc-1))
    _added_side_node_offsets[p+1] += _added_side_node_offsets[p];

  num_nodes = _added_side_node_offsets[n_proc-1];

  dof_id_type n_local_nodes = cast_int<dof_id_type>
    (std::distance(mesh.local_nodes_begin(),
                   mesh.local_nodes_end()));
  dof_id_type n_total_nodes = n_local_nodes;
  mesh.comm().sum(n_total_nodes);

  const dof_id_type max_nn   = mesh.max_node_id();
  const dof_id_type n_gaps = max_nn - n_total_nodes;
  const dof_id_type gaps_per_processor = n_gaps / n_proc;
  const dof_id_type remainder_gaps = n_gaps % n_proc;

  n_local_nodes = n_local_nodes +      // Actual nodes
                  gaps_per_processor + // Our even share of gaps
                  (mesh.processor_id() < remainder_gaps); // Leftovers

  mesh.comm().allgather(n_local_nodes, _true_node_offsets);
  for (auto p : make_range(n_proc-1))
    _true_node_offsets[p+1] += _true_node_offsets[p];
  libmesh_assert_equal_to(_true_node_offsets[n_proc-1], mesh.max_node_id());
}



dof_id_type ExodusII_IO_Helper::node_id_to_vec_id(dof_id_type n) const
{
  if (_added_side_node_offsets.empty())
    return n;

  // Find the processor id that has node_id in the parallel vec
  const auto lb = std::upper_bound(_true_node_offsets.begin(),
                                   _true_node_offsets.end(), n);
  libmesh_assert(lb != _true_node_offsets.end());
  const processor_id_type p = lb - _true_node_offsets.begin();

  return n + (p ? _added_side_node_offsets[p-1] : 0);
}



dof_id_type ExodusII_IO_Helper::added_node_offset_on(processor_id_type p) const
{
  libmesh_assert (p < _true_node_offsets.size());
  const dof_id_type added_node_offsets =
    (_added_side_node_offsets.empty() || !p) ? 0 :
    _added_side_node_offsets[p-1];
  return _true_node_offsets[p] + added_node_offsets;
}



int ExodusII_IO_Helper::Conversion::get_node_map(int i) const
{
  if (!node_map)
    return i;

  libmesh_assert_less (i, node_map->size());
  return (*node_map)[i];
}



int ExodusII_IO_Helper::Conversion::get_inverse_node_map(int i) const
{
  if (!inverse_node_map)
    return i;

  libmesh_assert_less (i, inverse_node_map->size());
  return (*inverse_node_map)[i];
}



int ExodusII_IO_Helper::Conversion::get_side_map(int i) const
{
  if (!side_map)
    return i;

  // If we asked for a side that doesn't exist, return an invalid_id
  // and allow higher-level code to handle it.
  if (static_cast<size_t>(i) >= side_map->size())
    return invalid_id;

  return (*side_map)[i];
}



int ExodusII_IO_Helper::Conversion::get_inverse_side_map(int i) const
{
  // For identity side mappings, we our convention is to return a 1-based index.
  if (!inverse_side_map)
    return i + 1;

  libmesh_assert_less (i, inverse_side_map->size());
  return (*inverse_side_map)[i];
}



/**
 * \returns The ith component of the shellface map for this element.
 * \note Nothing is currently using this.
 */
int ExodusII_IO_Helper::Conversion::get_shellface_map(int i) const
{
  if (!shellface_map)
    return i;

  libmesh_assert_less (i, shellface_map->size());
  return (*shellface_map)[i];
}



int ExodusII_IO_Helper::Conversion::get_inverse_shellface_map(int i) const
{
  if (!inverse_shellface_map)
    return i + 1;

  libmesh_assert_less (i, inverse_shellface_map->size());
  return (*inverse_shellface_map)[i];
}



ElemType ExodusII_IO_Helper::Conversion::libmesh_elem_type() const
{
  return libmesh_type;
}



std::string ExodusII_IO_Helper::Conversion::exodus_elem_type() const
{
  return exodus_type;
}



/**
 * \returns The shellface index offset.
 */
std::size_t ExodusII_IO_Helper::Conversion::get_shellface_index_offset() const
{
  return shellface_index_offset;
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

      // Properly terminate these C-style strings, just to be safe.
      data_table[i][0] = '\0';

      // Set pointer into the data_table
      data_table_pointers[i] = data_table[i].data();
    }
}



void ExodusII_IO_Helper::NamesData::push_back_entry(const std::string & name)
{
  libmesh_assert_less (counter, table_size);

  // 1.) Copy the C++ string into the vector<char>...
  size_t num_copied = name.copy(data_table[counter].data(), data_table[counter].size()-1);

  // 2.) ...And null-terminate it.
  data_table[counter][num_copied] = '\0';

  // Go to next row
  ++counter;
}



char ** ExodusII_IO_Helper::NamesData::get_char_star_star()
{
  return data_table_pointers.data();
}



char * ExodusII_IO_Helper::NamesData::get_char_star(int i)
{
  libmesh_error_msg_if(static_cast<unsigned>(i) >= table_size,
                       "Requested char * " << i << " but only have " << table_size << "!");

  return data_table[i].data();
}



} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_EXODUS_API
