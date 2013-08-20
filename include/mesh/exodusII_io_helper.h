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

#ifndef LIBMESH_EXODUSII_IO_HELPER_H
#define LIBMESH_EXODUSII_IO_HELPER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_EXODUS_API

// Local includes
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>

// Macro to simplify checking Exodus error codes
#define EX_CHECK_ERR(code, msg) \
  do { \
    if ((code) < 0) { \
      libMesh::err << (msg) << std::endl; \
      libmesh_error(); \
    } } while(0)


namespace libMesh
{


namespace exII {
  extern "C" {
#include "exodusII.h" // defines MAX_LINE_LENGTH, MAX_STR_LENGTH used later
  }
}

/**
 * This is the \p ExodusII_IO_Helper class.  This class hides the implementation
 * details of interfacing with the Exodus binary format.
 *
 * @author Johw W. Peterson, 2002.
 */
class ExodusII_IO_Helper : public ParallelObject
{
public:

  /**
   * Constructor. Automatically initializes all the private members of
   * the class.  Also allows you to set the verbosity level to v=true
   * (on) or v=false (off).  The second argument, if true, tells the class to only
   * perform its actions if running on processor zero.  If you initialize this
   * to false, the writing methods will run on all processors instead.
   */
  ExodusII_IO_Helper(const ParallelObject &parent,
                     bool v=false,
                     bool run_only_on_proc0=true);
  /**
   * Destructor.
   */
  virtual ~ExodusII_IO_Helper();

  /**
   * Returns true once create() has been successfully called, and
   * false otherwise.  create() is called to open a new file for
   * _writing_.
   */
  bool created();

  /**
   * Returns true once open() has been successfully called, and false
   * otherwise.  open() is called to open an existing for _reading_.
   */
  bool opened();

  /**
   * Get/set flag telling whether message printing is on or off.
   */
  void verbose (bool set_verbosity);

  /**
   * @returns the current element type.  Note: the default behavior is
   * for this value to be in all capital letters, e.g. \p HEX27.
   */
  const char* get_elem_type() const;

  /**
   * Opens an \p ExodusII mesh file named \p filename for reading.
   */
  void open(const char* filename);

  /**
   * Reads an \p ExodusII mesh file header.
   */
  void read_header();

  /**
   * Prints the \p ExodusII mesh file header, which includes the mesh
   * title, the number of nodes, number of elements, mesh dimension,
   * number of sidesets, and number of nodesets.
   */
  void print_header();

  /**
   * Reads the nodal data (x,y,z coordinates) from the \p ExodusII
   * mesh file.
   */
  void read_nodes();

  /**
   * Reads the optional \p node_num_map from the \p ExodusII mesh
   * file.
   */
  void read_node_num_map();

  /**
   * Prints the nodal information, by default to \p libMesh::out.
   */
  void print_nodes(std::ostream &out = libMesh::out);

  /**
   * Reads information for all of the blocks in the \p ExodusII mesh
   * file.
   */
  void read_block_info();

  /**
   * Get the block number for the given block index.
   */
  int get_block_id(int index);

  /**
   * Get the block name for the given block index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_block_name(int index);

  /**
   * Get the side set id for the given side set index.
   */
  int get_side_set_id(int index);

  /**
   * Get the side set name for the given side set index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_side_set_name(int index);

  /**
   * Get the node set id for the given node set index.
   */
  int get_node_set_id(int index);

  /**
   * Get the node set name for the given node set index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_node_set_name(int index);

  /**
   * Reads all of the element connectivity for block \p block in the
   * \p ExodusII mesh file.
   */
  void read_elem_in_block(int block);

  /**
   * Reads the optional \p node_num_map from the \p ExodusII mesh
   * file.
   */
  void read_elem_num_map();

  /**
   * Reads information about all of the sidesets in the \p ExodusII
   * mesh file.
   */
  void read_sideset_info();

  /**
   * Reads information about all of the nodesets in the \p ExodusII
   * mesh file.
   */
  void read_nodeset_info();

  /**
   * Reads information about sideset \p id and inserts it into the
   * global sideset array at the position \p offset.
   */
  void read_sideset(int id, int offset);

  /**
   * Reads information about nodeset \p id and inserts it into the
   * global nodeset array at the position \p offset.
   */
  void read_nodeset(int id);

  /**
   * Closes the \p ExodusII mesh file.
   */
  void close();

  /**
   * Generic inquiry, returns the value
   */
  int inquire(int req_info, std::string error_msg="");

  /**
   * Reads and stores the timesteps in the 'time_steps' array.
   */
  void read_time_steps();

  /**
   * Reads the nodal variable names and stores them in the 'nodal_var_names' array.
   */
  void read_nodal_var_names();

  /**
   * Reads the nodal values for the variable 'nodal_var_name' at the
   * specified time into the 'nodal_var_values' array.
   */
  void read_nodal_var_values(std::string nodal_var_name, int time_step);

  /**
   * Reads the elemental variable names and stores them in the 'elem_var_names' array.
   */
  void read_elemental_var_names();

  /**
   * Reads elemental values for the variable 'elemental_var_name' at the
   * specified timestep into the 'elem_var_values' array.
   */
  void read_elemental_var_values(std::string elemental_var_name, int time_step);

  /**
   * Opens an \p ExodusII mesh file named \p filename for writing.
   */
  virtual void create(std::string filename);

  /**
   * Initializes the Exodus file
   */
  virtual void initialize(std::string title, const MeshBase & mesh);

  /**
   * Initializes the Exodus file
   */
  void initialize_discontinuous(std::string title, const MeshBase & mesh);

  /**
   * Writes the nodal coordinates contained in "mesh"
   */
  virtual void write_nodal_coordinates(const MeshBase & mesh);

  /**
   * Writes the nodal coordinates contained in "mesh"
   */
  void write_nodal_coordinates_discontinuous(const MeshBase & mesh);

  /**
   * Writes the elements contained in "mesh". FIXME: This only works
   * for Meshes having a single type of element!
   */
  virtual void write_elements(const MeshBase & mesh);

  /**
   * Writes the elements contained in "mesh". FIXME: This only works
   * for Meshes having a single type of element!
   */
  void write_elements_discontinuous(const MeshBase & mesh);

  /**
   * Writes the sidesets contained in "mesh"
   */
  virtual void write_sidesets(const MeshBase & mesh);

  /**
   * Writes the nodesets contained in "mesh"
   */
  virtual void write_nodesets(const MeshBase & mesh);

  /**
   * Sets up the nodal variables
   */
  void initialize_element_variables(const MeshBase & mesh, std::vector<std::string> names);

  /**
   * Sets up the nodal variables
   */
  void initialize_nodal_variables(std::vector<std::string> names);

  /**
   * Sets up the global variables
   */
  void initialize_global_variables(const std::vector<std::string> & names);

  /**
   * Writes the time for the timestep
   */
  void write_timestep(int timestep, Real time);

  /**
   * Writes the vector of values to the element variables.
   */
  void write_element_values(const MeshBase & mesh, const std::vector<Number> & values, int timestep);

  /**
   * Writes the vector of values to a nodal variable.
   */
  void write_nodal_values(int var_id, const std::vector<Number> & values, int timestep);

  /**
   * Writes the vector of information records.
   */
  void write_information_records(const std::vector<std::string> & records);

  /**
   * Writes the vector of global variables.
   */
  void write_global_values(const std::vector<Number> & values, int timestep);

  /**
   * Sets the underlying value of the boolean flag
   * _use_mesh_dimension_instead_of_spatial_dimension.  By default,
   * the value of this flag is false.
   *
   * See the ExodusII_IO class documentation for a detailed
   * description of this flag.
   */
  void use_mesh_dimension_instead_of_spatial_dimension(bool val);

  /**
   * Allows you to set a vector that is added to the coordinates of all
   * of the nodes.  Effectively, this "moves" the mesh to a particular position
   */
  void set_coordinate_offset(Point p);

  /**
   * This is the \p ExodusII_IO_Helper Conversion class.  It provides
   * a data structure which contains \p ExodusII node/edge maps and
   * name conversions.  It's defined below.
   */
  class Conversion;

  /**
   * This is the \p ExodusII_IO_Helper ElementMap class.
   * It contains constant maps between the \p ExodusII naming/numbering
   * schemes and the canonical schemes used in this code.  It's defined
   * below.
   */
  class ElementMaps;

  /**
   * This is the \p ExodusII_IO_Helper NamesData class.
   * It manages the C data structure necessary for writing out named
   * entities to ExodusII files.
   */
  class NamesData;

  /**
   * Prints the message defined in \p msg. Can be turned off if
   * verbosity is set to 0.
   */
  void message(const std::string msg);

  /**
   * Prints the message defined in \p msg, and appends the number \p i
   * to the end of the message.  Useful for printing messages in
   * loops.  Can be turned off if verbosity is set to 0.
   */
  void message(const std::string msg, int i);

  // File identification flag
  int ex_id;

  // General error flag
  int ex_err;

  // Number of dimensions in the mesh
  int num_dim;

  // Number of global variables
  int num_globals;

  // Total number of nodes in the mesh
  int num_nodes;

  // Total number of elements in the mesh
  int num_elem;

  // Total number of element blocks
  int num_elem_blk;

  // Total number of node sets
  int num_node_sets;

  // Total number of element sets
  int num_side_sets;

  // Number of elements in this block
  int num_elem_this_blk;

  // Number of nodes in each element
  int num_nodes_per_elem;

  // Number of attributes for a given block
  int num_attr;

  // Total number of elements in all side sets
  int num_elem_all_sidesets;

  // Vector of the block identification numbers
  std::vector<int> block_ids;

  // Vector of nodes in an element
  std::vector<int> connect;

  // Vector of the sideset IDs
  std::vector<int> ss_ids;

  // Vector of the nodeset IDs
  std::vector<int> nodeset_ids;

  // Number of sides (edges/faces) in current set
  std::vector<int> num_sides_per_set;

  // Number of nodes in current set
  std::vector<int> num_nodes_per_set;

  // Number of distribution factors per set
  std::vector<int> num_df_per_set;

  // Number of distribution factors per set
  std::vector<int> num_node_df_per_set;

  // List of element numbers in all sidesets
  std::vector<int> elem_list;

  // Side (face/edge) number actually on the boundary
  std::vector<int> side_list;

  // Node number actually on the boundary
  std::vector<int> node_list;

  // Side (face/edge) id number
  std::vector<int> id_list;

  // Optional mapping from internal [0,num_nodes) to arbitrary indices
  std::vector<int> node_num_map;

  // Optional mapping from internal [0,num_elem) to arbitrary indices
  std::vector<int> elem_num_map;

  // x locations of node points
  std::vector<Real> x;

  // y locations of node points
  std::vector<Real> y;

  // z locations of node points
  std::vector<Real> z;

  //  Problem title (Use vector<char> to emulate a char*)
  std::vector<char> title;

  // Type of element in a given block
  std::vector<char> elem_type;

  // Maps libMesh element numbers to Exodus element numbers
  // gets filled in when write_elements gets called
  std::map<int, int> libmesh_elem_num_to_exodus;
  std::vector<int> exodus_elem_num_to_libmesh;

  // Map of all node numbers connected to local node numbers to their exodus numbering.
  // The exodus numbers are stored in here starting with 1
  std::map<int, int> libmesh_node_num_to_exodus;
  std::vector<int> exodus_node_num_to_libmesh;

  // The number of timesteps in the file, as returned by ex_inquire
  int num_time_steps;

  // The timesteps stored in the solution file, filled by read_time_steps()
  std::vector<Real> time_steps;

  // The number of nodal variables in the Exodus file
  int num_nodal_vars;

  // The names of the nodal variables stored in the Exodus file
  std::vector<std::string> nodal_var_names;

  // Holds the nodal variable values for a given variable, one value per node
  std::vector<Real> nodal_var_values;

  // The number of elemental variables in the Exodus file
  int num_elem_vars;

  // The names of the elemental variables stored in the Exodus file
  std::vector<std::string> elem_var_names;

  // Holds the elemental variable values for a given variable, one value per element
  std::vector<Real> elem_var_values;

  // Maps of Ids to named entities
  std::map<int, std::string> id_to_block_names;
  std::map<int, std::string> id_to_ss_names;
  std::map<int, std::string> id_to_ns_names;

 protected:
  // This flag gets set after the create() function has been successfully called.
  // We call create() to open an ExodusII file for writing.
  bool _created;

  // This flag gets set after the open() function has been successfully called.
  // We call open() to open an ExodusII file for reading.
  bool _opened;

  // On/Off message flag
  bool _verbose;

  // If true, whenever there is an I/O operation, only perform if if we are on processor 0.
  bool _run_only_on_proc0;

  // True once the elem vars are initialized
  bool _elem_vars_initialized;

  // True once the global vars are initialized
  bool _global_vars_initialized;

  /**
   * If true, use the Mesh's dimension (as determined by the dimension
   * of the elements comprising the mesh) instead of the mesh's
   * spatial dimension, when writing.  By default this is false.
   */
  bool _use_mesh_dimension_instead_of_spatial_dimension;

  Point _coordinate_offset;
};








class ExodusII_IO_Helper::Conversion
{
public:

  /**
   * Constructor.  Initializes the const private member
   * variables.
   */
  Conversion(const int* nm,       // node_map
	     size_t nm_size,
	     const int* inm,      // inverse_node_map
	     size_t inm_size,
	     const int* sm,       // side_map
	     size_t sm_size,
	     const int* ism,      // inverse_side_map
	     size_t ism_size,
	     const ElemType ct,   // "canonical" aka libmesh element type
	     std::string ex_type) // string representing the Exodus element type
    : node_map(nm),
      node_map_size(nm_size),
      inverse_node_map(inm),
      inverse_node_map_size(inm_size),
      side_map(sm),
      side_map_size(sm_size),
      inverse_side_map(ism),
      inverse_side_map_size(ism_size),
      canonical_type(ct),
      exodus_type(ex_type)
  {}

  /**
   * Returns the ith component of the node map for this
   * element.  The node map maps the exodusII node numbering
   * format to this library's format.
   */
  int get_node_map(int i) const
  {
    libmesh_assert_less (static_cast<size_t>(i), node_map_size);
    return node_map[i];
  }

  /**
   * Returns the ith component of the inverse node map for this
   * element.  The inverse node map maps the libmesh node numbering
   * to Exodus' node numbering.  Note that all elements except Hex27
   * currently have the same node numbering as libmesh elements.
   */
  int get_inverse_node_map(int i) const
  {
    libmesh_assert_less (static_cast<size_t>(i), inverse_node_map_size);
    return inverse_node_map[i];
  }

  /**
   * Returns the ith component of the side map for this
   * element.  The side map maps the exodusII side numbering
   * format to this library's format.
   */
  int get_side_map(int i) const
  {
    libmesh_assert_less (static_cast<size_t>(i), side_map_size);
    return side_map[i];
  }

  /**
   * Returns the ith component of the side map for this
   * element.  The side map maps the libMesh side numbering
   * format to this exodus's format.
   */
  int get_inverse_side_map(int i) const
  {
    libmesh_assert_less (static_cast<size_t>(i), inverse_side_map_size);
    return inverse_side_map[i];
  }

  /**
   * Returns the canonical element type for this
   * element.  The canonical element type is the standard
   * element type understood by this library.
   */
  ElemType get_canonical_type()    const { return canonical_type; }

  /**
   * Returns the string corresponding to the Exodus type for this element
   */
  std::string exodus_elem_type() const { return exodus_type; }


private:
  /**
   * Pointer to the node map for this element.
   */
  const int* node_map;

  /**
   * The size of the node map array, this helps with bounds checking...
   */
  size_t node_map_size;

  /**
   * Pointer to the inverse node map for this element.
   * For all elements except for the Hex27, this is the same
   * as the node map.
   */
  const int* inverse_node_map;

  /**
   * The size of the inverse node map array, this helps with bounds checking...
   */
  size_t inverse_node_map_size;

  /**
   * Pointer to the side map for this element.
   */
  const int* side_map;

  /**
   * The size of the side map array, this helps with bounds checking...
   */
  size_t side_map_size;

  /**
   * Pointer to the inverse side map for this element.
   */
  const int* inverse_side_map;

  /**
   * The size of the inverse side map array, this helps with bounds checking...
   */
  size_t inverse_side_map_size;

  /**
   * The canonical (i.e. standard for this library)
   * element type.
   */
  const ElemType canonical_type;

  /**
   * The string corresponding to the Exodus type for this element
   */
  const std::string exodus_type;
};






class ExodusII_IO_Helper::ElementMaps
{
public:

  /**
   * Constructor.
   */
  ElementMaps() {}

  /**
   * 1D node maps.  These define mappings from ExodusII-formatted
   * element numberings.
   */

  /**
   * The Edge2 node map.  Use this map for linear elements in 1D.
   */
  static const int edge2_node_map[2];

  /**
   * The Edge3 node map.  Use this map for quadratic elements in 1D.
   */
  static const int edge3_node_map[3];

  /**
   * 1D edge maps
   */
  // FIXME: This notion may or may not be defined in ExodusII

  /**
   * Maps the Exodus edge numbering for line elements.  Useful for
   * reading sideset information.
   */
  static const int edge_edge_map[2];

  /**
   * Maps the Exodus edge numbering for line elements.
   * Useful for writing sideset information.
   */
  static const int edge_inverse_edge_map[2];

  /**
   * 2D node maps.  These define mappings from ExodusII-formatted
   * element numberings.
   */

  /**
   * The Quad4 node map.  Use this map for bi-linear quadrilateral
   * elements in 2D.
   */
  static const int quad4_node_map[4];

  /**
   * The Quad8 node map.  Use this map for serendipity quadrilateral
   * elements in 2D.
   */
  static const int quad8_node_map[8];

  /**
   * The Quad9 node map.  Use this map for bi-quadratic quadrilateral
   * elements in 2D.
   */
  static const int quad9_node_map[9];

  /**
   * The Tri3 node map.  Use this map for linear triangles in 2D.
   */
  static const int tri3_node_map[3];

  /**
   * The Tri6 node map.  Use this map for quadratic triangular
   * elements in 2D.
   */
  static const int tri6_node_map[6];

  /**
   * 2D edge maps
   */

  /**
   * Maps the Exodus edge numbering for triangles.  Useful for reading
   * sideset information.
   */
  static const int tri_edge_map[3];

  /**
   * Maps the Exodus edge numbering for quadrilaterals.  Useful for
   * reading sideset information.
   */
  static const int quad_edge_map[4];

  /**
   * Maps the Exodus edge numbering for triangles.  Useful for writing
   * sideset information.
   */
  static const int tri_inverse_edge_map[3];

  /**
   * Maps the Exodus edge numbering for quadrilaterals.  Useful for
   * writing sideset information.
   */
  static const int quad_inverse_edge_map[4];

  /**
   * 3D maps.  These define mappings from ExodusII-formatted element
   * numberings.
   */

  /**
   * The Hex8 node map.  Use this map for bi-linear hexahedral
   * elements in 3D.
   */
  static const int hex8_node_map[8];

  /**
   * The Hex20 node map.  Use this map for serendipity hexahedral
   * elements in 3D.
   */
  static const int hex20_node_map[20];

  /**
   * The Hex27 node map.  Use this map for reading tri-quadratic
   * hexahedral elements in 3D.
   */
  static const int hex27_node_map[27];

  /**
   * The Hex27 inverse node map.  Use this map for writing
   * tri-quadratic hexahedral elements in 3D.
   */
  static const int hex27_inverse_node_map[27];

  /**
   * The Tet4 node map.  Use this map for linear tetrahedral elements
   * in 3D.
   */
  static const int tet4_node_map[4];

  /**
   * The Tet10 node map.  Use this map for quadratic tetrahedral
   * elements in 3D.
   */
  static const int tet10_node_map[10];

  /**
   * The Prism6 node map.
   */
  static const int prism6_node_map[6];

  /**
   * The Prism15 node map.  Use this map for "serendipity" prisms in
   * 3D.
   */
  static const int prism15_node_map[15];

  /**
   * The Prism18 node map.
   */
  static const int prism18_node_map[18];

  /**
   * The Pyramid5 node map.  Use this map for linear pyramid elements
   * in 3D.
   */
  static const int pyramid5_node_map[5];


  /**
   * 3D face maps.
   */

  /**
   * Maps the Exodus face numbering for general hexahedrals.
   * Useful for reading sideset information.
   */
  static const int hex_face_map[6];

  /**
   * Maps the Exodus face numbering for 27-noded hexahedrals.
   * Useful for reading sideset information.
   */
  static const int hex27_face_map[6];

  /**
   * Maps the Exodus face numbering for general tetrahedrals.
   * Useful for reading sideset information.
   */
  static const int tet_face_map[4];

  /**
   * Maps the Exodus face numbering for general prisms.
   * Useful for reading sideset information.
   */
  static const int prism_face_map[5];

  /**
   * Maps the Exodus face numbering for general pyramids.
   * Useful for reading sideset information.
   */
  static const int pyramid_face_map[5];

  /**
   * Maps the Exodus face numbering for general hexahedrals.
   * Useful for writing sideset information.
   */
  static const int hex_inverse_face_map[6];

  /**
   * Maps the Exodus face numbering for 27-noded hexahedrals.
   * Useful for writing sideset information.
   */
  static const int hex27_inverse_face_map[6];

  /**
   * Maps the Exodus face numbering for general tetrahedrals.
   * Useful for writing sideset information.
   */
  static const int tet_inverse_face_map[4];

  /**
   * Maps the Exodus face numbering for general prisms.
   * Useful for writing sideset information.
   */
  static const int prism_inverse_face_map[5];

  /**
   * Maps the Exodus face numbering for general pyramids.
   * Useful for writing sideset information.
   */
  static const int pyramid_inverse_face_map[5];

  /**
   * @returns a conversion object given an element type name.
   */
  ExodusII_IO_Helper::Conversion assign_conversion(std::string type_str);

  /**
   * @returns a conversion object given an element type.
   */
  ExodusII_IO_Helper::Conversion assign_conversion(const ElemType type);
};



/**
 * This class is useful for managing anything that requires a char**
 * input/output in ExodusII file.  You must know the number of strings
 * and the length of each one at the time you create it.
 */
class ExodusII_IO_Helper::NamesData
{
public:
  /**
   * Constructor.  Allocates enough storage to hold n_strings of
   * length string_length.  (Actually allocates string_length+1 characters
   * per string to account for the trailing NULL character.)
   */
  explicit
  NamesData(size_t n_strings, size_t string_length);

  /**
   * Adds another name to the current data table.
   */
  void push_back_entry(const std::string & name);

  /**
   * Provide access to the underlying C data table
   */
  char** get_char_star_star();

  /**
   * Provide access to the i'th underlying char*
   */
  char* get_char_star(int i);

private:
  // C++ data structures for managing string memory
  std::vector<std::vector<char> > data_table;
  std::vector<char*> data_table_pointers;

  size_t counter;
  size_t table_size;
};


} // namespace libMesh

#endif // LIBMESH_HAVE_EXODUS_API

#endif // LIBMESH_EXODUSII_IO_HELPER_H
