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

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>

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
		     bool v=false, bool run_only_on_proc0=true) :
    ParallelObject(parent),
    comp_ws(sizeof(Real)),
    io_ws(0),
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
    req_info(0),
    ret_int(0),
    num_elem_all_sidesets(0),
    ex_version(0.0),
    ret_float(0.0),
    ret_char(0),
    num_time_steps(0),
    _created(false),
    _verbose(v),
    _run_only_on_proc0(run_only_on_proc0),
    _elem_vars_initialized(false),
    _global_vars_initialized(false),
    _use_mesh_dimension_instead_of_spatial_dimension(false)
  {
    title.resize(MAX_LINE_LENGTH+1);
    elem_type.resize(MAX_STR_LENGTH);
  }

  /**
   * Destructor.  The only memory
   * allocated is for \p title and
   * \p elem_type.  This memory
   * is freed in the destructor.
   */
  virtual ~ExodusII_IO_Helper();

  /**
   * Returns true once create() has been successfully called, and
   * false otherwise.
   */
  bool created();

  /**
   * Get/set flag telling whether message printing is on or off.
   */
  void verbose (bool set_verbosity);

  /**
   * @returns the \p ExodusII
   * mesh dimension.
   */
  int get_num_dim()                const { return num_dim; }

  /**
   * @returns the total number of
   * global variables.
   */
  int get_num_globals()            const { return num_globals; }

  /**
   * @returns the total number of
   * nodes in the \p ExodusII mesh.
   */
  int get_num_nodes()              const { return num_nodes; }


  /**
   * @returns the total number of
   * elements in the \p ExodusII mesh.
   */
  int get_num_elem()               const { return num_elem; }

  /**
   * @returns the total number
   * of element blocks in
   * the \p ExodusII mesh.
   */
  int get_num_elem_blk()           const { return num_elem_blk; }

  /**
   * For a given block,
   * returns the total number
   * of elements.
   */
  int get_num_elem_this_blk()      const { return num_elem_this_blk; }

  /**
   * @returns the number of
   * nodes per element in
   * a given block. e.g.
   * for HEX27 it returns 27.
   */
  int get_num_nodes_per_elem()     const { return num_nodes_per_elem; }

  /**
   * @returns the total number
   * of sidesets in the \p ExodusII
   * mesh.  Each sideset contains
   * only one type of element.
   */
  int get_num_side_sets()          const { return num_side_sets; }

  /**
   * @returns the total number
   * of nodesets in the \p ExodusII
   * mesh.
   */
  int get_num_node_sets()          const { return num_node_sets; }

  //     /**
  //      * @returns the number of
  //      * elements in all the sidesets.
  //      * Effectively returns the
  //      * total number of elements
  //      * on the \p ExodusII mesh boundary.
  //      */
  //     int get_num_elem_all_sidesets()  const { return num_elem_all_sidesets; }

  /**
   * @returns the \f$ i^{th} \f$
   * node number in the
   * element connectivity
   * list for a given element.
   */
  int get_connect(int i)           const { return connect[i]; }

  /**
   * For a single sideset,
   * returns the total number of
   * elements in the sideset.
   */
  int get_num_sides_per_set(int i) const { return num_sides_per_set[i]; }

  /**
   * For a single nodeset,
   * returns the total number of
   * nodes in the nodeset.
   */
  int get_num_nodes_per_set(int i) const { return num_nodes_per_set[i]; }

  //     /**
  //      * @returns the \f$ i^{th} \f$ entry
  //      * in the element list.
  //      * The element list contains
  //      * the numbers of all elements
  //      * on the boundary.
  //      */
  //     int get_elem_list(int i)         const { return elem_list[i]; }

  /**
   * @return a constant reference to the \p elem_list.
   */
  const std::vector<int>& get_elem_list() const { return elem_list; }

  //     /**
  //      * @returns the \f$ i^{th} \f$ entry in
  //      * the side list.  This is
  //      * effectively the "side"
  //      * (face in 3D or edge in
  //      * 2D) number which lies
  //      * on the boundary.
  //      */
  //     int get_side_list(int i)         const { return side_list[i]; }

  /**
   * @return a constant reference to the \p side_list.
   */
  const std::vector<int>& get_side_list() const { return side_list; }

  /**
   * @return a constant reference to the \p node_list.
   */
  const std::vector<int>& get_node_list() const { return node_list; }

  /**
   * @return the nodeset id corresponding to the ith nodeset.
   */
  int get_nodeset_id(unsigned int i) const { return nodeset_ids[i]; }

  //     /**
  //      * @returns the \f$ i^{th} \f$ entry in
  //      * the id list.  This is the id
  //      * for the ith face on the boundary.
  //      */
  //     int get_id_list(int i)         const { return id_list[i]; }

  /**
   * @return a constant reference to the \p id_list.
   */
  const std::vector<int>& get_id_list() const { return id_list; }

  /**
   * @returns the current
   * element type.  Note:
   * the default behavior
   * is for this value
   * to be in all capital
   * letters, e.g. \p HEX27.
   */
  const char* get_elem_type()            const { return &elem_type[0]; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's x-coordinate.
   */
  Real get_x(int i) const { return x[i]; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's y-coordinate.
   */
  Real get_y(int i) const { return y[i]; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's z-coordinate.
   */
  Real get_z(int i) const { return z[i]; }

  /**
   * Opens an \p ExodusII mesh
   * file named \p filename
   * for reading.
   */
  void open(const char* filename);

  /**
   * Reads an \p ExodusII mesh
   * file header.
   */
  void read_header();

  /**
   * Prints the \p ExodusII
   * mesh file header,
   * which includes the
   * mesh title, the number
   * of nodes, number of
   * elements, mesh dimension,
   * number of sidesets, and
   * number of nodesets
   */
  void print_header();

  /**
   * Reads the nodal data
   * (x,y,z coordinates)
   * from the \p ExodusII mesh
   * file.
   */
  void read_nodes();

  /**
   * Reads the optional \p node_num_map
   * from the \p ExodusII mesh file.
   */
  void read_node_num_map();

  /**
   * Prints the nodal information,
   * by default to \p libMesh::out.
   */
  void print_nodes(std::ostream &out = libMesh::out);

  /**
   * Reads information for
   * all of the blocks in
   * the \p ExodusII mesh file.
   */
  void read_block_info();

  /**
   * Get the block number
   * for the given block index.
   */
  int get_block_id(int index);

  /**
   * Get the block name for the given block index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_block_name(int index);

  /**
   * Get the side set id
   * for the given side set index.
   */
  int get_side_set_id(int index);

  /**
   * Get the side set name for the given side set index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_side_set_name(int index);

  /**
   * Get the node set id
   * for the given node set index.
   */
  int get_node_set_id(int index);

  /**
   * Get the node set name for the given node set index if supplied in
   * the mesh file.  Otherwise an empty string is returned.
   */
  std::string get_node_set_name(int index);

  /**
   * Reads all of the element
   * connectivity for
   * block \p block in the
   * \p ExodusII mesh file.
   */
  void read_elem_in_block(int block);

  /**
   * Reads the optional \p node_num_map
   * from the \p ExodusII mesh file.
   */
  void read_elem_num_map();

  /**
   * Reads information about
   * all of the sidesets in
   * the \p ExodusII mesh file.
   */
  void read_sideset_info();

  /**
   * Reads information about
   * all of the nodesets in
   * the \p ExodusII mesh file.
   */
  void read_nodeset_info();

  /**
   * Reads information about
   * sideset \p id and
   * inserts it into the global
   * sideset array at the
   * position \p offset.
   */
  void read_sideset(int id, int offset);

  /**
   * Reads information about
   * nodeset \p id and
   * inserts it into the global
   * nodeset array at the
   * position \p offset.
   */
  void read_nodeset(int id);

  /**
   * Prints information
   * about all the sidesets.
   */
  void print_sideset_info();

  /**
   * Prints information
   * about all the nodesets.
   */
  void print_nodeset_info();

  /**
   * Closes the \p ExodusII
   * mesh file.
   */
  void close();

  /**
   * Generic inquiry, returns the value
   */
  int inquire(int req_info, std::string error_msg="");


  // For reading solutions:
  /*
   * Returns an array containing the timesteps in the file
   */
  const std::vector<Real>& get_time_steps();


  /*
   * Number of Nodal variables defined.
   */
  int get_num_nodal_vars(){ return num_nodal_vars; }


  /*
   * Returns an array containing the nodal var names in the file
   */
  const std::vector<std::string>& get_nodal_var_names();

  /*
   * Returns an array containing the nodal variable values
   * at the specified time
   */
  const std::vector<Real>& get_nodal_var_values(std::string nodal_var_name, int time_step);

  // For Writing Solutions
  /**
   * Opens an \p ExodusII mesh
   * file named \p filename
   * for writing.
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
   * Writes the elements contained in "mesh"
   * FIXME: This only works for Mesh's having a single type of element!
   */
  virtual void write_elements(const MeshBase & mesh);

  /**
   * Writes the elements contained in "mesh"
   * FIXME: This only works for Mesh's having a single type of element!
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
   * This is the \p ExodusII_IO_Helper Conversion class.
   * It provides a data structure which contains \p ExodusII node/edge
   * maps and name conversions.  It's defined below.
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
   * It manages the C datastructure necessary for writing out named
   * entities to ExodusII files.
   */
  class NamesData;

  //private:


  /**
   * All of the \p ExodusII
   * API functions return
   * an \p int error value.
   * This function checks
   * to see if the error has
   * been set, and if it has,
   * prints the error message
   * contained in \p msg.
   */
  void check_err(const int error, const std::string msg);

  /**
   * Prints the message defined
   * in \p msg. Can be turned off if
   * verbosity is set to 0.
   */
  void message(const std::string msg);

  /**
   * Prints the message defined
   * in \p msg, and appends the number
   * \p i to the end of the
   * message.  Useful for
   * printing messages in loops.
   * Can be turned off if
   * verbosity is set to 0.
   */
  void message(const std::string msg, int i);

  int   comp_ws;                       // ?
  int   io_ws;                         // ?
  int   ex_id;                         // File identification flag
  int   ex_err;                        // General error flag
  int   num_dim;                       // Number of dimensions in the mesh
  int   num_globals;                   // Number of global variables
  int   num_nodes;                     // Total number of nodes in the mesh
  int   num_elem;                      // Total number of elements in the mesh
  int   num_elem_blk;                  // Total number of element blocks
  int   num_node_sets;                 // Total number of node sets
  int   num_side_sets;                 // Total number of element sets
  int   num_elem_this_blk;             // Number of elements in this block
  int   num_nodes_per_elem;            // Number of nodes in each element
  int   num_attr;                      // Number of attributes for a given block
  int   req_info;                      // Generic required info tag
  int   ret_int;                       // Generic int returned by ex_inquire
  int   num_elem_all_sidesets;         // Total number of elements in all side sets
  std::vector<int> block_ids;          // Vector of the block identification numbers
  std::vector<int> connect;            // Vector of nodes in an element
  std::vector<int> ss_ids;             // Vector of the sideset IDs
  std::vector<int> nodeset_ids;        // Vector of the nodeset IDs
  std::vector<int> num_sides_per_set;  // Number of sides (edges/faces) in current set
  std::vector<int> num_nodes_per_set;  // Number of nodes in current set
  std::vector<int> num_df_per_set;     // Number of distribution factors per set
  std::vector<int> num_node_df_per_set;// Number of distribution factors per set
  std::vector<int> elem_list;          // List of element numbers in all sidesets
  std::vector<int> side_list;          // Side (face/edge) number actually on the boundary
  std::vector<int> node_list;          // Node number actually on the boundary
  std::vector<int> id_list;            // Side (face/edge) id number
  std::vector<int> node_num_map;       // Optional mapping from internal [0,num_nodes) to arbitrary indices
  std::vector<int> elem_num_map;       // Optional mapping from internal [0,num_elem) to arbitrary indices
  float ex_version;                    // Version of Exodus you are using
  float ret_float;                     // Generic float returned by ex_inquire
  std::vector<Real> x;                 // x locations of node points
  std::vector<Real> y;                 // y locations of node points
  std::vector<Real> z;                 // z locations of node points
  char    ret_char;                    // Generic char returned by ex_inquire
  // Use vectors of char to emulate char*'s
  std::vector<char> title;             //  Problem title
  std::vector<char> elem_type;         // Type of element in a given block

  // Maps libMesh element numbers to Exodus element numbers
  // gets filled in when write_elements gets called
  std::map<int, int> libmesh_elem_num_to_exodus;
  std::vector<int> exodus_elem_num_to_libmesh;

  /**
   * Map of all node numbers connected to local node numbers to their exodus numbering.
   *
   * The exodus numbers are stored in here starting with 1
   */
  std::map<int, int> libmesh_node_num_to_exodus;
  std::vector<int> exodus_node_num_to_libmesh;

  //Solution Data
  int num_time_steps;
  std::vector<Real> time_steps;
  int num_nodal_vars;
  std::vector<std::string> nodal_var_names;
  std::vector<Real> nodal_var_values;

  int num_elem_vars;

  // A pair of containers used to emulate a char** data
  // structure without having to worry about dynamic memory
  // allocation ourselves.
  std::vector<std::vector<char> > vvc;
  std::vector<char*> strings; // vector of pointers into vvc

  /**
   * Maps of Ids to named entities
   */
  std::map<int, std::string> id_to_block_names;
  std::map<int, std::string> id_to_ss_names;
  std::map<int, std::string> id_to_ns_names;

 protected:
  bool _created; // This flag gets set after the the create() function has been successfully called.
  bool _verbose; // On/Off message flag
  bool _run_only_on_proc0; // If true, whenever there is an I/O operation, only perform if if we are on processor 0.
  bool _elem_vars_initialized; // True once the elem vars are initialized
  bool _global_vars_initialized; // True once the global vars are initialized

  /**
   * If true, use the Mesh's dimension (as determined by the dimension
   * of the elements comprising the mesh) instead of the mesh's
   * spatial dimension, when writing.  By default this is false.
   */
  bool _use_mesh_dimension_instead_of_spatial_dimension;
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
   * 1D node maps.  These define
   * mappings from ExodusII-formatted
   * element numberings.
   */

  /**
   * The Edge2 node map.
   * Use this map for linear elements in 1D.
   */
  static const int edge2_node_map[2];

  /**
   * The Edge3 node map.
   * Use this map for quadratic elements in 1D.
   */
  static const int edge3_node_map[3];

  /**
   * 1D edge maps
   */
  // FIXME: This notion may or may not be defined in ExodusII

  /**
   * Maps the Exodus edge numbering for line elements.
   * Useful for reading sideset information.
   */
  static const int edge_edge_map[2];

  /**
   * Maps the Exodus edge numbering for line elements.
   * Useful for writing sideset information.
   */
  static const int edge_inverse_edge_map[2];

  /**
   * 2D node maps.  These define
   * mappings from ExodusII-formatted
   * element numberings.
   */

  /**
   * The Quad4 node map.
   * Use this map for bi-linear
   * quadrilateral elements in 2D.
   */
  static const int quad4_node_map[4];

  /**
   * The Quad8 node map.
   * Use this map for serendipity
   * quadrilateral elements in 2D.
   */
  static const int quad8_node_map[8];

  /**
   * The Quad9 node map.
   * Use this map for bi-quadratic
   * quadrilateral elements in 2D.
   */
  static const int quad9_node_map[9];

  /**
   * The Tri3 node map.
   * Use this map for linear
   * triangles in 2D.
   */
  static const int tri3_node_map[3];

  /**
   * The Tri6 node map.
   * Use this map for quadratic
   * triangular elements in 2D.
   */
  static const int tri6_node_map[6];

  /**
   * 2D edge maps
   */

  /**
   * Maps the Exodus edge numbering for triangles.
   * Useful for reading sideset information.
   */
  static const int tri_edge_map[3];

  /**
   * Maps the Exodus edge numbering for quadrilaterals.
   * Useful for reading sideset information.
   */
  static const int quad_edge_map[4];

  /**
   * Maps the Exodus edge numbering for triangles.
   * Useful for writing sideset information.
   */
  static const int tri_inverse_edge_map[3];

  /**
   * Maps the Exodus edge numbering for quadrilaterals.
   * Useful for writing sideset information.
   */
  static const int quad_inverse_edge_map[4];

  /**
   * 3D maps.  These define
   * mappings from ExodusII-formatted
   * element numberings.
   */

  /**
   * The Hex8 node map.
   * Use this map for bi-linear
   * hexahedral elements in 3D.
   */
  static const int hex8_node_map[8];

  /**
   * The Hex20 node map.
   * Use this map for serendipity
   * hexahedral elements in 3D.
   */
  static const int hex20_node_map[20];

  /**
   * The Hex27 node map.
   * Use this map for reading tri-quadratic
   * hexahedral elements in 3D.
   */
  static const int hex27_node_map[27];

  /**
   * The Hex27 inverse node map.
   * Use this map for writing tri-quadratic
   * hexahedral elements in 3D.
   */
  static const int hex27_inverse_node_map[27];

  /**
   * The Tet4 node map.
   * Use this map for linear
   * tetrahedral elements in 3D.
   */
  static const int tet4_node_map[4];

  /**
   * The Tet10 node map.
   * Use this map for quadratic
   * tetrahedral elements in 3D.
   */
  static const int tet10_node_map[10];

  /**
   * The Prism6 node map.
   */
  static const int prism6_node_map[6];

  /**
   * The Prism15 node map.
   * Use this map for "serendipity" prisms in 3D.
   */
  static const int prism15_node_map[15];

  /**
   * The Prism18 node map.
   */
  static const int prism18_node_map[18];

  /**
   * The Pyramid5 node map.
   * Use this map for linear
   * pyramid elements in 3D.
   */
  static const int pyramid5_node_map[5];


  /**
   * 3D face maps.  Are these ever used for anything?
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
 * This class is useful for managing the names for named entities in an ExodusII
 * fortmat
 */
class ExodusII_IO_Helper::NamesData
{
public:
  explicit
  NamesData(size_t size);
  ~NamesData();

  /**
   * Adds another name to the current data table
   */
  void push_back_entry(const std::string & name);

  /**
   * Writes the datastructure to the specified ExodusII file
   */
  int write_to_exodus(int ex_id, exII::ex_entity_type type);

private:
  char **data_table;
  size_t counter;
  size_t table_size;
};


} // namespace libMesh

#endif // LIBMESH_HAVE_EXODUS_API

#endif // LIBMESH_EXODUSII_IO_HELPER_H
