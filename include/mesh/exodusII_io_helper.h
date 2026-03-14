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

#ifndef LIBMESH_EXODUSII_IO_HELPER_H
#define LIBMESH_EXODUSII_IO_HELPER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_EXODUS_API

// Local includes
#include "libmesh/parallel_object.h"
#include "libmesh/point.h"
#include "libmesh/boundary_info.h" // BoundaryInfo::BCTuple
#include "libmesh/enum_elem_type.h" // INVALID_ELEM
#include "libmesh/exodus_header_info.h"

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>

// Macros to simplify checking Exodus error codes
#define EX_CHECK_ERR(code, msg)                 \
  do {                                          \
    if ((code) < 0) {                           \
      libmesh_error_msg((msg));                 \
    } } while (0)

#define EX_EXCEPTIONLESS_CHECK_ERR(code, msg)   \
  do {                                          \
    if ((code) < 0) {                           \
      libMesh::err << (msg) << std::endl;       \
      libmesh_exceptionless_error();            \
    } } while (0)

// Before we include a header wrapped in a namespace, we'd better make
// sure none of its dependencies end up in that namespace
#include <errno.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>

namespace libMesh
{

// Forward declarations
class MeshBase;
class DofObject;

/**
 * This is the \p ExodusII_IO_Helper class.  This class hides the
 * implementation details of interfacing with the Exodus binary
 * format.
 *
 * \author John W. Peterson
 * \date 2002
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
  ExodusII_IO_Helper(const ParallelObject & parent,
                     bool v=false,
                     bool run_only_on_proc0=true,
                     bool single_precision=false);
  /**
   * Special functions. This class does not manage any dynamically
   * allocated resources (file pointers, etc.) so it _should_ be
   * default copy/move constructable, but I don't know
   * if any existing code actually uses these operations.
   */
  ExodusII_IO_Helper (const ExodusII_IO_Helper &) = default;
  ExodusII_IO_Helper (ExodusII_IO_Helper &&) = default;
  virtual ~ExodusII_IO_Helper();

  /**
   * This class contains references so it can't be default
   * copy/move-assigned.
   */
  ExodusII_IO_Helper & operator= (const ExodusII_IO_Helper &) = delete;
  ExodusII_IO_Helper & operator= (ExodusII_IO_Helper &&) = delete;

  /**
   * \returns The ExodusII API version, in "nodot" format; e.g. 822
   * for 8.22
   */
  static int get_exodus_version();

  /**
   * \returns The current element type.
   *
   * \note The default behavior is for this value to be in all capital
   * letters, e.g. \p HEX27.
   */
  const char * get_elem_type() const;

  /**
   * Sets whether or not to write extra "side" elements.  This is useful for
   * plotting SIDE_DISCONTINUOUS data.
   */
  void set_add_sides(bool add_sides);

  bool get_add_sides();

  /**
   * Opens an \p ExodusII mesh file named \p filename.  If
   * read_only==true, the file will be opened with the EX_READ flag,
   * otherwise it will be opened with the EX_WRITE flag.
   */
  void open(const char * filename, bool read_only);

  /**
   * Reads an \p ExodusII mesh file header, leaving this object's
   * internal data structures unchanged.
   */
  ExodusHeaderInfo read_header() const;

  /**
   * Reads an \p ExodusII mesh file header, and stores required
   * information on this object.
   */
  void read_and_store_header_info();

  /**
   * Reads the QA records from an ExodusII file.  We can use this to
   * detect when e.g. CUBIT 14+ was used to generate a Mesh file, and
   * work around certain known bugs in that version.
   */
  void read_qa_records();

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
   * Reads the optional \p bex_cv_blocks from the \p ExodusII mesh
   * file.
   */
  void read_bex_cv_blocks();

  /**
   * Prints the nodal information, by default to \p libMesh::out.
   */
  void print_nodes(std::ostream & out_stream = libMesh::out);

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
   * Read in edge blocks, storing information in the BoundaryInfo object.
   */
  void read_edge_blocks(MeshBase & mesh);

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
   * Reads information about all of the elemsets in the \p ExodusII
   * mesh file.
   */
  void read_elemset_info();

  /**
   * Reads information about sideset \p id and inserts it into the
   * global sideset array at the position \p offset.
   */
  void read_sideset(int id, int offset);

  /**
   * Reads information about elemset \p id and inserts it into the
   * global elemset array at the position \p offset.
   */
  void read_elemset(int id, int offset);

  /**
   * New API that reads all nodesets simultaneously. This may be slightly
   * faster than reading them one at a time. Calls ex_get_concat_node_sets()
   * under the hood.
   */
  void read_all_nodesets();

  /**
   * Closes the \p ExodusII mesh file.
   *
   * This function is called from the ExodusII_IO destructor, so it should
   * not throw an exception.
   */
  void close() noexcept;

  /**
   * Reads and stores the timesteps in the 'time_steps' array.
   */
  void read_time_steps();

  /**
   * Reads the number of timesteps currently stored in the Exodus file
   * and stores it in the num_time_steps variable.
   */
  void read_num_time_steps();

  /**
   * Reads the nodal values for the variable 'nodal_var_name' at the
   * specified time into the 'nodal_var_values' array.
   */
  void read_nodal_var_values(std::string nodal_var_name, int time_step);

  /**
   * Reads elemental values for the variable 'elemental_var_name' at the
   * specified timestep into the 'elem_var_value_map' which is passed in.
   */
  void read_elemental_var_values(std::string elemental_var_name,
                                 int time_step,
                                 std::map<dof_id_type, Real> & elem_var_value_map);

  /**
   * Reads element "extra integer" values, using the variable names
   * specified in \p extra_integer_var_names, and converting the
   * stored values to integer format.  The returned value is a vector
   * (indexed in the same way as extra_integer_vars) of maps from
   * element id to extra integer value.
   *
   * Currently values from the last time step are used.
   */
  std::vector<std::map<dof_id_type, dof_id_type>> read_extra_integers(
    const std::vector<std::string> & extra_integer_var_names);

  /**
   * Helper function that takes a (1-based) Exodus node/elem id and
   * determines the corresponding libMesh Node/Elem id. Takes into account
   * whether the user has chosen to set the Node/Elem unique ids based on
   * the {node,elem}_num_map or to let libMesh set them.
   */
  dof_id_type get_libmesh_node_id(int exodus_node_id);
  dof_id_type get_libmesh_elem_id(int exodus_elem_id);

  /**
   * Helper function that conditionally sets the unique_id of the
   * passed-in Node/Elem.  Calling this function does nothing if
   * _set_unique_ids_from_maps == false, otherwise it sets the
   * unique_id based on the entries of the {node,elem_num_map}.  The
   * input index is assumed to be a zero-based index into the
   * {node,elem}_num_map array.
   */
  void conditionally_set_node_unique_id(
    MeshBase & mesh, Node * node, int zero_based_node_num_map_index);
  void conditionally_set_elem_unique_id(
    MeshBase & mesh, Elem * elem, int zero_based_elem_num_map_index);

private:

  /**
   * Internal implementation for the two sets of functions above.
   */
  dof_id_type get_libmesh_id(
    int exodus_id,
    const std::vector<int> & num_map);

  void set_dof_object_unique_id(
    MeshBase & mesh,
    DofObject * dof_object,
    int exodus_mapped_id);

public:

  /**
   * Opens an \p ExodusII mesh file named \p filename for writing.
   */
  virtual void create(std::string filename);

  /**
   * Initializes the Exodus file.
   */
  virtual void initialize(std::string title, const MeshBase & mesh, bool use_discontinuous=false);

  /**
   * Writes the nodal coordinates contained in "mesh"
   */
  virtual void write_nodal_coordinates(const MeshBase & mesh, bool use_discontinuous=false);

  /**
   * Writes the elements contained in "mesh". FIXME: This only works
   * for Meshes having a single type of element in each subdomain!
   *
   * If \p use_discontinuous is true, we break apart elements, so that
   * shared nodes on faces/edges/vertices can take different values
   * from different elements.  This is useful for plotting
   * discontinuous underlying variables
   *
   * If \p _add_sides is true, we also output side elements, so that
   * shared nodes on edges/vertices can take different values from
   * different elements.  This is useful for plotting
   * SIDE_DISCONTINUOUS representing e.g. inter-element fluxes.
   */
  virtual void write_elements(const MeshBase & mesh,
                              bool use_discontinuous=false);

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
  virtual void initialize_element_variables(std::vector<std::string> names,
                                            const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains);

  /**
   * Sets up the nodal variables
   */
  void initialize_nodal_variables(std::vector<std::string> names);

  /**
   * Sets up the global variables
   */
  void initialize_global_variables(std::vector<std::string> names);

  /**
   * Writes the time for the timestep
   */
  void write_timestep(int timestep, Real time);

  /**
   * Write elemsets stored on the Mesh to the exo file.
   */
  void write_elemsets(const MeshBase & mesh);

  /**
   * Write sideset data for the requested timestep.
   */
  void
  write_sideset_data (const MeshBase & mesh,
                      int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<boundary_id_type>> & side_ids,
                      const std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals);

  /**
   * Read sideset variables, if any, into the provided data structures.
   */
  void
  read_sideset_data (const MeshBase & mesh,
                     int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<boundary_id_type>> & side_ids,
                     std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals);

  /**
   * Similar to read_sideset_data(), but instead of creating one
   * std::map per sideset per variable, creates a single map of (elem,
   * side, boundary_id) tuples, and stores the exo file array indexing
   * for any/all sideset variables on that sideset (they are all the
   * same). This function does not actually call exII::ex_get_sset_var()
   * to get values, and can be useful if only the original ordering of
   * (elem, side) pairs in the exo file is required in cases where a
   * separate approach is used to read the sideset data arrays.
   */
  void
  get_sideset_data_indices (const MeshBase & mesh,
                            std::map<BoundaryInfo::BCTuple, unsigned int> & bc_array_indices);

  /**
   * Write nodeset data for the requested timestep.
   */
  void
  write_nodeset_data (int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                      const std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals);

  /**
   * Read nodeset variables, if any, into the provided data structures.
   */
  void
  read_nodeset_data (int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                     std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals);

  /**
   * Similar to read_nodeset_data(), but instead of creating one
   * std::map per nodeset per variable, creates a single map of
   * (node_id, boundary_id) tuples, and stores the exo file array
   * indexing for any/all nodeset variables on that nodeset (they are
   * all the same).
   */
  void
  get_nodeset_data_indices (std::map<BoundaryInfo::NodeBCTuple, unsigned int> & bc_array_indices);

  /**
   * Write elemset data for the requested timestep.
   */
  void
  write_elemset_data (int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                      const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals);

  /**
   * Read elemset variables, if any, into the provided data structures.
   */
  void
  read_elemset_data (int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                     std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals);

  /**
   * Similar to read_elemset_data(), but instead of creating one
   * std::map per elemset per variable, creates a single map of
   * (elem_id, elemset_id) tuples, and stores the exo file array
   * indexing for any/all elemset variables on that elemset (they are
   * all the same).
   */
  void
  get_elemset_data_indices (std::map<std::pair<dof_id_type, elemset_id_type>, unsigned int> & elemset_array_indices);

  /**
   * Writes the vector of values to the element variables.
   *
   * The 'values' vector is assumed to be in the order:
   * {(u1, u2, u3, ..., uN), (v1, v2, v3, ..., vN), (w1, w2, w3, ..., wN)}
   * where N is the number of elements.
   *
   * This ordering is produced by calls to ES::build_elemental_solution_vector().
   * ES::build_discontinuous_solution_vector(), on the other hand, produces an
   * element-major ordering. See the function below for that case.
   */
  void write_element_values
  (const MeshBase & mesh,
   const std::vector<Real> & values,
   int timestep,
   const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains);

  /**
   * Same as the function above, but assume the input 'values' vector is
   * in element-major order, i.e.
   * {(u1,v1,w1), (u2,v2,w2), ... (uN,vN,wN)}
   * This function is called by
   * ExodusII_IO::write_element_data_from_discontinuous_nodal_data()
   * because ES::build_discontinuous_solution_vector() builds the solution
   * vector in this order.
   *
   * \note If some variables are subdomain-restricted, then the tuples will
   * be of different lengths for each element, i.e.
   * {(u1,v1,w1), (u2,v2), ... (uN,vN,wN)}
   * if variable w is not active on element 2.
   */
  void write_element_values_element_major
  (const MeshBase & mesh,
   const std::vector<Real> & values,
   int timestep,
   const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains,
   const std::vector<std::string> & derived_var_names,
   const std::map<subdomain_id_type, std::vector<std::string>> & subdomain_to_var_names);

  /**
   * Writes the vector of values to a nodal variable.
   */
  void write_nodal_values(int var_id, const std::vector<Real> & values, int timestep);

  /**
   * Writes the vector of information records.
   */
  void write_information_records(const std::vector<std::string> & records);

  /**
   * Writes the vector of global variables.
   */
  void write_global_values(const std::vector<Real> & values, int timestep);

  /**
   * Uses ex_update() to flush buffers to file.
   */
  void update();

  /**
   * Reads the vector of global variables.
   */
  void read_global_values(std::vector<Real> & values, int timestep);

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
   * Set to true (the default) to write files in an HDF5-based file
   * format (when HDF5 is available), or to false to write files in
   * the old NetCDF3-based format.  If HDF5 is unavailable, this
   * setting does nothing.
   */
  void set_hdf5_writing(bool write_hdf5);

  /**
   * Sets the value of _write_as_dimension.
   *
   * This directly controls the num_dim which is written to the Exodus
   * file.  If non-zero, this value supersedes all other dimensions,
   * including:
   * 1.) MeshBase::spatial_dimension()
   * 2.) MeshBase::mesh_dimension()
   * 3.) Any value passed to use_mesh_dimension_instead_of_spatial_dimension()
   * This is useful/necessary for working around a bug in Paraview which
   * prevents the "Plot Over Line" filter from working on 1D meshes.
   */
  void write_as_dimension(unsigned dim);

  /**
   * Allows you to set a vector that is added to the coordinates of all
   * of the nodes.  Effectively, this "moves" the mesh to a particular position
   */
  void set_coordinate_offset(Point p);

  /**
   * Set how many characters to use in names when opening a file for
   * writing.
   */
  void set_max_name_length(unsigned int max_length);

  /**
   * \returns A vector with three copies of each element in the provided name vector,
   * starting with r_, i_ and a_ respectively. If the "write_complex_abs" parameter
   * is true (default), the complex modulus is written, otherwise only the real and
   * imaginary parts are written.
   */
  std::vector<std::string>
  get_complex_names(const std::vector<std::string> & names,
                    bool write_complex_abs) const;

  /**
   * returns a "tripled" copy of \p vars_active_subdomains, which is necessary in the
   * complex-valued case.
   */
  std::vector<std::set<subdomain_id_type>>
  get_complex_vars_active_subdomains
  (const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains,
   bool write_complex_abs) const;

  /**
   * Takes a map from subdomain id -> vector of active variable names
   * as input and returns a corresponding map where the original
   * variable names have been replaced by their complex counterparts.
   * Used by the ExodusII_IO::write_element_data_from_discontinuous_nodal_data()
   * function.
   */
  std::map<subdomain_id_type, std::vector<std::string>>
  get_complex_subdomain_to_var_names
  (const std::map<subdomain_id_type, std::vector<std::string>> & subdomain_to_var_names,
   bool write_complex_abs) const;

  /**
   * This is the \p ExodusII_IO_Helper Conversion class.  It provides
   * a data structure which contains \p ExodusII node/edge maps and
   * name conversions.  It's defined below.
   */
  class Conversion;

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
  void message(std::string_view msg);

  /**
   * Prints the message defined in \p msg, and appends the number \p i
   * to the end of the message.  Useful for printing messages in
   * loops.  Can be turned off if verbosity is set to 0.
   */
  void message(std::string_view msg, int i);

  // File identification flag
  int ex_id;

  // General error flag
  int ex_err;

  // struct which contains data fields from the Exodus file header
  ExodusHeaderInfo header_info;

  // Problem title (Use vector<char> to emulate a char *)
  std::vector<char> & title;

  // Number of dimensions in the mesh
  int & num_dim;

  // Total number of nodes in the mesh
  int & num_nodes;

  // Total number of elements in the mesh
  int & num_elem;

  // Smallest element id which exceeds every element id in the mesh.
  // (this may exceed num_elem due to mapping)
  int end_elem_id() const;

  // Total number of element blocks
  int & num_elem_blk;

  // Total number of edges
  int & num_edge;

  // Total number of edge blocks. The sum of the number of edges in
  // each block must equal num_edge.
  int & num_edge_blk;

  // Total number of node sets
  int & num_node_sets;

  // Total number of side sets
  int & num_side_sets;

  // Total number of element sets
  int & num_elem_sets;

  // Number of global variables
  int num_global_vars;

  // Number of sideset variables
  int num_sideset_vars;

  // Number of nodeset variables
  int num_nodeset_vars;

  // Number of elemset variables
  int num_elemset_vars;

  // Number of elements in this block
  int num_elem_this_blk;

  // Number of nodes in each element
  int num_nodes_per_elem;

  // Number of attributes for a given block
  int num_attr;

  // Total number of elements in all side sets
  int num_elem_all_sidesets;

  // Total number of elements in all elem sets
  int num_elem_all_elemsets;

  // Vector of element block identification numbers
  std::vector<int> block_ids;

  // Vector of edge block identification numbers
  std::vector<int> edge_block_ids;

  // Vector of nodes in an element
  std::vector<int> connect;

  // Vector of the sideset IDs
  std::vector<int> ss_ids;

  // Vector of the nodeset IDs
  std::vector<int> nodeset_ids;

  // Vector of the elemset IDs
  std::vector<int> elemset_ids;

  // Number of sides in each sideset
  std::vector<int> num_sides_per_set;

  // Number of nodes in each nodeset
  std::vector<int> num_nodes_per_set;

  // Number of elems in each elemset
  std::vector<int> num_elems_per_set;

  // Number of distribution factors per sideset
  std::vector<int> num_df_per_set;

  // Number of distribution factors per nodeset
  std::vector<int> num_node_df_per_set;

  // Number of distribution factors per elemset
  std::vector<int> num_elem_df_per_set;

  // Starting indices for each nodeset in the node_sets_node_list vector.
  // Used in the calls to ex_{put,get}_concat_node_sets().
  std::vector<int> node_sets_node_index;

  // Starting indices for each nodeset in the node_sets_dist_fact vector.
  // Used in the calls to ex_{put,get}_concat_node_sets().
  std::vector<int> node_sets_dist_index;

  // Node ids for all nodes in nodesets, concatenated together.
  // Used in the calls to ex_{put,get}_concat_node_sets().
  std::vector<int> node_sets_node_list;

  // Distribution factors for all nodes in all nodesets, concatenated together.
  // Used in the calls to ex_{put,get}_concat_node_sets().
  std::vector<Real> node_sets_dist_fact;

  // List of element numbers in all sidesets
  std::vector<int> elem_list;

  // Side (face/edge) number actually on the boundary
  std::vector<int> side_list;

  // Side (face/edge) id number
  std::vector<int> id_list;

  // List of element numbers in all elemsets
  std::vector<int> elemset_list;

  // List of elemset ids for all elements in elemsets
  std::vector<int> elemset_id_list;

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

  // Spline weights associated with node points, in IGA meshes
  std::vector<Real> w;

  // Number of Bezier Extraction coefficient vectors in a block
  unsigned int bex_num_elem_cvs;

  // Bezier Extraction connectivity indices, in IGA meshes
  std::vector<std::vector<long unsigned int>> bex_cv_conn;

  // Bezier Extraction coefficient vectors, in IGA meshes
  // bex_dense_constraint_vecs[block_num][vec_num][column_num] = coef
  std::vector<std::vector<std::vector<Real>>> bex_dense_constraint_vecs;

  // Type of element in a given block
  std::vector<char> elem_type;

  // Maps libMesh element numbers to Exodus element numbers
  // gets filled in when write_elements gets called
  std::map<dof_id_type, dof_id_type> libmesh_elem_num_to_exodus;
  std::vector<int> exodus_elem_num_to_libmesh;

  // Map of all node numbers connected to local node numbers to their exodus numbering.
  // The exodus numbers are stored in here starting with 1
  std::map<dof_id_type, dof_id_type> libmesh_node_num_to_exodus;
  std::vector<int> exodus_node_num_to_libmesh;

  // The number of timesteps in the file, as returned by ex_inquire
  int num_time_steps;

  // The timesteps stored in the solution file, filled by read_time_steps()
  std::vector<Real> time_steps;

  // The number of nodal variables in the Exodus file
  int num_nodal_vars;

  // The names of the nodal variables stored in the Exodus file
  std::vector<std::string> nodal_var_names;

  // Holds the nodal variable values for a given variable, one value
  // per node, indexed by libMesh node id.
  // This is a map so it can handle Nemesis files as well.
  std::map<dof_id_type, Real> nodal_var_values;

  // The number of elemental variables in the Exodus file
  int num_elem_vars;

  // The names of the elemental variables stored in the Exodus file
  std::vector<std::string> elem_var_names;

  // Holds the elemental variable values for a given variable, one value per element
  std::vector<Real> elem_var_values;

  // The names of the global variables stored in the Exodus file
  std::vector<std::string> global_var_names;

  // The names of the sideset variables stored in the Exodus file
  std::vector<std::string> sideset_var_names;

  // The names of the nodeset variables stored in the Exodus file
  std::vector<std::string> nodeset_var_names;

  // The names of the elemset variables stored in the Exodus file
  std::vector<std::string> elemset_var_names;

  // Maps of Ids to named entities
  std::map<int, std::string> id_to_block_names;
  std::map<int, std::string> id_to_edge_block_names;
  std::map<int, std::string> id_to_ss_names;
  std::map<int, std::string> id_to_ns_names;
  std::map<int, std::string> id_to_elemset_names;

  // On/Off message flag
  bool verbose;

  // Same as the ExodusII_IO flag by the same name. This flag is
  // also set whenever ExodusII_IO::set_unique_ids_from_maps() is called.
  bool set_unique_ids_from_maps;

  // This flag gets set after the Exodus file has been successfully opened for writing.
  // Both the create() and open() (if called with EX_WRITE) functions may set this flag.
  bool opened_for_writing;

  // This flag gets set after the open() function has been successfully called.
  // We call open() to open an ExodusII file for reading.
  bool opened_for_reading;

  // When either create() or open() is called, the Helper stores the
  // name of the opened file as current_filename.  This way, the
  // ExodusII_IO object can check to see if, on subsequent writes, the
  // user is asking to write to a *different* filename from the one
  // that is currently open, and signal an error.  The current
  // ExodusII_IO implementation is designed to work with a single file
  // only, so if you want to write to multiple Exodus files, use a
  // different ExodusII_IO object for each one.
  std::string current_filename;

  /**
   * Wraps calls to exII::ex_get_var_names() and exII::ex_get_var_param().
   * The enumeration controls whether nodal, elemental, global, etc.
   * variable names are read and which class members are filled in.
   * NODAL:     num_nodal_vars   nodal_var_names
   * ELEMENTAL: num_elem_vars    elem_var_names
   * GLOBAL:    num_global_vars  global_var_names
   * SIDESET:   num_sideset_vars sideset_var_names
   * NODESET:   num_nodeset_vars nodeset_var_names
   */
  enum ExodusVarType {NODAL=0, ELEMENTAL=1, GLOBAL=2, SIDESET=3, NODESET=4, ELEMSET=5};
  void read_var_names(ExodusVarType type);

  const ExodusII_IO_Helper::Conversion &
  get_conversion(const ElemType type) const;

  const ExodusII_IO_Helper::Conversion &
  get_conversion(std::string type_str) const;

  /*
   * Returns the sum of node "offsets" that are to be expected from a
   * parallel nodal solution vector that has had "fake" nodes added on
   * each processor.  This plus a node id gives a valid nodal solution
   * vector index.
   */
  dof_id_type node_id_to_vec_id(dof_id_type n) const;

  /*
   * Returns the sum of both added node "offsets" on processors 0
   * through p-1 and real nodes added on processors 0 to p.
   * This is the starting index for added nodes' data.
   */
  dof_id_type added_node_offset_on(processor_id_type p) const;

protected:
  /**
   * Calculate _added_side_node_offsets needed to add "fake" side
   * elements to the given mesh
   */
  void calculate_added_side_node_offsets(const MeshBase & mesh);

  /**
   * When appending: during initialization, check that variable names
   * in the file match those you attempt to initialize with.
   */
  void check_existing_vars(ExodusVarType type, std::vector<std::string> & names, std::vector<std::string> & names_from_file);

  /**
   * Wraps calls to exII::ex_put_var_names() and exII::ex_put_var_param().
   * The enumeration controls whether nodal, elemental, or global
   * variable names are read and which class members are filled in.
   */
  void write_var_names(ExodusVarType type, const std::vector<std::string> & names);

  // If true, whenever there is an I/O operation, only perform if if we are on processor 0.
  bool _run_only_on_proc0;

  // This flag gets set after the create() function has been successfully called.
  bool _opened_by_create;

  // True once the elem vars are initialized
  bool _elem_vars_initialized;

  // True once the global vars are initialized
  bool _global_vars_initialized;

  // True once the nodal vars are initialized
  bool _nodal_vars_initialized;

  // If true, use the Mesh's dimension (as determined by the dimension
  // of the elements comprising the mesh) instead of the mesh's
  // spatial dimension, when writing.  By default this is false.
  bool _use_mesh_dimension_instead_of_spatial_dimension;

  // If true, write an HDF5 file when available.  If false, write the
  // old format.
  bool _write_hdf5;

  // The maximum name length to use when writing
  unsigned int _max_name_length;

  // Set once the elem num map has been read
  int _end_elem_id;

  // Use this for num_dim when writing the Exodus file.  If non-zero, supersedes
  // any value set in _use_mesh_dimension_instead_of_spatial_dimension.
  unsigned _write_as_dimension;

  // On output, shift every point by _coordinate_offset
  Point _coordinate_offset;

  // If true, forces single precision I/O
  bool _single_precision;

  /**
   * If we're adding "fake" sides to visualize SIDE_DISCONTINUOUS
   * variables, _added_side_node_offsets[p] gives us the total
   * solution vector offset to use on processor p+1 from the nodes on
   * those previous ranks' sides.
   */
  std::vector<dof_id_type> _added_side_node_offsets;

  /**
   * If we're adding "fake" sides to visualize SIDE_DISCONTINUOUS
   * variables, we also need to know how many real nodes from previous
   * ranks are taking up space in a solution vector.
   */
  std::vector<dof_id_type> _true_node_offsets;

  /**
   * This class facilitates inline conversion of an input data vector
   * to a different precision level, depending on the underlying type
   * of Real and whether or not the single_precision flag is set. This
   * should be used whenever floating point data is being written to
   * the Exodus file. Note that if no precision conversion has to take
   * place, there should be very little overhead involved in using
   * this object.
   */
  struct MappedOutputVector
  {
    // If necessary, allocates space to store a version of vec_in in a
    // different precision than it was input with.
    MappedOutputVector(const std::vector<Real> & vec_in,
                       bool single_precision_in);

    ~MappedOutputVector() = default;

    // Returns void * pointer to either the mapped data or the
    // original data, as necessary.
    void * data();

  private:
    const std::vector<Real> & our_data;
    bool single_precision;
    std::vector<double> double_vec;
    std::vector<float> float_vec;
  };

  /**
   * This class facilitates reading in vectors from Exodus file that
   * may be of a different floating point type than Real. It employs
   * basically the same approach as the MappedOutputVector, just going
   * in the opposite direction. For more information, see the
   * MappedOutputVector class docs.
   */
  struct MappedInputVector
  {
    MappedInputVector(std::vector<Real> & vec_in,
                      bool single_precision_in);
    ~MappedInputVector();

    // Returns void * pointer to either the mapped data or the
    // original data, as necessary.
    void * data();

  private:
    std::vector<Real> & our_data;
    bool single_precision;
    std::vector<double> double_vec;
    std::vector<float> float_vec;
  };


protected:

  /**
   * read_var_names() dispatches to this function.  We need to
   * override it slightly for Nemesis.
   */
  virtual void read_var_names_impl(const char * var_type,
                                   int & count,
                                   std::vector<std::string> & result);

  /**
   * Set to true iff we want to write separate "side" elements too.
   */
  bool _add_sides = false;

  /**
   * Map of subdomains to element numbers.
   */
  std::map<subdomain_id_type, std::vector<dof_id_type>> _subdomain_map;

  /**
   * One beyond the last "real" subdomain id, in cases where we have
   * "fake" subdomain ids for added sides
   */
  subdomain_id_type _subdomain_id_end;

  /**
   * Method for constructing _subdomain_map
   */
  void build_subdomain_map(const MeshBase & mesh, bool local);

private:
  /**
   * write_var_names() dispatches to this function.
   */
  void write_var_names_impl(const char * var_type,
                            int & count,
                            const std::vector<std::string> & names);

  /**
   * Defines equivalence classes of Exodus element types that map to
   * libmesh ElemTypes.
   */
  std::map<std::string, ElemType> element_equivalence_map;
  void init_element_equivalence_map();

  /**
   * Associates libMesh ElemTypes with node/face/edge/etc. mappings
   * of the corresponding Exodus element types.
   *
   * We have to map based on both ElemType and mesh dimension, because
   * Exodus treats "TRI" side numbering in two different ways
   * depending on whether a triangle is embedded in a 2D or a 3D mesh.
   */
  std::map<int, std::map<ElemType, ExodusII_IO_Helper::Conversion>> conversion_map;
  void init_conversion_map();
};



class ExodusII_IO_Helper::Conversion
{
public:

  /**
   * Constructor. Zero initializes all variables.
   */
  Conversion()
    : node_map(nullptr),
      inverse_node_map(nullptr),
      side_map(nullptr),
      inverse_side_map(nullptr),
      shellface_map(nullptr),
      inverse_shellface_map(nullptr),
      shellface_index_offset(0),
      libmesh_type(INVALID_ELEM),
      dim(0),
      n_nodes(0),
      exodus_type("")
  {}

  /**
   * \returns The ith component of the node map for this element.
   *
   * The node map maps the exodusII node numbering format to this
   * library's format.
   */
  int get_node_map(int i) const;

  /**
   * \returns The ith component of the inverse node map for this
   * element.
   *
   * The inverse node map maps the libmesh node numbering to Exodus'
   * node numbering.
   *
   * \note All elements except Hex27 currently have the same node
   * numbering as libmesh elements.
   */
  int get_inverse_node_map(int i) const;

  /**
   * \returns The ith component of the side map for this element.
   *
   * The side map maps the exodusII side numbering format to this
   * library's format.
   */
  int get_side_map(int i) const;

  /**
   * \returns The ith component of the side map for this element.
   *
   * The side map maps the libMesh side numbering format to this
   * exodus's format.
   */
  int get_inverse_side_map(int i) const;

  /**
   * \returns The ith component of the shellface map for this element.
   * \note Nothing is currently using this.
   */
  int get_shellface_map(int i) const;

  /**
   * \returns The ith component of the inverse shellface map for this element.
   */
  int get_inverse_shellface_map(int i) const;

  /**
   * \returns The canonical element type for this element.
   *
   * The canonical element type is the standard element type
   * understood by this library.
   */
  ElemType libmesh_elem_type() const;

  /**
   * \returns The string corresponding to the Exodus type for this element.
   */
  std::string exodus_elem_type() const;

  /**
   * \returns The shellface index offset.
   */
  std::size_t get_shellface_index_offset() const;

  /**
   * An invalid_id that can be returned to signal failure in case
   * something goes wrong.
   */
  static const int invalid_id;

  /**
   * Pointer to the node map for this element.
   */
  const std::vector<int> * node_map;

  /**
   * Pointer to the inverse node map for this element.
   * For all elements except for the Hex27, this is the same
   * as the node map.
   */
  const std::vector<int> * inverse_node_map;

  /**
   * Pointer to the side map for this element.
   */
  const std::vector<int> * side_map;

  /**
   * Pointer to the inverse side map for this element.
   */
  const std::vector<int> * inverse_side_map;

  /**
   * Pointer to the shellface map for this element. Only the inverse
   * is actually used currently, this one is provided for completeness
   * and libmesh_ingore()d to avoid warnings.
   */
  const std::vector<int> * shellface_map;

  /**
   * Pointer to the inverse shellface map for this element.
   */
  const std::vector<int> * inverse_shellface_map;

  /**
   * The shellface index offset defines the offset due to a difference between libMesh
   * and Exodus in indexing sidesets. This is relevant for shell elements, for
   * example, since Exodus includes extra "shell face" sides in that case.
   */
  size_t shellface_index_offset;

  /**
   * The canonical (i.e. standard for this library)
   * element type.
   */
  ElemType libmesh_type;

  /**
   * The element dimension; useful since we don't seem to have a cheap
   * way to look this up from ElemType
   */
  int dim;

  /**
   * The number of nodes per element; useful likewise
   */
  int n_nodes;

  /**
   * The string corresponding to the Exodus type for this element
   */
  std::string exodus_type;
};



/**
 * This class is useful for managing anything that requires a char **
 * input/output in ExodusII file.  You must know the number of strings
 * and the length of each one at the time you create it.
 */
class ExodusII_IO_Helper::NamesData
{
public:
  /**
   * Constructor.  Allocates enough storage to hold n_strings of
   * length string_length.  (Actually allocates string_length+1 characters
   * per string to account for the trailing '\0' character.)
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
  char ** get_char_star_star();

  /**
   * Provide access to the i'th underlying char *
   */
  char * get_char_star(int i);

private:
  // C++ data structures for managing string memory
  std::vector<std::vector<char>> data_table;
  std::vector<char *> data_table_pointers;

  size_t counter;
  size_t table_size;
};


inline void ExodusII_IO_Helper::set_add_sides(bool add_sides)
{
  _add_sides = add_sides;
}


inline bool ExodusII_IO_Helper::get_add_sides()
{
  return _add_sides;
}


inline int ExodusII_IO_Helper::end_elem_id() const
{
  libmesh_assert_equal_to(std::size_t(num_elem), elem_num_map.size());
  return _end_elem_id;
}


} // namespace libMesh

#endif // LIBMESH_HAVE_EXODUS_API

#endif // LIBMESH_EXODUSII_IO_HELPER_H
