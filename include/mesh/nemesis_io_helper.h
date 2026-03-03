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

#ifndef LIBMESH_NEMESIS_IO_HELPER_H
#define LIBMESH_NEMESIS_IO_HELPER_H

#include "libmesh/libmesh_config.h"

#if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)

// Local headers
#include "libmesh/exodusII_io_helper.h"

// C++ headers
#include <set>
#include <vector>

namespace libMesh
{

// Forward declarations
class EquationSystems;
template <typename T> class NumericVector;


/**
 * This is the \p Nemesis_IO_Helper class.  Think of it as
 * a big struct with storage for all the stuff one might
 * want to pull from a Nemesis file.  Derived from
 * ExodusII_IO_Helper object, since Nemesis is based on
 * the same file format.
 *
 * \author John W. Peterson
 * \author Roy Stogner
 * \date 2008
 * \date 2020
 */
class Nemesis_IO_Helper : public ExodusII_IO_Helper
{
public:
  /**
   * Constructor.
   */
  explicit
  Nemesis_IO_Helper(const ParallelObject & parent,
                    bool verbose=false, bool single_precision=false);

  /**
   * Destructor.
   */
  virtual ~Nemesis_IO_Helper();

  /**
   * Set the flag indicating whether the complex modulus should be
   * written when complex numbers are enabled. By default this flag
   * is set to true.
   */
  void write_complex_magnitude (bool val);

  /**
   * Reading functions.  These just allocate memory for you and call the Nemesis
   * routines of the same name.  They also handle error checking for the Nemesis
   * return value.  Be careful calling these at random, some depend on others
   * being called first...
   */

  /**
   * Reads the node ids of nodeset \p id and stores them in the \p
   * node_list member of this class.
   *
   * \note This used to be an ExodusII_IO_Helper function but it use
   * was completely replaced by read_all_nodesets(). For now, it is
   * still used by the Nemesis reader so we have moved it here.
   */
  void read_nodeset(int id);

  /**
   * Fills: num_nodes_global, num_elems_global, num_elem_blks_global,
   * num_node_sets_global, num_side_sets_global
   * Call after: read_and_store_header_info()
   * Call before: Any other get_* function from this class
   */
  void get_init_global();

  /**
   * Fills: global_sideset_ids, num_global_side_counts, num_global_side_df_counts
   * Call after: get_init_global()
   */
  void get_ss_param_global();
  void get_ns_param_global();
  void get_eb_info_global();
  void get_init_info();
  void get_loadbal_param();
  void get_elem_map();
  void get_node_map();
  void get_cmap_params();
  void get_node_cmap();
  void get_elem_cmap();

  /**
   * Writing functions.
   */

  /**
   * Writes basic info about the partitioning to file
   * .) num_proc - number of processors
   * .) num_proc_in_file - number of processors in the current file - generally equal to 1
   * .) ftype = "s" for scalar load-balance file, "p" for parallel file
   */
  void put_init_info(unsigned num_proc,
                     unsigned num_proc_in_file,
                     const char * ftype);

  /**
   * Writes global information including:
   * .) global number of nodes
   * .) global number of elems
   * .) global number of element blocks
   * .) global number of node sets
   * .) global number of side sets
   */
  void put_init_global(dof_id_type num_nodes_global,
                       dof_id_type num_elems_global,
                       unsigned num_elem_blks_global,
                       unsigned num_node_sets_global,
                       unsigned num_side_sets_global);

  /**
   * Writes global block information to the file
   * .) global_elem_blk_ids - list of block IDs for all blocks present in the mesh
   * .) global_elem_blk_cnts - number of elements in each block for the global mesh
   *
   * Must be called after put_init_global().
   */
  void put_eb_info_global(std::vector<int> & global_elem_blk_ids,
                          std::vector<int> & global_elem_blk_cnts);

  /**
   * This function writes information about global node sets.
   * .) global_nodeset_ids - vector of global node set IDs
   * .) num_global_node_counts - vector of global node counts contained in each global node set
   * .) num_global_df_count - vector of global distribution factors in each global node set
   *
   * Must be called after put_init_global()
   */
  void put_ns_param_global(std::vector<int> & global_nodeset_ids,
                           std::vector<int> & num_global_node_counts,
                           std::vector<int> & num_global_node_df_counts);

  /**
   * This function writes information about global side sets.
   * .) global_sideset_ids - vector of global side set IDs
   * .) num_global_side_counts - vector of global side counts contained in each global side set
   * .) num_global_df_count - vector of global distribution factors in each global side set
   *
   * Must be called after put_init_global()
   */
  void put_ss_param_global(std::vector<int> & global_sideset_ids,
                           std::vector<int> & num_global_side_counts,
                           std::vector<int> & num_global_side_df_counts);



  /**
   * Writes load balance parameters, some of which are described below:
   * .) num_internal_nodes - nodes "wholly" owned by the current processor
   * .) num_border_nodes - nodes local to a processor but residing in an element
   *                       which also has nodes on other processors
   * .) num_external_nodes - nodes that reside on other processors but whose element
   *                         "partially" resides on the current processor --
   *                          we assert this should be zero on reading!
   * .) num_border_elems - elements local to this processor but whose nodes reside
   *                       on other processors as well.
   * .) processor - ID of the processor for which information is to be written
   */
  void put_loadbal_param(unsigned num_internal_nodes,
                         unsigned num_border_nodes,
                         unsigned num_external_nodes,
                         unsigned num_internal_elems,
                         unsigned num_border_elems,
                         unsigned num_node_cmaps,
                         unsigned num_elem_cmaps);

  /**
   * Outputs initial information for communication maps.
   *
   * \note The order of the arguments specified in the Nemesis User's
   * Manual is \e wrong.  The correct order is (ids, counts, ids,
   * counts).  Must be called after put_loadbal_param().
   */
  void put_cmap_params(std::vector<int> & node_cmap_ids,
                       std::vector<int> & node_cmap_node_cnts,
                       std::vector<int> & elem_cmap_ids,
                       std::vector<int> & elem_cmap_elem_cnts);

  /**
   * Outputs *all* of the nodal communication maps for this processor.  Internally,
   * this function loops over all communication maps and calls
   * Nemesis::ne_put_node_cmap() for each one.
   *
   * .) node_cmap_node_ids = Nodal IDs of the FEM nodes in this communication map
   * .) node_cmap_proc_ids = processor IDs associated with each of the nodes in node_ids
   *
   * In the Nemesis file, these all appear to be written to the same chunks of data:
   * n_comm_nids and n_comm_proc, but don't rely on these names...
   *
   * \note This class contains \p node_cmap_node_ids and \p
   * node_cmap_proc_ids which can be used when calling this function.
   *
   * Must be called after put_cmap_params().
   */
  void put_node_cmap(std::vector<std::vector<int>> & node_cmap_node_ids,
                     std::vector<std::vector<int>> & node_cmap_proc_ids);

  /**
   * Outputs IDs of internal, border, and external nodes.
   * LibMesh asserts that the number of external nodes is zero in the
   * Nemesis files it reads
   */
  void put_node_map(std::vector<int> & node_mapi,
                    std::vector<int> & node_mapb,
                    std::vector<int> & node_mape);

  /**
   * Writes information about elemental communication map.
   *
   * \note This class contains \p elem_cmap_elem_ids, \p
   * elem_cmap_side_ids, abd \p elem_cmap_proc_ids which can be used
   * when calling this function.
   *
   * Must be called after put_cmap_params().
   */
  void put_elem_cmap(std::vector<std::vector<int>> & elem_cmap_elem_ids,
                     std::vector<std::vector<int>> & elem_cmap_side_ids,
                     std::vector<std::vector<int>> & elem_cmap_proc_ids);

  /**
   * Outputs IDs of internal and border elements.
   *
   * Must be called after ne_put_loadbal_param().
   */
  void put_elem_map(std::vector<int> & elem_mapi,
                    std::vector<int> & elem_mapb);

  /**
   * This function is specialized from ExodusII_IO_Helper to write only the
   * nodal coordinates stored on the local piece of the Mesh.
   */
  virtual void write_nodal_coordinates(const MeshBase & mesh, bool use_discontinuous=false) override;

  /**
   * This function is specialized to write the connectivity.
   */
  virtual void write_elements(const MeshBase & mesh, bool use_discontinuous=false) override;

  /**
   * Writes the sidesets for this processor.
   */
  virtual void write_sidesets(const MeshBase & mesh) override;

  /**
   * Writes the nodesets for this processor.
   */
  virtual void write_nodesets(const MeshBase & mesh) override;

  /**
   * Specialization of the initialize function from ExodusII_IO_Helper that
   * also writes global initial data to file.
   */
  virtual void initialize(std::string title, const MeshBase & mesh, bool use_discontinuous=false) override;

  /**
   * This function uses global communication routines to determine the
   * number of element blocks across the entire mesh.
   */
  void compute_num_global_elem_blocks(const MeshBase & pmesh);

  /**
   * This function builds the libmesh -> exodus and exodus -> libmesh
   * node and element maps.  These maps allow us to have a consistent
   * numbering scheme within an Exodus file, given an existing globally
   * consistent numbering scheme from LibMesh.
   */
  void build_element_and_node_maps(const MeshBase & pmesh);

  /**
   * Takes a parallel solution vector containing the node-major
   * solution vector for all variables and outputs it to the files.
   * \param parallel_soln
   * \param names A vector containing the names of _all_ variables in parallel_soln.
   * \param timestep To be passed to the ExodusII_IO_Helper::write_nodal_values() function.
   * \param output_names A vector containing the names of variables in parallel_soln that should actually be written (whitelist).
   *
   * \note This version of write_nodal_solution() is called by the
   * parallel version of Nemesis_IO::write_nodal_data(), which is
   * called by MeshOutput::write_equation_systems() for parallel I/O
   * formats like Nemesis.  The other version is still available to
   * continue supporting things like NamebasedIO::write_nodal_data(),
   * but this version should be preferred when running in parallel.
   */
  void write_nodal_solution(const NumericVector<Number> & parallel_soln,
                            const std::vector<std::string> & names,
                            int timestep,
                            const std::vector<std::string> & output_names);

  /**
   * Outputs EquationSystems current_local_solution nodal values.
   */
  void write_nodal_solution(const EquationSystems & es,
                            const std::vector<std::pair<unsigned int, unsigned int>> & var_nums,
                            int timestep,
                            const std::vector<std::string> & output_names);

  /**
   * Takes a solution vector containing the solution for all variables and outputs it to the files
   */
  void write_nodal_solution(const std::vector<Number> & values,
                            const std::vector<std::string> & names,
                            int timestep);

  /**
   * Override the Exodus Helper's implementation of this function so
   * that it works correctly in parallel.
   */
  virtual
  void initialize_element_variables(std::vector<std::string> names,
                                    const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains) override;
  /**
   * Writes the vector of elemental variable values, one variable and
   * one subdomain at a time.
   */
  void write_element_values(const MeshBase & mesh,
                            const EquationSystems & es,
                            const std::vector<std::pair<unsigned int, unsigned int>> &var_nums,
                            int timestep,
                            const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains);

  /**
   * Given base_filename, foo.e, constructs the Nemesis filename
   * foo.e.X.Y, where X=n. CPUs and Y=processor ID
   */
  std::string construct_nemesis_filename(std::string_view base_filename);

  /**
   * Member data
   */

  /**
   * All (?) Nemesis functions return an int.  If it's negative that signals an error!
   * Internally, we use the ExodusII_IO_Helper::check_err() function to check for errors.
   */
  int nemesis_err_flag;

  /**
   * Global initial information.  The names are self-explanatory
   * for the most part.  Used with Nemesis::ne_get_init_global().
   */
  int num_nodes_global;
  int num_elems_global;
  int num_elem_blks_global;
  int num_node_sets_global;
  int num_side_sets_global;

  /**
   * The number of processors for which the NEMESIS I file was created.
   * To be used with Nemesis::ne_get_init_info().
   */
  int num_proc;

  /**
   * The number of processors for which the NEMESIS I file stores information.
   * This is generally equal to 1 (1 CPU/file) at least for the splitting Derek gave us.
   * To be used with Nemesis::ne_get_init_info().
   */
  int num_proc_in_file;

  /**
   * The type of file to be written. Either 's', for a scalar
   * load-balance file, or 'p' for a parallel file.
   * To be used with Nemesis::ne_get_init_info().
   */
  char ftype;

  // Stores node ids read in by the read_nodeset() function
  std::vector<int> node_list;

  /**
   * Containers for reading global sideset (boundary conditions) information.  Each vector will
   * eventually have num_side_sets_global entries, and be used in calls to
   * Nemesis::ne_get_ss_param_global().
   *
   * It's an error to call ne_get_ss_param_global when num_side_sets_global==0
   */
  std::vector<int> global_sideset_ids;
  std::vector<int> num_global_side_counts;
  std::vector<int> num_global_side_df_counts;


  /**
   * Containers for reading global nodeset information.  One vector entry per nodeset.
   * Each vector will eventually have num_node_sets_global entries, and
   * will be used in calls to Nemesis::ne_get_ns_param_global().
   *
   * It's an error to call ne_get_ns_param_global when num_node_sets_global==0
   */
  std::vector<int> global_nodeset_ids;
  std::vector<int> num_global_node_counts;
  std::vector<int> num_global_node_df_counts;


  /**
   * Read the global element block IDs and counts.  These vectors will
   * eventually have num_elem_blks_global entries.  To be used with
   * Nemesis::ne_get_eb_info_global().
   */
  std::vector<int> global_elem_blk_ids;
  std::vector<int> global_elem_blk_cnts;

  /**
   * libMesh numbered node ids attached to local elems.
   */
  std::set<int> nodes_attached_to_local_elems;

  /**
   * This is the block connectivity, i.e. for each subdomain (block) there
   * is an element connectivity list. This map associates the block ID to that vector.
   */
  std::map<int, std::vector<int>> block_id_to_elem_connectivity;

  /**
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */

  /**
   * The number of FEM nodes contained in FEM elements wholly owned by the current processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_internal_nodes;

  /**
   * The number of FEM nodes local to a processor but residing in an
   * element which also has FEM nodes on other processors.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_border_nodes;

  /**
   * The number of FEM nodes that reside on another processor but
   * whose element partially resides on the current processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_external_nodes;

  /**
   * The number of internal FEM elements. Elements local to this processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_internal_elems;

  /**
   * The number of border FEM elements. Elements local to this
   * processor but whose FEM nodes reside on other processors as well.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_border_elems;

  /**
   * The number of nodal communication maps for this processor. (One
   * per neighboring proc?)
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_node_cmaps;

  /**
   * The number of elemental communication maps for this
   * processor. (One per neighboring proc?)
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_elem_cmaps;



  /**
   * Vector which stores internal element IDs.  Will have length
   * num_internal_elems.
   * To be used with Nemesis::ne_get_elem_map().
   */
  std::vector<int> elem_mapi;

  /**
   * Vector which stores border element IDs.  Will have length
   * num_border_elems.
   * To be used with Nemesis::ne_get_elem_map().
   */
  std::vector<int> elem_mapb;



  /**
   * Vector which stores internal node IDs.  Will have length
   * num_internal_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mapi;

  /**
   * Vector which stores border node IDs.  Will have length
   * num_border_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mapb;

  /**
   * Vector which stores external node IDs.  Will have length
   * num_external_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mape;


  /**
   * Vectors for storing the communication map parameters.
   * Each will eventually have length num_node_cmaps OR
   * num_elem_cmaps as appropriate.
   * For use with Nemesis::ne_get_cmap_params().
   */
  std::vector<int> node_cmap_ids;
  std::vector<int> node_cmap_node_cnts;
  std::vector<int> elem_cmap_ids;
  std::vector<int> elem_cmap_elem_cnts;


  /**
   * 2 vectors of vectors for storing the node communication IDs for this processor.
   * There will be num_node_cmaps rows, row i will have node_cmap_node_cnts[i] entries.
   * To be used with Nemesis::ne_get_node_cmap().
   *
   * Remark: node_cmap_proc_ids is a vector, all entries of which are = node_cmap_ids[i]
   * Not sure what the point of that is...
   */
  std::vector<std::vector<int>> node_cmap_node_ids;
  std::vector<std::vector<int>> node_cmap_proc_ids;


  /**
   * 3 vectors of vectors for storing element communication IDs for this processor.
   * There will be num_elem_cmaps rows, row i will have elem_cmap_elem_cnts[i] entries.
   * To be used with Nemesis::ne_get_elem_cmap().
   */
  std::vector<std::vector<int>> elem_cmap_elem_ids;
  std::vector<std::vector<int>> elem_cmap_side_ids;
  std::vector<std::vector<int>> elem_cmap_proc_ids;

  /**
   * By default, when complex numbers are enabled, for each variable
   * we write out three values: the real part, "r_u" the imaginary
   * part, "i_u", and the complex modulus, a_u := sqrt(r_u*r_u +
   * i_u*i_u), which is also the value returned by
   * std::abs(std::complex).  Since the modulus is not an independent
   * quantity, we can set this flag to false and save some file space
   * by not writing out.
   */
  bool write_complex_abs;

protected:
  /**
   * read_var_names() dispatches to this function.  We need to
   * override it slightly for Nemesis.
   */
  virtual void read_var_names_impl(const char * var_type,
                                   int & count,
                                   std::vector<std::string> & result) override;

private:
  /**
   * This map keeps track of the number of elements in each subdomain
   * (block) for *this* processor.
   */
  std::map<subdomain_id_type, unsigned> local_subdomain_counts;

  /**
   * The set which will eventually contain the IDs of "border nodes".  These are nodes
   * that lie on the boundary between one or more processors.
   */
  std::set<unsigned> border_node_ids;

  /**
   * Another map to store sets of intersections with each other processor
   * (other than ourself, of course).  A node which appears in one of these
   * vectors belongs to element owned by at least this processor and one other.
   */
  std::map<unsigned, std::set<unsigned>> proc_nodes_touched_intersections;

  /**
   * Typedef for an iterator into the data structure above.
   */
  typedef std::map<unsigned, std::set<unsigned>>::iterator proc_nodes_touched_iterator;

  /**
   * Map between processor ID and (element,side) pairs bordering that processor ID.
   */
  std::map<unsigned, std::set<std::pair<unsigned,unsigned>>> proc_border_elem_sets;

  /**
   * Typedef for an iterator into the data structure above.
   */
  typedef std::map<unsigned, std::set<std::pair<unsigned,unsigned>>>::iterator proc_border_elem_sets_iterator;

  /**
   * A set of internal node IDs for this processor.
   */
  std::set<unsigned> internal_node_ids;

  /**
   * A set of internal elem IDs for this processor.
   */
  std::set<unsigned> internal_elem_ids;

  /**
   * A set of border elem IDs for this processor.
   */
  std::set<unsigned> border_elem_ids;

  /**
   * This function uses global communication routines to determine the
   * number of nodesets across the entire mesh.
   */
  void compute_num_global_nodesets(const MeshBase & pmesh);

  /**
   * This function uses global communication routines to determine the
   * number of sidesets across the entire mesh.
   */
  void compute_num_global_sidesets(const MeshBase & pmesh);

  /**
   * This function constructs the set of border node IDs present
   * on the current mesh.  These are nodes which live on the "border"
   * between elements which live on different processors.
   */
  void compute_border_node_ids(const MeshBase & pmesh);

  /**
   * This function constructs the set of border and internal element IDs
   * and internal node IDs present on the current mesh.
   */
  void compute_internal_and_border_elems_and_internal_nodes(const MeshBase & pmesh);

  /**
   * This function determines the communication map parameters
   * which will eventually be written to file
   */
  void compute_communication_map_parameters();

  /**
   * Compute the node communication maps (really just pack vectors)
   * in preparation for writing them to file.
   */
  void compute_node_communication_maps();

  /**
   * Compute the node maps (really just pack vectors) which
   * map the nodes to internal, border, and external nodes in
   * the file.
   */
  void compute_node_maps();

  /**
   * This function computes element communication maps (really
   * just packs vectors) in preparation for writing them to file.
   */
  void compute_elem_communication_maps();

  /**
   * This function computes element maps (really just packs vectors)
   * which map the elements to internal and border elements.
   */
  void compute_element_maps();

  /**
   * This function writes exodus-specific initialization information.
   * This information is slightly different when you are working with
   * Nemesis, as it depends on some global information being known.
   */
  void write_exodus_initialization_info(const MeshBase & pmesh,
                                        const std::string & title);
};

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)
#endif // LIBMESH_NEMESIS_IO_HELPER_H
