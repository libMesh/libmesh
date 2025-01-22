// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EXODUSII_IO_H
#define LIBMESH_EXODUSII_IO_H


// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"
#include "libmesh/parallel_object.h"
#include "libmesh/boundary_info.h" // BoundaryInfo::BCTuple
#include "libmesh/exodus_header_info.h"

namespace libMesh
{

// Forward declarations
class EquationSystems;
class ExodusII_IO_Helper;
class MeshBase;
class System;

/**
 * The \p ExodusII_IO class implements reading meshes in the
 * \p ExodusII file format from Sandia National Labs.  By
 * default, LibMesh expects ExodusII files to have a ".exd"
 * or ".e" file extension.
 *
 * \author Benjamin Kirk
 * \author John Peterson
 * \date 2004
 * \brief Handles reading and writing of Exodus binary files.
 */
class ExodusII_IO : public MeshInput<MeshBase>,
                    public MeshOutput<MeshBase>,
                    public ParallelObject
{
public:

  /**
   * Constructor.  Takes a writable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  explicit
  ExodusII_IO (MeshBase & mesh,
               bool single_precision=false);

  /**
   * ExodusII_IO special functions:
   * - Can't be (default) copy constructed or assigned since they
   *   contain a unique_ptr member.
   * - Can't be (default) move assigned because the helper class
       contains references.
   * - The destructor is responsible for closing the file
   */
  ExodusII_IO (ExodusII_IO &&) = default;
  ExodusII_IO (const ExodusII_IO &) = delete;
  ExodusII_IO & operator= (const ExodusII_IO &) = delete;
  ExodusII_IO & operator= (ExodusII_IO &&) = delete;
  virtual ~ExodusII_IO ();

  /**
   * \returns The ExodusII API version, in "nodot" format (e.g. 822
   * for 8.22), or 0 if ExodusII support is not available.
   */
  static int get_exodus_version();

  /**
   * This method implements reading a mesh from a specified file.
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.
   */
  virtual void read (const std::string & name) override;

  /**
   * Read only the header information, instead of the entire
   * mesh. After the header is read, the file is closed and the
   * MeshBase object remains unchanged. This capability is useful if
   * you only need to know the mesh "metadata" and should be faster
   * than reading in all the nodes and elems.
   *
   * The header information currently includes:
   * .) Title string
   * .) Mesh spatial dimension
   * .) Number of nodes
   * .) Number of elements
   * .) Number of element blocks
   * .) Number of node sets
   * .) Number of side sets
   * .) Number of edge blocks/edges
   */
  ExodusHeaderInfo read_header (const std::string & name);

  /**
   * This method implements writing a mesh to a specified file.
   *
   * Note that writes may be buffered for efficiency, and so may not
   * reach disk until after the file has been closed, which happens
   * when the \p ExodusII_IO object is destructed.
   */
  virtual void write (const std::string & fname) override;

  /**
   * Set the flag indicating if we should be verbose.
   */
  void verbose (bool set_verbosity);

  /**
   * Set the flag indicating whether the complex modulus should be
   * written when complex numbers are enabled. By default this flag
   * is set to true.
   */
  void write_complex_magnitude (bool val);

  /**
   * By default, we only write out the elements physically stored in
   * the mesh.  If we have any SIDE_DISCONTINUOUS variables, however,
   * we cannot easily output them on elements with sides that share
   * vertices (and in 3D, edges).  We can set this flag to instead
   * create extra "side elements" on which to visualize such
   * variables.
   *
   * By default this flag is set to false.
   */
  void write_added_sides (bool val);

  virtual bool get_add_sides () override;

  /**
   * \returns An array containing the timesteps in the file.
   */
  const std::vector<Real> & get_time_steps();

  /**
   * \returns The number of timesteps currently stored in the Exodus
   * file.
   *
   * Knowing the number of time steps currently stored in the file is
   * sometimes necessary when appending, so we can know where to start
   * writing new data.  Throws an error if the file is not currently
   * open for reading or writing.
   */
  int get_num_time_steps();

  /**
   * If we read in a nodal solution while reading in a mesh, we can attempt
   * to copy that nodal solution into an EquationSystems object.
   */
  void copy_nodal_solution(System & system,
                           std::string system_var_name,
                           std::string exodus_var_name,
                           unsigned int timestep=1);

  /**
   * If we read in a elemental solution while reading in a mesh, we can attempt
   * to copy that elemental solution into an EquationSystems object.
   */
  void copy_elemental_solution(System & system,
                               std::string system_var_name,
                               std::string exodus_var_name,
                               unsigned int timestep=1);

  /**
   * Copy global variables into scalar variables of a System object.
   */
  void copy_scalar_solution(System & system,
                            std::vector<std::string> system_var_names,
                            std::vector<std::string> exodus_var_names,
                            unsigned int timestep=1);

  /**
   * Given an elemental variable and a time step, returns a mapping from the
   * elements (top parent) unique IDs to the value of the elemental variable at
   * the corresponding time step index.
   * Note that this function MUST only be called before renumbering!
   * This function is essentially a wrapper for read_elemental_var_values from
   * the exodus helper (which is not accessible outside this class).
   * \param elemental_var_name Name of an elemental variable
   * \param timestep The corresponding time step index
   * \param unique_id_to_value_map The map to be filled
   */
  void read_elemental_variable(std::string elemental_var_name,
                               unsigned int timestep,
                               std::map<unsigned int, Real> & unique_id_to_value_map);

  /**
   * Given a vector of global variables and a time step, returns the values
   * of the global variable at the corresponding time step index.
   * \param global_var_names Vector of names of global variables
   * \param timestep The corresponding time step index
   * \param global_values The vector to be filled
   */
  void read_global_variable(std::vector<std::string> global_var_names,
                            unsigned int timestep,
                            std::vector<Real> & global_values);

  /**
   * Writes a exodusII file with discontinuous data
   */
  void write_discontinuous_exodusII (const std::string & name,
                                     const EquationSystems & es,
                                     const std::set<std::string> * system_names=nullptr);

  /**
   * Writes a discontinuous solution at a specific timestep
   * \param fname Name of the file to be written
   * \param es EquationSystems object which contains the solution vector
   * \param timestep The timestep to write out. (should be _1_ indexed)
   * \param time The current simulation time
   * \param system_names Optional list of systems to write solutions for.
   */
  void write_timestep_discontinuous (const std::string &fname,
                                     const EquationSystems &es,
                                     const int timestep,
                                     const Real time,
                                     const std::set<std::string> * system_names=nullptr);

  /**
   * Write out element solution.
   */
  void write_element_data (const EquationSystems & es);

  /**
   * Similar to the function above, but instead of only handling
   * (CONSTANT, MONOMIAL) data, writes out a general discontinuous
   * solution field, e.g. (FIRST, L2_LAGRANGE) or (SECOND, MONOMIAL)
   * as a number of elemental fields equal to the number of vertices in
   * each element. For example, if you have a (FIRST, L2_LAGRANGE)
   * variable "u" defined on HEX8 elements, calling this function
   * would by default write 8 elemental fields named u_elem_node_0,
   * u_elem_node_1, u_elem_node_2, etc.
   *
   * This may be useful if you have a viz tool which is capable of
   * interpreting this element data as a discontinuous solution field.
   * Note that (CONSTANT, MONOMIAL) data is still written as a single
   * value per element, as it makes no sense to write n_vertices
   * copies of the same value.
   *
   * The 'var_suffix' parameter, which defaults to "_elem_node_", is
   * used to generate the elemental variable names, and is inserted
   * between the base variable name and the node id which the variable
   * applies to, e.g. "u_elem_node_0", "u_elem_node_1", etc.
   */
  void
  write_element_data_from_discontinuous_nodal_data
  (const EquationSystems & es,
   const std::set<std::string> * system_names = nullptr,
   const std::string & var_suffix = "_elem_node_");

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;
  using MeshOutput<MeshBase>::write_nodal_data_discontinuous;

  /**
   * Write out a nodal solution.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) override;

  /**
   * Write out a discontinuous nodal solution.
   */
  void write_nodal_data_discontinuous (const std::string &,
                                       const std::vector<Number> &,
                                       const std::vector<std::string> &) override;

  /**
   * Write out global variables.
   */
  void write_global_data (const std::vector<Number> &,
                          const std::vector<std::string> &);

  /**
   * Write out information records.
   */
  void write_information_records (const std::vector<std::string> &);

  /**
   * Writes out the solution at a specific timestep.
   * \param fname Name of the file to write to
   * \param es EquationSystems object which contains the solution vector.
   * \param timestep The timestep to write out, should be _1_ indexed.
   * \param time The current simulation time.
   * \param system_names Optional list of systems to write solutions for.
   */
  void write_timestep (const std::string & fname,
                       const EquationSystems & es,
                       const int timestep,
                       const Real time,
                       const std::set<std::string> * system_names=nullptr);

  /**
   * Writes out the solution for no specific time or timestep.
   * \param fname Name of the file to write to
   * \param es EquationSystems object which contains the solution vector.
   * \param system_names Optional list of systems to write solutions for.
   */
  virtual void write_equation_systems (const std::string & fname,
                                       const EquationSystems & es,
                                       const std::set<std::string> * system_names=nullptr) override;

  /**
   * Write elemsets stored on the Mesh to file.
   *
   * \note An elemset is a concept which is related to (but distinct
   * from) both elem blocks and sidesets/nodesets. Elements in an elemset
   * can be from multiple different blocks. In addition, one can
   * define an elemset variable, which is like an elemental variable
   * but exists only on the elements defined in the set.
   */
  void write_elemsets();

  /**
   * The Exodus format can also store values on sidesets. This can be
   * thought of as an alternative to defining an elemental variable
   * field on lower-dimensional elements making up a part of the
   * boundary. The inputs to the function are:
   * .) var_names[i] is the name of the ith sideset variable to be written to file.
   * .) side_ids[i] is a set of side_ids where var_names[i] is active.
   * .) bc_vals[i] is a map from (elem,side,id) BCTuple objects to the
   *    corresponding real-valued data.
   *
   * \note You must have already written the mesh by calling
   * e.g. write() before calling this function, because it uses the
   * existing ordering of the Exodus sidesets.
   */
  void
  write_sideset_data (int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<boundary_id_type>> & side_ids,
                      const std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals);

  /**
   * Similar to write_sideset_data(), this function is used to read
   * the data at a particular timestep. TODO: currently _all_ the
   * sideset variables are read, but we might want to change this to
   * only read the requested ones.
   */
  void
  read_sideset_data (int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<boundary_id_type>> & side_ids,
                     std::vector<std::map<BoundaryInfo::BCTuple, Real>> & bc_vals);

  /**
   * Similar to read_sideset_data(), but instead of creating one
   * std::map per sideset per variable, creates a single map of (elem,
   * side, boundary_id) tuples, and stores the exo file array indices
   * for any/all sideset variables on that sideset (they are all the
   * same). In cases where there are hundreds of sideset variables on
   * a single sideset, it is more efficient to store the array indices
   * in a quickly searchable data structure than to repeat the
   * indexing once per variable as is done in the read_sideset_data()
   * case.
   */
  void
  get_sideset_data_indices (std::map<BoundaryInfo::BCTuple, unsigned int> & bc_array_indices);

  /**
   * The Exodus format can also store values on nodesets. This can be
   * thought of as an alternative to defining a nodal variable
   * field on lower-dimensional elements making up a part of the
   * boundary. The inputs to the function are:
   * .) var_names[i] is the name of the ith sideset variable to be written to file.
   * .) node_boundary_ids[i] is a set of node_ids where var_names[i] is active.
   * .) bc_vals[i] is a map from (node-id, boundary-id) NodeBCTuple objects to the
   *    corresponding real-valued data.
   *
   * \note You must have already written the mesh by calling
   * e.g. write() before calling this function, because it uses the
   * existing ordering of the Exodus nodesets.
   */
  void
  write_nodeset_data (int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                      const std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals);

  /**
   * Read all the nodeset data at a particular timestep. TODO:
   * currently _all_ the nodeset variables are read, but we might want
   * to change this to only read the requested ones.
   */
  void
  read_nodeset_data (int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<boundary_id_type>> & node_boundary_ids,
                     std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> & bc_vals);

  /**
   * Similar to read_nodeset_data(), but instead of creating one
   * std::map per nodeset per variable, creates a single map of
   * (node_id, boundary_id) tuples, and stores the exo file array indices
   * for any/all nodeset variables on that nodeset (they are all the
   * same). In cases where there are hundreds of nodeset variables on
   * a single nodeset, it is more efficient to store the array indices
   * in a quickly searchable data structure than to repeat the
   * indexing once per variable as is done in the read_nodeset_data()
   * case.
   */
  void
  get_nodeset_data_indices (std::map<BoundaryInfo::NodeBCTuple, unsigned int> & bc_array_indices);

  /**
   * The Exodus format can also store values on elemsets. The inputs to the function are:
   * .) var_names[i] is the name of the ith elemset variable to be written to file.
   * .) elemset_ids_in[i] is a set of elemset ids where var_names[i] is active.
   * .) elemset_vals[i] is a map from (elem-id, elemset-id) pairs to the
   *    corresponding real-valued data.
   *
   * \note You must have already written the mesh by calling
   * e.g. write() (and write_elemsets()) before calling this function,
   * because it uses the ordering of the elemsets that have already
   * been written to the Exodus file.
   */
  void
  write_elemset_data (int timestep,
                      const std::vector<std::string> & var_names,
                      const std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                      const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals);

  /**
   * Read all the elemset data at a particular timestep.
   */
  void
  read_elemset_data (int timestep,
                     std::vector<std::string> & var_names,
                     std::vector<std::set<elemset_id_type>> & elemset_ids_in,
                     std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & elemset_vals);

  /**
   * Similar to read_elemset_data(), but instead of creating one
   * std::map per nodeset per variable, creates a single map of
   * (elem_id, elemset_id) pairs -> exo file array indices
   * for any/all variables on that elemset (they are all the
   * same). In cases where there are hundreds of elemset variables on
   * a single elemset, it is more efficient to store the array indices
   * in a quickly searchable data structure than to repeat the
   * indexing once per variable as is done in the read_elemset_data()
   * case.
   */
  void
  get_elemset_data_indices (std::map<std::pair<dof_id_type, elemset_id_type>, unsigned int> & elemset_array_indices);

  /**
   * Set the elemental variables in the Exodus file to be read into extra
   * element integers. The names of these elemental variables will be used to
   * name the extra element integers.
   */
  void set_extra_integer_vars(const std::vector<std::string> & extra_integer_vars);

  /**
   * Sets the list of variable names to be included in the output.
   * This is _optional_.  If this is never called then all variables
   * will be present. If this is called and an empty vector is supplied
   * no variables will be output. Setting the allow_empty = false will
   * result in empty vectors supplied here to also be populated with all
   * variables.
   */
  void set_output_variables(const std::vector<std::string> & output_variables,
                            bool allow_empty = true);

  /**
   * In the general case, meshes containing 2D elements can be
   * manifolds living in 3D space, thus by default we write all
   * meshes with the Exodus dimension set to LIBMESH_DIM =
   * mesh.spatial_dimension().
   *
   * In certain cases, however, the user may know his 2D mesh actually
   * lives in the z=0 plane, and therefore wants to write a truly 2D
   * Exodus mesh.  In such a case, he should call this function with
   * val=true.
   */
  void use_mesh_dimension_instead_of_spatial_dimension(bool val);

  /**
   * Directly control the num_dim which is written to the Exodus file.
   * If non-zero, this value supersedes all other dimensions, including:
   * 1.) MeshBase::spatial_dimension()
   * 2.) MeshBase::mesh_dimension()
   * 3.) Any value passed to use_mesh_dimension_instead_of_spatial_dimension()
   * This is useful/necessary for working around a bug in Paraview which
   * prevents the "Plot Over Line" filter from working on 1D meshes.
   */
  void write_as_dimension(unsigned dim);

  /**
   * Allows you to set a vector that is added to the coordinates of
   * all of the nodes. Effectively, this "moves" the mesh to a
   * particular position.
   *
   * \deprecated As requested by Roy in libmesh PR #90, this function
   * was "deprecated on arrival". There is not really a suitable
   * replacement for it in the works, however. The same effect *could*
   * be achieved by calling MeshTools::Modification::translate()
   * twice, but that approach seems inefficient in the case of very
   * large problems with millions of nodes. That said, this should
   * probably become a base class API so that it works for all the
   * different IO subclasses.
   */
  void set_coordinate_offset(Point p);

  /**
   * If true, this flag will cause the ExodusII_IO object to attempt to
   * open an existing file for writing, rather than creating a new file.
   * Obviously this will only work if the file already exists.
   */
  void append(bool val);

  /**
   * Return list of the elemental variable names
   */
  const std::vector<std::string> & get_elem_var_names();

  /**
   * Return list of the nodal variable names
   */
  const std::vector<std::string> & get_nodal_var_names();

  /**
   * Return list of the global variable names
   */
  const std::vector<std::string> & get_global_var_names();

  /**
   * Returns a const reference to the elem_num_map, which is a vector
   * that is created when a Mesh is read from file.  LibMesh will
   * number its mesh elements consistently with the elem_num_map
   * array, except that the indices in this array are 1-based, and
   * libmesh always uses a 0-based numbering. For example, given:
   * elem_num_map = [4,2,3,1]
   * libmesh will assign the first element it reads from the Exodus
   * file elem->id() == 3, the second will get elem->id() == 1, and so
   * on. We note that not all Exodus files contain an elem_num_map,
   * and in that case, calling this function will return a reference
   * to a vector containing the 1-based identity array, [1,2,3,...]
   */
  const std::vector<int> & get_elem_num_map() const;

  /**
   * Identical to the behavior of get_elem_num_map(), but for the
   * node_num_map instead.
   */
  const std::vector<int> & get_node_num_map() const;

#ifdef LIBMESH_HAVE_EXODUS_API
  /**
   * Return a reference to the ExodusII_IO_Helper object.
   */
  ExodusII_IO_Helper & get_exio_helper();
#endif

  /**
   * Set to true (the default) to write files in an HDF5-based file
   * format (when HDF5 is available), or to false to write files in
   * the old NetCDF3-based format.  If HDF5 is unavailable, this
   * setting does nothing.
   */
  void set_hdf5_writing(bool write_hdf5);

  /**
   * Set to true (false is the default) to generate independent nodes
   * for every Bezier Extraction element in an input file containing
   * them, causing those elements to be disconnected from each other.
   * (Their vertices will be multiple distinct nodes at overlapping
   * points.)  If this is false, Bezier Extraction elements will be
   * connected where possible, thereby using fewer redundancies and
   * less memory, but the input mesh will need to be conforming.
   */
  void set_discontinuous_bex(bool disc_bex);

  /**
   * This function factors out a bunch of code which is common to the
   * write_nodal_data() and write_nodal_data_discontinuous() functions
   */
  void write_nodal_data_common(std::string fname,
                               const std::vector<std::string> & names,
                               bool continuous=true);

private:
  /**
   * Only attempt to instantiate an ExodusII helper class
   * if the Exodus API is defined.  This class will have no
   * functionality when LIBMESH_HAVE_EXODUS_API is not defined.
   */
#ifdef LIBMESH_HAVE_EXODUS_API
  std::unique_ptr<ExodusII_IO_Helper> exio_helper;

  /**
   * Stores the current value of the timestep when calling
   * ExodusII_IO::write_timestep().
   */
  int _timestep;

  /**
   * should we be verbose?
   */
  bool _verbose;

  /**
   * Default false.  If true, files will be opened with EX_WRITE
   * rather than created from scratch when writing.
   */
  bool _append;
#endif

  /**
   * An optional list of variables in the EXODUS file that are
   * to be used to set extra integers when loading the file into a mesh.
   * The variable names will be used to name the extra integers.
   */
  std::vector<std::string> _extra_integer_vars;

  /**
   * The names of the variables to be output.
   * If this is empty then all variables are output.
   */
  std::vector<std::string> _output_variables;

  /**
   * Flag which controls the behavior of _output_variables:
   * .) If true, _output_variables is allowed to remain empty.
   * .) If false, if _output_variables is empty it will be populated
   *    with a complete list of all variables.
   * .) By default, calling set_output_variables() sets this flag to
   *    true, but it provides an override.
   */
  bool _allow_empty_variables;

  /**
   * By default, when complex numbers are enabled, for each variable
   * we write out three values: the real part, "r_u" the imaginary
   * part, "i_u", and the complex modulus, a_u := sqrt(r_u*r_u +
   * i_u*i_u), which is also the value returned by
   * std::abs(std::complex).  Since the modulus is not an independent
   * quantity, we can set this flag to false and save some file space
   * by not writing out.
   */
  bool _write_complex_abs;

  /**
   * Set to true (false is the default) to generate independent nodes
   * for every Bezier Extraction element.
   */
  bool _disc_bex;
};


} // namespace libMesh


#endif // LIBMESH_EXODUSII_IO_H
