// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes

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
   * Destructor.
   */
  virtual ~ExodusII_IO ();

  /**
   * This method implements reading a mesh from a specified file.
   * Open the file named \p name and read the mesh in Sandia National Lab's
   * ExodusII format. This is the method to use for reading in meshes generated
   * by cubit.  Works in 2D for \p TRIs, \p TRI6s, \p QUAD s, and \p QUAD9s.
   * Works in 3D for \p TET4s, \p TET10s, \p HEX8s, and \p HEX27s.
   */
  virtual void read (const std::string & name) libmesh_override;

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string & fname) libmesh_override;

  /**
   * Set the flag indicating if we should be verbose.
   */
  void verbose (bool set_verbosity);

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
   * Backward compatibility version of function that takes a single variable name.
   *
   * \deprecated Use the version of copy_nodal_solution() that takes two names.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void copy_nodal_solution(System & system,
                           std::string var_name,
                           unsigned int timestep=1);
#endif

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
                                     const std::set<std::string> * system_names=libmesh_nullptr);

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
                                     const std::set<std::string> * system_names=libmesh_nullptr);

  /**
   * Write out element solution.
   */
  void write_element_data (const EquationSystems & es);

  /**
   * Bring in base class functionality for name resolution and to
   * avoid warnings about hidden overloaded virtual functions.
   */
  using MeshOutput<MeshBase>::write_nodal_data;

  /**
   * Write out a nodal solution.
   */
  virtual void write_nodal_data (const std::string &,
                                 const std::vector<Number> &,
                                 const std::vector<std::string> &) libmesh_override;

  /**
   * Write out a discontinuous nodal solution.
   */
  void write_nodal_data_discontinuous (const std::string &,
                                       const std::vector<Number> &,
                                       const std::vector<std::string> &) libmesh_override;

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
                       const std::set<std::string> * system_names=libmesh_nullptr);

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
   * The names of the variables to be output.
   * If this is empty then all variables are output.
   */
  std::vector<std::string> _output_variables;

  /**
   * This function factors out a bunch of code which is common to the
   * write_nodal_data() and write_nodal_data_discontinuous() functions
   */
  void write_nodal_data_common(std::string fname,
                               const std::vector<std::string> & names,
                               bool continuous=true);

  /**
   * If true, _output_variables is allowed to remain empty.
   * If false, if _output_variables is empty it will be populated with a complete list of all variables
   * By default, calling set_output_variables() sets this flag to true, but it provides an override.
   */
  bool _allow_empty_variables;

};


} // namespace libMesh


#endif // LIBMESH_EXODUSII_IO_H
