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
 * @author Benjamin Kirk, John Peterson, 2004.
 */

// ------------------------------------------------------------
// ExodusII_IO class definition
class ExodusII_IO : public MeshInput<MeshBase>,
		    public MeshOutput<MeshBase>,
		    public ParallelObject
{
 public:

  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  explicit
  ExodusII_IO (MeshBase& mesh);

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
  virtual void read (const std::string& name);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& fname);

  /**
   * Set the flag indicating if we should be verbose.
   */
  void verbose (bool set_verbosity);

  /**
   * Returns an array containing the timesteps in the file
   */
  const std::vector<Real>& get_time_steps();

  /**
   * Backward compatibility version of function that takes a single variable name
   */
  void copy_nodal_solution(System& system, std::string var_name, unsigned int timestep=1);

  /**
   * If we read in a nodal solution while reading in a mesh, we can attempt
   * to copy that nodal solution into an EquationSystems object.
   */
  void copy_nodal_solution(System& es, std::string system_var_name, std::string exodus_var_name, unsigned int timestep=1);

  /**
   * If we read in a elemental solution while reading in a mesh, we can attempt
   * to copy that elemental solution into an EquationSystems object.
   */
  void copy_elemental_solution(System& es, std::string system_var_name, std::string exodus_var_name, unsigned int timestep=1);

  /**
   * Writes a exodusII file with discontinuous data
   */
  void write_discontinuous_exodusII (const std::string& name, const EquationSystems& es);

  /**
   * Write out element solution.
   */
  void write_element_data (const EquationSystems& es);

  /**
   * Write out a nodal solution.
   */
  void write_nodal_data (const std::string&,
			 const std::vector<Number>&,
			 const std::vector<std::string>&);

  /**
   * Write out a discontinuous nodal solution.
   */
  void write_nodal_data_discontinuous (const std::string&,
                                       const std::vector<Number>&,
                                       const std::vector<std::string>&);

  /**
   * Write out global variables.
   */
  void write_global_data (const std::vector<Number>&,
                          const std::vector<std::string>&);

  /**
   * Write out information records.
   */
  void write_information_records (const std::vector<std::string>&);

  /**
   * Writes out the solution at a specific timestep.
   * @param timestep The timestep to write out, should be _1_ indexed.
   */
  void write_timestep (const std::string& fname,
		       const EquationSystems& es,
		       const int timestep,
		       const Real time);

  /**
   * Sets the list of variable names to be included in the output.
   * This is _optional_.  If this is never called then all variables
   * will be present.
   */
  void set_output_variables(const std::vector<std::string> & output_variables);

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
   * Allows you to set a vector that is added to the coordinates of all
   * of the nodes.  Effectively, this "moves" the mesh to a particular position
   */
  void set_coordinate_offset(Point p);

 private:
  /**
   * Only attempt to instantiate an ExodusII helper class
   * if the Exodus API is defined.  This class will have no
   * functionality when LIBMESH_HAVE_EXODUS_API is not defined.
   */
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO_Helper *exio_helper;
#endif

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
   * The names of the variables to be output.
   * If this is empty then all variables are output.
   */
  std::vector<std::string> _output_variables;
};


} // namespace libMesh


#endif // LIBMESH_EXODUSII_IO_H
