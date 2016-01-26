// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_ENSIGHT_IO_H
#define LIBMESH_ENSIGHT_IO_H

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <map>
#include <string>
#include <vector>

namespace libMesh
{

// Forward declarations
class EquationSystems;

/**
 * This class implements writing meshes and solutions in Ensight's Gold format.
 *
 * \author Camata
 * \date 2009
 * \author J. W. Peterson (refactoring and iostreams implementation)
 * \date 2016
 */
class EnsightIO : public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.
   */
  EnsightIO (const std::string & filename,
             const EquationSystems & eq);

  /**
   * Empty destructor.
   */
  ~EnsightIO () {}

  /**
   * Tell the EnsightIO interface to output the finite element (not
   * SCALAR) variable named "s".  Note: you must call add_scalar() or
   * add_vector() (see below) at least once, otherwise only the Mesh
   * will be written out.
   */
  void add_scalar (const std::string & system,
                   const std::string & scalar_description,
                   const std::string & s);

  /**
   * Tell the EnsightIO interface that the variables (u,v) constitute
   * a vector.  Note: u and v must have the same FEType, and be
   * defined in the same system.
   */
  void add_vector (const std::string & system,
                   const std::string & vec_description,
                   const std::string & u,
                   const std::string & v);

  /**
   * Tell the EnsightIO interface that the variables (u, v, w)
   * constitute a vector.  Note: Requires a 3D mesh, u, v, and must
   * have the same FEType, and must be defined in the same system.
   */
  void add_vector (const std::string & system,
                   const std::string & vec_description,
                   const std::string & u,
                   const std::string & v,
                   const std::string & w);
  /**
   * Calls write_ascii() and write_case().
   * Writes case, mesh, and solution files named:
   * name.case             (contains a description of other files)
   * name.geo000           (mesh)
   * name_{varname}.scl000 (one file per scalar variable)
   * name_{vecname}.vec000 (one file per vector variable)
   */
  void write (Real time = 0);

  /**
   * Calls this->write(0);
   */
  virtual void write (const std::string & name) libmesh_override;

private:
  // Represents the vectors that are used by the EnsightIO
  struct Vectors
  {
    std::string description;
    std::vector<std::string> components;
  };

  // Represents the scalars
  struct Scalars
  {
    std::string scalar_name;
    std::string description;
  };

  // Store the variables of system
  struct SystemVars
  {
    std::vector<Vectors> EnsightVectors;
    std::vector<Scalars> EnsightScalars;
  };

  // private methods
  // write solution in ascii format file
  void write_ascii (Real time = 0);
  void write_scalar_ascii (const std::string & sys, const std::string & var);
  void write_vector_ascii (const std::string & sys, const std::vector<std::string> & vec, const std::string & var_name);
  void write_solution_ascii ();
  void write_geometry_ascii ();
  void write_case();

  // private Attributes
  std::string _ensight_file_name;
  std::vector<Real> _time_steps;

  // mapping from system names to variable names+descriptions
  typedef std::map <std::string, SystemVars> system_vars_map_t;
  system_vars_map_t _system_vars_map;

  // Reference to the EquationSystems we were constructed with
  const EquationSystems & _equation_systems;

  // static mapping between libmesh ElemTypes and Ensight element strings.
  static std::map<ElemType, std::string> _element_map;

  // Static function used to build the _element_map.
  static std::map<ElemType, std::string> build_element_map();
};


} // namespace libMesh


#endif // LIBMESH_ENSIGHT_IO_H
