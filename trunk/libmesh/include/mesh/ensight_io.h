// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef __ensight_io_h__
#define __ensight_io_h__

// The c++ include
#include <string>
#include <map>
#include <vector>

// libMesh includes
#include "libmesh.h"
#include "enum_elem_type.h"
#include "mesh_base.h"
#include "mesh_output.h"

namespace libMesh
{

// Forward declarations
class EquationSystems;

/**
 * This class implements writing meshes and solutions in Ensight's Gold format.
 * \p author Camata
 */
class EnsightIO : public MeshOutput<MeshBase>
{
 public:

  /**
   * Constructor.
   */
  EnsightIO (const std::string &filename, const EquationSystems &eq);

  ~EnsightIO ();

  /**
   * add 2D vector: Tell the EnsightIO interface that the variables u and v are a vector.
   * Note that u and v should be the same variables defined in the system.
   */
  void add_vector (const std::string &system, const std::string &vec_description,
		   const std::string &u, const std::string &v);

  /**
   * add 3D vector: tell the EnsightIO interface that the variables u, v and w are vector components
   */
  void add_vector (const std::string &system, const std::string &vec_description,
		   const std::string &u, const std::string &v, const std::string &w);

  /**
   * add scalar: tell the EnsightIO interface that the variable s is a scalar
   */
  void add_scalar (const std::string &system, const std::string &scalar_description,
		   const std::string &s);

  /**
   * write solution
   */
  virtual void write (const std::string &name);

  /**
   * write solution
   */
  void write (const double time = 0);

  bool& has_mesh_refinement();

private:

  // Define aux. structures

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


  typedef std::map <std::string, SystemVars>           SystemsVarsMap;
  typedef std::map <std::string, SystemVars>::iterator SystemsVarsMapIterator;
  typedef std::pair<std::string, SystemVars>           SystemsVarsValue;

  // private methods
  // write solution in ascii format file
  void write_ascii (const double time = 0);
  void write_scalar_ascii (const std::string &sys, const std::string &var);
  void write_vector_ascii (const std::string &sys, const std::vector<std::string> &vec, const std::string &var_name);
  void write_solution_ascii ();
  void write_geometry_ascii ();


  void write_case();
  void elem_type_to_string (ElemType, char*);

  // private Attributes
  std::string     _ensight_file_name;
  std::vector<double>   _time_steps;
  SystemsVarsMap   _systems_vars_map;
  const EquationSystems &_equation_systems;
};


} // namespace libMesh


#endif // #define __ensight_io_h__
