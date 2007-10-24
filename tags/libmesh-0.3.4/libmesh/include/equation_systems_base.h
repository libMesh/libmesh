// $Id: equation_systems_base.h,v 1.3 2003-04-09 02:30:14 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __equation_systems_base_h__
#define __equation_systems_base_h__

// C++ includes
#include <set>
#include <map>
#include <vector>
#include <string>

// Local Includes
#include "mesh_common.h"
#include "enum_solver_package.h"


// Forward Declarations
class Mesh;



/**
 * This is the base class for the \p EquationSystems<T_sys>,
 * providing rudimentary functionality concerning flags,
 * parameters etc.  The interesting things are handled in
 * the derived class.
 *
 * @author Benjamin S. Kirk, 2002-2003
 */

// ------------------------------------------------------------
// EquationSystemsBase class definition
class EquationSystemsBase
{
protected:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  Derived classes decide whether
   * to use PETSc or LASPACK.
   */
  EquationSystemsBase (Mesh& mesh,
		       const SolverPackage sp);

  /**
   * Destructor.
   */
  virtual ~EquationSystemsBase ();
 
  /**
   * Returns tha data structure to a pristine state.
   */
  void clear ();
  
  /**
   * Initialize all the systems
   */
  void init ();
  

public:

  /**
   * @returns \p true if the flag \p fl is set, returns
   * \p false otherwise.
   */
  bool flag (const std::string& fl) const;

  /**
   * Defines the flag \p fl as \p true. This flag will
   * be used for all systems.
   */ 
  void set_flag (const std::string& fl);

  /**
   * Undefines the flag \p fl.
   */
  void unset_flag (const std::string& fl);

  /**
   * @returns the number of flags. 
   */
  unsigned int n_flags () const;
  
  /**
   * @returns the parameter value assoicated with \p id.
   */
  Real parameter (const std::string& id) const;

  /**
   * Defines the value of parameter \p id as \p value.
   * This parameter will be used for all systems, so
   * this method makes sure the parameter isn't already
   * set to avoid accidental overwriting.
   */
  Real & set_parameter (const std::string& id);
  
  /**
   * Undefines the value of parameter \p id. 
   */
  void unset_parameter (const std::string& id);

  /**
   * @returns the number of parameters. 
   */
  unsigned int n_parameters () const;

  /**
   * Fill the input vector \p var_names with the names
   * of the variables for each system.
   */
  virtual void build_variable_names (std::vector<std::string>& var_names) = 0;

  /**
   * Fill the input vector \p soln with the solution values for the
   * system named \p name.  Note that the input
   * vector \p soln will only be assembled on processor 0, so this
   * method is only applicable to outputting plot files from processor 0.
   */
  virtual void build_solution_vector (std::vector<Number>& soln,
				      std::string& system_name,
				      std::string& variable_name) = 0;
  
  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).  Note that the input
   * vector \p soln will only be assembled on processor 0, so this
   * method is only applicable to outputting plot files from processor 0.
   */
  virtual void build_solution_vector (std::vector<Number>& soln) = 0;

  /**
   * @returns a constant reference to the mesh
   */
  const Mesh & get_mesh() const;

  /**
   * @returns a reference to the mesh
   */
  Mesh & get_mesh();

  /**
   * @returns the solver package type currently in use
   */
  SolverPackage get_solver_package() const 
  { return _solver_package; }



protected:

  /**
   * @returns a string containing information about the
   * flags and parameters.
   */
  std::string get_info() const;
    
  /**
   * The mesh data structure
   */ 
  Mesh& _mesh;
  
  /**
   * Flag indicating what linear solver package to use
   */
  const SolverPackage _solver_package;

  /**
   * Data structure to hold user-specified flags.
   */
  std::set<std::string> _flags;

  /**
   * Data structore to hold user-specified parameters 
   */
  std::map<std::string, Real> _parameters;
};



// ------------------------------------------------------------
// EquationSystems inline methods
inline
const Mesh & EquationSystemsBase::get_mesh () const
{
  return _mesh;
}



inline
Mesh & EquationSystemsBase::get_mesh ()
{
  return _mesh;
}


inline
unsigned int EquationSystemsBase::n_flags () const
{
  return _flags.size();
}


inline
unsigned int EquationSystemsBase::n_parameters () const
{
  return _parameters.size();
}

#endif
