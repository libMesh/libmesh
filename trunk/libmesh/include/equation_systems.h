// $Id: equation_systems.h,v 1.9 2003-02-10 22:03:22 benkirk Exp $

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



#ifndef __equation_systems_h__
#define __equation_systems_h__

// C++ includes
#include <set>
#include <map>

// Local Includes
#include "mesh_common.h"
#include "xdr_cxx.h"
#include "enum_solver_package.h"


// Forward Declarations
#define SystemData GeneralSystem
class SystemData;
class Mesh;



/**
 * This contains one or more equation systems that are
 * to be solved in a simulation.  These equation systems
 * are identified by a user-specified name and are solved
 * in the order that they are declared.
 *
 * @author Benjamin S. Kirk, 2002-2003
 */

// ------------------------------------------------------------
// EquationSystems class definition
class EquationSystems
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  By default Petsc data structures
   * will be used, but this is optional via the \p sp
   * flag.
   */
  EquationSystems (const Mesh& mesh,

#if defined(HAVE_PETSC)
		   
		   // Default to PETSC solvers if they are there
		   const SolverPackage sp = PETSC_SOLVERS
		   
#elif defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)
		   
		   // Try LASPACK if PETSC is not available
		   const SolverPackage sp = LASPACK_SOLVERS
		   
#else
		   
		   // No linear solvers for you!
		   const SolverPackage sp = INVALID_SOLVER_PACKAGE
#endif
		   
		   );
  /**
   * Destructor.
   */
  ~EquationSystems ();
 
  /**
   * Returns tha data structure to a pristine state.
   */
  void clear ();
  
  /**
   * Initialize all the systems
   */
  void init ();

  /**
   * @returns the number of equation systems.
   */
  unsigned int n_systems() const;

  /**
   * Add the system named \p name to the systems array.
   */
  void add_system (const std::string& name);
  
  /**
   * Remove the system named \p name from the systems array.
   */
  void delete_system (const std::string& name);

  /**
   * @returns the total number of variables in all
   * systems.
   */
  unsigned int n_vars () const;
  
  /**
   * @returns the total number of degrees of freedom
   * in all systems.
   */
  unsigned int n_dofs () const;

  /**
   * @returns a reference to the system named \p name.
   */
  SystemData & operator () (const std::string& name);

  /**
   * @returns a constant reference to the system name
   */
  const SystemData & operator () (const std::string& name) const;

  /**
   * @returns a reference to system number \p num.
   */
  SystemData & operator () (const unsigned int num);

  /**
   * @returns a constant reference to system number \p num.
   */
  const SystemData & operator () (const unsigned int num) const;

  /**
   * @returns the name of the system number num.
   */
  const std::string & name (const unsigned int num) const;
  
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
   * Fill the input vector \p var_names with the names
   * of the variables for each system.
   */
  void build_variable_names (std::vector<std::string>& var_names);

  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names())
   */
  void build_solution_vector (std::vector<Complex>& soln);

  /**
   * @returns a constant reference to the mesh
   */
  const Mesh & get_mesh() const;
  
  /**
   * Read & initialize the systems from disk using the XDR data format. 
   * This format allows for machine-independent binary output.
   *
   * Note that the equation system can be defined without initializing
   * the data vectors to any solution values.  This can be done
   * by calling the routine with the read_data flag set to false.
   */
  void read(const std::string& name,
	    const Xdr::XdrMODE,
	    const bool read_header=true,
	    const bool read_data=true);

  /**
   * Write the systems to disk using the XDR data format. 
   * This format allows for machine-independent binary output.
   *
   * Note that the solution data can be omitted by calling
   * this routine with the write_data flag set to false.
   */
  void write(const std::string& name,
	     const Xdr::XdrMODE,
	     const bool write_data=true);

  /**
   * Prints information about the equation systems.
   */
  void print_info () const {std::cout << get_info() << std::endl; };

  /**
   * @returns a string containing information about the
   * equation systems.
   */
  std::string get_info() const;
  
  
 protected:

  
  /**
   * The mesh data structure
   */ 
  const Mesh& _mesh;
  
  /**
   * Flag indicating what linear solver package to use
   */
  const SolverPackage _solver_package;

  /**     
   * Data structure that holds the systems.
   */
  std::map<std::string, SystemData*> _systems;
  
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
unsigned int EquationSystems::n_systems () const
{
  return _systems.size();
};



inline
const Mesh & EquationSystems::get_mesh () const
{
  return _mesh;
};


#endif
