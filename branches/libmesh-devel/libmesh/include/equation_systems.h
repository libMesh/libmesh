// $Id: equation_systems.h,v 1.23.2.3 2003-05-06 17:13:32 benkirk Exp $

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
#include <vector>
#include <string>

// Local Includes
#include "mesh_common.h"
#include "enum_xdr_mode.h"

// Forward Declarations
class Mesh;
class SystemBase;


/**
 * This is the base class for the \p EquationSystems<T_sys>,
 * providing rudimentary functionality concerning flags,
 * parameters etc.  The interesting things are handled in
 * the derived class.
 *
 * @author Benjamin S. Kirk, 2002-2003
 */

// ------------------------------------------------------------
// EquationSystems class definition
class EquationSystems
{
public:
  
  /**
   * Constructor.
   */
  EquationSystems (const Mesh& mesh);

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
   * Reinitialize all the systems
   */
  void reinit ();

  /**
   * @returns the number of equation systems.
   */
  unsigned int n_systems() const;

  /**
   * @returns a constant reference to the system named \p name.
   */
  const SystemBase & get_system(const std::string& name) const;
 
  /**
   * @returns a constant reference to the system named \p name.
   */
  SystemBase& get_system(const std::string& name);

  /**
   * @returns a constant reference to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys& get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys& get_system (const std::string& name);
  
  /**
   * @returns a reference to the system named \p name.
   */
  SystemBase& operator () (const std::string& name);
 
  /**
   * @returns a constant reference to the system name
   */
  const SystemBase& operator () (const std::string& name) const;
 
  /**
   * @returns a reference to system number \p num.
   */
  SystemBase& operator () (const unsigned int num);
 
  /**
   * @returns a constant reference to system number \p num.
   */
  const SystemBase& operator () (const unsigned int num) const;
  
  /**
   * Add the system of type \p system_type named \p name to the
   * systems array. Must be overloaded in the derived classes.
   */
  void add_system (const std::string& system_type,
		   const std::string& name);
  
  /**
   * Add the system named \p name to the systems array.
   * Must be overloaded in the derived classes.
   */
  template <typename T_sys>
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
   * Fill the input vector \p var_names with the names
   * of the variables for each system.
   */
  void build_variable_names (std::vector<std::string>& var_names) const;

  /**
   * Fill the input vector \p soln with the solution values for the
   * system named \p name.  Note that the input
   * vector \p soln will only be assembled on processor 0, so this
   * method is only applicable to outputting plot files from processor 0.
   */
  void build_solution_vector (std::vector<Number>& soln,
                              std::string& system_name,
                              std::string& variable_name) const;
  
  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).  Note that the input
   * vector \p soln will only be assembled on processor 0, so this
   * method is only applicable to outputting plot files from processor 0.
   */
  void build_solution_vector (std::vector<Number>& soln) const;
  
  /**
   * Read & initialize the systems from disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Note that the equation system can be defined without initializing
   * the data vectors to any solution values.  This can be done
   * by calling the routine with the read_data flag set to false.
   */
  void read(const std::string& name,
            const libMeshEnums::XdrMODE,
            const bool read_header=true,
            const bool read_data=true,
            const bool read_additional_data=true);

  /**
   * Write the systems to disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Note that the solution data can be omitted by calling
   * this routine with the write_data flag set to false.
   */
  void write(const std::string& name,
             const libMeshEnums::XdrMODE,
             const bool write_data=true,
             const bool write_additional_data=true) const;

  /**
   * @returns \p true when this equation system contains
   * identical data, up to the given threshold.  Delegates
   * most of the comparisons to perform to the responsible
   * systems
   */
  bool compare (const EquationSystems& other_es, 
                const Real threshold,
                const bool verbose) const;

  /**
   * Prints information about the equation systems.
   */
  void print_info () const;



  
  //-------------------------------------------------
  // Flags & parameters
  // 
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
   * @returns a constant reference to the mesh
   */
  const Mesh & get_mesh() const;

  
protected:

  
  /**
   * @returns a string containing information about the
   * flags and parameters.
   */
  std::string get_info() const;
    
  /**
   * The mesh data structure
   */ 
  const Mesh& _mesh;

  /**
   * Data structure holding the systems.
   */
  std::map<std::string, SystemBase*> _systems;
  
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
const Mesh & EquationSystems::get_mesh () const
{
  return _mesh;
}



inline
unsigned int EquationSystems::n_systems () const
{
  return _systems.size();
}



inline
unsigned int EquationSystems::n_flags () const
{
  return _flags.size();
}



inline
unsigned int EquationSystems::n_parameters () const
{
  return _parameters.size();
}



template <typename T_sys>
inline
const T_sys& EquationSystems::get_system (const std::string& name) const
{
  std::map<std::string, SystemBase*>::const_iterator
    pos = _systems.find(name);
  
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named " << name << " found!"
		<< std::endl;
      error();
    }

  return dynamic_cast<const T_sys&>(*(pos->second));
}



template <typename T_sys>
inline
T_sys& EquationSystems::get_system (const std::string& name)
{
  std::map<std::string, SystemBase*>::const_iterator pos = _systems.find(name);
  
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named " << name << " found!"
		<< std::endl;
      error();
    }

  return dynamic_cast<T_sys&>(*(pos->second));
}



#endif
