// $Id: equation_systems.h,v 1.16 2005-06-06 14:53:18 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_common.h"
#include "parameters.h"
#include "system.h"
#include "enum_xdr_mode.h"

// HP aCC needs these for some reason
#ifdef __HP_aCC
# include "frequency_system.h"
# include "transient_system.h"
# include "newmark_system.h"
# include "steady_system.h"
#endif

// Forward Declarations
class MeshData;
//class System;
class Elem;
class Mesh;

/**
 * This is the \p EquationSystems class.  It is in charge
 * of handling all the various equation systems defined
 * for a \p Mesh.  It may have multiple systems, which may
 * be active or inactive, so that at different solution
 * stages only a sub-set may be solved for.  Also, through
 * the templated access, @e different types of systems
 * may be handled.  Also other features, like flags, 
 * parameters, I/O etc are provided.
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
  EquationSystems (Mesh& mesh, MeshData* mesh_data=NULL);

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
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const unsigned int num);

  /**
   * @returns a constant reference to the system named \p name.
   */
  const System & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   */
  System & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   */
  const System & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   */
  System & get_system (const unsigned int num);
  
  /**
   * @returns a reference to the system named \p name.
   */
  System & operator () (const std::string& name);
 
  /**
   * @returns a constant reference to the system name
   */
  const System & operator () (const std::string& name) const;
 
  /**
   * @returns a reference to system number \p num.
   */
  System & operator () (const unsigned int num);
 
  /**
   * @returns a constant reference to system number \p num.
   */
  const System & operator () (const unsigned int num) const;
  
  /**
   * Add the system of type \p system_type named \p name to the
   * systems array.
   */
  System & add_system (const std::string& system_type,
		       const std::string& name);
  
  /**
   * Add the system named \p name to the systems array.
   */
  template <typename T_sys>
  T_sys & add_system (const std::string& name);
  
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
   * Returns the number of active degrees of freedom
   * for the EquationSystems object.
   */
  unsigned int n_active_dofs() const;
  
  /**
   * Call \p solve on all the individual equation systems.
   */
  void solve ();
  
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
                              const std::string& system_name,
                              const std::string& variable_name = "all_vars") const;
  
  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).
   */
  void build_solution_vector (std::vector<Number>& soln) const;
  
  /**
   * Fill the input vector \p soln with solution values.  The
   * entries will be in variable-major format (corresponding to
   * the names from \p build_variable_names()).
   */
  void build_discontinuous_solution_vector (std::vector<Number>& soln) const;
  
  /**
   * Read & initialize the systems from disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Note that the equation system can be defined without initializing
   * the data vectors to any solution values.  This can be done
   * by calling the routine with the read_data flag set to false.
   */
  void read (const std::string& name,
	     const libMeshEnums::XdrMODE,
	     const bool read_header=true,
	     const bool read_data=true,
	     const bool read_additional_data=false);

  /**
   * Write the systems to disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Note that the solution data can be omitted by calling
   * this routine with the write_data flag set to false.
   */
  void write (const std::string& name,
	      const libMeshEnums::XdrMODE,
	      const bool write_data=true,
	      const bool write_additional_data=false) const;

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
   * @returns a string containing information about the
   * systems, flags, and parameters.
   */
  std::string get_info() const;
    
  /**
   * Prints information about the equation systems.
   */
  void print_info (std::ostream& os=std::cout) const;

  /**
   * Same as above, but allows you to also use stream syntax.
   */
  friend std::ostream& operator << (std::ostream& os, const EquationSystems& es);

  /**
   * @returns a constant reference to the mesh
   */
  const Mesh & get_mesh() const;

  /**
   * @returns a reference to the mesh
   */
  Mesh & get_mesh();

  /**
   * @returns a constant reference to the mesh_data
   */
  const MeshData & get_mesh_data() const;

  /**
   * @returns a reference to the mesh_data
   */
  MeshData & get_mesh_data();


  /**
   * Data structure holding arbitrary parameters.
   */
  Parameters parameters;

  
protected:

  
  /**
   * The mesh data structure
   */
  Mesh& _mesh;

  /**
   * A pointer to the MeshData object you would like to use.
   * Can be NULL.
   */
  MeshData* _mesh_data;

  /**
   * Data structure holding the systems.
   */
  std::map<std::string, System*> _systems;

  /**
   * Typedef for system iterators
   */
  typedef std::map<std::string, System*>::iterator       system_iterator;

  /**
   * Typedef for constatnt system iterators
   */
  typedef std::map<std::string, System*>::const_iterator const_system_iterator;

private:
  /**
   * This function is used in the implementation of add_system,
   * it loops over the nodes and elements of the Mesh, adding the
   * system to each one.  The main reason to separate this part
   * is to avoid coupling this header file to mesh.h, and elem.h.
   */
  void _add_system_to_nodes_and_elems();
};



// ------------------------------------------------------------
// EquationSystems inline methods
inline
const Mesh & EquationSystems::get_mesh () const
{
  return _mesh;
}



inline
Mesh & EquationSystems::get_mesh ()
{
  return _mesh;
}


inline
const MeshData & EquationSystems::get_mesh_data () const
{
  assert (_mesh_data != NULL);
  return *_mesh_data;
}


inline
MeshData & EquationSystems::get_mesh_data ()
{
  assert (_mesh_data != NULL);
  return *_mesh_data;
}



inline
unsigned int EquationSystems::n_systems () const
{
  return _systems.size();
}




template <typename T_sys>
inline
T_sys & EquationSystems::add_system (const std::string& name)
{
  T_sys* ptr = NULL;
  
  if (!_systems.count(name))
    {
      ptr = new T_sys(*this, name, this->n_systems());

      _systems.insert (std::make_pair(name, ptr));   
    }
  else
    {
      std::cerr << "ERROR: There was already a system"
		<< " named " << name
		<< std::endl;

      error();
    }

  assert (ptr != NULL);
  
  // Tell all the \p DofObject entities to add a system.
  this->_add_system_to_nodes_and_elems();

  // Return a dynamically casted reference to the newly added System.
  return *(dynamic_cast<T_sys*>(ptr));
}




template <typename T_sys>
inline
const T_sys & EquationSystems::get_system (const unsigned int num) const
{
  assert (num < this->n_systems());


  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
    {
      std::cerr << "ERROR: no system number " << num << " found!"
		<< std::endl;
      error();
    }

  // Attempt dynamic cast
  T_sys* ptr = dynamic_cast<T_sys*>(pos->second);

  // Check for failure of dynamic cast
  if (ptr == NULL)
    {
      std::cerr << "ERROR: cannot convert system "
		<< num << " to requested type!"
		<< std::endl;
      error();
    }
  
  return *ptr;
}




template <typename T_sys>
inline
T_sys & EquationSystems::get_system (const unsigned int num)
{
  assert (num < this->n_systems());
  
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
    {
      std::cerr << "ERROR: no system number " << num << " found!"
		<< std::endl;
      error();
    }

  // Attempt dynamic cast
  T_sys* ptr = dynamic_cast<T_sys*>(pos->second);

  // Check for failure of dynamic cast
  if (ptr == NULL)
    {
      std::cerr << "ERROR: cannot convert system "
		<< num << " to requested type!"
		<< std::endl;
      error();
    }

  return *ptr; 
}






template <typename T_sys>
inline
const T_sys & EquationSystems::get_system (const std::string& name) const
{
  const_system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named \"" << name << "\" found!"
		<< std::endl;
      error();
    }

  // Attempt dynamic cast
  T_sys* ptr = dynamic_cast<T_sys*>(pos->second);

  // Check for failure of dynamic cast
  if (ptr == NULL)
    {
      std::cerr << "ERROR: cannot convert system \""
		<< name << "\" to requested type!"
		<< std::endl;
      error();
    }

  return *ptr; 
}






template <typename T_sys>
inline
T_sys & EquationSystems::get_system (const std::string& name)
{
  system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named " << name << " found!"
		<< std::endl;
      error();
    }

  // Attempt dynamic cast
  T_sys* ptr = dynamic_cast<T_sys*>(pos->second);

  // Check for failure of dynamic cast
  if (ptr == NULL)
    {
      std::cerr << "ERROR: cannot convert system \""
		<< name << "\" to requested type!"
		<< std::endl;
      error();
    }

  return *ptr; 
}







inline
const System & EquationSystems::get_system (const std::string& name) const
{
  return this->get_system<System>(name);
}



inline
System & EquationSystems::get_system (const std::string& name)
{
  return this->get_system<System>(name);
}



inline
const System & EquationSystems::get_system (const unsigned int num) const
{
  return this->get_system<System>(num);
}



inline
System & EquationSystems::get_system (const unsigned int num)
{
  return this->get_system<System>(num);
}



inline
System & EquationSystems::operator () (const std::string& name)
{
  deprecated(); // Use the get_system() interface directly instead.
  return this->get_system (name);
}
  
 
inline
const System & EquationSystems::operator () (const std::string& name) const
{
  deprecated(); // Use the get_system() interface directly instead.
  return this->get_system (name);
}
  
 
 
inline
System & EquationSystems::operator () (const unsigned int num)
{
  deprecated(); // Use the get_system() interface directly instead.
  return this->get_system (num);
}



inline
const System & EquationSystems::operator () (const unsigned int num) const
{
  deprecated(); // Use the get_system() interface directly instead.
  return this->get_system (num);
}



#endif
