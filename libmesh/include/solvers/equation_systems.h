// $Id: equation_systems.h,v 1.24 2007-10-03 20:31:42 roystgnr Exp $

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
   * Define enumeration to set properties in EquationSystems::read()
   */
  enum ReadFlags { READ_HEADER          = 1,
                   READ_DATA            = 2,
                   READ_ADDITIONAL_DATA = 4 };

  /**
   * Define enumeration to set properties in EquationSystems::write()
   */
  enum WriteFlags { WRITE_DATA            = 1,
                    WRITE_ADDITIONAL_DATA = 2 };
  
  /**
   * Constructor.
   */
  EquationSystems (Mesh& mesh, MeshData* mesh_data=NULL);

  /**
   * Destructor.  Should be virtual, since the user may want to derive
   * subclasses of EquationSystems.
   */
  virtual ~EquationSystems ();
 
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
   * Updates local values for all the systems
   */
  void update ();

  /**
   * @returns the number of equation systems.
   */
  unsigned int n_systems() const;

  /**
   * @returns true if the system named \p name exists within
   * this EquationSystems object.
   */
  bool has_system (const std::string& name) const;

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
  
//   /**
//    * @returns a reference to the system named \p name.
//    */
//   System & operator () (const std::string& name);
 
//   /**
//    * @returns a constant reference to the system name
//    */
//   const System & operator () (const std::string& name) const;
 
//   /**
//    * @returns a reference to system number \p num.
//    */
//   System & operator () (const unsigned int num);
 
//   /**
//    * @returns a constant reference to system number \p num.
//    */
//   const System & operator () (const unsigned int num) const;
  
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
   * This function is now deprecated - write the
   * libmesh-devel mailing list if you need it reimplemented.
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
   *
   * By default this function solves each equation system once,
   * in the order they were added.  For more sophisticated decoupled
   * problems the user may with to override this behavior in a derived
   * class.
   */
  virtual void solve ();
  
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
   * Set which sections of the file to read by bitwise OR'ing the 
   * EquationSystems::ReadFlags enumeration together. For example, to 
   * read all sections of the file, set read_flags to:
   * (READ_HEADER | READ_DATA | READ_ADDITIONAL_DATA)
   *
   * Note that the equation system can be defined without initializing
   * the data vectors to any solution values.  This can be done
   * by omitting READ_DATA in the read_flags parameter.
   */
  void read (const std::string& name,
	     const libMeshEnums::XdrMODE,
             const unsigned int read_flags=(READ_HEADER | READ_DATA));

  /**
   * Write the systems to disk using the XDR data format.
   * This format allows for machine-independent binary output.
   *
   * Set the writing properties using the EquationSystems::WriteFlags
   * enumeration. Set which sections to write out by bitwise OR'ing
   * the enumeration values. Write everything by setting write_flags to:
   * (WRITE_DATA | WRITE_ADDITIONAL_DATA)
   *
   * Note that the solution data can be omitted by calling
   * this routine with WRITE_DATA omitted in the write_flags argument.
   */
  void write (const std::string& name,
	      const libMeshEnums::XdrMODE,
              const unsigned int write_flags=WRITE_DATA) const;

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
   * @returns true when the _mesh_data pointer is not NULL.
   * This is needed because get_mesh_data will fail if it is NULL
   */
  const bool has_mesh_data() const;

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
const bool EquationSystems::has_mesh_data () const 
{
  return (_mesh_data!=NULL);
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

      // Tell all the \p DofObject entities to add a system.
      this->_add_system_to_nodes_and_elems();
    }
  else
    {
      // We now allow redundant add_system calls, to make it
      // easier to load data from files for user-derived system
      // subclasses
//      std::cerr << "ERROR: There was already a system"
//		<< " named " << name
//		<< std::endl;

//      error();

      ptr = &(this->get_system<T_sys>(name));
    }

  // Return a dynamically casted reference to the newly added System.
  return *ptr;
}



inline
bool EquationSystems::has_system (const std::string& name) const
{
  if (_systems.find(name) == _systems.end())
    return false;
  return true;
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



// inline
// System & EquationSystems::operator () (const std::string& name)
// {
//   deprecated(); // Use the get_system() interface directly instead.
//   return this->get_system (name);
// }
  
 
// inline
// const System & EquationSystems::operator () (const std::string& name) const
// {
//   deprecated(); // Use the get_system() interface directly instead.
//   return this->get_system (name);
// }
  
 
 
// inline
// System & EquationSystems::operator () (const unsigned int num)
// {
//   deprecated(); // Use the get_system() interface directly instead.
//   return this->get_system (num);
// }



// inline
// const System & EquationSystems::operator () (const unsigned int num) const
// {
//   deprecated(); // Use the get_system() interface directly instead.
//   return this->get_system (num);
// }



#endif
