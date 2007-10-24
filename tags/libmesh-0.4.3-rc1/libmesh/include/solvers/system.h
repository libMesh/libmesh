// $Id: system.h,v 1.3 2004-06-02 15:08:41 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __system_h__
#define __system_h__

// C++ includes

// Local Includes
#include "libmesh_common.h"
#include "fe_type.h"
#include "dof_map.h"
#include "auto_ptr.h"
#include "numeric_vector.h"
#include "reference_counted_object.h"



// Forward Declarations
class System;
class EquationSystems;
class Mesh;
class Xdr;
template <typename T> class SparseMatrix;
template <typename T> class LinearSolverInterface;

/**
 * This is the base class for classes which contain 
 * information related to any physical process that might be simulated.  
 * Such information may range from the actual solution values to 
 * algorithmic flags that may be used to control the numerical methods
 * employed.  In general, use an \p EqnSystems<T_sys> object to handle
 * one or more of the children of this class.
 * Note that templating \p EqnSystems relaxes the use of virtual members.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// System class definition
class System : public ReferenceCountedObject<System>
{
protected:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  System (EquationSystems& es,
	  const std::string& name,
	  const unsigned int number);
  
public:

  
  /**
   * Destructor.
   */
  virtual ~System ();

  /**
   * The type of system.
   */
  typedef System sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();
  
  /** 
   * Initializes degrees of freedom on the current mesh.
   * Sets the 
   */
  void init ();
  
  /**
   * Reinitializes degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications. @e Should be overloaded in derived classes.
   */
  virtual void reinit ();
  
  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  virtual void update ();
  
  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * @e Should be overloaded in derived classes.
   */
  virtual void assemble ();

  /**
   * Solves the system.  Must be overloaded in derived systems.
   */
  virtual void solve () = 0;
  
  /**
   * @returns \p true when the other system contains
   * identical data, up to the given threshold.  Outputs
   * some diagnostic info when \p verbose is set.
   */
  virtual bool compare (const System& other_system, 
			const Real threshold,
			const bool verbose) const;

  /**
   * @returns the system name.
   */
  const std::string & name () const;

  /**
   * @returns the type of system, helpful in identifying
   * which system type to use when reading equation system
   * data from file.  Should be overloaded in derived classes.
   */
  virtual std::string system_type () const { return "BasicSystem"; }

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh. 
   */
  void project_vector (NumericVector<Number>&) const;

  /**
   * Projects the vector defined on the old mesh onto the
   * new mesh. The original vector is unchanged and the new vector
   * is passed through the second argument.
   */
  void project_vector (const NumericVector<Number>&,
		       NumericVector<Number>&) const;
  
  /**
   * @returns the system number.   
   */
  unsigned int number () const;
  
  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on all processors.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number>& global_soln) const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on processor \p dest_proc.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Number>& global_soln,
			       const unsigned int dest_proc) const;

  /**
   * @returns a constant reference to this systems's \p _mesh.
   */
  const Mesh & get_mesh() const;
  
  /**
   * @returns a reference to this systems's \p _mesh.
   */
  Mesh & get_mesh();
  
  /**
   * @returns a constant reference to this system's \p _dof_map.
   */
  const DofMap & get_dof_map() const;
  
  /**
   * @returns a writeable reference to this system's \p _dof_map.
   */
  DofMap & get_dof_map();

  /**
   * @returns \p true if the system is active, \p false otherwise.
   * An active system will be solved.
   */
  bool active () const;

  /**
   * Activates the system.  Only active systems are solved.
   */
  void activate ();

  /**
   * Deactivates the system.  Only active systems are solved.
   */
  void deactivate ();
    
  /**
   * Adds the additional vector \p vec_name to this system.  Only
   * allowed @e prior to \p init().  All the additional vectors
   * are similarly distributed, like the \p solution,
   * and inititialized to zero.
   */
  NumericVector<Number> & add_vector (const std::string& vec_name);

  /**
   * @returns a const reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  const NumericVector<Number> & get_vector(const std::string& vec_name) const;

  /**
   * @returns a writeable reference to this system's @e additional vector
   * named \p vec_name.  Access is only granted when the vector is already
   * properly initialized.
   */
  NumericVector<Number> & get_vector(const std::string& vec_name);

  /**
   * @returns the number of vectors (in addition to the solution) handled by this system
   * This is the size of the \p _vectors map
   */
  unsigned int n_vectors () const;

  /**
   * @returns the number of variables in the system
   */
  unsigned int n_vars() const;

  /**
   * @returns the number of degrees of freedom in the system
   */
  unsigned int n_dofs() const;

  /**
   * Returns the number of active degrees of freedom
   * for this System.
   */
  unsigned int n_active_dofs() const;

  /**
   * @returns the number of constrained degrees of freedom
   * in the system.
   */
  unsigned int n_constrained_dofs() const;
  
  /**
   * @returns the number of degrees of freedom local
   * to this processor
   */
  unsigned int n_local_dofs() const;
  
  /**
   * Adds the variable \p var to the list of variables
   * for this system.
   */
  void add_variable (const std::string& var,
		     const FEType& type);

  /**
   * Adds the variable \p var to the list of variables
   * for this system.  Same as before, but assumes \p LAGRANGE
   * as default value for \p FEType.family.
   */
  void add_variable (const std::string& var,
		     const Order order = FIRST,
		     const FEFamily = LAGRANGE);

  /**
   * @returns the name of variable \p i.
   */
  const std::string & variable_name(const unsigned int i) const;
  
  /**
   * @returns the variable number assoicated with
   * the user-specified variable named \p var.
   */
  unsigned short int variable_number (const std::string& var) const;

  /**
   * @returns the finite element type variable number \p i.
   */
  const FEType & variable_type (const unsigned int i) const;

  /**
   * @returns the finite element type for variable \p var.
   */
  const FEType & variable_type (const std::string& var) const;

  /**
   * Reads the basic data for this System.
   */
  void read (Xdr& io, 
	     const bool read_header=true,
	     const bool read_additional_data=true);

  /**
   * Reads additional data, namely vectors, for this System.
   */
  void read_data (Xdr& io,
		  const bool read_additional_data=true);

  /**
   * Writes the basic data for this System.
   */
  void write (Xdr& io) const;

  /**
   * Writes additional data, namely vectors, for this System.
   */
  void write_data (Xdr& io,
		   const bool write_additional_data = true) const;

  /**
   * @returns a string containing information about the
   * system.
   */
  std::string get_info () const;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function (void fptr(EquationSystems& es,
				       const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function (void fptr(EquationSystems& es,
					   const std::string& name));
  
  /**
   * Function that initializes the system.
   */
  void (* init_system) (EquationSystems& es,
			const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_system) (EquationSystems& es,
			    const std::string& name);



  //--------------------------------------------------
  // The solution and solution access members
  
  /**
   * @returns the current solution for the specified global
   * DOF.
   */
  Number current_solution (const unsigned int global_dof_number) const;

  /**
   * Data structure to hold solution values.
   */
  AutoPtr<NumericVector<Number> > solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.  This vector is necessarily larger
   * than the \p solution vector in the case of a parallel
   * simulation.  The \p update() member is used to synchronize
   * the contents of the \p solution and \p current_local_solution
   * vectors.
   */
  AutoPtr<NumericVector<Number> > current_local_solution;


protected:

  
  /**
   * Initializes the data for the system.  Note that this is called
   * before any user-supplied intitialization function so that all
   * required storage will be available.
   */
  virtual void init_data ();

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update ();
  
  /**
   * A name associated with this system.
   */
  const std::string _sys_name;

  /**
   * The number associated with this system
   */
  const unsigned int _sys_number;

  /**
   * The names of the variables associated with this
   * system.
   */
  std::vector<std::string> _var_names;
  
  /**
   * The finite element type for the variables associated
   * with this system.
   */
  std::map<std::string, FEType> _var_type;
  
  /**
   * The variable numbers corresponding to user-specified
   * names.
   */
  std::map<std::string, unsigned short int> _var_num;

  /**
   * Flag stating if the system is active or not.
   */
  bool _active;
  
  /**
   * Some systems need an arbitrary number of vectors.
   * This map allows names to be associated with arbitrary
   * vectors.  All the vectors in this map will be distributed
   * in the same way as the solution vector.
   */
  std::map<std::string, NumericVector<Number>* > _vectors;

  /**
   * \p true when additional vectors may still be added, \p false otherwise.
   */
  bool _can_add_vectors;



  // -------------------------------------------------
  // Necessary classes
  //
  
  /**
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  DofMap _dof_map;

  /**
   * Constant reference to the \p EquationSystems object
   * used for the simulation.
   */
  EquationSystems& _equation_systems;
  
  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.   
   */
  Mesh& _mesh;
};



// ------------------------------------------------------------
// System inline methods
inline
const std::string & System::name() const
{
  return _sys_name;
}



inline
unsigned int System::number() const
{
  return _sys_number;
}



inline
const Mesh & System::get_mesh() const
{
  return _mesh;
}



inline
Mesh & System::get_mesh()
{
  return _mesh;
}



inline
const DofMap & System::get_dof_map() const
{
  return _dof_map;
}



inline
DofMap & System::get_dof_map()
{
  return _dof_map;
}



inline
bool System::active() const
{
  return _active;  
}



inline
void System::activate ()
{
  _active = true;
}



inline
void System::deactivate ()
{
  _active = false;
}



inline
unsigned int System::n_vars() const
{
  return _var_names.size();
}



inline
const std::string & System::variable_name (const unsigned int i) const
{
  assert (i < n_vars());

  return _var_names[i];
}



inline
const FEType & System::variable_type (const unsigned int i) const
{
  return variable_type(_var_names[i]);
}



inline
const FEType & System::variable_type (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = _var_type.find(var);

  assert (pos != _var_type.end());
  
  return pos->second;
}



inline
unsigned int System::n_dofs() const
{
  return _dof_map.n_dofs();
}


inline
unsigned int System::n_active_dofs() const
{
  return this->n_dofs() - this->n_constrained_dofs();
}


inline
unsigned int System::n_constrained_dofs() const
{
#ifdef ENABLE_AMR

  return _dof_map.n_constrained_dofs();

#else

  return 0;

#endif
}



inline
unsigned int System::n_local_dofs() const
{
  return _dof_map.n_dofs_on_processor (libMesh::processor_id());
}



inline
unsigned int System::n_vectors () const
{
 return _vectors.size(); 
}



inline
Number System::current_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < _dof_map.n_dofs());
  assert (global_dof_number < current_local_solution->size());
   
  return (*current_local_solution)(global_dof_number);
}


#endif
