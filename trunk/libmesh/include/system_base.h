// $Id: system_base.h,v 1.5 2003-02-17 01:23:01 benkirk Exp $

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



#ifndef __system_base_h__
#define __system_base_h__

// C++ includes
#include <set>

// Local Includes
#include "mesh_common.h"
#include "mesh.h"
#include "dof_map.h"
#include "fe_type.h"
#include "auto_ptr.h"
#include "enum_solver_package.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "linear_solver_interface.h"
#include "reference_counted_object.h"


// Forward Declarations
class SystemBase;


/**
 * This is the base class for classes which contain 
 * information related to any physical process that might be simulated.  
 * Such information may range from the actual solution values to 
 * algorithmic flags that may be used to control the numerical methods
 * employed.  In general, use an \p EqnSystems<T_sys> object to handle
 * one or more of the children of this class.
 * Note that templating \p EqnSystems relaxes the use of virtual members.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// SystemBase class definition
class SystemBase : public ReferenceCountedObject<SystemBase>
{
protected:

  /**
   * Constructor.  Optionally initializes required
   * data structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  SystemBase (Mesh&               mesh,
	      const std::string&  name,
	      const unsigned int  number,
	      const SolverPackage solver_package);
  
  /**
   * Clear all the data structures associated with
   * the system.  Protected so that only children can
   * use it.
   */
  void clear ();
  
  /**
   * Initializes degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications. @e Should be overloaded in derived classes.
   * Protected so that only children can use it.
   */
  void init ();
  
  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * @e Should be overloaded in derived classes.
   * Protected so that only children can use it.
   */
  void assemble ();

public:

  /**
   * Destructor.
   */
  virtual ~SystemBase ();

  /**
   * @returns the system name.
   */
  const std::string & name () const;

  /**
   * @returns the system number.   
   */
  unsigned int number () const;
  
  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on all processors.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Complex>& global_soln) const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on processor \p dest_proc.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<Complex>& global_soln,
			       const unsigned int dest_proc) const;
  
  /**
   * @returns a constant reference to this system's \p _dof_map.
   */
  const DofMap & get_dof_map() const { return _dof_map; }
  
  /**
   * @returns a writeable reference to this system's \p _dof_map.
   */
  DofMap & get_dof_map() { return _dof_map; }

  /**
   * @returns the number of variables in the system
   */
  unsigned int n_vars() const;

  /**
   * @returns the number of degrees of freedom in the system
   */
  unsigned int n_dofs() const;

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
   * for this system.  Same as before, but assumes LAGRANGE
   * as default value for FEType.family.
   */
  void add_variable (const std::string& var,
		     const Order order);

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
   * Data structure to hold solution values.
   */
  AutoPtr<NumericVector> solution;

  /**
   * Data structure to hold the system right-hand-side.
   */
  AutoPtr<NumericVector> rhs;
  
  /**
   * Data structure to hold the system matrix.
   */
  AutoPtr<SparseMatrix> matrix;

  /**
   * Interface to the Petsc solvers.
   */
  AutoPtr<LinearSolverInterface> linear_solver_interface;


protected:

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
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  DofMap _dof_map;

  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.   
   */
  Mesh& _mesh;
};



// ------------------------------------------------------------
// SystemBase inline methods
inline
const std::string & SystemBase::name() const
{
  return _sys_name;
}



inline
unsigned int SystemBase::number() const
{
  return _sys_number;
}



inline
unsigned int SystemBase::n_vars() const
{
  assert (!_var_names.empty());

  return _var_names.size();
}



inline
const std::string & SystemBase::variable_name (const unsigned int i) const
{
  assert (i < n_vars());

  return _var_names[i];
}



inline
const FEType & SystemBase::variable_type (const unsigned int i) const
{
  return variable_type(_var_names[i]);
}



inline
const FEType & SystemBase::variable_type (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = _var_type.find(var);

  assert (pos != _var_type.end());
  
  return pos->second;
}



inline
unsigned int SystemBase::n_dofs() const
{
  return _dof_map.n_dofs();
}





inline
unsigned int SystemBase::n_constrained_dofs() const
{
#ifdef ENABLE_AMR

  return _dof_map.n_constrained_dofs();

#else

  std::cerr << "This feature requires AMR!" << std::endl;
  error();
  return 0;

#endif
}




inline
unsigned int SystemBase::n_local_dofs() const
{
  return _dof_map.n_dofs_on_processor(_mesh.processor_id());
}



#endif
