// $Id: system_data.h,v 1.4 2003-01-21 19:24:35 benkirk Exp $

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



#ifndef __system_data_h__
#define __system_data_h__

// C++ includes
#include <set>

// Local Includes
#include "mesh_common.h"
#include "mesh.h"
#include "dof_map.h"
#include "fe_type.h"
#include "petsc_interface.h"
#include "equation_systems.h"

// Forward Declarations



/**
 * This class contains information related to any physical
 * process that might be simulated.  Such information may
 * range from the actual solution values to algorithmic
 * flags that may be used to control the numerical methods
 * employed.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// SystemData class definition

class SystemData
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  SystemData (EquationSystems& es,
	      const std::string& name);

  /**
   * Destructor.
   */
  ~SystemData ();

  /**
   * @returns the system name.
   */
  std::string name () const;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();

  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  void update ();

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on all processors.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<number>& global_soln) const;

  /**
   * Fill the input vector \p global_soln so that it contains
   * the global solution on processor \p dest_proc.
   * Requires communication with all other processors.
   */
  void update_global_solution (std::vector<number>& global_soln,
			       const unsigned int dest_proc) const;

  /**
   * @returns a constant reference to this system's \p dof_map.
   */
  const DofMap & get_dof_map() const { return dof_map; };
  
  /**
   * @returns a writeable reference to this system's \p dof_map.
   */
  DofMap & get_dof_map() { return dof_map; };
  
  /**
   * Initialize the degrees of freedom and other 
   * required data on the current mesh.  Note that the matrix
   * is not initialized at this time since it may not be required
   * for all applications.
   */
  void init ();
  
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
  std::string variable_name(const unsigned int i) const;
  
  /**
   * @returns the variable number assoicated with
   * the user-specified variable named \p var.
   */
  unsigned short int variable_number (const std::string& var) const;

  /**
   * @returns the approximation order of variable number \p i.
   */
  Order variable_order (const unsigned int i) const;

  /**
   * @returns the approximation order of variable \p var.
   */
  Order variable_order (const std::string& var) const;

  /**
   * @returns the finite element type variable number \p i.
   */
  FEType variable_type (const unsigned int i) const;

  /**
   * @returns the finite element type for variable \p var.
   */
  FEType variable_type (const std::string& var) const;
  
  /**
   * @returns the current solution for the specified global
   * DOF.
   */
  number current_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the old solution for the specified global
   * DOF.
   */
  number old_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the older solution for the specified global
   * DOF.
   */
  number older_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the current solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number current_nodal_solution (const unsigned int node,
				 const std::string& var,
				 const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the current
   * solution.
   */
  number current_nodal_solution (const unsigned int node,
				 const unsigned short int var=0,
				 const unsigned int index=0) const;
  
  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the old solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number old_nodal_solution (const unsigned int node,
			     const std::string& var,
			     const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the old
   * solution.
   */
  number old_nodal_solution (const unsigned int node,
			     const unsigned short int var=0,
			     const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the older solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number older_nodal_solution (const unsigned int node,
			       const std::string& var,
			       const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the older
   * solution.
   */
  number older_nodal_solution (const unsigned int node,
			       const unsigned short int var=0,
			       const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the current solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number current_elem_solution (const unsigned int elem,
				const std::string& var,
				const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the current
   * solution.
   */
  number current_elem_solution (const unsigned int elem,
				const unsigned short int var=0,
				const unsigned int index=0) const;
  
  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the old solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number old_elem_solution (const unsigned int elem,
			    const std::string& var,
			    const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the old
   * solution.
   */
  number old_elem_solution (const unsigned int elem,
			    const unsigned short int var=0,
			    const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the older solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  number older_elem_solution (const unsigned int elem,
			      const std::string& var,
			      const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the older
   * solution.
   */
  number older_elem_solution (const unsigned int elem,
			      const unsigned short int var=0,
			      const unsigned int index=0) const;
  
  /**
   * Data structure describing the relationship between
   * nodes, variables, etc... and degrees of freedom.
   */
  DofMap dof_map;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function(void fptr(EquationSystems& es,
				      const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function(void fptr(EquationSystems& es,
					  const std::string& name));
  
  /**
   * Function that initializes the system.
   */
  void (* init_system) (EquationSystems& es,
			const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble) (EquationSystems& es,
		     const std::string& name);
  
#ifdef HAVE_PETSC

  /**
   * Solve the linear system.
   */
  std::pair<unsigned int, real> solve ();
  
  /**
   * Data structure to hold solution values.
   */
  PetscVector solution;

  /**
   * Data structure to hold the system right-hand-side.
   */
  PetscVector rhs;
  
  /**
   * Data structure to hold the system matrix.
   */
  PetscMatrix matrix;

  /**
   * Interface to the Petsc solvers.
   */
  PetscInterface petsc_interface;
  
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  PetscVector current_local_solution;
  
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand for the previous time step.
   */
  PetscVector old_local_solution;
  
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand two time steps ago.
   */
  PetscVector older_local_solution;
  
#endif

private:


  /**
   * A name associated with this system.
   */
  const std::string sys_name;
  
  /**
   * Reference to the \p equation_systems data structure used
   * for the simulation.   
   */
  EquationSystems& equation_systems;

  /**
   * Constant reference to the \p mesh data structure used
   * for the simulation.   
   */
  const Mesh& mesh;

  /**
   * The names of the variables associated with this
   * system.
   */
  std::vector<std::string> var_names;
  
  /**
   * The finite element type for the variables associated
   * with this system.
   */
  std::map<std::string, FEType> var_type;
  
  /**
   * The variable numbers correspoinding to user-specified
   * names.
   */
  std::map<std::string, unsigned short int> var_num;
};



// ------------------------------------------------------------
// SystemData inline methods
inline
std::string SystemData::name() const
{
  return sys_name;
};



inline
unsigned int SystemData::n_vars() const
{
  assert (!var_names.empty());

  return var_names.size();
};



inline
std::string SystemData::variable_name (unsigned int i) const
{
  assert (i < n_vars());

  return var_names[i];
};



inline
Order SystemData::variable_order (const unsigned int i) const
{
  return variable_type(var_names[i]).order;
};



inline
Order SystemData::variable_order (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = var_type.find(var);

  assert (pos != var_type.end());
  
  return pos->second.order;
};



inline
FEType SystemData::variable_type (const unsigned int i) const
{
  return variable_type(var_names[i]);
};



inline
FEType SystemData::variable_type (const std::string& var) const
{
  std::map<std::string, FEType>::const_iterator
    pos = var_type.find(var);

  assert (pos != var_type.end());
  
  return pos->second;
};



inline
unsigned int SystemData::n_dofs() const
{
  return dof_map.n_dofs();
};





inline
unsigned int SystemData::n_constrained_dofs() const
{
#ifdef ENABLE_AMR

  return dof_map.n_constrained_dofs();

#else

  std::cerr << "This feature requires AMR!" << std::endl;
  error();
  return 0;

#endif
};




inline
unsigned int SystemData::n_local_dofs() const
{
  return dof_map.n_dofs_on_processor(mesh.processor_id());
};



inline
number SystemData::current_solution (const unsigned int global_dof_number) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  assert (global_dof_number < dof_map.n_dofs());
  
  return current_local_solution(global_dof_number);

#endif
};



inline
number SystemData::old_solution (const unsigned int global_dof_number) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  assert (global_dof_number < dof_map.n_dofs());
  
  return old_local_solution(global_dof_number);

#endif
};



inline
number SystemData::older_solution (const unsigned int global_dof_number) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  assert (global_dof_number < dof_map.n_dofs());
  
  return older_local_solution(global_dof_number);

#endif
};



inline
number SystemData::current_nodal_solution (const unsigned int node,
					   const std::string& var,
					   const unsigned int index) const
{
  return current_nodal_solution (node, variable_number(var), index);
};



inline
number SystemData::current_nodal_solution (const unsigned int node,
					   const unsigned short int var,
					   const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return current_local_solution(dof_map.node_dof_number(node, var, index));

#endif
};



inline
number SystemData::old_nodal_solution (const unsigned int node,
				       const std::string& var,
				       const unsigned int index) const
{
  return old_nodal_solution (node, variable_number(var), index);
};



inline
number SystemData::old_nodal_solution (const unsigned int node,
				       const unsigned short int var,
				       const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return old_local_solution(dof_map.node_dof_number(node, var, index));

#endif
};



inline
number SystemData::older_nodal_solution (const unsigned int node,
					 const std::string& var,
					 const unsigned int index) const
{
  return older_nodal_solution (node, variable_number(var), index);
};



inline
number SystemData::older_nodal_solution (const unsigned int node,
					 const unsigned short int var,
					 const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return older_local_solution(dof_map.node_dof_number(node, var, index));

#endif
};



inline
number SystemData::current_elem_solution (const unsigned int elem,
					  const std::string& var,
					  const unsigned int index) const
{
  return current_elem_solution (elem, variable_number(var), index);
};



inline
number SystemData::current_elem_solution (const unsigned int elem,
					  const unsigned short int var,
					  const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return current_elem_solution(dof_map.elem_dof_number(elem, var, index));

#endif
};



inline
number SystemData::old_elem_solution (const unsigned int elem,
				      const std::string& var,
				      const unsigned int index) const
{
  return old_elem_solution (elem, variable_number(var), index);
};



inline
number SystemData::old_elem_solution (const unsigned int elem,
				      const unsigned short int var,
				      const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return old_local_solution(dof_map.elem_dof_number(elem, var, index));

#endif
};



inline
number SystemData::older_elem_solution (const unsigned int elem,
					const std::string& var,
					const unsigned int index) const
{
  return older_elem_solution (elem, variable_number(var), index);
};



inline
number SystemData::older_elem_solution (const unsigned int elem,
					const unsigned short int var,
					const unsigned int index) const
{
#ifndef HAVE_PETSC

  std::cerr << "ERROR: this feature requires Petsc support!"
	    << std::endl;

  error();

  return 0.;
  
#else

  return older_local_solution(dof_map.elem_dof_number(elem, var, index));

#endif
};


#endif
