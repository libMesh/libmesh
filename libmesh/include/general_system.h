// $Id: general_system.h,v 1.7 2003-03-20 11:51:23 ddreyer Exp $

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



#ifndef __general_system_h__
#define __general_system_h__

// C++ includes

// Local Includes
#include "equation_systems.h"
#include "system_base.h"


// Forward Declarations
class GeneralSystem;


/**
 * This class provides a specific system class.  It aims
 * at general systems, offering some backward memory,
 * likely suitable for nonlinear problems.
 */

// ------------------------------------------------------------
// GeneralSystem class definition

class GeneralSystem : public SystemBase
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  GeneralSystem (EquationSystems<GeneralSystem>& es,
		 const std::string&              name,
		 const unsigned int              number,
		 const SolverPackage             solver_package);

  /**
   * Destructor.
   */
  ~GeneralSystem ();
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  void init ();

  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  void update ();
 
  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  void assemble ();
  
  /**
   * Assemble & solve the linear system.
   */
  std::pair<unsigned int, Real> solve ();

  /*   /\** */
  /*    * @returns the approximation order of variable number \p i. */
  /*    *\/ */
  /*   Order variable_order (const unsigned int i) const; */

  /*   /\** */
  /*    * @returns the approximation order of variable \p var. */
  /*    *\/ */
  /*   Order variable_order (const std::string& var) const; */
  
  /**
   * @returns "General System".  Helps in identifying
   * the system type in an equation system file.
   */
  static const std::string system_type () { return "General"; }

  //-----------------------------------------------------------------
  // access to the solution data fields
  /**
   * @returns the current solution for the specified global
   * DOF.
   */
  Number current_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the current solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number current_nodal_solution (const unsigned int node,
				 const std::string& var,
				 const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the current
   * solution.
   */
  Number current_nodal_solution (const unsigned int node,
				 const unsigned short int var=0,
				 const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the current solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number current_elem_solution (const unsigned int elem,
				const std::string& var,
				const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the current
   * solution.
   */
  Number current_elem_solution (const unsigned int elem,
				const unsigned short int var=0,
				const unsigned int index=0) const;
  
  /**
   * @returns the old solution for the specified global
   * DOF.
   */
  Number old_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the older solution for the specified global
   * DOF.
   */
  Number older_solution (const unsigned int global_dof_number) const;
  
  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the old solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number old_nodal_solution (const unsigned int node,
			     const std::string& var,
			     const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the old
   * solution.
   */
  Number old_nodal_solution (const unsigned int node,
			     const unsigned short int var=0,
			     const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at node point
   * \p node in the mesh for the older solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number older_nodal_solution (const unsigned int node,
			       const std::string& var,
			       const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * node point \p node in the mesh for the older
   * solution.
   */
  Number older_nodal_solution (const unsigned int node,
			       const unsigned short int var=0,
			       const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the old solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number old_elem_solution (const unsigned int elem,
			    const std::string& var,
			    const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the old
   * solution.
   */
  Number old_elem_solution (const unsigned int elem,
			    const unsigned short int var=0,
			    const unsigned int index=0) const;

  /**
   * @returns the value of variable \p var at elem point
   * \p elem in the mesh for the older solution.
   * Slower than using the integer id for the variable
   * since this requires a lookup.
   */
  Number older_elem_solution (const unsigned int elem,
			      const std::string& var,
			      const unsigned int index=0) const;
  
  /**
   * @returns the value of variable number \p var at
   * elem point \p elem in the mesh for the older
   * solution.
   */
  Number older_elem_solution (const unsigned int elem,
			      const unsigned short int var=0,
			      const unsigned int index=0) const;
  
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  AutoPtr<NumericVector<Number> > current_local_solution;
      
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand for the previous time step.
   */
  AutoPtr<NumericVector<Number> > old_local_solution;
  
  /**
   * All the values I need to compute my contribution
   * to the simulation at hand two time steps ago.
   */
  AutoPtr<NumericVector<Number> > older_local_solution;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function(void fptr(EquationSystems<GeneralSystem>& es,
				      const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function(void fptr(EquationSystems<GeneralSystem>& es,
					  const std::string& name));
  
  /**
   * Function that initializes the system.
   */
  void (* init_system_fptr) (EquationSystems<GeneralSystem>& es,
			     const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_fptr) (EquationSystems<GeneralSystem>& es,
			  const std::string& name);


protected:

  /**
   * Reference to the \p equation_systems data structure 
   * that handles us.   
   */
  EquationSystems<GeneralSystem>& equation_systems;

};



// ------------------------------------------------------------
// GeneralSystem inline methods
/* inline */
/* Order GeneralSystem::variable_order (const unsigned int i) const */
/* { */
/*   return variable_type(var_names[i]).order; */
/* } */



/* inline */
/* Order GeneralSystem::variable_order (const std::string& var) const */
/* { */
/*   std::map<std::string, FEType>::const_iterator */
/*     pos = var_type.find(var); */

/*   assert (pos != var_type.end()); */
  
/*   return pos->second.order; */
/* } */



inline
Number GeneralSystem::current_solution (const unsigned int global_dof_number) const
{
  
  assert (global_dof_number < _dof_map.n_dofs());
  
  return (*current_local_solution)(global_dof_number);
}



inline
Number GeneralSystem::current_nodal_solution (const unsigned int node,
					      const std::string& var,
					      const unsigned int index) const
{
  return current_nodal_solution (node, variable_number(var), index);
}



inline
Number GeneralSystem::current_nodal_solution (const unsigned int node,
					      const unsigned short int var,
					      const unsigned int index) const
{
  return (*current_local_solution)(_mesh.node(node).dof_number(number(), var, index));
}




inline
Number GeneralSystem::current_elem_solution (const unsigned int elem,
					     const std::string& var,
					     const unsigned int index) const
{
  return current_elem_solution(elem, variable_number(var), index);
}



inline
Number GeneralSystem::current_elem_solution (const unsigned int elem,
					     const unsigned short int var,
					     const unsigned int index) const
{
  return (*current_local_solution)(_mesh.elem(elem)->dof_number(number(), var, index));
}




inline
Number GeneralSystem::old_solution (const unsigned int global_dof_number) const
{
  assert (global_dof_number < _dof_map.n_dofs());
  
  return (*old_local_solution)(global_dof_number);
}



inline
Number GeneralSystem::older_solution (const unsigned int global_dof_number) const
{
  assert (global_dof_number < _dof_map.n_dofs());
  
  return (*older_local_solution)(global_dof_number);
}



inline
Number GeneralSystem::old_nodal_solution (const unsigned int node,
					  const std::string& var,
					  const unsigned int index) const
{
  return old_nodal_solution (node, variable_number(var), index);
}



inline
Number GeneralSystem::old_nodal_solution (const unsigned int node,
					  const unsigned short int var,
					  const unsigned int index) const
{
  return (*old_local_solution)(_mesh.node(node).dof_number(number(), var, index));
}



inline
Number GeneralSystem::older_nodal_solution (const unsigned int node,
					    const std::string& var,
					    const unsigned int index) const
{
  return older_nodal_solution(node, variable_number(var), index);
}



inline
Number GeneralSystem::older_nodal_solution (const unsigned int node,
					    const unsigned short int var,
					    const unsigned int index) const
{
  return (*older_local_solution)(_mesh.node(node).dof_number(number(), var, index));
}



inline
Number GeneralSystem::old_elem_solution (const unsigned int elem,
					 const std::string& var,
					 const unsigned int index) const
{
  return old_elem_solution (elem, variable_number(var), index);
}



inline
Number GeneralSystem::old_elem_solution (const unsigned int elem,
					 const unsigned short int var,
					 const unsigned int index) const
{
  return (*old_local_solution)(_mesh.elem(elem)->dof_number(number(), var, index));
}



inline
Number GeneralSystem::older_elem_solution (const unsigned int elem,
					   const std::string& var,
					   const unsigned int index) const
{
  return older_elem_solution (elem, variable_number(var), index);
}



inline
Number GeneralSystem::older_elem_solution (const unsigned int elem,
					   const unsigned short int var,
					   const unsigned int index) const
{
  return (*older_local_solution)(_mesh.elem(elem)->dof_number(number(), var, index));
}






#endif
