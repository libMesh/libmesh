// $Id: transient_system.h,v 1.1 2003-11-05 22:26:42 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __transient_system_h__
#define __transient_system_h__

// C++ includes

// Local Includes
#include "steady_system.h"


// Forward Declarations


/**
 * This class provides a specific system class.  It aims
 * at transient systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p SystemBase.
 */

// ------------------------------------------------------------
// TransientSystem class definition

class TransientSystem : public SteadySystem
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  TransientSystem (EquationSystems& es,
		   const std::string& name,
		   const unsigned int number);

  /**
   * Destructor.
   */
  ~TransientSystem ();
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  void reinit ();
   
  /**
   * Update the local values to reflect the solution
   * on neighboring processors.
   */
  void update ();

//   /**
//    * Assemble the linear system.  Does not
//    * actually call the solver.
//    */
//   void assemble ();
  
//   /**
//    * Assemble & solve the linear system.
//    */
//   std::pair<unsigned int, Real> solve ();
  
  /**
   * @returns \p "Transient".  Helps in identifying
   * the system type in an equation system file.
   */
  std::string system_type () const { return "Transient"; }

//   /**
//    * Register a user function to use in initializing the system.
//    */
//   void attach_init_function (void fptr(EquationSystems& es,
// 				       const std::string& name));
  
//   /**
//    * Register a user function to use in assembling the system
//    * matrix and RHS.
//    */
//   void attach_assemble_function (void fptr(EquationSystems& es,
// 					   const std::string& name));
  
//   /**
//    * Function that initializes the system.
//    */
//   void (* init_system) (EquationSystems& es,
// 			const std::string& name);
  
//   /**
//    * Function that assembles the system.
//    */
//   void (* assemble_system) (EquationSystems& es,
// 			    const std::string& name);


  
  //-----------------------------------------------------------------
  // access to the solution data fields
  
  /**
   * @returns the old solution (at the previous timestep)
   * for the specified global DOF.
   */
  Number old_solution (const unsigned int global_dof_number) const;

  /**
   * @returns the older solution (two timesteps ago)
   * for the specified global DOF.
   */
  Number older_solution (const unsigned int global_dof_number) const;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  AutoPtr<NumericVector<Number> > old_local_solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  AutoPtr<NumericVector<Number> > older_local_solution;


protected:
  

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  void init_data ();

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  void re_update ();
};



// ------------------------------------------------------------
// TransientSystem inline methods
inline
Number TransientSystem::old_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < _dof_map.n_dofs());
  assert (global_dof_number < old_local_solution->size());
   
  return (*old_local_solution)(global_dof_number);
}



inline
Number TransientSystem::older_solution (const unsigned int global_dof_number) const
{
  // Check the sizes
  assert (global_dof_number < _dof_map.n_dofs());
  assert (global_dof_number < older_local_solution->size());
   
  return (*older_local_solution)(global_dof_number);
}



#endif
