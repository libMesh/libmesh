// $Id: newmark_system.h,v 1.2 2003-04-09 19:26:57 ddreyer Exp $

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



#ifndef __newmark_system_h__
#define __newmark_system_h__

// C++ includes

// Local Includes
#include "system_base.h"


// Forward Declarations
class NewmarkSystem;
template <class T_sys> class EquationSystems;


/**
 * This class contains a specific system class.
 * It provides an implicit time integration scheme
 * known as the Newmark method.
 *
 * In the algorithm implemented here the system is solved for
 * displacements.
 * Curently the Newmark scheme is implemented for constant
 * time step sizes only. This time step is stored in the
 * \p EquationSystems parameter named \p "Newmark \p time \p step".
 * For the case of constant time steps the matrix only has to be
 * assembled once, whereas the rhs has to be updated in each timestep.
 * Default values of the Newmark parameters \p alpha and \p delta
 * used for time integration are provided.
 * For details refer to the examples section.
 */

// ------------------------------------------------------------
// NewmarkSystem class definition

class NewmarkSystem : public SystemBase
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  NewmarkSystem (EquationSystems<NewmarkSystem>& es,
		 const std::string&              name,
		 const unsigned int              number,
		 const SolverPackage             solver_package);

  /**
   * Destructor.
   */
  ~NewmarkSystem ();
  
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
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  void assemble ();
  
  /**
   * Apply initial conditions.
   */
  void initial_conditions ();

  /**
   * Compute the global matrix by adding up scaled
   * mass damping and stiffness matrix.
   */
  void compute_matrix ();


  /**
   * Update the rhs.
   */
  void update_rhs ();


  /**
   * Update displacement, velocity and acceleration.  
   */
  void update_u_v_a ();


  /**
   * Solve the linear system. Does not call assemble.
   */
  std::pair<unsigned int, Real> solve ();

  
  /**
   * @returns \p "Newmark".  Helps in identifying
   * the system type in an equation system file.
   */
  static const std::string system_type () { return "Newmark"; }

  /**
   * @returns \p true when the other system contains
   * identical data, up to the given threshold.
   */
  bool compare (const NewmarkSystem& other_system, 
		const Real threshold,
		const bool verbose) const;


  /**
   * Set the time step size and the newmark parameter alpha and
   * delta and calculate the constant parameters used for
   * time integration.
   */
  void set_newmark_parameters (const Real delta_T,
			       const Real alpha=.25,
			       const Real delta=.5);


  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_cond_function(void fptr(EquationSystems<NewmarkSystem>& es,
					   const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function(void fptr(EquationSystems<NewmarkSystem>& es,
					  const std::string& name));

  /**
   * Register a required user function to use in assembling/solving the system.  
   * Right before solving for a new frequency or time step, this function is called.
   */
  void attach_solve_function(void fptr(EquationSystems<NewmarkSystem>& es,
				       const std::string& name));

  /**
   * Function that applies initial conditions.
   */
  void (* init_cond_fptr) (EquationSystems<NewmarkSystem>& es,
			     const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_fptr) (EquationSystems<NewmarkSystem>& es,
			  const std::string& name);


 /**
   * Function that computes time or frequency dependent data of the system.
   */
  void (* solve_fptr) (EquationSystems<NewmarkSystem>& es,
			const std::string& name);


protected:

  /**
   * Reference to the \p EquationSystems<NewmarkSystem> data structure 
   * that handles us.   
   */
  EquationSystems<NewmarkSystem>& _equation_systems;


private:

    /**
     * Constants used for the time integration.
     */
    Real _a_0;
    Real _a_1;
    Real _a_2;
    Real _a_3;
    Real _a_4;
    Real _a_5;
    Real _a_6;
    Real _a_7;

    /**
     * Returns true if the matrix assambly is finished.
     */
    bool _finished_assemble;
  
};



// ------------------------------------------------------------
// NewmarkSystem inline methods



#endif
