// $Id: frequency_system.h,v 1.1 2003-02-12 02:03:47 ddreyer Exp $

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



#ifndef __frequency_system_h__
#define __frequency_system_h__

// C++ includes

// Local Includes
#include "mesh_config.h"
#include "equation_systems.h"
#include "system_base.h"

// Forward Declarations


/*
 * For the moment, only PETSc provides complex support
 */
#if defined(USE_COMPLEX_NUMBERS) && defined(HAVE_PETSC)



/**
 * \p FrequencySystem provides a specific system class
 * for frequency-dependent (linear) systems, providing
 * the typical setting of mass, damping and stiffness
 * matrices.  This class only appears when complex numbers
 * are enabled.
 */

// ------------------------------------------------------------
// FrequencySystem class definition

class FrequencySystem : public SystemBase
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FrequencySystem (EquationSystems& es,
		   const std::string& name,
		   const SolverPackage solver_package);

  /**
   * Destructor.
   */
  ~FrequencySystem ();
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();

  /**
   * Clears only the frequency-related data.
   */
  void clear_frequencies ();

  /**
   * Calls \p SystemBase::init() (check its documentation
   * concerning matrix initialization) and initializes the 
   * member data fields.  
   */
  void init ();
 
  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  void assemble ();
  
  /**
   * Set the frequency range for which the
   * system should be solved, must be called
   * prior to solution.  Also initializes
   * the \p _solutions vectors.
   */
  void set_frequencies (const Real base_freq,
			const Real freq_step,
			const unsigned int n_steps);


  /**
   * @returns the number of frequencies to solve
   */
  unsigned int n_frequencies () const 
    { assert (_have_freq); return _n_freq; };

  /**
   * Solve the linear system.  Prior to solving
   * the system, will check whether the mass, damping
   * stiffness matrices have already been assembled.
   */
  std::vector< std::pair<unsigned int, Real> > solve ();

  /**
   * @returns a const reference to the frequencies to solve
   */
  const std::vector<Real>& get_frequencies () const
    { assert (_have_freq); return _frequencies; };

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
  void (* init_system_fptr) (EquationSystems& es,
			     const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_fptr) (EquationSystems& es,
			  const std::string& name);
  
  /**
   * Data structure to hold the mass matrix
   */
  AutoPtr<SparseMatrix> mass;
  
  /**
   * Data structure to hold the damping matrix
   */
  AutoPtr<SparseMatrix> damping;

  /**
   * Data structure to hold the stiffness matrix.
   */
  AutoPtr<SparseMatrix> stiffness;


protected:

  /**
   * Reference to the \p equation_systems data structure 
   * that handles us.   
   */
  EquationSystems& equation_systems;


  /**
   * vector containing the solutions for the given
   * frequencies
   */
  std::vector< NumericVector* > _solutions;

  /**
   * In this linear, frequency-dependent setting,
   * the overall system matrices \p mass, \p damping
   * and \p stiffness only have to be assembled once,
   * before multiple solutions may be obtained for
   * different frequencies.  When the matrices already 
   * have been asssembled, this \p bool is \p true,
   * otherwise \p false.
   */
  bool _is_assembled;

  /**
   * is \p true when the frequency vector is set,
   * otherwise false.
   */
  bool _have_freq;

  /**
   * the frequencies at which to solve the system
   */
  std::vector<Real> _frequencies;

  /**
   * the number of frequencies to solve
   */
  unsigned int _n_freq;

};



// ------------------------------------------------------------
// FrequencySystem inline methods






#endif // if defined(USE_COMPLEX_NUMBERS) && defined(HAVE_PETSC)


#endif // ifndef __frequency_system_h__
