// $Id: frequency_system.h,v 1.10 2003-05-15 23:34:33 benkirk Exp $

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
#include <string>
#include <vector>


// Local Includes
#include "mesh_config.h"

/*
 * Frequency domain solutions only possible with complex arithmetic
 */
#if defined(USE_COMPLEX_NUMBERS)

#include "steady_system.h"

// Forward Declarations




/**
 * \p FrequencySystem provides a specific system class
 * for frequency-dependent (linear) systems. 
 * Generally two solution flavors are possible:
 *
 * - @e fast solution of moderately-sized systems:
 *   For moderate numbers of dof, it is possible to keep the mass, 
 *   damping and stiffness contribution of all elements in memory.
 *   These values may be stored in additional matrices.  The user-provided
 *   functions \p _assemble_fptr and \p _solve_fptr then should compute
 *   the element contributions, and simply add these contributions
 *   to give the frequency-dependent overall matrix, respectively.
 *   For details refer to the examples section.
 *
 * - solution of @e large systems with some support for multiple frequencies:
 *   When there is not enough space to keep the frequency-independent
 *   contributions in memory, the user need only provide a function
 *   \p _solve_fptr which assembles the overall, frequency-dependent
 *   matrix for the current frequency given in the parameter section
 *   of \p EquationSystems<FrequencySystem> named \p current \p frequency.
 *
 * \author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// FrequencySystem class definition

class FrequencySystem : public SteadySystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FrequencySystem (EquationSystems& es,
		   const std::string& name,
		   const unsigned int number);
  /**
   * Destructor.
   */
  ~FrequencySystem ();
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();
 
//   /**
//    * Reinitializes the member data fields associated with
//    * the system, so that, e.g., \p assemble() may be used.
//    */
//   void reinit ();

  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  void assemble ();

  
  /**
   * Solves the linear system for the
   * \f$ [ \texttt{n\_start\_in, n\_stop\_in} ]^{th} \f$ 
   * frequencies.  When neither \p n_start nor \p n_stop are given, solves for
   * all frequencies.  The solution vectors are stored in automatically
   * allocated vectors named \p solution_nnnn.  For access to these vectors,
   * see \p SystemBase. When calling this, the frequency range should better
   * be already set.
   */
  std::vector< std::pair<unsigned int, Real> > solve (const unsigned int n_start_in = static_cast<unsigned int>(-1),
						      const unsigned int n_stop_in  = static_cast<unsigned int>(-1));
  
  /**
   * @returns \p "Frequency".  Helps in identifying
   * the system type in an equation system file.
   */
  std::string system_type () const { return "Frequency"; }


  //--------------------------------------------------------
  // Methods specific to the FrequencySystem
  //
  
  /**
   * Set the frequency range for which the
   * system should be solved.  \p n_freq frequencies
   * are created, where the lowest and highest frequencies
   * are \p base_freq, and \p base_freq+freq_step*(n_freq-1),
   * respectively.  Calls to this of the form
   * \p set_frequencies_by_steps(30.) lets this object
   * solve the system for only this single frequency.
   */
  void set_frequencies_by_steps (const Real base_freq,
				 const Real freq_step=0.,
				 const unsigned int n_freq=0);

  /**
   * Set the frequency range for which the system should 
   * be solved.  \p n_freq frequencies are equally 
   * distributed in the interval 
   * \f$ [ \texttt{min\_freq, max\_freq} ] \f$ .
   */
  void set_frequencies_by_range (const Real min_freq,
				 const Real max_freq,
				 const unsigned int n_freq);

  /**
   * @returns the number of frequencies to solve
   */
  unsigned int n_frequencies () const { return _frequencies.size(); }

  /**
   * @returns a const reference to the frequencies to solve
   */
  const std::vector<Real>& get_frequencies () const
    { assert (_finished_set_frequencies); return _frequencies; }
  
  /**
   * Register a required user function to use in assembling/solving the system.
   * It is intended to compute @e frequency-dependent data.  For proper
   * work of \p FrequencySystem, at least @e this function has to be provided 
   * by the user.
   */
  void attach_solve_function(void fptr(EquationSystems& es,
				       const std::string& name));
  
  /**
   * Function that computes frequency-dependent data of the system.
   */
  void (* solve_system) (EquationSystems& es,
			 const std::string& name);

  /**
   * @returns a string of the form \p "frequency_x", where \p x is
   * the integer \p n.  Useful for identifying frequencies and 
   * solution vectors in the parameters set of \p _equation_systems.
   */
  std::string form_freq_param_name(const unsigned int n) const;

  /**
   * @returns a string of the form \p "solution_x", where \p x is
   * the integer \p n.  Useful for identifying frequencies and 
   * solution vectors in the vectors map of \p SystemBase.
   */
  std::string form_solu_vec_name(const unsigned int n) const;


protected:


  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.  
   * The frequenices have to be set @e prior to calling 
   * \p init().
   */
  void init_data ();
  
  /**
   * Sets the current frequency to the \p n-th entry in the vector
   * \p _frequencies.
   */
  void set_current_frequency(unsigned int n);

  /**
   * true when we have frequencies to solve for. 
   * Setting the frequencies is the first step after
   * creating/adding a \p FrequencySystem.
   */
  bool _finished_set_frequencies;

  /**
   * true when we have finished the \p init() phase.
   * This is the second step, and requires that
   * frequencies have already been set.
   */
  bool _finished_init;

  /**
   * true when we have finished the \p assemble() phase.
   * This is the third step and requires that
   * a) frequencies have already been set, and 
   * b) the system has been initialized.
   *
   * In this linear, frequency-dependent setting,
   * the overall system matrices \p mass, \p damping
   * and \p stiffness only have to be assembled once,
   * before multiple solutions may be obtained for
   * different frequencies. 
   */
  bool _finished_assemble;

  /**
   * The frequencies for which to solve the frequency-dependent
   * system.
   */
  std::vector<Real> _frequencies;
};



// ------------------------------------------------------------
// FrequencySystem inline methods



#endif // if defined(USE_COMPLEX_NUMBERS)

#endif // ifndef __frequency_system_h__
