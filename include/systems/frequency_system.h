// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FREQUENCY_SYSTEM_H
#define LIBMESH_FREQUENCY_SYSTEM_H

#include "libmesh/libmesh_config.h"

// Frequency domain solutions only possible with complex arithmetic
#if defined(LIBMESH_USE_COMPLEX_NUMBERS)

// Local Includes
#include "libmesh/linear_implicit_system.h"

// C++ includes
#include <string>
#include <vector>

namespace libMesh
{

/**
 * \p FrequencySystem provides a specific system class
 * for frequency-dependent (linear) systems.
 * Generally two solution flavors are possible:
 *
 * - @e fast solution of moderately-sized systems:
 *   For moderate numbers of dof, it is possible to keep
 *   frequency-independent matrices in memory.  For this,
 *   simply provide multiple frequencies prior to \p ini().
 *   Also provide functions \p _assemble_fptr and \p _solve_fptr
 *   that should compute the element contributions, and simply add
 *   these contributions to give the frequency-dependent overall
 *   matrix, respectively. For details see the examples section.
 *
 * - solution of @e large systems:
 *   When there is not enough space to keep the frequency-independent
 *   contributions in memory, the user need only provide a function
 *   \p _solve_fptr which assembles the overall, frequency-dependent
 *   matrix for the current frequency given in the parameter section
 *   of \p EquationSystems<FrequencySystem> named \p current \p frequency.
 *   For this to work, only provide @e one frequency.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class FrequencySystem : public LinearImplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FrequencySystem (EquationSystems & es,
                   const std::string & name_in,
                   const unsigned int number_in);
  /**
   * Destructor.
   */
  ~FrequencySystem ();

  /**
   * Clear all the data structures associated with
   * the system, but leave the frequencies untouched.
   * The frequencies belong to the \p EquationSystems
   * object.
   */
  virtual void clear () libmesh_override;

  /**
   * The full clear method also clears the frequencies
   * (stored as parameters of the \p EquationSystems
   * object).
   */
  void clear_all ();

  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  virtual void assemble () libmesh_override;

  /**
   * Solves the system for all frequencies.
   */
  virtual void solve () libmesh_override;

  /**
   * Solves the linear system for the
   * \f$ [ \texttt{n\_start, n\_stop} ]^{th} \f$
   * frequencies. The solution vectors are stored in automatically
   * allocated vectors named \p solution_nnnn.  For access to these vectors,
   * see \p System. When calling this, the frequency range should better
   * be already set.
   */
  void solve (const unsigned int n_start,
              const unsigned int n_stop);

  /**
   * @returns \p "Frequency".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const libmesh_override { return "Frequency"; }


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
   * By default, the solution for each frequency is copied
   * to a separate \p NumericVector.  This feature
   * can be disabled with \p allocate_solution_duplicates=false.
   */
  void set_frequencies_by_steps (const Real base_freq,
                                 const Real freq_step=0.,
                                 const unsigned int n_freq=1,
                                 const bool allocate_solution_duplicates=true);

  /**
   * Set the frequency range for which the system should
   * be solved.  \p n_freq frequencies are equally
   * distributed in the interval
   * \f$ [ \texttt{min\_freq, max\_freq} ] \f$ .
   * By default, the solution for each frequency is copied
   * to a separate \p NumericVector.  This feature
   * can be disabled with \p allocate_solution_duplicates=false.
   */
  void set_frequencies_by_range (const Real min_freq,
                                 const Real max_freq,
                                 const unsigned int n_freq,
                                 const bool allocate_solution_duplicates=true);

  /**
   * Set the frequency range by simply copying the values
   * from \p frequencies.
   * By default, the solution for each frequency is copied
   * to a separate \p NumericVector.  This feature
   * can be disabled with \p allocate_solution_duplicates=false.
   */
  void set_frequencies (const std::vector<Real> & frequencies,
                        const bool allocate_solution_duplicates=true);

  /**
   * @returns the number of frequencies to solve
   */
  unsigned int n_frequencies () const;

  /**
   * Register a required user function to use in assembling/solving the system.
   * It is intended to compute @e frequency-dependent data.  For proper
   * work of \p FrequencySystem, at least @e this function has to be provided
   * by the user.
   */
  void attach_solve_function(void fptr(EquationSystems & es,
                                       const std::string & name));

  /**
   * Function that computes frequency-dependent data of the system.
   */
  void (* solve_system) (EquationSystems & es,
                         const std::string & name);

  /**
   * @returns the number of iterations and the final residual.
   */
  std::pair<unsigned int, Real> get_rval (unsigned int n) const;

  /**
   * @returns a string of the form \p "frequency_x", where \p x is
   * the integer \p n.  Useful for identifying frequencies and
   * solution vectors in the parameters set of \p _equation_systems.
   */
  std::string form_freq_param_name(const unsigned int n) const;

  /**
   * @returns a string of the form \p "solution_x", where \p x is
   * the integer \p n.  Useful for identifying frequencies and
   * solution vectors in the vectors map of \p System.
   */
  std::string form_solu_vec_name(const unsigned int n) const;


protected:


  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   * The frequenices have to be set @e prior to calling
   * \p init().
   */
  virtual void init_data () libmesh_override;

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
   * when the solution for each frequency should be
   * stored in an additional vector, then this \p bool
   * is \p true, otherwise \p false.
   * \p _keep_solution_duplicates is implicitly set
   * through the \p set_frequencies methods.
   */
  bool _keep_solution_duplicates;

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
   * The number of iterations and the final residual
   * when the Ax=b is solved for multiple frequencies.
   */
  std::vector<std::pair<unsigned int, Real> > vec_rval;

};



// ------------------------------------------------------------
// FrequencySystem inline methods
inline
std::pair<unsigned int, Real> FrequencySystem::get_rval (unsigned int n) const
{
  libmesh_assert_less (n, vec_rval.size());

  return vec_rval[n];
}


} // namespace libMesh

#endif // if defined(LIBMESH_USE_COMPLEX_NUMBERS)

#endif // LIBMESH_FREQUENCY_SYSTEM_H
