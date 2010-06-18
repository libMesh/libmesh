// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef __qn_transient_scm_system_h__
#define __qn_transient_scm_system_h__

// Configuration data
#include "libmesh_config.h"

// This class derives from RBSCMSystem, which requires SLEPC
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "rb_scm_system.h"

/**
 * This class is part of the rbOOmit framework.
 *
 * QNTransientSCMSystem extends RBSCMSystem in order to
 * implement the Succesive Constraint Method (SCM) for
 * transient quadratically non-linear problems.
 *
 * @author David J. Knezevic 2009
 */

// ------------------------------------------------------------
// QNTransientSCMSystem class definition

class QNTransientSCMSystem : public RBSCMSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  QNTransientSCMSystem (EquationSystems& es,
                        const std::string& name,
                        const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~QNTransientSCMSystem ();

  /**
   * The type of system.
   */
  typedef QNTransientSCMSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef RBSCMSystem Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * If the number of basis functions in the associated RBSystem
   * changes then call this function to bring this SCM system
   * up to date.
   */
  void resize_to_new_n_bfs();

  /**
   * Clear and resize the SCM data vectors. Overloaded
   * to also clear C_J_RB_coeffs.
   */
  virtual void resize_SCM_vectors ();

  /**
   * Perform the SCM greedy algorithm to develop a lower
   * bound over the training set. Overloaded to do some
   * extra steps required for the transient quadratically
   * nonlinear case.
   */
  virtual void perform_SCM_greedy();

  /**
   * Attach the function defining theta_c, i.e. the parameter dependent
   * function corresponding to the trilinear form.
   */
  void attach_theta_c(theta_q_fptr theta_c_in);

  /**
   * Evaluate theta_c at the current parameter.
   */
  Real eval_theta_c();

  /**
   * Overload eval_theta_q_a because the first n_bfs functions
   * have to be given a special definition, i.e. they now
   * come from the RB coefficients.
   */
  virtual Real eval_theta_q_a(unsigned int q);

  /**
   * Overload get_Q_a to increment by n_bfs.
   */
  virtual unsigned int get_Q_a() { return Parent::get_Q_a() + get_n_basis_functions(); }

  /**
   * Get the specified entry of RB_solution_data, the RB coefficients
   * corresponding to the training parameters.
   */
  std::vector<Real> get_training_RB_coeffs(unsigned int index);

  /**
   * Set the specified entry of training_RB_coeffs, the RB coefficients
   * corresponding to the training parameters.
   */
  void set_training_RB_coeffs(unsigned int index, std::vector<Real> coeffs);

  /**
   * Set current_RB_coeffs.
   */
  void set_current_RB_coeffs(std::vector<Real> coeffs) { current_RB_coeffs = coeffs; }

  /**
   * Get/set the n_bfs parameter; this determines the number of
   * affine forms we require.
   */
  void set_n_basis_functions(unsigned int n_bfs_in) { this->n_bfs = n_bfs_in; }
  unsigned int get_n_basis_functions()              { return n_bfs; }

  /**
   * Get/get the number of time-steps.
   */
  unsigned int get_K() { return _K; }
  void set_K(unsigned int K) { this->_K = K; }

  /**
   * Get/set the time-step size.
   */
  Number get_dt() { return dt; }
  void set_dt(Number dt) { this->dt = dt; }

  /**
   * Useful function that helps us explore the magnitude of the
   * stability constants. Here we compute \rho based on the
   * truth solution (stored in this->current_newton_iterate) and hence
   * doesn't rely on RB approximation. Must be implemented in
   * a derived class.
   */
  virtual Number compute_truth_stability_constant() { libmesh_not_implemented(); return 0.; }

  /**
   * Overload training parameter initialization to also include time
   * as a parameter.
   */
  virtual void initialize_training_parameters(const std::vector<Number>& mu_min_vector,
                                              const std::vector<Number>& mu_max_vector,
                                              const unsigned int n_param_samples,
                                              const std::vector<bool> log_param_scale,
                                              const bool deterministic=true);

  /**
   * Overwrite the training parameters with new_training_set. Overloaded to also
   * handle the extra temporal parameters in the QNTransientSCM case.
   */
  virtual void load_training_set(std::vector< std::vector<Real> >& new_training_set);

  /**
   * Write out all the data to text files in order to segregate the
   * Offline stage from the Online stage.
   */
  virtual void write_offline_data_to_files(const std::string& directory_name = "offline_data");

  /**
   * Read in the saved Offline reduced basis data
   * to initialize the system for Online solves.
   */
  virtual void read_offline_data_from_files(const std::string& directory_name = "offline_data");

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Load the n_bfs trilinear form operators from the
   * associated RB system.
   */
  virtual void add_scaled_symm_Aq(unsigned int q_a, Real scalar);

  /**
   * Load the specified training parameter into current_params.
   * Local to this processor. Overloaded to load RB_coeffs as well.
   */
  virtual void load_training_parameter_locally(unsigned int index);

  /**
   * Load the specified training parameter into current_params.
   * Loads the parameter on all processors. Overloaded to load RB_coeffs as well.
   */
  virtual void load_training_parameter_globally(unsigned int index);

  /**
   * Helper function to set current_parameters to the
   * specified parameters saved in C_J. Overloaded to set RB_coeffs as well.
   */
  virtual void set_current_parameters_from_C_J(unsigned int C_J_index);

  /**
   * Helper function to save current_parameters in
   * saved_parameters. Overloaded to save RB_coeffs as well.
   */
  virtual void save_current_parameters();

  /**
   * Helper functiont to (re)load current_parameters
   * from saved_parameters. Overloaded to save RB_coeffs as well.
   */
  virtual void reload_current_parameters();

  /**
   * Static utility function for rounding to the nearest
   * integer.
   */
  static int round(Number x)
  {
    return static_cast<int>(x > 0.0 ? x + 0.5 : x - 0.5);
  }

  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * The RB solution coefficients that we are currently considering!
   */
  std::vector<Number> current_RB_coeffs;

  /**
   * Vector storing the RB coefficients corresponding to the i^th entry of C_J.
   */
  std::vector< std::vector<Number> > C_J_RB_coeffs;

  /**
   * The RB solution data for each training parameter and all
   * time-levels.
   */
  std::vector< NumericVector<Real>* > training_RB_coeffs;

  /**
   * A vector for saving RB coefficients; useful in get_SCM_LB
   * for swapping out current_RB_coeffs.
   */
  std::vector<Number> saved_RB_coeffs;

  /**
   * The dimension of the reduced basis.
   */
  unsigned int n_bfs;

private:

  //----------- PRIVATE DATA MEMBERS -----------//


  /**
   * The current time-level.
   */
  unsigned int _k;

  /**
   * The total number of time-steps.
   */
  unsigned int _K;

  /**
   * The time-step size.
   */
  Number dt;

  /**
   * The number of training samples in time for each parameter
   * training sample in the training set. Needs to divide
   * _K+1 with no remainder.
   */
  unsigned int n_time_samples;

  /**
   * The function pointer for theta_c (for the trilinear form).
   */
  theta_q_fptr theta_c;

};

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
