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

#ifndef __rb_scm_system_h__
#define __rb_scm_system_h__

// Configuration data
#include "libmesh_config.h"

// Currently, the RBSCMSystem is only usable
// if SLEPc is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "eigen_system.h"
#include "rb_base.h"

/**
 * This class is part of the rbOOmit framework.
 *
 * RBSCMSystem implements the the Successive Constraint Method (SCM)
 * for computing rigorous lower bounds for stability constants.
 *
 * @author David J. Knezevic 2009
 */

// ------------------------------------------------------------
// RBSCMSystem class definition

class RBSCMSystem : public RBBase<EigenSystem>
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBSCMSystem (EquationSystems& es,
               const std::string& name,
               const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBSCMSystem ();

  /**
   * The type of system.
   */
  typedef RBSCMSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef RBBase<EigenSystem> Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Clear and resize the SCM data vectors. Overload
   * in subclass as necessary.
   */
  virtual void resize_SCM_vectors ();

  /**
   * Overload the solve function to solve the condensed
   * eigenproblem with constrained DOFs stripped out.
   */
  virtual void solve();

  /**
   * Overload get_eigenpair to retrieve the eigenpair for
   * the condensed eigensolve.
   */
  virtual std::pair<Real, Real> get_eigenpair(unsigned int i);

  /**
   * This function is called before truth eigensolves in
   * compute_SCM_bounding_box and evaluate_stability_constant.
   * Overload it to set specific properties to optimize
   * eigensolver performance. The argument refers to
   * the operator index in compute_SCM_bounding_box;
   * a negative value of the argument indicates we are
   * not performing a bounding box solve.
   */
  virtual void set_eigensolver_properties(int ) { }

  /**
   * Set the name of the associated RB system --- we need
   * this to load the (symmetrized) affine operators.
   */
  void set_RB_system_name(const std::string& name) { RB_system_name = name; }

  /**
   * Get/set SCM_eps: tolerance for SCM greedy.
   */
  Real get_SCM_eps() const          { return SCM_eps; }
  void set_SCM_eps(Real SCM_eps_in) { this->SCM_eps = SCM_eps_in; }

  /**
   * Get/set SCM_M: the number of nearby points to a given parameter
   * that we use to construct the SCM LP.
   */
  unsigned int get_SCM_M() const     { return SCM_M; }
  void set_SCM_M(unsigned int SCM_M_in) { this->SCM_M = SCM_M_in; }

  /**
   * Perform the SCM greedy algorithm to develop a lower bound
   * over the training set.
   */
  virtual void perform_SCM_greedy();

  /**
   * Evaluate single SCM lower bound.
   */
  virtual Number get_SCM_LB();

  /**
   * Evaluate single SCM upper bound.
   */
  virtual Number get_SCM_UB();

  /**
   * Get stability constraints (i.e. the values of coercivity/
   * inf-sup/stability constants at the parameter values chosen
   * during the greedy); we store one constraint for each element
   * of C_J.
   */
  Real get_C_J_stability_constraint(unsigned int j) const;

  /**
   * Get entries of SCM_UB_vector, which stores the
   * vector y, corresponding to the minimizing eigenvectors
   * for the elements of C_J.
   */
  Real get_SCM_UB_vector(unsigned int j, unsigned int q);

  /**
   * Get size of the set C_J.
   */
  unsigned int get_C_J_size() { return C_J.size(); }

  /**
   * Get entry of C_J.
   */
  std::vector<Real> get_C_J_entry(unsigned int j);

  /**
   * Get entry of C_J_stability_vector.
   */
  Real get_C_J_stability_value(unsigned int j) { return C_J_stability_vector[j]; }

  /**
   * Get B_min and B_max.
   */
  Real get_B_min(unsigned int i) const;
  Real get_B_max(unsigned int i) const;

  /**
   * Attach the deflation space defined by the specified vector, can
   * be useful in solving constrained eigenvalue problems.
   *
   * This function is called at the start of perform_SCM_greedy and by
   * default is does nothing. Overload in subclass to attach a specific
   * vector.
   */
  virtual void attach_deflation_space() {}

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

  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The filename of the text file from which we read in the
   * problem parameters. We use getpot.h to perform the reading.
   */
  std::string parameters_filename;

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Add the scaled symmetrized affine matrix from the
   * associated RBSystem to matrix_A.
   */
  virtual void add_scaled_symm_Aq(unsigned int q_a, Real scalar);

  /**
   * Copy over the matrix to store in matrix_B,
   * usually this is the mass or inner-product
   * matrix, but needs to be implemented in subclass.
   */
  virtual void load_matrix_B() { libmesh_not_implemented(); }

  /**
   * Compute the SCM bounding box.
   */
  virtual void compute_SCM_bounding_box();

  /**
   * Compute the stability constant for current_parameters
   * by solving a generalized eigenvalue problem over the
   * truth space.
   */
  virtual void evaluate_stability_constant();

  /**
   * Helper function to set current_parameters to the
   * specified parameters saved in C_J.
   */
  virtual void set_current_parameters_from_C_J(unsigned int C_J_index);

  /**
   * Compute upper and lower bounds for each SCM training point.
   * Return a pair containing the maximum SCM error, and the
   * index of the parameter in the training set at which the max
   * error is achieved.
   */
  virtual std::pair<unsigned int,Real> compute_SCM_bounds_on_training_set();

  /**
   * Enrich C_J by adding the element of SCM_training_samples
   * that has the largest gap between alpha_LB and alpha_LB.
   */
  virtual void enrich_C_J(unsigned int new_C_J_index);

  /**
   * Set B_min and B_max.
   */
  void set_B_min(unsigned int i, Real B_min_val);
  void set_B_max(unsigned int i, Real B_max_val);

  /**
   * Set stability constraints (i.e. the values of coercivity/
   * inf-sup/stability constants at the parameter values chosen
   * during the greedy); we store one constraint for each element
   * of C_J.
   */
  void set_C_J_stability_constraint(unsigned int j, Real stability_constraint_in);

  /**
   * Set entries of SCM_UB_vector, which stores the
   * vector y, corresponding to the minimizing eigenvectors
   * for the elements of C_J.
   */
  void set_SCM_UB_vector(unsigned int j, unsigned int q, Real y_q);

  /**
   * Compute the inner product between two vectors using the system's
   * matrix_B.
   */
  Real B_inner_product(const NumericVector<Number>& v, const NumericVector<Number>& w) const;

  /**
   * Compute the inner product between two vectors using
   * matrix Aq.
   */
  Real Aq_inner_product(unsigned int q,
                        const NumericVector<Number>& v,
                        const NumericVector<Number>& w);

  /**
   * Helper function to compute the distance between
   * two parameters, using the Euclidean norm.
   *
   * We make this virtual because in some cases (e.g. Burgers)
   * we use a weighted distance.
   */
  virtual Real param_dist(const std::vector<Real>& mu_1, const std::vector<Real>& mu_2);

  /**
   * Helper function which provides an error
   * indicator to be used in the SCM greedy.
   * Overload in subclasses to specialize behavior.
   */
  virtual Real SCM_greedy_error_indicator(Real LB, Real UB) { return fabs(UB-LB)/fabs(UB); }

  /**
   * Helper function to save current_parameters in
   * saved_parameters.
   */
  virtual void save_current_parameters();

  /**
   * Helper functiont to (re)load current_parameters
   * from saved_parameters.
   */
  virtual void reload_current_parameters();


  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * The condensed matrices, with Dirichlet dofs removed.
   */
  AutoPtr< SparseMatrix<Number> > condensed_matrix_A;
  AutoPtr< SparseMatrix<Number> > condensed_matrix_B;

  /**
   * SCM tolerance, where SCM_eps \in (0,1).
   */
  Real SCM_eps;

  /**
   * Parameter used in SCM algorithm.
   */
  unsigned int SCM_M;

  /**
   * B_min, B_max define the bounding box.
   */
  std::vector<Real> B_min;
  std::vector<Real> B_max;

  /**
   * Vector storing the greedily selected parameters
   * during SCM training.
   */
  std::vector< std::vector<Real> > C_J;

  /**
   * Vector storing the (truth) stability values
   * at the parameters in C_J.
   */
  std::vector<Real> C_J_stability_vector;

  /**
   * This matrix stores the infimizing vectors
   * y_1(\mu),...,y_Q_a(\mu), for each \mu in
   * C_J, which are used in computing the SCM
   * upper bounds.
   */
  std::vector< std::vector<Real> > SCM_UB_vectors;

  /**
   * Vector in which to save a parameter set. Useful
   * in get_SCM_LB, for example.
   */
  std::vector<Real> saved_parameters;

  /**
   * The name of the associated RB system.
   */
  std::string RB_system_name;

private:

};

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
