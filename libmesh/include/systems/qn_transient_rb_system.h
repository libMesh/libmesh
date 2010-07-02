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

#ifndef __qn_transient_rb_system_h__
#define __qn_transient_rb_system_h__

// Configuration data
#include "libmesh_config.h"

// This class requires QNTransientSCMSystem, which is not
// defined if SLEPc is not present.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "transient_rb_system.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * QNTransientRBSystem extends TransientRBSystem
 * in order to implement the reduced basis
 * method for time-dependent quadratically
 * non-linear problems.
 *
 * @author David J. Knezevic 2009
 */

// ------------------------------------------------------------
// QNTransientRBSystem class definition

class QNTransientRBSystem : public TransientRBSystem
{
public:


  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  QNTransientRBSystem (EquationSystems& es,
                       const std::string& name,
                       const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~QNTransientRBSystem ();

  /**
   * The type of system.
   */
  typedef QNTransientRBSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef TransientRBSystem Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Clear the basis functions and all basis-function-dependent data.
   * Overload in subclasses to clear any extra data.
   * Overloaded to also clear the C-dependent data.
   */
  virtual void clear_basis_function_dependent_data();

  /**
   * Attach user-defined assembly routine
   * for the trilinear form operators.
   */
  void attach_C(theta_q_fptr theta_q_c, affine_assembly_fptr C_assembly);

  /**
   * Attach user-defined assembly routine
   * for the truth element interior assembly.
   */
  void attach_truth_assembly(affine_assembly_fptr truth_intrr_assembly,
                             affine_assembly_fptr truth_bndry_assembly);

  /**
   * Evaluate theta_q_c at the current parameter.
   */
  Real eval_theta_c();

  /**
   * Perform a truth solve at the current parameter.
   * Overload to solve the nonlinear truth system
   * using Newton's method.
   */
  virtual Real truth_solve(int write_interval);

  /**
   * Perform online solve for current_params
   * with the N basis functions. Overload this
   * to solve the nonlinear RB system using
   * Newton's method.
   */
  virtual Real RB_solve(unsigned int N);

  /**
   * Set the nonlinear tolerance for Newton's
   * method for both the truth and RB solves.
   */
  void set_nonlinear_tolerance(Number nonlinear_tolerance)
    { this->nonlinear_tolerance = nonlinear_tolerance; }
    
  /**
   * Get the nonlinear tolerance for Newton's
   * method for both the truth and RB solves.
   */
  Real get_nonlinear_tolerance()
    { return nonlinear_tolerance; }

  /**
   * Set the maximum number of Newton steps
   * for both the truth and RB solves.
   */
  void set_n_newton_steps(unsigned int n_newton_steps_in)
    { this->n_newton_steps = n_newton_steps_in; }
    
  /**
   * Get the maximum number of Newton steps
   * for both the truth and RB solves.
   */
  unsigned int get_n_newton_steps()
    { return n_newton_steps; }

  /**
   * Get a pointer to C_n, the n^th affine operator
   * (where 1 <= n <= Nmax) that arises from the
   * quadratic nonlinearity.
   */
  SparseMatrix<Number>* get_C_n(unsigned int n);

  /**
   * Get the nominal stability factor lower bound.
   * This will typically be overloaded in a subclass
   * to be a function of the current parameters.
   */
  virtual Number get_nominal_rho_LB() { return 1.; }

  /**
   * Assemble matrix C_n, i.e. c(\xi_n, \phi_j, \phi_i), where
   * \xi_n is the n^th RB basis function.
   */
  void assemble_C_n_matrix(unsigned int n, SparseMatrix<Number>* input_matrix);

  /**
   * Add a scaled version of matrix C_n to input_matrix. If symmetrize==true,
   * then we add symmetrized C_n.
   */
  void add_scaled_Cn(Real scalar, unsigned int n, SparseMatrix<Number>* input_matrix, bool symmetrize);

  /**
   * Assemble matrix C (the assembly depends on the value of current_newton_iterate)
   * and add to input_matrix. If symmetrize==true,
   * then we add symmetrized C.
   */
  void add_scaled_current_C(Real scalar, SparseMatrix<Number>* input_matrix, bool symmetrize);

  /**
   * Overload this function to set current_newton_iterate = vec.
   */
  virtual void set_context_solution_vec(NumericVector<Number>& vec);

  /**
   * Compute the trilinear form operators for _all_ basis functions.
   * This function can be used to initialize the system after we have
   * read in a set of basis functions.
   */
  void update_all_trilinear_form_operators();

  /**
   * Overload this function to write out the extra
   * Offline data specific to quadratically nonlinear
   * problems.
   */
  virtual void write_offline_data_to_files(const std::string& directory_name = "offline_data");

  /**
   * Overload this function to read in the extra
   * Offline data specific to quadratically nonlienar
   * problems.
   */
  virtual void read_offline_data_from_files(const std::string& directory_name = "offline_data");


  //----------- PUBLIC DATA MEMBERS -----------//


  /**
   * Vector to store the current Newton iterate for truth_solve().
   */
  AutoPtr< NumericVector<Number> > current_newton_iterate;

  /**
   * Storage for the reduced basis trilinear form that arises
   * from the quadratic nonlinearity.
   */
  std::vector< std::vector< std::vector<Number> > > RB_trilinear_form;

  /**
   * Boolean flag to indicate whether we should use a nominal value
   * of rho_LB (i.e. lower bound for the stability constant) or
   * compute a lower bound using the SCM.
   */
  bool use_nominal_rho_LB;

  /**
   * These vectors store the upper and lower bounds for the stability
   * constant from the most recent RB_solve.
   */
  std::vector<Real> rho_LB_vector;
  std::vector<Real> rho_UB_vector;

protected:

  /**
   * Initializes the member data fields associated with
   * the system.
   */
  virtual void init_data ();

  /**
   * Builds a RBContext object with enough information to do
   * evaluations on each element. Overload to create a
   * QNTransientRBContext.
   */
  virtual AutoPtr<RBContext> build_context();

  /**
   * Update the system after enriching the RB space.
   */
  virtual void update_system();

  /**
   * Compute the reduced basis matrices for the current basis.
   */
  virtual void update_RB_system_matrices();

  /**
   * Compute the terms that are combined `online'
   * to determine the dual norm of the residual.
   */
  virtual void update_residual_terms(bool compute_inner_products);

  /**
   * Compute the operators arising from the trilinear form
   * for the recently added basis functions. This assembly loop
   * calls an RBContext-based assembly function.
   */
  void update_trilinear_form_operators();

  /**
   * Imposes truth initial condition.
   */
  virtual void initialize_truth();

  /**
   * Assemble the truth matrix and right-hand side
   * for current_parameters. Overload to add
   * the term due to the quadratic nonlinearity.
   */
  virtual void truth_assembly();

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution. Overloaded to handle the
   * quadratic nonlinearity as well.
   */
  virtual Number compute_residual_dual_norm(const unsigned int N);

  /**
   * Get the SCM lower bound for the current parameters
   * from the associated SCM system.
   */
  virtual Real get_SCM_lower_bound();

  /**
   * Get the SCM upper bound for the current parameters
   * from the associated SCM system.
   */
  virtual Real get_SCM_upper_bound();

private:

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * The nonlinear tolerance for the Newton iteration.
   */
  Number nonlinear_tolerance;

  /**
   * The maximum number of Newton iterations.
   */
  unsigned int n_newton_steps;

  /**
   * Parameters for the Eisenstat linear tolerance formula.
   */
  bool use_eisenstat;   /* if false, a fixed tolerance of eisenstat_eta0 will be used */
  Real eisenstat_eta0;  /* 0.1, 0.01, 1e-3 */
  Real eisenstat_gamma; /* .9, 1 */
  Real eisenstat_alpha; /* 2, 0.5*(1.+std::sqrt(5.)) */

  /**
   * We store Nmax trilinear form operators to speed up RB matrix assembly.
   */
  std::vector< SparseMatrix<Real>* > C_n_vector;

  /**
   * Vectors storing the representors for the residual terms arising
   * from the quadratic nonlinearity.
   */
  std::vector< std::vector< NumericVector<Number>* > > C_representor;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   */
  std::vector< std::vector< std::vector<Number> > > Fq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Mq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > Aq_C_representor_norms;
  std::vector< std::vector< std::vector< std::vector<Number> > > > C_C_representor_norms;

  /**
   * Function pointer for assembling the trilinear
   * form operators.
   */
  affine_assembly_fptr C_assembly;

  /**
   * Function pointer for assembling the truth
   * matrix and right-hand side: element interior terms.
   */
  affine_assembly_fptr truth_intrr_assembly;

  /**
   * Function pointer for assembling the truth
   * matrix and right-hand side: boundary terms.
   */
  affine_assembly_fptr truth_bndry_assembly;

  /**
   * The function pointer for theta_c (for the trilinear form).
   */
  theta_q_fptr theta_c;

};

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK

#endif
