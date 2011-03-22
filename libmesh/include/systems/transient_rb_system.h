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

#ifndef __transient_rb_system_h__
#define __transient_rb_system_h__

#include "transient_system.h"
#include "rb_system.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * TransientRBSystem extends RBSystem to add
 * functionality relevant in the time-dependent
 * case.
 *
 * @author David J. Knezevic 2009
 */

// ------------------------------------------------------------
// TransientRBSystem class definition

class TransientRBSystem : public TransientSystem<RBSystem>
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  TransientRBSystem (EquationSystems& es,
            const std::string& name,
            const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~TransientRBSystem ();

  /**
   * The type of system.
   */
  typedef TransientRBSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef TransientSystem<RBSystem> Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Perform a truth solve at the current parameter.
   */
  virtual Real truth_solve(int write_interval);

  /**
   * Train the reduced basis. Overloaded so that we can set the
   * flag compute_truth_projection_error to true so that the calls
   * to truth_solve during the basis construction will compute the
   * projection error. Other calls to truth_solve generally do not
   * need to perform these projection calculations.
   */
  virtual Real train_reduced_basis(const std::string& directory_name = "offline_data");

  /**
   * Function that indicates when to terminate the Greedy
   * basis training.
   */
  virtual bool greedy_termination_test(Real training_greedy_error, int count);

  /**
   * Clear the basis functions and all basis-function-dependent data.
   * Overloaded to also clear the M-representors.
   */
  virtual void clear_basis_function_dependent_data();

  /**
   * Assemble and store all the affine operators.
   * Overload to assemble the mass matrix operators.
   */
  virtual void assemble_all_affine_operators();
  
  /**
   * Override to assemble the L2 matrix as well.
   */
  virtual void assemble_misc_matrices();

  /**
   * Assemble the L2 matrix.
   */
  void assemble_L2_matrix(SparseMatrix<Number>* input_matrix);

  /**
   * Assemble the mass matrix at the current parameter 
   * and store it in input_matrix.
   */
  void assemble_mass_matrix(SparseMatrix<Number>* input_matrix);
  
  /**
   * Add the scaled mass matrix (assembled for the current parameter)
   * to input_matrix.
   */
  void add_scaled_mass_matrix(Number scalar,
                              SparseMatrix<Number>* input_matrix);

  /**
   * Perform a matrix-vector multiplication with the current mass matrix
   * and store the result in dest.
   */
  void mass_matrix_scaled_matvec(Number scalar,
                                 NumericVector<Number>& dest,
                                 NumericVector<Number>& arg);

  /**
   * Build a new TransientRBEvaluation object and add
   * it to the rb_evaluation_objects vector.
   */
  virtual void add_new_rb_evaluation_object();
  
  /**
   * Get Q_m, the number of terms in the affine
   * expansion for the mass operator.
   */
  virtual unsigned int get_Q_m() { return theta_q_m_vector.size(); }
  
  /**
   * Attach user-defined assembly routine
   * for the L2 matrix.
   */
  void attach_L2_assembly(affine_assembly_fptr L2_assembly);

  /**
   * Attach parameter-dependent function and user-defined assembly routine
   * for affine operator. Boundary assembly defaults to NULL here since
   * we typically only need interior assembly for the mass operator.
   */
  void attach_M_q(theta_q_fptr theta_q_m,
                  affine_assembly_fptr M_q_intrr_assembly,
                  affine_assembly_fptr M_q_bndry_assembly=NULL);

  /**
   * Assemble the q^th affine term of the mass matrix and store it in input_matrix.
   */
  void assemble_Mq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix);
  
  /**
   * Get a pointer to M_q.
   */
  SparseMatrix<Number>* get_M_q(unsigned int q);
  
  /**
   * Evaluate theta_q_m at the current parameter.
   */
  Number eval_theta_q_m(unsigned int q);

  /**
   * Override initialize_RB_system to also initialize
   * RB_ic_proj_rhs_all_N, if necessary.
   */
  virtual void initialize_RB_system(bool do_not_assemble);

  /**
   * Assemble the truth system in the transient linear case.
   */
  virtual void truth_assembly();

  /**
   * Get/set max_truth_solves, the maximum number of RB
   * truth solves we are willing to compute in the transient
   * case. Note in the steady state case max_truth_solves is
   * not needed since equivalent to Nmax.
   */
  int get_max_truth_solves() const                   { return max_truth_solves; }
  void set_max_truth_solves(int max_truth_solves_in) { this->max_truth_solves = max_truth_solves_in; }

  /**
   * Get/set dt, the time-step size.
   */
  Real get_dt() const        { return dt; }
  void set_dt(const Real dt) { this->dt = dt; }

  /**
   * Get/set POD_tol
   */
  Real get_POD_tol() const                { return POD_tol; }
  void set_POD_tol(const Real POD_tol_in) { this->POD_tol = POD_tol_in; }

  /**
   * Get/set euler_theta, parameter that determines
   * the temporal discretization.
   * euler_theta = 0   ---> Forward Euler
   * euler_theta = 0.5 ---> Crank-Nicolson
   * euler_theta = 1   ---> Backward Euler
   */
  Real get_euler_theta() const                 { return euler_theta; }
  void set_euler_theta(const Real euler_theta) { libmesh_assert((0. <= euler_theta ) && (euler_theta <= 1.));
                                                 this->euler_theta = euler_theta; }

  /**
   * Get/set the current time-level.
   */
  unsigned int get_time_level() const       { return _k; }
  void set_time_level(const unsigned int k) { libmesh_assert(_k <= _K); this->_k = k; }

  /**
   * Get/set K, the total number of time-steps.
   */
  unsigned int get_K() const       { return _K; }
  void set_K(const unsigned int K) { this->_K = K; }

  /**
   * Get/set M, the number of basis functions we add to the
   * RB space from each POD
   */
  unsigned int get_delta_N() const       { return delta_N; }
  void set_delta_N(const unsigned int delta_N) { this->delta_N = delta_N; }

  /**
   * Load the RB solution from the current time-level
   * into the libMesh solution vector.
   */
  virtual void load_RB_solution();

  /**
   * Get column k (i.e. the current time level) of temporal_data.
   * This gives access to the truth projection error data. If
   * the RB basis is empty, then this corresponds to the truth
   * solution data itself.
   */
  const NumericVector<Number>& get_error_temporal_data();

  /**
   * Compute the L2 projection of the initial condition
   * onto the RB space for 1 <= N <= RB_size and store
   * each projection in RB_initial_condition_matrix.
   */
  void update_RB_initial_condition_all_N();

  /**
   * Specifies the residual scaling on the numerator to
   * be used in the a posteriori error bound. Overload
   * in subclass in order to obtain the desired error bound.
   */
  virtual Real residual_scaling_numer(Real alpha_LB);

  /**
   * Specifies the residual scaling on the denominator to
   * be used in the a posteriori error bound. Overload
   * in subclass in order to obtain the desired error bound.
   */
  virtual Real residual_scaling_denom(Real alpha_LB);

  /**
   * Overload write_offline_data_to_files in order to
   * write out the mass matrix and initial condition
   * data as well.
   */
  virtual void write_offline_data_to_files(const std::string& directory_name = "offline_data",
                                           const RBDataIO io_flag = ALL_DATA);

  /**
   * Overload read_offline_data_from_files in order to
   * read in the mass matrix and initial condition
   * data as well.
   */
  virtual void read_offline_data_from_files(const std::string& directory_name = "offline_data",
                                            const RBDataIO io_flag = ALL_DATA);


  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The L2 matrix.
   */
  AutoPtr< SparseMatrix<Number> > L2_matrix;

  /**
   * Vector storing the Q_m matrices from the mass operator
   */
  std::vector< SparseMatrix<Number>* > M_q_vector;

  /**
   * The true error data for all time-levels from the
   * most recent RB_solve.
   */
//   std::vector< Number > true_error_all_k;

  /**
   * The truth outputs for all time-levels from the
   * most recent truth_solve.
   */
  std::vector< std::vector<Number> > truth_outputs_all_k;

  /**
   * Vectors storing the mass matrix representors.
   */
  std::vector< std::vector< NumericVector<Number>* > > M_q_representor;

  /**
   * Boolean flag to indicate whether we are using a non-zero initialization.
   * If we are, then an initialization function must be attached to the system.
   */
  bool nonzero_initialization;
  
  /**
   * Boolean flag that indicates whether we will compute the projection error
   * for the truth solution into the RB space (at every time level).
   * This typically only needs to true during a call to train_reduced_basis.
   */
  bool compute_truth_projection_error;

  /**
   * The filename of the file containing the initial
   * condition projected onto the truth mesh.
   */
  std::string init_filename;

protected:

  /**
   * Initializes the mesh-dependent data.
   */
  virtual void init_data ();

  /**
   * Read in the parameters from file and set up the system
   * accordingly.
   */
  virtual void process_parameters_file ();

  /**
   * Helper function that actually allocates all the data
   * structures required by this class.
   */
  virtual void allocate_data_structures();

  /**
   * This function imposes a truth initial condition,
   * defaults to zero initial condition if the flag
   * nonzero_initialization is true.
   */
  virtual void initialize_truth();
  
  /**
   * Override to use the L2 product matrix for output
   * dual norm solves for transient state problems.
   */
  virtual void assemble_matrix_for_output_dual_solves();

  /**
   * Initialize RB space by adding the truth initial condition
   * as the first RB basis function.
   */
  void add_IC_to_RB_space();

  /**
   * Add a new basis functions to the RB space. In the transient
   * case we first perform a POD of the time-dependent "truth"
   * and then add a certain number of POD modes to the reduced basis.
   */
  virtual void enrich_RB_space();

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
   * Set column k (i.e. the current time level) of temporal_data to the
   * difference between the current solution and the orthogonal
   * projection of the current solution onto the current RB space.
   */
  Number set_error_temporal_data();

  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * Time step size.
   */
  Real dt;

  /**
   * Parameter that defines the temporal discretization:
   * euler_theta = 0   ---> Forward Euler
   * euler_theta = 0.5 ---> Crank-Nicolson
   * euler_theta = 1   ---> Backward Euler
   */
  Real euler_theta;

  /**
   * The current time-level, 0 <= _k <= _K.
   */
  unsigned int _k;

  /**
   * Total number of time-steps.
   */
  unsigned int _K;

  /**
   * If positive, this tolerance determines the number of POD modes we
   * add to the space on a call to enrich_RB_space(). If negative, we add
   * delta_N POD modes.
   */
  Real POD_tol;

  /**
   * Maximum number of truth solves in the POD-Greedy. This can be
   * different from Nmax in the transient case since we may add
   * more than one basis function per truth solve.
   * If negative, it's ignored.
   */
  int max_truth_solves;
  
  /**
   * Function pointer for assembling the L2 matrix.
   */
  affine_assembly_fptr L2_assembly;

  /**
   * Vector of parameter-dependent functions for assembling
   * the mass operator.
   */
  std::vector<theta_q_fptr> theta_q_m_vector;
  
  /**
   * Vectors of function pointers for assembling the mass operator.
   */
  std::vector<affine_assembly_fptr> M_q_intrr_assembly_vector;
  std::vector<affine_assembly_fptr> M_q_bndry_assembly_vector;
  
  /**
   * The vector that stores the right-hand side for the initial
   * condition projections.
   */
  DenseVector<Number> RB_ic_proj_rhs_all_N;

private:

  /**
   * Private non-virtual helper function to encapsulate
   * the code to clear the basis-function-related data.
   */
  void clear_basis_helper();

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * Dense matrix to store the data that we use for the temporal POD.
   */
  std::vector< NumericVector<Number>* > temporal_data;

};

} // namespace libMesh

#endif
