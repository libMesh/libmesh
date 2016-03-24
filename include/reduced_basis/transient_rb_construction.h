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

#ifndef LIBMESH_TRANSIENT_RB_CONSTRUCTION_H
#define LIBMESH_TRANSIENT_RB_CONSTRUCTION_H

// rbOOmit includes
#include "libmesh/rb_construction.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/rb_temporal_discretization.h"

// libMesh includes
#include "libmesh/transient_system.h"

// C++ includes

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * TransientRBConstruction extends RBConstruction to add
 * functionality relevant in the time-dependent case.
 *
 * We can handle time controls on the RHS as h(t)*f(x,\f$ \mu \f$).
 * See Martin Grepl's thesis for more details.
 *
 * \author David J. Knezevic
 * \date 2009
 */
class TransientRBConstruction : public TransientSystem<RBConstruction>, public RBTemporalDiscretization
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  TransientRBConstruction (EquationSystems & es,
                           const std::string & name,
                           const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~TransientRBConstruction ();

  /**
   * The type of system.
   */
  typedef TransientRBConstruction sys_type;

  /**
   * The type of the parent.
   */
  typedef TransientSystem<RBConstruction> Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () libmesh_override;

  /**
   * Allocate all the data structures necessary for the construction
   * stage of the RB method. This function also performs
   * matrix and vector assembly of the "truth" affine expansion.
   *
   * Override to check that theta and assembly expansions are consistently
   * sized.
   */
  virtual void initialize_rb_construction(bool skip_matrix_assembly=false,
                                          bool skip_vector_assembly=false) libmesh_override;

  /**
   * Perform a truth solve at the current parameter.
   */
  virtual Real truth_solve(int write_interval) libmesh_override;

  /**
   * Train the reduced basis. Overloaded so that we can set the
   * flag compute_truth_projection_error to true so that the calls
   * to truth_solve during the basis construction will compute the
   * projection error. Other calls to truth_solve generally do not
   * need to perform these projection calculations.
   */
  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true) libmesh_override;

  /**
   * Read in the parameters from file and set up the system
   * accordingly.
   */
  virtual void process_parameters_file (const std::string & parameters_filename) libmesh_override;

  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info() libmesh_override;

  /**
   * Function that indicates when to terminate the Greedy
   * basis training.
   */
  virtual bool greedy_termination_test(
    Real abs_greedy_error, Real initial_greedy_error, int count) libmesh_override;

  /**
   * Assemble and store all the affine operators.
   * Overload to assemble the mass matrix operators.
   */
  virtual void assemble_all_affine_operators() libmesh_override;

  /**
   * Override to assemble the L2 matrix as well.
   */
  virtual void assemble_misc_matrices() libmesh_override;

  /**
   * Assemble the L2 matrix.
   */
  void assemble_L2_matrix(SparseMatrix<Number> * input_matrix,
                          bool apply_dirichlet_bc=true);

  /**
   * Assemble the mass matrix at the current parameter
   * and store it in input_matrix.
   */
  void assemble_mass_matrix(SparseMatrix<Number> * input_matrix);

  /**
   * Add the scaled mass matrix (assembled for the current parameter)
   * to input_matrix.
   */
  void add_scaled_mass_matrix(Number scalar,
                              SparseMatrix<Number> * input_matrix);

  /**
   * Perform a matrix-vector multiplication with the current mass matrix
   * and store the result in dest.
   */
  void mass_matrix_scaled_matvec(Number scalar,
                                 NumericVector<Number> & dest,
                                 NumericVector<Number> & arg);

  /**
   * Set the L2 object.
   */
  void set_L2_assembly(ElemAssembly & L2_assembly_in);

  /**
   * @return a reference to the L2 assembly object
   */
  ElemAssembly & get_L2_assembly();

  /**
   * Assemble the q^th affine term of the mass matrix and store it in input_matrix.
   */
  void assemble_Mq_matrix(unsigned int q,
                          SparseMatrix<Number> * input_matrix,
                          bool apply_dirichlet_bc=true);

  /**
   * Get a pointer to M_q.
   */
  SparseMatrix<Number> * get_M_q(unsigned int q);

  /**
   * Get a pointer to non_dirichlet_M_q.
   */
  SparseMatrix<Number> * get_non_dirichlet_M_q(unsigned int q);

  /**
   * Get a map that stores pointers to all of the matrices.
   */
  virtual void get_all_matrices(std::map<std::string, SparseMatrix<Number> *> & all_matrices) libmesh_override;

  /**
   * Assemble the truth system in the transient linear case.
   */
  virtual void truth_assembly() libmesh_override;

  /**
   * Get/set max_truth_solves, the maximum number of RB
   * truth solves we are willing to compute in the transient
   * case. Note in the steady state case max_truth_solves is
   * not needed since equivalent to Nmax.
   */
  int get_max_truth_solves() const                   { return max_truth_solves; }
  void set_max_truth_solves(int max_truth_solves_in) { this->max_truth_solves = max_truth_solves_in; }

  /**
   * Get/set POD_tol
   */
  Real get_POD_tol() const                { return POD_tol; }
  void set_POD_tol(const Real POD_tol_in) { this->POD_tol = POD_tol_in; }

  /**
   * Set delta_N, the number of basis functions we add to the
   * RB space from each POD
   */
  void set_delta_N(const unsigned int new_delta_N) { this->delta_N = new_delta_N; }

  /**
   * Load the RB solution from the current time-level
   * into the libMesh solution vector.
   */
  virtual void load_rb_solution() libmesh_override;

  /**
   * Get the column of temporal_data corresponding to the current time level.
   * This gives access to the truth projection error data. If
   * the RB basis is empty, then this corresponds to the truth
   * solution data itself.
   */
  const NumericVector<Number> & get_error_temporal_data();

  /**
   * Compute the L2 projection of the initial condition
   * onto the RB space for 1 <= N <= RB_size and store
   * each projection in RB_initial_condition_matrix.
   */
  void update_RB_initial_condition_all_N();

  /**
   * Write out all the Riesz representor data to files. Override
   * to write out transient data too.
   */
  virtual void write_riesz_representors_to_files(const std::string & riesz_representors_dir,
                                                 const bool write_binary_residual_representors) libmesh_override;

  /**
   * Write out all the Riesz representor data to files. Override
   * to read in transient data too.
   */
  virtual void read_riesz_representors_from_files(const std::string & riesz_representors_dir,
                                                  const bool write_binary_residual_representors) libmesh_override;


  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The L2 matrix.
   */
  UniquePtr< SparseMatrix<Number> > L2_matrix;

  /**
   * The L2 matrix without Dirichlet conditions enforced.
   * (This is only computed if store_non_dirichlet_operators == true.)
   */
  UniquePtr< SparseMatrix<Number> > non_dirichlet_L2_matrix;

  /**
   * Vector storing the Q_m matrices from the mass operator
   */
  std::vector<SparseMatrix<Number> *> M_q_vector;

  /**
   * We sometimes also need a second set of M_q matrices
   * that do not have the Dirichlet boundary conditions
   * enforced.
   */
  std::vector<SparseMatrix<Number> *> non_dirichlet_M_q_vector;

  /**
   * The truth outputs for all time-levels from the
   * most recent truth_solve.
   */
  std::vector<std::vector<Number> > truth_outputs_all_k;

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
   * Helper function that actually allocates all the data
   * structures required by this class.
   */
  virtual void allocate_data_structures() libmesh_override;

  /**
   * Override assemble_affine_expansion to also initialize
   * RB_ic_proj_rhs_all_N, if necessary.
   */
  virtual void assemble_affine_expansion(bool skip_matrix_assembly,
                                         bool skip_vector_assembly) libmesh_override;

  /**
   * This function imposes a truth initial condition,
   * defaults to zero initial condition if the flag
   * nonzero_initialization is true.
   */
  virtual void initialize_truth();

  /**
   * Override to return the L2 product matrix for output
   * dual norm solves for transient state problems.
   */
  virtual SparseMatrix<Number> & get_matrix_for_output_dual_solves() libmesh_override;

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
  virtual void enrich_RB_space() libmesh_override;

  /**
   * Update the system after enriching the RB space.
   */
  virtual void update_system() libmesh_override;

  /**
   * Compute the reduced basis matrices for the current basis.
   */
  virtual void update_RB_system_matrices() libmesh_override;

  /**
   * Compute the terms that are combined `online'
   * to determine the dual norm of the residual.
   */
  virtual void update_residual_terms(bool compute_inner_products) libmesh_override;

  /**
   * Set column k (i.e. the current time level) of temporal_data to the
   * difference between the current solution and the orthogonal
   * projection of the current solution onto the current RB space.
   */
  Number set_error_temporal_data();

  //----------- PROTECTED DATA MEMBERS -----------//

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
  ElemAssembly * L2_assembly;

  /**
   * The vector that stores the right-hand side for the initial
   * condition projections.
   */
  DenseVector<Number> RB_ic_proj_rhs_all_N;

private:

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * Dense matrix to store the data that we use for the temporal POD.
   */
  std::vector< NumericVector<Number> * > temporal_data;
};

} // namespace libMesh

#endif // LIBMESH_TRANSIENT_RB_CONSTRUCTION_H
