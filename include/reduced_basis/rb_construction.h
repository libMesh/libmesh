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

#ifndef LIBMESH_RB_CONSTRUCTION_H
#define LIBMESH_RB_CONSTRUCTION_H

// rbOOmit includes
#include "libmesh/rb_construction_base.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dg_fem_context.h"
#include "libmesh/dirichlet_boundaries.h"

// C++ includes

namespace libMesh
{

class RBThetaExpansion;
class RBAssemblyExpansion;
class RBEvaluation;
class ElemAssembly;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBConstruction implements the Construction stage
 * of the certified reduced basis method for
 * steady-state elliptic parametrized PDEs.
 *
 * \author David J. Knezevic
 * \date 2009
 */
class RBConstruction : public RBConstructionBase<LinearImplicitSystem>
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBConstruction (EquationSystems & es,
                  const std::string & name,
                  const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - Destructor is defaulted out-of-line
   */
  RBConstruction (RBConstruction &&) = default;
  RBConstruction (const RBConstruction &) = delete;
  RBConstruction & operator= (const RBConstruction &) = delete;
  RBConstruction & operator= (RBConstruction &&) = delete;
  virtual ~RBConstruction ();

  /**
   * The type of system.
   */
  typedef RBConstruction sys_type;

  /**
   * Assembles & solves the linear system A*x=b for the specified
   * matrix \p input_matrix and right-hand side \p rhs.
   */
  virtual void solve_for_matrix_and_rhs (LinearSolver<Number> & input_solver,
                                         SparseMatrix<Number> & input_matrix,
                                         NumericVector<Number> & input_rhs);

  /**
   * Set the RBEvaluation object.
   */
  void set_rb_evaluation(RBEvaluation & rb_eval_in);

  /**
   * Get a reference to the RBEvaluation object.
   */
  RBEvaluation & get_rb_evaluation();

  /**
   * \returns \p true if rb_eval is initialized. False, otherwise.
   */
  bool is_rb_eval_initialized() const;

  /**
   * Get a reference to the RBThetaExpansion object that
   * that belongs to rb_eval.
   */
  RBThetaExpansion & get_rb_theta_expansion();

  /**
   * Set the rb_assembly_expansion object.
   */
  void set_rb_assembly_expansion(RBAssemblyExpansion & rb_assembly_expansion_in);

  /**
   * \returns A reference to the rb_assembly_expansion object
   */
  RBAssemblyExpansion & get_rb_assembly_expansion();

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * The type of the parent.
   */
  typedef RBConstructionBase<LinearImplicitSystem> Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () override;

  /**
   * \returns A string indicating the type of the system.
   */
  virtual std::string system_type () const override;

  /**
   * Perform a "truth" solve, i.e. solve the finite element system at
   * at the parameters currently set in the system. This is used
   * extensively in training the reduced basis, since "truth snapshots"
   * are employed as basis functions.
   */
  virtual Real truth_solve(int plot_solution);

  /**
   * Train the reduced basis. This can use different approaches, e.g. Greedy
   * or POD, which are chosen using the RB_training_type member variable.
   *
   * In the case that we use Greedy training, this function returns the
   * final maximum a posteriori error bound on the training set.
   */
  virtual Real train_reduced_basis(const bool resize_rb_eval_data=true);

  /**
   * Train the reduced basis using the "Greedy algorithm."
   *
   * Each stage of the Greedy algorithm involves solving the reduced basis
   * over a large training set and selecting the parameter at which the
   * reduced basis error bound is largest, then performing a truth_solve
   * at that parameter and enriching the reduced basis with the corresponding
   * snapshot.
   *
   * \p resize_rb_eval_data is a boolean flag to indicate whether or not we
   * call rb_eval->resize_data_structures(Nmax). True by default, but we may
   * set it to false if, for example, we are continuing from a previous
   * training run and don't want to clobber the existing rb_eval data.
   *
   * \returns The final maximum a posteriori error bound on the training set.
   */
  Real train_reduced_basis_with_greedy(const bool resize_rb_eval_data);

  /**
   * This function computes one basis function for each rhs term. This is
   * useful in some cases since we can avoid doing a full greedy if we know
   * that we do not have any "left-hand side" parameters, for example.
   */
  void enrich_basis_from_rhs_terms(const bool resize_rb_eval_data=true);

  /**
   * Train the reduced basis using Proper Orthogonal Decomposition (POD).
   * This is an alternative to train_reduced_basis(), which uses the RB greedy
   * algorithm. In contrast to the RB greedy algorithm, POD requires us to
   * perform truth solves at all training samples, which can be computationally
   * intensive.
   *
   * The main advantage of using POD is that it does not rely on the RB error
   * indicator. The RB error indicator typically stagnates due to rounding
   * error at approximately square-root of machine precision, since it involves
   * taking the square-root of a sum of terms that cancel. This error indicator
   * stagnation puts a limit on the accuracy level that can be achieved with
   * the RB greedy algorithm, so for cases where we need higher accuracy, the
   * POD approach is a good alternative.
   */
  void train_reduced_basis_with_POD();

  /**
   * (i) Compute the a posteriori error bound for each set of parameters
   * in the training set, (ii) set current_parameters to the parameters that
   * maximize the error bound, and (iii) return the maximum error bound.
   */
  virtual Real compute_max_error_bound();

  /**
   * Return the parameters chosen during the i^th step of
   * the Greedy algorithm.
   */
  const RBParameters & get_greedy_parameter(unsigned int i);

  /**
   * Get/set the relative tolerance for the basis training.
   */
  void set_rel_training_tolerance(Real new_training_tolerance)
  {this->rel_training_tolerance = new_training_tolerance; }
  Real get_rel_training_tolerance() { return rel_training_tolerance; }

  /**
   * Get/set the absolute tolerance for the basis training.
   */
  void set_abs_training_tolerance(Real new_training_tolerance)
  {this->abs_training_tolerance = new_training_tolerance; }
  Real get_abs_training_tolerance() { return abs_training_tolerance; }

  /**
   * Get/set the boolean to indicate if we normalize the RB error in the greedy.
   */
  void set_normalize_rb_bound_in_greedy(bool normalize_rb_bound_in_greedy_in)
  {this->normalize_rb_bound_in_greedy = normalize_rb_bound_in_greedy_in; }
  bool get_normalize_rb_bound_in_greedy() { return normalize_rb_bound_in_greedy; }

  /**
   * Get/set the string that determines the training type.
   */
  void set_RB_training_type(const std::string & RB_training_type_in);
  const std::string & get_RB_training_type() const;

  /**
   * Get/set Nmax, the maximum number of RB
   * functions we are willing to compute.
   */
  unsigned int get_Nmax() const    { return Nmax; }
  virtual void set_Nmax(unsigned int Nmax);

  /**
   * Load the i^th RB function into the RBConstruction
   * solution vector.
   */
  virtual void load_basis_function(unsigned int i);

  /**
   * Load the RB solution from the most recent solve with rb_eval
   * into this system's solution vector.
   */
  virtual void load_rb_solution();

  /**
   * The slow (but simple, non-error prone) way to compute the residual dual norm.
   * Useful for error checking.
   */
  Real compute_residual_dual_norm_slow(const unsigned int N);

  /**
   * Get a pointer to inner_product_matrix. Accessing via this
   * function, rather than directly through the class member allows
   * us to do error checking (e.g. inner_product_matrix is not
   * defined in low-memory mode).
   */
  SparseMatrix<Number> * get_inner_product_matrix();

  /**
   * Get the non-Dirichlet (or more generally no-constraints) version
   * of the inner-product matrix. This is useful for performing multiplications
   * on vectors that already have constraints enforced.
   */
  SparseMatrix<Number> * get_non_dirichlet_inner_product_matrix();

  /**
   * Get the non-Dirichlet inner-product matrix if it's available,
   * otherwise get the inner-product matrix with constraints.
   */
  SparseMatrix<Number> * get_non_dirichlet_inner_product_matrix_if_avail();

  /**
   * Get a pointer to Aq.
   */
  SparseMatrix<Number> * get_Aq(unsigned int q);

  /**
   * Get a pointer to non_dirichlet_Aq.
   */
  SparseMatrix<Number> * get_non_dirichlet_Aq(unsigned int q);

  /**
   * Get a pointer to non_dirichlet_Aq if it's available, otherwise
   * get Aq.
   */
  SparseMatrix<Number> * get_non_dirichlet_Aq_if_avail(unsigned int q);

  /**
   * Allocate all the data structures necessary for the construction
   * stage of the RB method. This function also performs
   * matrix and vector assembly of the "truth" affine expansion.
   *
   * We can optionally skip the matrix or vector assembly steps by setting
   * skip_matrix_assembly = true, or skip_vector_assembly = true,
   * respectively.
   */
  virtual void initialize_rb_construction(bool skip_matrix_assembly=false,
                                          bool skip_vector_assembly=false);

  /**
   * Get a pointer to Fq.
   */
  NumericVector<Number> * get_Fq(unsigned int q);

  /**
   * Get a pointer to non-Dirichlet Fq.
   */
  NumericVector<Number> * get_non_dirichlet_Fq(unsigned int q);

  /**
   * Get a pointer to non_dirichlet_Fq if it's available, otherwise
   * get Fq.
   */
  NumericVector<Number> * get_non_dirichlet_Fq_if_avail(unsigned int q);

  /**
   * Get a pointer to the n^th output.
   */
  NumericVector<Number> * get_output_vector(unsigned int n, unsigned int q_l);

  /**
   * Get a pointer to non-Dirichlet output vector.
   */
  NumericVector<Number> * get_non_dirichlet_output_vector(unsigned int n, unsigned int q_l);

  /**
   * Get a map that stores pointers to all of the matrices.
   */
  virtual void get_all_matrices(std::map<std::string, SparseMatrix<Number> *> & all_matrices);

  /**
   * Get a map that stores pointers to all of the vectors.
   */
  virtual void get_all_vectors(std::map<std::string, NumericVector<Number> *> & all_vectors);

  /**
   * Get a map that stores pointers to all of the vectors.
   */
  virtual void get_output_vectors(std::map<std::string, NumericVector<Number> *> & all_vectors);

  /**
   * Assemble the matrices and vectors for this system.
   * Optionally skip matrix or vector assembly (e.g. we may want to
   * read data in from disk instead).
   */
  virtual void assemble_affine_expansion(bool skip_matrix_assembly,
                                         bool skip_vector_assembly);

  /**
   * Assemble the inner product matrix and store it in input_matrix.
   */
  void assemble_inner_product_matrix(SparseMatrix<Number> * input_matrix, bool apply_dof_constraints=true);

  /**
   * Assemble the q^th affine matrix and store it in input_matrix.
   */
  void assemble_Aq_matrix(unsigned int q, SparseMatrix<Number> * input_matrix, bool apply_dof_constraints=true);

  /**
   * Assemble the q^th affine vector and store it in input_matrix.
   */
  void assemble_Fq_vector(unsigned int q, NumericVector<Number> * input_vector, bool apply_dof_constraints=true);

  /**
   * Add the scaled q^th affine matrix to input_matrix. If symmetrize==true, then
   * we symmetrize Aq before adding it.
   */
  void add_scaled_Aq(Number scalar, unsigned int q_a,
                     SparseMatrix<Number> * input_matrix,
                     bool symmetrize);

  /**
   * Write out all the Riesz representor data to files.
   * \p residual_representors_dir specifies which directory to write to.
   * \p write_binary_residual_representors specifies whether we write binary or ASCII data.
   */
  virtual void write_riesz_representors_to_files(const std::string & riesz_representors_dir,
                                                 const bool write_binary_residual_representors);

  /**
   * Read in all the Riesz representor data from files.
   * \p directory_name specifies which directory to read from.
   * \p io_flag specifies whether we read in all data, or only
   * a basis (in)dependent subset.
   */
  virtual void read_riesz_representors_from_files(const std::string & riesz_representors_dir,
                                                  const bool write_binary_residual_representors);


  /**
   * This function computes all of the residual representors, can be useful
   * when restarting a basis training computation.
   * If \p compute_inner_products is false, we just compute the residual Riesz
   * representors, whereas if true, we also compute all the corresponding inner
   * product terms.
   */
  virtual void recompute_all_residual_terms(const bool compute_inner_products=true);

  /**
   * Read in from the file specified by \p parameters_filename
   * and set the this system's member variables accordingly.
   */
  virtual void process_parameters_file(const std::string & parameters_filename);

  /**
   * Set the state of this RBConstruction object based on the arguments
   * to this function.
   */
  void set_rb_construction_parameters(unsigned int n_training_samples_in,
                                      bool deterministic_training_in,
                                      unsigned int training_parameters_random_seed_in,
                                      bool quiet_mode_in,
                                      unsigned int Nmax_in,
                                      Real rel_training_tolerance_in,
                                      Real abs_training_tolerance_in,
                                      bool normalize_rb_error_bound_in_greedy_in,
                                      const std::string & RB_training_type_in,
                                      RBParameters mu_min_in,
                                      RBParameters mu_max_in,
                                      std::map<std::string, std::vector<Real>> discrete_parameter_values_in,
                                      std::map<std::string,bool> log_scaling,
                                      std::map<std::string, std::vector<Number>> * training_sample_list=nullptr);

  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info();

  /**
   * Print out a matrix that shows the orthogonality of the RB basis functions.
   * This is a helpful debugging tool, e.g. orthogonality can be degraded
   * due to finite precision arithmetic.
   */
  void print_basis_function_orthogonality();

  /**
   * Get delta_N, the number of basis functions we
   * add to the RB space per iteration of the greedy
   * algorithm. For steady-state systems, this should
   * be 1, but can be more than 1 for time-dependent
   * systems.
   */
  unsigned int get_delta_N() const { return delta_N; }

  /**
   * Set the rb_assembly_expansion object.
   */
  void set_inner_product_assembly(ElemAssembly & inner_product_assembly_in);

  /**
   * \returns A reference to the inner product assembly object
   */
  ElemAssembly & get_inner_product_assembly();

  /**
   * Specify the coefficients of the A_q operators to be used in the
   * energy inner-product.
   */
  void set_energy_inner_product(const std::vector<Number> & energy_inner_product_coeffs_in);

  /**
   * It is sometimes useful to be able to zero vector entries
   * that correspond to constrained dofs.
   */
  void zero_constrained_dofs_on_vector(NumericVector<Number> & vector);

  /**
   * @return true if the most recent truth solve gave a zero solution.
   */
  virtual bool check_if_zero_truth_solve();

#ifdef LIBMESH_ENABLE_DIRICHLET
  /**
   * It's helpful to be able to generate a DirichletBoundary that stores a ZeroFunction in order
   * to impose Dirichlet boundary conditions.
   */
  static std::unique_ptr<DirichletBoundary> build_zero_dirichlet_boundary_object();
#endif

  /**
   * Setter for the flag determining if convergence should be
   * checked after each solve.
   */
  void set_convergence_assertion_flag(bool flag);

  /**
   * Get/set flag to pre-evaluate the theta functions
   */
  bool get_preevaluate_thetas_flag() const;
  void set_preevaluate_thetas_flag(bool flag);


  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * Vector storing the values of the error bound
   * for each parameter in the training set --- the
   * parameter giving the largest error bound is
   * chosen for the next snapshot in the Greedy basis
   * training.
   */
  std::vector<Real> training_error_bounds;

  /**
   * We store an extra linear solver object which we can optionally
   * use for solving all systems in which the system matrix is set
   * to inner_product_matrix.
   */
  std::unique_ptr<LinearSolver<Number>> inner_product_solver;

  /**
   * Also, we store a pointer to an extra linear solver. This can be
   * useful if we want to pass in the linear solver from somewhere
   * else. For example, if a solver is already primed elsewhere
   * then it can be more efficient to use that solver.
   */
  LinearSolver<Number> * extra_linear_solver;

  /**
   * The inner product matrix.
   */
  std::unique_ptr<SparseMatrix<Number>> inner_product_matrix;

  /**
   * Vector storing the truth output values from the most
   * recent truth solve.
   */
  std::vector<Number > truth_outputs;

  /**
   * The vector storing the dual norm inner product terms
   * for each output.
   */
  std::vector<std::vector<Number >> output_dual_innerprods;

  /**
   * Vector storing the residual representors associated with the
   * right-hand side.
   * These are basis independent and hence stored here, whereas
   * the Aq_representors are stored in RBEvaluation
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> Fq_representor;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online. We store
   * the Fq representor norms here because they are independent
   * of a reduced basis space. The basis dependent representors
   * are stored in RBEvaluation.
   */
  std::vector<Number> Fq_representor_innerprods;

  /**
   * Boolean flag to indicate if we skip residual calculations
   * in train_reduced_basis. This should only be used in
   * special cases, e.g. when we know a priori that we want
   * exactly one basis function and hence we do not need the
   * residual based error indicator.
   */
  bool skip_residual_in_train_reduced_basis;

  /**
   * Boolean flag to indicate whether we exit the greedy if
   * we select the same parameters twice in a row. In some
   * problems this indicates that the greedy has "saturated"
   * typically due to numerical rounding effects.
   */
  bool exit_on_repeated_greedy_parameters;

  /**
   * Boolean flag to indicate whether we impose "fluxes"
   * (i.e. element boundary contributions to the weak form)
   * on internal element boundaries in the assembly routines.
   */
  bool impose_internal_fluxes;

  /**
   * In some cases meshes are intentionally created with degenerate sides
   * as a way to represent, say, triangles using a hex-only mesh. In this
   * situation we should detect and skip any degenerate sides in order to
   * prevent zero or negative element Jacobian errors.
   */
  bool skip_degenerate_sides;

  /**
   * Boolean flag to indicate whether we compute the RB_inner_product_matrix.
   * This is false by default in RBConstruction since (in the default implementation)
   * the RB inner-product matrix will just be the identity. But we may need the
   * inner-product matrix subclasses.
   */
  bool compute_RB_inner_product;

  /**
   * Boolean flag to indicate whether we store a second copy of each
   * affine operator and vector which does not have Dirichlet bcs
   * enforced.
   */
  bool store_non_dirichlet_operators;

  /**
   * Boolean flag to indicate whether we store a second copy of the
   * basis without constraints or dof transformations applied to it.
   * This is necessary when we have dof transformations and need
   * to calculate the residual R(U) = C^T F - C^T A C U, since we
   * need to evaluate R(U) using the untransformed basis U rather
   * than C U to avoid "double applying" dof transformations in C.
   */
  bool store_untransformed_basis;

  /**
   * A boolean flag to indicate whether or not we initialize the
   * Greedy algorithm by performing rb_solves on the training set
   * with an "empty" (i.e. N=0) reduced basis space.
   */
  bool use_empty_rb_solve_in_greedy;

  /**
   * A boolean flag to indicate whether or not the Fq representor norms
   * have already been computed --- used to make sure that we don't
   * recompute them unnecessarily.
   */
  bool Fq_representor_innerprods_computed;

protected:

  /**
   * Helper function that actually allocates all the data
   * structures required by this class.
   */
  virtual void allocate_data_structures();

  /**
   * Assemble the truth matrix and right-hand side
   * for current_parameters.
   */
  virtual void truth_assembly();

  /**
   * Builds a DGFEMContext object with enough information to do
   * evaluations on each element. We use DGFEMContext since it
   * allows for both DG and continuous Galerkin formulations.
   */
  virtual std::unique_ptr<DGFEMContext> build_context();

  /**
   * Return the matrix for the output residual dual
   * norm solves. By default we use the inner product matrix
   * for steady state problems.
   */
  virtual SparseMatrix<Number> & get_matrix_for_output_dual_solves();

  /**
   * Function that indicates when to terminate the Greedy
   * basis training. Override in subclasses to specialize.
   */
  virtual bool greedy_termination_test(Real abs_greedy_error,
                                       Real initial_greedy_error,
                                       int count);

  /**
   * Update the list of Greedily chosen parameters with
   * current_parameters.
   */
  void update_greedy_param_list();

  /**
   * This function loops over the mesh and applies the specified
   * interior and/or boundary assembly routines, then adds the
   * scaled result to input_matrix and/or input_vector.
   * If symmetrize==true then we assemble the symmetric part
   * of the matrix, 0.5*(A + A^T)
   */
  void add_scaled_matrix_and_vector(Number scalar,
                                    ElemAssembly * elem_assembly,
                                    SparseMatrix<Number> * input_matrix,
                                    NumericVector<Number> * input_vector,
                                    bool symmetrize=false,
                                    bool apply_dof_constraints=true);

  /**
   * This function is called from add_scaled_matrix_and_vector()
   * before each element matrix and vector are assembled into their
   * global counterparts. By default it is a no-op, but it could be
   * used to apply any user-defined transformations immediately prior
   * to assembly. We use DGFEMContext since it allows for both DG and
   * continuous Galerkin formulations.
   */
  virtual void post_process_elem_matrix_and_vector
  (DGFEMContext & /*context*/) {}

  /**
   * Similarly, provide an opportunity to post-process the truth
   * solution after the solve is complete. By default this is a no-op,
   * but it could be used to apply any required user-defined post
   * processing to the solution vector. Note: the truth solution is
   * stored in the "solution" member of this class, which is inherited
   * from the parent System class several levels up.
   */
  virtual void post_process_truth_solution() {}

  /**
   * Set current_local_solution = vec so that we can access vec
   * from FEMContext during assembly. Override in subclasses if
   * different behavior is required.
   */
  virtual void set_context_solution_vec(NumericVector<Number> & vec);

  /**
   * Assemble and store all the inner-product
   * matrix, the constraint matrix (for constrained
   * problems) and the mass matrix (for time-dependent
   * problems).
   */
  virtual void assemble_misc_matrices();

  /**
   * Assemble and store all Q_a affine operators
   * as well as the inner-product matrix.
   */
  virtual void assemble_all_affine_operators();

  /**
   * Assemble and store the affine RHS vectors.
   */
  virtual void assemble_all_affine_vectors();

  /**
   * Assemble and store the output vectors.
   */
  virtual void assemble_all_output_vectors();

  /**
   * Compute and store the dual norm of each output functional.
   */
  virtual void compute_output_dual_innerprods();

  /**
   * Compute the terms that are combined `online'
   * to determine the dual norm of the residual. Here we
   * compute the terms associated with the right-hand side.
   * These terms are basis independent, hence we separate
   * them from the rest of the calculations that are done in
   * update_residual_terms.
   * By default,
   * inner product terms are also computed, but you can turn
   * this feature off e.g. if you are already reading in that
   * data from files.
   */
  virtual void compute_Fq_representor_innerprods(bool compute_inner_products=true);

  /**
   * Add a new basis function to the RB space. This is called
   * during train_reduced_basis.
   */
  virtual void enrich_RB_space();

  /**
   * Update the system after enriching the RB space; this calls
   * a series of functions to update the system properly.
   */
  virtual void update_system();

  /**
   * \returns The RB error bound for the current parameters.
   *
   * Used in the Greedy algorithm to select the next parameter.
   */
  virtual Real get_RB_error_bound();

  /**
   * Compute the reduced basis matrices for the current basis.
   */
  virtual void update_RB_system_matrices();

  /**
   * Compute the terms that are combined `online'
   * to determine the dual norm of the residual.  By default,
   * inner product terms are also computed, but you can turn
   * this feature off e.g. if you are already reading in that
   * data from files.
   */
  virtual void update_residual_terms(bool compute_inner_products=true);

  /**
   * Initialize the FEMContext prior to performing
   * an element loop.
   * Reimplement this in derived classes in order to
   * call FE::get_*() as the particular physics requires.
   */
  virtual void init_context(FEMContext &) {
    // Failing to rederive init_context() means your FE objects don't
    // know what to compute.
    libmesh_deprecated();
  }

  /**
   * Getter for the flag determining if convergence should be
   * checked after each solve.
   */
  bool get_convergence_assertion_flag() const;

  /**
   * Check if the linear solver reports convergence.
   * Throw an error when that is not the case.
   */
  void check_convergence(LinearSolver<Number> & input_solver);

  /**
   * Get/set the current training parameter index
   */
  unsigned int get_current_training_parameter_index() const;
  void set_current_training_parameter_index(unsigned int index);

  /**
   * Return the evaluated theta functions at the given training parameter index.
   */
  const std::vector<Number> & get_evaluated_thetas(unsigned int training_parameter_index) const;

  /*
   * Pre-evaluate the theta functions on the entire (local) training parameter set.
   */
  void preevaluate_thetas();

  //----------- PROTECTED DATA MEMBERS -----------//

  /**
   * Maximum number of reduced basis
   * functions we are willing to use.
   */
  unsigned int Nmax;

  /**
   * The number of basis functions that we add at each greedy step.
   * This defaults to 1 in the steady case.
   */
  unsigned int delta_N;

  /**
   * A boolean flag to indicate whether or not the output dual norms
   * have already been computed --- used to make sure that we don't
   * recompute them unnecessarily.
   */
  bool output_dual_innerprods_computed;

  /**
   * A boolean flag to indicate whether to check for proper convergence
   * after each solve.
   */
  bool assert_convergence;

private:

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * The current RBEvaluation object we are using to
   * perform the Evaluation stage of the reduced basis
   * method.
   */
  RBEvaluation * rb_eval;

  /**
   * This member holds the (parameter independent) assembly functors
   * that define the "affine expansion" of the PDE that we are solving.
   */
  RBAssemblyExpansion * rb_assembly_expansion;

  /**
   * Pointer to inner product assembly.
   */
  ElemAssembly * inner_product_assembly;

  /**
   * Boolean to indicate whether we're using the energy inner-product.
   * If this is false then we use inner_product_assembly instead.
   */
  bool use_energy_inner_product;

  /**
   * We may optionally want to use the "energy inner-product" rather
   * than the inner-product assembly specified in inner_product_assembly.
   * In this case the inner-product will be defined by sum_q^Q k_q * A_q.
   * Here we provide the k_q values that will be used.
   * (Note that a true "energy-inner product" would obtain the k_q from
   * the theta_q's, but this is different for each parameter choice so
   * we just provide a fixed set of k_q's here to ensure that the
   * inner-product is parameter independent)
   */
  std::vector<Number> energy_inner_product_coeffs;

  /**
   * Vector storing the Q_a matrices from the affine expansion
   */
  std::vector<std::unique_ptr<SparseMatrix<Number>>> Aq_vector;

  /**
   * Vector storing the Q_f vectors in the affine decomposition
   * of the right-hand side.
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> Fq_vector;

  /**
   * The libMesh vectors that define the output functionals.
   * Each row corresponds to the affine expansion of an output.
   */
  std::vector<std::vector<std::unique_ptr<NumericVector<Number>>>> outputs_vector;

  /**
   * We may also need a second set of matrices/vectors
   * that do not have the Dirichlet boundary conditions
   * enforced.
   */
  std::vector<std::unique_ptr<SparseMatrix<Number>>> non_dirichlet_Aq_vector;
  std::vector<std::unique_ptr<NumericVector<Number>>> non_dirichlet_Fq_vector;
  std::vector<std::vector<std::unique_ptr<NumericVector<Number>>>> non_dirichlet_outputs_vector;
  std::unique_ptr<SparseMatrix<Number>> non_dirichlet_inner_product_matrix;

  /**
   * Relative and absolute tolerances for training reduced basis
   * using the Greedy scheme.
   */
  Real rel_training_tolerance;
  Real abs_training_tolerance;

  /**
   * This boolean indicates if we normalize the RB error in the greedy using
   * RBEvaluation::get_error_bound_normalization().
   */
  bool normalize_rb_bound_in_greedy;

  /**
   * This string indicates the type of training that we will use.
   * Options are:
   *  - Greedy: Reduced basis greedy algorithm
   *  - POD: Proper Orthogonal Decomposition
   */
  std::string RB_training_type;

  /**
   * In cases where we have dof transformations such as a change of
   * coordinates at some nodes we need to store an extra set of basis
   * functions which have not had dof transformations applied to them.
   * These vectors are required in order to compute the residual in
   * the error indicator.
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> _untransformed_basis_functions;

  /**
   * We also store a copy of the untransformed solution in order to
   * create _untransformed_basis_functions.
   */
  std::unique_ptr<NumericVector<Number>> _untransformed_solution;

  /**
   * Flag to indicate if we preevaluate the theta functions
   */
  bool _preevaluate_thetas_flag;

  /**
   * The current training parameter index during reduced basis training.
   */
  unsigned int _current_training_parameter_index;

  /**
   * Storage of evaluated theta functions at a set of parameters. This
   * can be used to store all of our theta functions at training samples
   * instead of re-evaluating the same values repeatedly during training.
   */
  std::vector<std::vector<Number>> _evaluated_thetas;
};

} // namespace libMesh

#endif // LIBMESH_RB_CONSTRUCTION_H
