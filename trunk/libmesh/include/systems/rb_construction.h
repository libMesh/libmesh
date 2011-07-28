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

#ifndef __rb_construction_h__
#define __rb_construction_h__

// rbOOmit includes
#include "rb_construction_base.h"
#include "rb_evaluation.h"

// libMesh includes
#include "linear_implicit_system.h"
#include "dense_vector.h"
#include "dense_matrix.h"
#include "fem_context.h"
#include "elem_assembly.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBConstruction implements the Construction stage
 * of the certified reduced basis method for
 * steady-state elliptic parametrized PDEs.
 *
 * @author David J. Knezevic, 2009
 */


/**
 * First, define a class that allows us to assemble a list of Dirichlet
 * dofs. This list is used by RBConstruction to impose homoegenous Dirichlet
 * boundary conditions (in the RB method, non-homogeneous Dirichlet
 * boundary conditions should be imposed via lifting functions).
 */
class DirichletDofAssembly : public ElemAssembly
{
public:
  std::set<unsigned int> dirichlet_dofs_set;
};

// ------------------------------------------------------------
// RBConstruction class definition

class RBConstruction : public RBConstructionBase<LinearImplicitSystem>
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBConstruction (EquationSystems& es,
                  const std::string& name,
                  const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBConstruction ();

  /**
   * The type of system.
   */
  typedef RBConstruction sys_type;

  /**
   * @returns a clever pointer to the system.
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
  virtual void clear ();

  /**
   * @returns a string indicating the type of the system.
   */
  virtual std::string system_type () const;

  /**
   * Perform a "truth" solve, i.e. solve the finite element system at
   * at the parameters currently set in the system. This is used
   * extensively in training the reduced basis, since "truth snapshots"
   * are employed as basis functions.
   */
  virtual Real truth_solve(int plot_solution);

  /**
   * Train the reduced basis. This is the crucial function in the Offline
   * stage: it generates the reduced basis using the "Greedy algorithm."
   * Each stage of the Greedy algorithm involves solving the reduced basis
   * over a large training set and selecting the parameter at which the
   * reduced basis error bound is largest, then performing a truth_solve
   * at that parameter and enriching the reduced basis with the corresponding
   * snapshot.
   *
   * \p resize_rb_eval_data is a boolean flag to indicate whether or not we
   * call rb_eval->resize_data_structures(Nmax). True by default, by we may
   * set it to false if, for example, we are continuing from a previous
   * training run and don't want to clobber the existing rb_eval data.
   *
   * @returns the final maximum a posteriori error bound on the training set.
   */
  virtual Real train_reduced_basis(const std::string& directory_name = "offline_data",
                                   const bool resize_rb_eval_data=true);

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
  std::vector<Real> get_greedy_parameter(unsigned int i);

  /**
   * Get/set the tolerance for the basis training.
   */
  void set_training_tolerance(Real training_tolerance)
    {this->training_tolerance = training_tolerance; }
  Real get_training_tolerance() { return training_tolerance; }

  /**
   * Get/set Nmax, the maximum number of RB
   * functions we are willing to compute.
   */
  unsigned int get_Nmax() const    { return Nmax; }
  virtual void set_Nmax(unsigned int Nmax);

  /**
   * Set the quiet_mode flag. If quiet == false then
   * we print out a lot of extra information
   * during the Offline stage.
   */
  void set_quiet_mode(bool quiet_mode_in)
    { this->quiet_mode = quiet_mode_in; }

  /**
   * Is the system in quiet mode?
   */
  bool is_quiet() const
    { return this->quiet_mode; }

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
   * Register a user-defined object to be used in initializing the
   * lists of Dirichlet dofs.
   */
  void attach_dirichlet_dof_initialization (DirichletDofAssembly* dirichlet_init);

  /**
   * Call the user-defined Dirichlet dof initialization function.
   */
  void initialize_dirichlet_dofs();

  /**
   * Attach parameter-dependent function and user-defined assembly routine
   * for affine operator (both interior and boundary assembly).
   */
  virtual void attach_A_q(RBTheta* theta_q_a,
                          ElemAssembly* A_q_assembly);

  /**
   * Attach parameter-dependent function and user-defined assembly routine
   * for affine vector. (Interior assembly and boundary assembly).
   */
  virtual void attach_F_q(RBTheta* theta_q_f,
                          ElemAssembly* F_q_assembly);

  /**
   * Attach user-defined assembly routine
   * for the inner-product matrix.
   */
  void attach_inner_prod_assembly(ElemAssembly* IP_assembly);

  /**
   * Attach user-defined assembly routine
   * for the constraint matrix.
   */
  void attach_constraint_assembly(ElemAssembly* constraint_assembly_in);

  /**
   * Attach user-defined assembly routine for output. 
   * (Interior assembly and boundary assembly).
   * In this case we pass in vector arguments to allow for Q_l > 1.
   */
  virtual void attach_output(std::vector<RBTheta*> theta_q_l,
                             std::vector<ElemAssembly*> output_assembly);

  /**
   * Attach user-defined assembly routine for output. 
   * (Interior assembly and boundary assembly).
   * This function provides simpler syntax in the case that Q_l = 1; we
   * do not need to use a vector in this case.
   */
  virtual void attach_output(RBTheta* theta_q_l,
                             ElemAssembly* output_assembly);

  /**
   * Attach the full affine expansion associated with this system.
   * We assert that the assembly vectors have sizes consistent
   * with \p rb_theta_expansion_in. Also, this is not used for attaching
   * outputs (use attach_output for that), hence we assert that
   * rb_theta_expansion_in->n_outputs() == 0.
   */
  void attach_affine_expansion(RBThetaExpansion& rb_theta_expansion_in,
                               std::vector<ElemAssembly*> A_q_assembly_vector_in,
                               std::vector<ElemAssembly*> F_q_assembly_vector_in);

  /**
   * Get a pointer to inner_product_matrix. Accessing via this
   * function, rather than directly through the class member allows
   * us to do error checking (e.g. inner_product_matrix is not
   * defined in low-memory mode).
   */
  SparseMatrix<Number>* get_inner_product_matrix();

  /**
   * Get a pointer to non_dirichlet_inner_product_matrix.
   * Accessing via this function, rather than directly through
   * the class member allows us to do error checking (e.g.
   * non_dirichlet_inner_product_matrix is not
   * defined in low-memory mode, and we need
   * store_non_dirichlet_operators==true).
   */
  SparseMatrix<Number>* get_non_dirichlet_inner_product_matrix();

  /**
   * Get a pointer to A_q.
   */
  SparseMatrix<Number>* get_A_q(unsigned int q);

  /**
   * Get a pointer to non_dirichlet_A_q.
   */
  SparseMatrix<Number>* get_non_dirichlet_A_q(unsigned int q);

  /**
   * Allocate all the data structures necessary for the construction
   * stage of the RB method. This function also performs
   * matrix and vector assembly of the "truth" affine expansion.
   */
  virtual void initialize_rb_construction();

  /**
   * Get a pointer to F_q.
   */
  NumericVector<Number>* get_F_q(unsigned int q);

  /**
   * Get a pointer to non-Dirichlet F_q.
   */
  NumericVector<Number>* get_non_dirichlet_F_q(unsigned int q);

  /**
   * Get a pointer to the n^th output.
   */
  NumericVector<Number>* get_output_vector(unsigned int n, unsigned int q_l);

  /**
   * Assemble the inner product matrix and store it in input_matrix.
   */
  void assemble_inner_product_matrix(SparseMatrix<Number>* input_matrix, bool apply_dirichlet_bc=true);

  /**
   * Assemble the constraint matrix and store it in input_matrix.
   */
  void assemble_constraint_matrix(SparseMatrix<Number>* input_matrix);

  /**
   * Assemble the constraint matrix and add it to input_matrix.
   */
  void assemble_and_add_constraint_matrix(SparseMatrix<Number>* input_matrix);

  /**
   * Assemble the q^th affine matrix and store it in input_matrix.
   */
  void assemble_Aq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix, bool apply_dirichlet_bc=true);

  /**
   * Assemble the q^th affine vector and store it in input_matrix.
   */
  void assemble_Fq_vector(unsigned int q, NumericVector<Number>* input_vector, bool apply_dirichlet_bc=true);

  /**
   * Add the scaled q^th affine matrix to input_matrix. If symmetrize==true, then
   * we symmetrize Aq before adding it.
   */
  void add_scaled_Aq(Number scalar, unsigned int q_a,
                     SparseMatrix<Number>* input_matrix,
                     bool symmetrize);

  /**
   * Write out all the Riesz representor data to files.
   * \p residual_representors_dir specifies which directory to write to.
   * \p write_binary_residual_representors specifies whether we write binary or ASCII data.
   */
  virtual void write_riesz_representors_to_files(const std::string& riesz_representors_dir,
                                                 const bool write_binary_residual_representors);

  /**
   * Read in all the Riesz representor data from files.
   * \p directory_name specifies which directory to read from.
   * \p io_flag specifies whether we read in all data, or only
   * a basis (in)dependent subset.
   */
  virtual void read_riesz_representors_from_files(const std::string& riesz_representors_dir,
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
   * Read in the parameters from file specified by the member variable
   * "parameters_filename" and set the this system's member variables
   * accordingly.
   */
  virtual void process_parameters_file(const std::string& parameters_filename);
  
  /**
   * Print out info that describes the current setup of this RBConstruction.
   */
  virtual void print_info();

  /**
   * Build a new RBEvaluation object.
   * @return an AutoPtr to the new RBEvaluation object.
   */
  virtual AutoPtr<RBEvaluation> build_rb_evaluation();

  /**
   * Get delta_N, the number of basis functions we
   * add to the RB space per iteration of the greedy
   * algorithm. For steady-state systems, this should
   * be 1, but can be more than 1 for time-dependent
   * systems.
   */
  unsigned int get_delta_N() const { return delta_N; }


  //----------- PUBLIC DATA MEMBERS -----------//

  /**
   * The current RBEvaluation object we are using to
   * perform the Evaluation stage of the reduced basis
   * method.
   */
  RBEvaluation* rb_eval;

  /**
   * Vector storing the values of the error bound
   * for each parameter in the training set --- the
   * parameter giving the largest error bound is
   * chosen for the next snapshot in the Greedy basis
   * training.
   */
  std::vector<Real> training_error_bounds;

  /**
   * The inner product matrix.
   */
  AutoPtr< SparseMatrix<Number> > inner_product_matrix;

  /**
   * The inner product matrix without Dirichlet conditions enforced.
   * (This is only computed if store_non_dirichlet_operators == true.)
   */
  AutoPtr< SparseMatrix<Number> > non_dirichlet_inner_product_matrix;

  /**
   * The constraint matrix, e.g. the pressure matrix entries
   * in a Stokes problem.
   */
  AutoPtr< SparseMatrix<Number> > constraint_matrix;

  /**
   * Vector storing the truth output values from the most
   * recent truth solve.
   */
  std::vector< Number > truth_outputs;

  /**
   * The vector storing the dual norm inner product terms
   * for each output.
   */
  std::vector< std::vector< Number > > output_dual_norms;

  /**
   * Vector storing the residual representors associated with the
   * right-hand side.
   * These are basis independent and hence stored here, whereas
   * the A_q_representors are stored in RBEvaluation
   */
  std::vector< NumericVector<Number>* > F_q_representor;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online. We store
   * the Fq representor norms here because they are independent
   * of a reduced basis space. The basis dependent representors
   * are stored in RBEvaluation.
   */
  std::vector<Number> Fq_representor_norms;

  /**
   * Set storing the global Dirichlet dof indices.
   * We update this set in initialize_dirichlet_dofs().
   */
  std::set<unsigned int> global_dirichlet_dofs_set;

  /**
   * Boolean flag to indicate whether this is a constrained problem
   * (e.g. Stokes) or not.
   */
  bool constrained_problem;

  /**
   * Boolean flag to indicate whether we store just one matrix in order
   * to conserve memory. If single_matrix_mode=false then we store
   * each matrix that is required for the affine decomposition of the PDE.
   */
  bool single_matrix_mode;

  /**
   * Boolean flag to indicate whether we reuse the preconditioner
   * on consecutive Offline solves for updating residual data.
   */
  bool reuse_preconditioner;

  /**
   * Boolean flag to indicate whether we use an absolute or
   * relative error bound in the Greedy algorithm for training
   * a Reduced Basis.
   */
  bool use_relative_bound_in_greedy;

  /**
   * Boolean flag to indicate whether train_reduced_basis writes
   * out offline data after each truth solve (to allow continuing
   * in case the code crashes or something).
   */
  bool write_data_during_training;

  /**
   * Boolean flag to indicate whether we require a Dirichlet boundary
   * condition on internal mesh dofs. This is useful for example in problems
   * with subdomain-only-variables.
   */
  bool impose_internal_dirichlet_BCs;

  /**
   * Boolean flag to indicate whether we impose "fluxes"
   * (i.e. element boundary contributions to the weak form)
   * on internal element boundaries in the assembly routines.
   */
  bool impose_internal_fluxes;
  
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
   * Public member variable which we use to determine whether or
   * not we enforce hanging-dof and/or periodic constraints exactly.
   * This is primarily important in nonlinear problems where we may
   * "undersolve" Newton iterates for the sake of efficiency.
   */
  bool enforce_constraints_exactly;

  /**
   * A boolean flag to indicate whether or not we initialize the
   * Greedy algorithm by performing rb_solves on the training set
   * with an "empty" (i.e. N=0) reduced basis space.
   */
  bool use_empty_rb_solve_in_greedy;

protected:
  
  /**
   * Helper function that actually allocates all the data
   * structures required by this class.
   */
  virtual void allocate_data_structures();

  /**
   * Assemble and store the Dirichlet dof lists, the
   * affine and output vectors.
   * Optionally assemble and store all the affine matrices if we
   * are not in low-memory mode.
   */
  virtual void assemble_affine_expansion();

  /**
   * Assemble the truth matrix and right-hand side
   * for current_parameters.
   */
  virtual void truth_assembly();

  /**
   * Builds a FEMContext object with enough information to do
   * evaluations on each element.
   */
  virtual AutoPtr<FEMContext> build_context();
  
  /**
   * Define the matrix assembly for the output residual dual
   * norm solves. By default we use the inner product matrix
   * for steady state problems.
   */
  virtual void assemble_matrix_for_output_dual_solves();

  /**
   * Function that indicates when to terminate the Greedy
   * basis training. Overload in subclasses to specialize.
   */
  virtual bool greedy_termination_test(Real training_greedy_error, int count);

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
                                    ElemAssembly* elem_assembly,
                                    SparseMatrix<Number>* input_matrix,
                                    NumericVector<Number>* input_vector,
                                    bool symmetrize=false,
                                    bool apply_dirichlet_bc=true);

  /**
   * Set current_local_solution = vec so that we can access vec
   * from FEMContext during assembly. Override in subclasses if
   * different behavior is required.
   */
  virtual void set_context_solution_vec(NumericVector<Number>& vec);

  /**
   * This function loops over the mesh and assembles the
   * matrix-vector product and stores the scaled result
   * in dest.
   */
  void assemble_scaled_matvec(Number scalar,
                              ElemAssembly* elem_assembly,
                              NumericVector<Number>& dest,
                              NumericVector<Number>& arg);

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
  virtual void compute_output_dual_norms();

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
  virtual void compute_Fq_representor_norms(bool compute_inner_products=true);

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
   * This function returns the RB error bound for the current parameters and
   * is used in the Greedy algorithm to select the next parameter.
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
  virtual void init_context(FEMContext& ) {}

  /**
   * Set the dofs on the Dirichlet boundary (stored in
   * the set global_dirichlet_dofs) to zero
   * in the right-hand side vector.
   */
  void zero_dirichlet_dofs_on_rhs();

  /**
   * Set the dofs on the Dirichlet boundary (stored in
   * the set global_dirichlet_dofs) to zero
   * in the vector temp.
   */
  void zero_dirichlet_dofs_on_vector(NumericVector<Number>& temp);

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
   * Flag to indicate whether we print out extra information during
   * the Offline stage.
   */
  bool quiet_mode;

  /**
   * Pointer to inner product assembly.
   */
  ElemAssembly* inner_prod_assembly;

  /**
   * Function pointer for assembling the constraint
   * matrix.
   */
  ElemAssembly* constraint_assembly;

  /**
   * Vectors storing the function pointers to the assembly
   * routines for the affine operators, both interior and boundary
   * assembly.
   */
  std::vector<ElemAssembly*> A_q_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the rhs affine vectors.
   */
  std::vector<ElemAssembly*> F_q_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the outputs. Element interior part.
   */
  std::vector< std::vector<ElemAssembly*> > output_assembly_vector;

  /**
   * A boolean flag to indicate whether or not the output dual norms
   * have already been computed --- used to make sure that we don't
   * recompute them unnecessarily.
   */
  bool output_dual_norms_computed;

  /**
   * A boolean flag to indicate whether or not the Fq representor norms
   * have already been computed --- used to make sure that we don't
   * recompute them unnecessarily.
   */
  bool Fq_representor_norms_computed;

private:

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * Vector storing the Q_a matrices from the affine expansion
   */
  std::vector< SparseMatrix<Number>* > A_q_vector;

  /**
   * Vector storing the Q_f vectors in the affine decomposition
   * of the right-hand side.
   */
  std::vector< NumericVector<Number>* > F_q_vector;

  /**
   * We sometimes also need a second set of matrices/vectors
   * that do not have the Dirichlet boundary conditions
   * enforced.
   */
  std::vector< SparseMatrix<Number>* > non_dirichlet_A_q_vector;
  std::vector< NumericVector<Number>* > non_dirichlet_F_q_vector;

  /**
   * The libMesh vectors that define the output functionals.
   * Each row corresponds to the affine expansion of an output.
   */
  std::vector< std::vector< NumericVector<Number>* > > outputs_vector;

  /**
   * Tolerance for training reduced basis using the Greedy scheme.
   */
  Real training_tolerance;

  /**
   * Assembly object that is used to initialize the list of Dirichlet dofs.
   */
  DirichletDofAssembly* _dirichlet_list_init;

};

} // namespace libMesh

#endif
