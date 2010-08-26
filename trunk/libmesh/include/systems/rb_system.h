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

#ifndef __rb_system_h__
#define __rb_system_h__

#include "linear_implicit_system.h"
#include "dense_vector.h"
#include "dense_matrix.h"
#include "rb_base.h"
#include "fem_context.h"

namespace libMesh
{

/**
 * This class is part of the rbOOmit framework.
 *
 * RBSystem implements the certified reduced
 * basis framework in the steady-state case.
 *
 * @author David J. Knezevic, 2009
 */

/**
 * The function pointer typedef for
 * assembly of affine operators.
 */
typedef void (*affine_assembly_fptr)(FEMContext&, System&);

/**
 * The function pointer typedef for
 * assembly of Dirichlet and non-Dirichlet
 * dof lists.
 */
typedef void (*dirichlet_list_fptr)(FEMContext&, System&, std::set<unsigned int>&);

// ------------------------------------------------------------
// RBSystem class definition

class RBSystem : public RBBase<LinearImplicitSystem>
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  RBSystem (EquationSystems& es,
            const std::string& name,
            const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~RBSystem ();

  /**
   * The type of system.
   */
  typedef RBSystem sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * The type of the parent.
   */
  typedef RBBase<LinearImplicitSystem> Parent;

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
   * Perform online solve with the N RB basis functions, for the
   * set of parameters in current_params, where 1 <= N <= RB_size.
   */
  virtual Real RB_solve(unsigned int N);

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
   * @returns the final maximum a posteriori error bound on the training set.
   */
  Real train_reduced_basis(const std::string& directory_name = "offline_data");

  /**
   * Clear the basis functions and all basis-function-dependent data.
   * Overload in subclasses to clear any extra data.
   */
  virtual void clear_basis_function_dependent_data();

  /**
   * Return the parameters chosen during the i^th step of
   * the Greedy algorithm.
   */
  std::vector<Real> get_greedy_parameter(unsigned int i);

  /**
   * Set the name of the eigen_system that performs the SCM.
   */
  void set_eigen_system_name(const std::string& name);

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
  void set_Nmax(unsigned int Nmax);

  /**
   * Get Q_f, the number of terms in the affine
   * expansion for the right-hand side.
   */
  unsigned int get_Q_f() const { return theta_q_f_vector.size(); }

  /**
   * Get n_outputs, the number output functionals.
   */
  unsigned int get_n_outputs() const { return theta_q_l_vector.size(); }
  
  /**
   * Get the number of affine terms associated with the specified output.
   */
  unsigned int get_Q_l(unsigned int output_index) const;

  /**
   * Set the quiet flag. If quiet == false then
   * we print out a lot of extra information
   * during the Offline stage.
   */
  void set_quiet(bool quiet)
    { this->quiet = quiet; }

  /**
   * Load the i^th RB function into the RBSystem
   * solution vector.
   */
  void load_basis_function(unsigned int i);

  /**
   * Load the i^th RB function into vec.
   */
  void load_basis_function(unsigned int i, NumericVector<Number>& vec);

  /**
   * Get the current number of basis functions.
   */
  unsigned int get_n_basis_functions() const { return basis_functions.size(); }

  /**
   * Get a const reference to basis function i.
   */
  const NumericVector<Number> & get_bf(unsigned int i) const
    { libmesh_assert(i<=basis_functions.size());
      return *basis_functions[i]; };

  /**
   * Load the RB solution from the most recent solve
   * into the libMesh solution vector.
   */
  virtual void load_RB_solution();

  /**
   * Register a user function to use in initializing the
   * lists of Dirichlet and non-Dirichlet dofs.
   */
  void attach_dirichlet_dof_initialization (dirichlet_list_fptr dirichlet_init);

  /**
   * Call the user-defined Dirichlet dof initialization function.
   */
  void initialize_dirichlet_dofs();

  /**
   * Override attach_theta_q_a to just throw an error. Should
   * use attach_A_q in RBSystem and its subclasses.
   */
  virtual void attach_theta_q_a(theta_q_fptr )
  {
    std::cout << "Error: Cannot use attach_theta_q_a in RBSystem. Use attach_A_q instead." << std::endl;
    libmesh_error();
  }

  /**
   * Attach parameter-dependent function and user-defined assembly routine
   * for affine operator (both interior and boundary assembly).
   */
  void attach_A_q(theta_q_fptr theta_q_a,
                  affine_assembly_fptr A_q_intrr_assembly,
                  affine_assembly_fptr A_q_bndry_assembly);

  /**
   * Attach parameter-dependent function and user-defined assembly routine
   * for affine vector. (Interior assembly and boundary assembly).
   */
  void attach_F_q(theta_q_fptr theta_q_f,
                  affine_assembly_fptr F_q_intrr_assembly,
                  affine_assembly_fptr F_q_bdnry_assembly);

  /**
   * Attach user-defined assembly routine
   * for the inner-product matrix.
   */
  void attach_inner_prod_assembly(affine_assembly_fptr IP_assembly);

  /**
   * Attach user-defined assembly routine
   * for the constraint matrix.
   */
  void attach_constraint_assembly(affine_assembly_fptr constraint_assembly_in);

  /**
   * Attach user-defined assembly routine for output. 
   * (Interior assembly and boundary assembly).
   * In this case we pass in vector arguments to allow for Q_l > 1.
   */
  void attach_output(std::vector<theta_q_fptr> theta_q_l,
                     std::vector<affine_assembly_fptr> output_intrr_assembly,
                     std::vector<affine_assembly_fptr> output_bndry_assembly);

  /**
   * Get a pointer to A_q.
   */
  SparseMatrix<Number>* get_A_q(unsigned int q);

  /**
   * Evaluate theta_q_f at the current parameter.
   */
  Number eval_theta_q_f(unsigned int q);

  /**
   * Evaluate theta_q_l at the current parameter.
   */
  Number eval_theta_q_l(unsigned int output_index, unsigned int q_l);

  /**
   * Assemble and store the Dirichlet dof lists, the
   * affine and output vectors.
   * Assemble and store all the affine matrices if we
   * are not in low-memory mode.
   */
  virtual void perform_initial_assembly();

  /**
   * Get a pointer to F_q.
   */
  NumericVector<Number>* get_F_q(unsigned int q);

  /**
   * Get a pointer to the n^th output.
   */
  NumericVector<Number>* get_output_vector(unsigned int n, unsigned int q_l);

  /**
   * Assemble the inner product matrix and store it in input_matrix.
   */
  void assemble_inner_product_matrix(SparseMatrix<Number>* input_matrix);

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
  void assemble_Aq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix);

  /**
   * Add the scaled q^th affine matrix to input_matrix. If symmetrize==true, then
   * we symmetrize Aq before adding it.
   */
  void add_scaled_Aq(Number scalar, unsigned int q_a,
                     SparseMatrix<Number>* input_matrix,
                     bool symmetrize);

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
  
  /**
   * This function computes all of the residual terms, can be useful
   * when restarting a basis training computation.
   */
  virtual void recompute_all_residual_terms();


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
   * The inner product matrix.
   */
  AutoPtr< SparseMatrix<Number> > inner_product_matrix;

  /**
   * The constraint matrix, e.g. the pressure matrix entries
   * in a Stokes problem.
   */
  AutoPtr< SparseMatrix<Number> > constraint_matrix;

  /**
   * The libMesh vectors storing the finite element coefficients
   * of the RB basis functions.
   */
  std::vector< NumericVector<Number>* > basis_functions;

  /**
   * Dense matrices for the RB computations.
   */
  std::vector< DenseMatrix<Number> > RB_A_q_vector;

  /**
   * Dense vector for the RHS.
   */
  std::vector< DenseVector<Number> > RB_F_q_vector;

  /**
   * The RB solution vector.
   */
  DenseVector<Number> RB_solution;

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
   * The vectors storing the RB output vectors.
   */
  std::vector< std::vector< DenseVector<Number> > > RB_output_vectors;

  /**
   * The vectors storing the RB output values and
   * corresponding error bounds.
   */
  std::vector< Number > RB_outputs;
  std::vector< Real > RB_output_error_bounds;

  /**
   * The list of parameters selected by the Greedy algorithm.
   */
  std::vector< std::vector<Real> > greedy_param_list;

  /**
   * Boolean flag to indicate whether this is a constrained problem
   * (e.g. Stokes) or not.
   */
  bool constrained_problem;

  /**
   * Boolean flag to indicate whether the basis functions are written
   * out from the Offline stage or read in during the Online stage.
   */
  bool store_basis_functions;

  /**
   * Boolean flag to indicate whether the residual representors are
   * written out from the offline stage.  For large problems, reading
   * them back in during a restart can be much faster than recomputing
   * them.
   */
  bool store_representors;

  /**
   * Boolean flag to indicate whether or not we are in "low-memory" mode.
   * In low-memory mode, we do not store any extra sparse matrices.
   */
  bool low_memory_mode;

  /**
   * Boolean flag to indicate whether we reuse the preconditioner
   * on consecutive Offline solves for updating residual data.
   */
  bool reuse_preconditioner;

  /**
   * Boolean flag to indicate whether RB_solve returns an absolute
   * or relative error bound. True => relative, false => absolute.
   */
  bool return_rel_error_bound;

  /**
   * Boolean flag to indicate whether train_reduced_basis writes
   * out offline data after each truth solve (to allow continuing
   * in case the code crashes or something).
   */
  bool write_data_during_training;

  /**
   * Boolean flag to indicate whether we require a Dirichlet boundary
   * condition on internal mesh dofs, for example in a problem in which
   * parameters are not defined on some subdomains.
   */
  bool impose_internal_dirichlet_BCs;

  /**
   * Boolean flag to indicate whether we impose flux on internal element
   * boundaries.
   */
  bool impose_internal_fluxes;

  /**
   * The filename of the text file from which we read in the
   * problem parameters. We use getpot.h to perform the reading.
   */
  std::string parameters_filename;

  /**
   * Defines extra quadrature order (can be positive or negative).
   * Default value is 0.
   */
  int extra_quadrature_order;
  
  /**
   * Public member variable which we use to determine whether or
   * not we enforce hanging-dof and/or periodic constraints exactly.
   * This is primarily important in nonlinear problems where we may
   * "undersolve" Newton iterates for the sake of efficiency.
   */
  bool enforce_constraints_exactly;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

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
                                    affine_assembly_fptr intrr_assembly,
                                    affine_assembly_fptr bndry_assembly,
                                    SparseMatrix<Number>* input_matrix,
                                    NumericVector<Number>* input_vector,
                                    bool symmetrize=false);

  /**
   * Set current_local_solution = vec so that we can access vec
   * from FEMContext during assembly. Overload in subclasses if
   * different behavior is required (e.g. in QNTransientRBSystem)
   */
  virtual void set_context_solution_vec(NumericVector<Number>& vec);

  /**
   * This function loops over the mesh and assembles the
   * matrix-vector product and stores the scaled result
   * in dest.
   */
  void assemble_scaled_matvec(Number scalar,
                              affine_assembly_fptr intrr_assembly,
                              affine_assembly_fptr bndry_assembly,
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
  void assemble_all_affine_vectors();

  /**
   * Assemble and store the output vectors.
   */
  void assemble_all_output_vectors();

  /**
   * Compute the dual norm of the residual for the solution
   * saved in RB_solution_vector.
   */
  virtual Real compute_residual_dual_norm(const unsigned int N);

  /**
   * Compute and store the dual norm of each output functional.
   */
  virtual void compute_output_dual_norms();

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
   * Compute the a posteriori error bound for each set of parameters
   * in training_parameters and return the pair containing the max
   * error and the index of the parameter that induces that error.
   */
  virtual Real compute_a_posteriori_bounds();

  /**
   * Get the SCM lower bound at the current parameter value.
   */
  virtual Real get_SCM_lower_bound();

  /**
   * Get the SCM upper bound at the current parameter value.
   */
  virtual Real get_SCM_upper_bound();

  /**
   * Specifies the residual scaling on the denominator to
   * be used in the a posteriori error bound. Overload
   * in subclass in order to obtain the desired error bound.
   */
  virtual Real residual_scaling_denom(Real alpha_LB);

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
   * Controls whether or not XDR (binary) files are written out for
   * the basis functions.  The binary file size can be as small
   * as 1/3 the size of an ASCII file.
   */
  bool write_binary_basis_functions;

  /**
   * Controls wether XDR (binary) files are read for
   * the basis functions.  Note: if you wrote ASCII basis functions
   * during a previous run but want to start writing XDR, set
   * read_binary_basis_functions = false and
   * write_binary_basis_functions = true.
   */
  bool read_binary_basis_functions;

  /**
   * Controls whether or not XDR (binary) files are written out for
   * the residual respresentors.  The binary file size can be as small
   * as 1/3 the size of an ASCII file.
   */
  bool write_binary_residual_representors;

  /**
   * Controls wether XDR (binary) files are read for
   * the residual representors.  Note: if you wrote ASCII representors
   * during a previous run but want to start writing XDR, set
   * read_binary_residual_representors = false and
   * write_binary_residual_representors = true.
   */
  bool read_binary_residual_representors;

  /**
   * Flag to indicate whether we print out extra information during
   * the Offline stage.
   */
  bool quiet;

  /**
   * Vectors storing the residual representors.
   */
  std::vector< NumericVector<Number>* > F_q_representor;
  std::vector< std::vector< NumericVector<Number>* > > A_q_representor;

  /**
   * Vectors storing the residual representor inner products
   * to be used in computing the residuals online.
   */
  std::vector<Number> Fq_representor_norms;
  std::vector< std::vector< std::vector<Number> > > Fq_Aq_representor_norms;
  std::vector< std::vector< std::vector<Number> > > Aq_Aq_representor_norms;

  /**
   * The name of the RBSCMSystem system that performs
   * the SCM.
   */
  std::string eigen_system_name;

  /**
   * Function pointer for assembling the inner product
   * matrix.
   */
  affine_assembly_fptr inner_prod_assembly;

  /**
   * Function pointer for assembling the constraint
   * matrix.
   */
  affine_assembly_fptr constraint_assembly;

  /**
   * Vectors storing the function pointers to the assembly
   * routines for the affine operators, both interior and boundary
   * assembly.
   */
  std::vector<affine_assembly_fptr> A_q_intrr_assembly_vector;
  std::vector<affine_assembly_fptr> A_q_bndry_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the affine vectors. Element interior part.
   */
  std::vector<affine_assembly_fptr> F_q_intrr_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the affine vectors. Boundary part.
   */
  std::vector<affine_assembly_fptr> F_q_bndry_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the outputs. Element interior part.
   */
  std::vector< std::vector<affine_assembly_fptr> > output_intrr_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the outputs. Boundary part.
   */
  std::vector< std::vector<affine_assembly_fptr> > output_bndry_assembly_vector;

private:

  /**
   * Private non-virtual helper function to encapsulate
   * the code to clear the basis-function-related data.
   */
  void clear_basis_helper();

  //----------- PRIVATE DATA MEMBERS -----------//

  /**
   * Vector storing the function pointers to the theta_q_f (affine rhs vectors).
   */
  std::vector<theta_q_fptr> theta_q_f_vector;

  /**
   * Vector storing the function pointers to the theta_q_l (for the outputs).
   */
  std::vector< std::vector<theta_q_fptr> > theta_q_l_vector;

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
   * The libMesh vectors that define the output functionals.
   * Each row corresponds to the affine expansion of an output.
   */
  std::vector< std::vector< NumericVector<Number>* > > outputs_vector;

  /**
   * Tolerance for training reduced basis using the Greedy scheme.
   */
  Real training_tolerance;

  /**
   * A boolean flag to indicate whether or not update_residual_terms has
   * been called before --- used to make sure that the representors for
   * the Fq are only computed on the first call to update_residual_terms.
   */
  bool update_residual_terms_called;

  /**
   * Function that initializes the lists of Dirichlet and
   * non-Dirichlet dofs.
   */
  dirichlet_list_fptr _dirichlet_list_init;

  /**
   * Member variable to store the value of Nmax that was read in from
   * file and was used in init_data to initialize the RBSystem.
   */
  unsigned int initial_Nmax;

  /**
   * Boolean flag to indicate whether the RBSystem has been initialized.
   */
  bool RB_system_initialized;
};

} // namespace libMesh

#endif
