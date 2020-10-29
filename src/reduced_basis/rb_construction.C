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

// rbOOmit includes
#include "libmesh/rb_construction.h"
#include "libmesh/rb_assembly_expansion.h"
#include "libmesh/rb_evaluation.h"
#include "libmesh/elem_assembly.h"

// LibMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/linear_solver.h"
#include "libmesh/getpot.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/timestamp.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/dg_fem_context.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/coupling_matrix.h"
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

// C++ includes
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdlib.h> // mkstemps on Linux
#include <unistd.h> // mkstemps on MacOS

namespace libMesh
{

RBConstruction::RBConstruction (EquationSystems & es,
                                const std::string & name_in,
                                const unsigned int number_in)
  : Parent(es, name_in, number_in),
    inner_product_solver(LinearSolver<Number>::build(es.comm())),
    extra_linear_solver(nullptr),
    inner_product_matrix(SparseMatrix<Number>::build(es.comm())),
    skip_residual_in_train_reduced_basis(false),
    exit_on_repeated_greedy_parameters(true),
    impose_internal_fluxes(false),
    skip_degenerate_sides(true),
    compute_RB_inner_product(false),
    store_non_dirichlet_operators(false),
    store_untransformed_basis(false),
    use_empty_rb_solve_in_greedy(true),
    Fq_representor_innerprods_computed(false),
    Nmax(0),
    delta_N(1),
    output_dual_innerprods_computed(false),
    assert_convergence(true),
    rb_eval(nullptr),
    inner_product_assembly(nullptr),
    use_energy_inner_product(false),
    rel_training_tolerance(1.e-4),
    abs_training_tolerance(1.e-12),
    normalize_rb_bound_in_greedy(false),
    RB_training_type("Greedy"),
    _preevaluate_thetas_flag(false)
{
  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;
}

RBConstruction::~RBConstruction () = default;

void RBConstruction::clear()
{
  LOG_SCOPE("clear()", "RBConstruction");

  Parent::clear();

  Aq_vector.clear();
  Fq_vector.clear();
  outputs_vector.clear();

  if (store_non_dirichlet_operators)
    {
      non_dirichlet_Aq_vector.clear();
      non_dirichlet_Fq_vector.clear();
      non_dirichlet_outputs_vector.clear();
    }

  // Also delete the Fq representors
  Fq_representor.clear();

  // Set Fq_representor_innerprods_computed flag to false now
  // that we've cleared the Fq representors
  Fq_representor_innerprods_computed = false;
}

std::string RBConstruction::system_type () const
{
  return "RBConstruction";
}

void RBConstruction::solve_for_matrix_and_rhs(LinearSolver<Number> & input_solver,
                                              SparseMatrix<Number> & input_matrix,
                                              NumericVector<Number> & input_rhs)
{
  // This is similar to LinearImplicitSysmte::solve()

  // Get a reference to the EquationSystems
  const EquationSystems & es =
    this->get_equation_systems();

  // If the linear solver hasn't been initialized, we do so here.
  input_solver.init();

  // Get the user-specifiied linear solver tolerance
  const double tol  =
    double(es.parameters.get<Real>("linear solver tolerance"));

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);

  // It's good practice to clear the solution vector first since it can
  // affect convergence of iterative solvers
  solution->zero();
  rval = input_solver.solve (input_matrix, *solution, input_rhs, tol, maxits);

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;

  get_dof_map().enforce_constraints_exactly(*this);

  // Update the system after the solve
  this->update();
}

void RBConstruction::set_rb_evaluation(RBEvaluation & rb_eval_in)
{
  rb_eval = &rb_eval_in;
}

RBEvaluation & RBConstruction::get_rb_evaluation()
{
  libmesh_error_msg_if(!rb_eval, "Error: RBEvaluation object hasn't been initialized yet");

  return *rb_eval;
}

bool RBConstruction::is_rb_eval_initialized() const
{
  return (rb_eval != nullptr);
}

RBThetaExpansion & RBConstruction::get_rb_theta_expansion()
{
  return get_rb_evaluation().get_rb_theta_expansion();
}

void RBConstruction::process_parameters_file (const std::string & parameters_filename)
{
  // First read in data from input_filename
  GetPot infile(parameters_filename);

  const unsigned int n_training_samples = infile("n_training_samples",0);
  const bool deterministic_training = infile("deterministic_training",false);
  unsigned int training_parameters_random_seed_in =
    static_cast<unsigned int>(-1);
  training_parameters_random_seed_in = infile("training_parameters_random_seed",
                                              training_parameters_random_seed_in);
  const bool quiet_mode_in = infile("quiet_mode", quiet_mode);
  const unsigned int Nmax_in = infile("Nmax", Nmax);
  const Real rel_training_tolerance_in = infile("rel_training_tolerance",
                                                rel_training_tolerance);
  const Real abs_training_tolerance_in = infile("abs_training_tolerance",
                                                abs_training_tolerance);

  // Initialize value to false, let the input file value override.
  const bool normalize_rb_bound_in_greedy_in = infile("normalize_rb_bound_in_greedy",
                                                      false);

  const std::string RB_training_type_in = infile("RB_training_type", "Greedy");

  // Read in the parameters from the input file too
  unsigned int n_continuous_parameters = infile.vector_variable_size("parameter_names");
  RBParameters mu_min_in;
  RBParameters mu_max_in;
  for (unsigned int i=0; i<n_continuous_parameters; i++)
    {
      // Read in the parameter names
      std::string param_name = infile("parameter_names", "NONE", i);

      {
        Real min_val = infile(param_name, 0., 0);
        mu_min_in.set_value(param_name, min_val);
      }

      {
        Real max_val = infile(param_name, 0., 1);
        mu_max_in.set_value(param_name, max_val);
      }
    }

  std::map<std::string, std::vector<Real>> discrete_parameter_values_in;

  unsigned int n_discrete_parameters = infile.vector_variable_size("discrete_parameter_names");
  for (unsigned int i=0; i<n_discrete_parameters; i++)
    {
      std::string param_name = infile("discrete_parameter_names", "NONE", i);

      unsigned int n_vals_for_param = infile.vector_variable_size(param_name);
      std::vector<Real> vals_for_param(n_vals_for_param);
      for (auto j : make_range(vals_for_param.size()))
        vals_for_param[j] = infile(param_name, 0., j);

      discrete_parameter_values_in[param_name] = vals_for_param;
    }

  std::map<std::string,bool> log_scaling_in;
  // For now, just set all entries to false.
  // TODO: Implement a decent way to specify log-scaling true/false
  // in the input text file
  for (const auto & pr : mu_min_in)
    log_scaling_in[pr.first] = false;

  // Set the parameters that have been read in
  set_rb_construction_parameters(n_training_samples,
                                 deterministic_training,
                                 training_parameters_random_seed_in,
                                 quiet_mode_in,
                                 Nmax_in,
                                 rel_training_tolerance_in,
                                 abs_training_tolerance_in,
                                 normalize_rb_bound_in_greedy_in,
                                 RB_training_type_in,
                                 mu_min_in,
                                 mu_max_in,
                                 discrete_parameter_values_in,
                                 log_scaling_in);
}

void RBConstruction::set_rb_construction_parameters(
                                                    unsigned int n_training_samples_in,
                                                    bool deterministic_training_in,
                                                    unsigned int training_parameters_random_seed_in,
                                                    bool quiet_mode_in,
                                                    unsigned int Nmax_in,
                                                    Real rel_training_tolerance_in,
                                                    Real abs_training_tolerance_in,
                                                    bool normalize_rb_bound_in_greedy_in,
                                                    const std::string & RB_training_type_in,
                                                    RBParameters mu_min_in,
                                                    RBParameters mu_max_in,
                                                    std::map<std::string, std::vector<Real>> discrete_parameter_values_in,
                                                    std::map<std::string,bool> log_scaling_in,
                                                    std::map<std::string, std::vector<Number>> * training_sample_list)
{
  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  set_training_random_seed(training_parameters_random_seed_in);

  // Set quiet mode
  set_quiet_mode(quiet_mode_in);

  // Initialize RB parameters
  set_Nmax(Nmax_in);

  set_rel_training_tolerance(rel_training_tolerance_in);
  set_abs_training_tolerance(abs_training_tolerance_in);

  set_normalize_rb_bound_in_greedy(normalize_rb_bound_in_greedy_in);

  set_RB_training_type(RB_training_type_in);

  // Initialize the parameter ranges and the parameters themselves
  initialize_parameters(mu_min_in, mu_max_in, discrete_parameter_values_in);

  initialize_training_parameters(this->get_parameters_min(),
                                 this->get_parameters_max(),
                                 n_training_samples_in,
                                 log_scaling_in,
                                 deterministic_training_in);   // use deterministic parameters

  if (training_sample_list)
    {
      // Note that we must call initialize_training_parameters() before
      // load_training_set() in order to initialize the parameter vectors.
      load_training_set(*training_sample_list);
    }
}

void RBConstruction::print_info()
{
  // Print out info that describes the current setup
  libMesh::out << std::endl << "RBConstruction parameters:" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "Nmax: " << Nmax << std::endl;
  libMesh::out << "Greedy relative error tolerance: " << get_rel_training_tolerance() << std::endl;
  libMesh::out << "Greedy absolute error tolerance: " << get_abs_training_tolerance() << std::endl;
  libMesh::out << "Do we normalize RB error bound in greedy? " << get_normalize_rb_bound_in_greedy() << std::endl;
  libMesh::out << "RB training type: " << get_RB_training_type() << std::endl;
  if (is_rb_eval_initialized())
    {
      libMesh::out << "Aq operators attached: " << get_rb_theta_expansion().get_n_A_terms() << std::endl;
      libMesh::out << "Fq functions attached: " << get_rb_theta_expansion().get_n_F_terms() << std::endl;
      libMesh::out << "n_outputs: " << get_rb_theta_expansion().get_n_outputs() << std::endl;
      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        libMesh::out << "output " << n << ", Q_l = " << get_rb_theta_expansion().get_n_output_terms(n) << std::endl;
    }
  else
    {
      libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
    }
  libMesh::out << "Number of parameters: " << get_n_params() << std::endl;
  for (const auto & pr : get_parameters())
    if (!is_discrete_parameter(pr.first))
      {
        libMesh::out <<   "Parameter " << pr.first
                     << ": Min = " << get_parameter_min(pr.first)
                     << ", Max = " << get_parameter_max(pr.first) << std::endl;
      }

  print_discrete_parameter_values();
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << "quiet mode? " << is_quiet() << std::endl;
  libMesh::out << std::endl;
}

void RBConstruction::print_basis_function_orthogonality()
{
  std::unique_ptr<NumericVector<Number>> temp = solution->clone();

  for (unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
    {
      for (unsigned int j=0; j<get_rb_evaluation().get_n_basis_functions(); j++)
        {
          get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));
          Number value = temp->dot( get_rb_evaluation().get_basis_function(i) );

          libMesh::out << value << " ";
        }
      libMesh::out << std::endl;
    }
  libMesh::out << std::endl;
}

void RBConstruction::set_rb_assembly_expansion(RBAssemblyExpansion & rb_assembly_expansion_in)
{
  rb_assembly_expansion = &rb_assembly_expansion_in;
}

RBAssemblyExpansion & RBConstruction::get_rb_assembly_expansion()
{
  libmesh_error_msg_if(!rb_assembly_expansion, "Error: RBAssemblyExpansion object hasn't been initialized yet");

  return *rb_assembly_expansion;
}

void RBConstruction::set_inner_product_assembly(ElemAssembly & inner_product_assembly_in)
{
  use_energy_inner_product = false;
  inner_product_assembly = &inner_product_assembly_in;
}

ElemAssembly & RBConstruction::get_inner_product_assembly()
{
  libmesh_error_msg_if(use_energy_inner_product,
                       "Error: inner_product_assembly not available since we're using energy inner-product");

  libmesh_error_msg_if(!inner_product_assembly,
                       "Error: inner_product_assembly hasn't been initialized yet");

  return *inner_product_assembly;
}

void RBConstruction::set_energy_inner_product(const std::vector<Number> & energy_inner_product_coeffs_in)
{
  use_energy_inner_product = true;
  energy_inner_product_coeffs = energy_inner_product_coeffs_in;
}

void RBConstruction::zero_constrained_dofs_on_vector(NumericVector<Number> & vector)
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  const DofMap & dof_map = get_dof_map();

  for (dof_id_type i=dof_map.first_dof(); i<dof_map.end_dof(); i++)
    {
      if (get_dof_map().is_constrained_dof(i))
        {
          vector.set(i, 0.);
        }
    }
#endif

  vector.close();
}

bool RBConstruction::check_if_zero_truth_solve()
{
  return (solution->l2_norm() == 0.);
}

void RBConstruction::initialize_rb_construction(bool skip_matrix_assembly,
                                                bool skip_vector_assembly)
{
  if (!skip_matrix_assembly && !skip_vector_assembly)
    {
      // Check that the theta and assembly objects are consistently sized
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_A_terms(), get_rb_assembly_expansion().get_n_A_terms());
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_F_terms(), get_rb_assembly_expansion().get_n_F_terms());
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_outputs(), get_rb_assembly_expansion().get_n_outputs());
      for (unsigned int i=0; i<get_rb_theta_expansion().get_n_outputs(); i++)
        {
          libmesh_assert_equal_to (get_rb_theta_expansion().get_n_output_terms(i),
                                   get_rb_assembly_expansion().get_n_output_terms(i));
        }
    }

  // Perform the initialization
  allocate_data_structures();
  assemble_affine_expansion(skip_matrix_assembly, skip_vector_assembly);

  // inner_product_solver performs solves with the same matrix every time
  // hence we can set reuse_preconditioner(true).
  inner_product_solver->reuse_preconditioner(true);

  // The primary solver is used for truth solves and other solves that
  // require different matrices, so set reuse_preconditioner(false).
  get_linear_solver()->reuse_preconditioner(false);

}

void RBConstruction::assemble_affine_expansion(bool skip_matrix_assembly,
                                               bool skip_vector_assembly)
{
  if (!skip_matrix_assembly)
    {
      // Assemble and store all of the matrices
      this->assemble_misc_matrices();
      this->assemble_all_affine_operators();
    }

  if (!skip_vector_assembly)
    {
      // Assemble and store all of the vectors
      this->assemble_all_affine_vectors();
      this->assemble_all_output_vectors();
    }
}

void RBConstruction::allocate_data_structures()
{
  // Resize vectors for storing mesh-dependent data
  Aq_vector.resize(get_rb_theta_expansion().get_n_A_terms());
  Fq_vector.resize(get_rb_theta_expansion().get_n_F_terms());

  // Resize the Fq_representors and initialize each to nullptr.
  // These are basis independent and hence stored here, whereas
  // the Aq_representors are stored in RBEvaluation
  Fq_representor.resize(get_rb_theta_expansion().get_n_F_terms());

  // Initialize vectors for the inner products of the Fq representors
  // These are basis independent and therefore stored here.
  unsigned int Q_f_hat = get_rb_theta_expansion().get_n_F_terms()*(get_rb_theta_expansion().get_n_F_terms()+1)/2;
  Fq_representor_innerprods.resize(Q_f_hat);

  // Resize the output vectors
  outputs_vector.resize(get_rb_theta_expansion().get_n_outputs());
  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    outputs_vector[n].resize( get_rb_theta_expansion().get_n_output_terms(n) );

  // Resize the output dual norm vectors
  output_dual_innerprods.resize(get_rb_theta_expansion().get_n_outputs());
  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      unsigned int Q_l_hat = get_rb_theta_expansion().get_n_output_terms(n)*(get_rb_theta_expansion().get_n_output_terms(n)+1)/2;
      output_dual_innerprods[n].resize(Q_l_hat);
    }

  {
    DofMap & dof_map = this->get_dof_map();

    dof_map.attach_matrix(*inner_product_matrix);
    inner_product_matrix->init();
    inner_product_matrix->zero();

    if(store_non_dirichlet_operators)
      {
        // We also need a non-Dirichlet inner-product matrix
        non_dirichlet_inner_product_matrix = SparseMatrix<Number>::build(this->comm());
        dof_map.attach_matrix(*non_dirichlet_inner_product_matrix);
        non_dirichlet_inner_product_matrix->init();
        non_dirichlet_inner_product_matrix->zero();
      }

    for (unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
      {
        // Initialize the memory for the matrices
        Aq_vector[q] = SparseMatrix<Number>::build(this->comm());
        dof_map.attach_matrix(*Aq_vector[q]);
        Aq_vector[q]->init();
        Aq_vector[q]->zero();
      }

    // We also need to initialize a second set of non-Dirichlet operators
    if (store_non_dirichlet_operators)
      {
        non_dirichlet_Aq_vector.resize(get_rb_theta_expansion().get_n_A_terms());
        for (unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
          {
            // Initialize the memory for the matrices
            non_dirichlet_Aq_vector[q] = SparseMatrix<Number>::build(this->comm());
            dof_map.attach_matrix(*non_dirichlet_Aq_vector[q]);
            non_dirichlet_Aq_vector[q]->init();
            non_dirichlet_Aq_vector[q]->zero();
          }
      }
  }

  // Initialize the vectors
  for (unsigned int q=0; q<get_rb_theta_expansion().get_n_F_terms(); q++)
    {
      // Initialize the memory for the vectors
      Fq_vector[q] = NumericVector<Number>::build(this->comm());
      Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    }

  // We also need to initialize a second set of non-Dirichlet operators
  if (store_non_dirichlet_operators)
    {
      non_dirichlet_Fq_vector.resize(get_rb_theta_expansion().get_n_F_terms());
      for (unsigned int q=0; q<get_rb_theta_expansion().get_n_F_terms(); q++)
        {
          // Initialize the memory for the vectors
          non_dirichlet_Fq_vector[q] = NumericVector<Number>::build(this->comm());
          non_dirichlet_Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
        }
    }

  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        // Initialize the memory for the truth output vectors
        outputs_vector[n][q_l] = NumericVector<Number>::build(this->comm());
        outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      }

  if (store_non_dirichlet_operators)
    {
      non_dirichlet_outputs_vector.resize(get_rb_theta_expansion().get_n_outputs());
      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        {
          non_dirichlet_outputs_vector[n].resize( get_rb_theta_expansion().get_n_output_terms(n) );
          for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
            {
              // Initialize the memory for the truth output vectors
              non_dirichlet_outputs_vector[n][q_l] = NumericVector<Number>::build(this->comm());
              non_dirichlet_outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }
        }
    }

  // Resize truth_outputs vector
  truth_outputs.resize(this->get_rb_theta_expansion().get_n_outputs());
}

std::unique_ptr<DGFEMContext> RBConstruction::build_context ()
{
  return libmesh_make_unique<DGFEMContext>(*this);
}

void RBConstruction::add_scaled_matrix_and_vector(Number scalar,
                                                  ElemAssembly * elem_assembly,
                                                  SparseMatrix<Number> * input_matrix,
                                                  NumericVector<Number> * input_vector,
                                                  bool symmetrize,
                                                  bool apply_dof_constraints)
{
  LOG_SCOPE("add_scaled_matrix_and_vector()", "RBConstruction");

  bool assemble_matrix = (input_matrix != nullptr);
  bool assemble_vector = (input_vector != nullptr);

  if (!assemble_matrix && !assemble_vector)
    return;

  const MeshBase & mesh = this->get_mesh();

  // First add any node-based terms (e.g. point loads)

  // Make a std::set of all the nodes that are in 1 or more
  // nodesets. We only want to call get_nodal_values() once per Node
  // per ElemAssembly object, regardless of how many nodesets it
  // appears in.
  std::set<dof_id_type> nodes_with_nodesets;
  for (const auto & t : mesh.get_boundary_info().build_node_list())
    nodes_with_nodesets.insert(std::get<0>(t));

  for (const auto & id : nodes_with_nodesets)
    {
      const Node & node = mesh.node_ref(id);

      // If node is on this processor, then all dofs on node are too
      // so we can do the add below safely
      if (node.processor_id() == this->comm().rank())
        {
          // Get the values to add to the rhs vector
          std::vector<dof_id_type> nodal_dof_indices;
          DenseMatrix<Number> nodal_matrix;
          DenseVector<Number> nodal_rhs;
          elem_assembly->get_nodal_values(nodal_dof_indices,
                                          nodal_matrix,
                                          nodal_rhs,
                                          *this,
                                          node);

          // Perform any required user-defined postprocessing on
          // the matrix and rhs.
          //
          // TODO: We need to postprocess node matrices and vectors
          // in some cases (e.g. when rotations are applied to
          // nodes), but since we don't have a FEMContext at this
          // point we would need to have a different interface
          // taking the DenseMatrix, DenseVector, and probably the
          // current node that we are on...
          // this->post_process_elem_matrix_and_vector(nodal_matrix, nodal_rhs);

          if (!nodal_dof_indices.empty())
            {
              if (assemble_vector)
                input_vector->add_vector(nodal_rhs, nodal_dof_indices);

              if (assemble_matrix)
                input_matrix->add_matrix(nodal_matrix, nodal_dof_indices);
            }
        }
    }

  std::unique_ptr<DGFEMContext> c = this->build_context();
  DGFEMContext & context  = cast_ref<DGFEMContext &>(*c);

  this->init_context(context);

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      if(elem->type() == NODEELEM)
        {
          // We skip NODEELEMs here since we assume that we
          // do not perform any assembly directly on NODEELEMs
          continue;
        }

      // Subdivision elements need special care:
      // - skip ghost elements
      // - init special quadrature rule
      std::unique_ptr<QBase> qrule;
      if (elem->type() == TRI3SUBDIVISION)
        {
          const Tri3Subdivision * gh_elem = static_cast<const Tri3Subdivision *> (elem);
          if (gh_elem->is_ghost())
            continue ;
          // A Gauss quadrature rule for numerical integration.
          // For subdivision shell elements, a single Gauss point per
          // element is sufficient, hence we use extraorder = 0.
          const int extraorder = 0;
          FEBase * elem_fe = nullptr;
          context.get_element_fe( 0, elem_fe );

          qrule = elem_fe->get_fe_type().default_quadrature_rule (2, extraorder);

          // Tell the finite element object to use our quadrature rule.
          elem_fe->attach_quadrature_rule (qrule.get());
        }

      context.pre_fe_reinit(*this, elem);
      context.elem_fe_reinit();
      elem_assembly->interior_assembly(context);

      const unsigned char n_sides = context.get_elem().n_sides();
      for (context.side = 0; context.side != n_sides; ++context.side)
        {
          // May not need to apply fluxes on non-boundary elements
          if ((context.get_elem().neighbor_ptr(context.get_side()) != nullptr) && !impose_internal_fluxes)
            continue;

          // skip degenerate sides with zero area
          if( (context.get_elem().side_ptr(context.get_side())->volume() <= 0.) && skip_degenerate_sides)
            continue;

          context.side_fe_reinit();
          elem_assembly->boundary_assembly(context);

          if (context.dg_terms_are_active())
            {
              input_matrix->add_matrix (context.get_elem_elem_jacobian(),
                                        context.get_dof_indices(),
                                        context.get_dof_indices());

              input_matrix->add_matrix (context.get_elem_neighbor_jacobian(),
                                        context.get_dof_indices(),
                                        context.get_neighbor_dof_indices());

              input_matrix->add_matrix (context.get_neighbor_elem_jacobian(),
                                        context.get_neighbor_dof_indices(),
                                        context.get_dof_indices());

              input_matrix->add_matrix (context.get_neighbor_neighbor_jacobian(),
                                        context.get_neighbor_dof_indices(),
                                        context.get_neighbor_dof_indices());
            }
        }

      // Do any required user post-processing before symmetrizing and/or applying
      // constraints.
      //
      // We only do this if apply_dof_constraints is true because we want to be
      // able to set apply_dof_constraints=false in order to obtain a matrix
      // A with no dof constraints or dof transformations, as opposed to C^T A C,
      // which includes constraints and/or dof transformations. Here C refers to
      // the matrix that imposes dof constraints and transformations on the
      // solution u.
      //
      // Matrices such as A are what we store in our "non_dirichlet" operators, and
      // they are useful for computing terms such as (C u_i)^T A (C u_j) (e.g. see
      // update_RB_system_matrices()),  where C u is the result of a "truth_solve",
      // which includes calls to both enforce_constraints_exactly() and
      // post_process_truth_solution(). If we use C^T A C to compute these terms then
      // we would "double apply" the matrix C, which can give incorrect results.
      if (apply_dof_constraints)
        this->post_process_elem_matrix_and_vector(context);

      // Need to symmetrize before imposing
      // periodic constraints
      if (assemble_matrix && symmetrize)
        {
          DenseMatrix<Number> Ke_transpose;
          context.get_elem_jacobian().get_transpose(Ke_transpose);
          context.get_elem_jacobian() += Ke_transpose;
          context.get_elem_jacobian() *= 0.5;
        }

      // As discussed above, we can set apply_dof_constraints=false to
      // get A instead of C^T A C
      if (apply_dof_constraints)
        {
          // Apply constraints, e.g. Dirichlet and periodic constraints
          this->get_dof_map().constrain_element_matrix_and_vector
            (context.get_elem_jacobian(),
             context.get_elem_residual(),
             context.get_dof_indices(),
             /*asymmetric_constraint_rows*/ false );
        }

      // Scale and add to global matrix and/or vector
      context.get_elem_jacobian() *= scalar;
      context.get_elem_residual() *= scalar;

      if (assemble_matrix)
        {

          CouplingMatrix * coupling_matrix = get_dof_map()._dof_coupling;
          if (!coupling_matrix)
            {
              // If we haven't defined a _dof_coupling matrix then just add
              // the whole matrix
              input_matrix->add_matrix (context.get_elem_jacobian(),
                                        context.get_dof_indices() );
            }
          else
            {
              // Otherwise we should only add the relevant submatrices
              for (unsigned int var1=0; var1<n_vars(); var1++)
                {
                  ConstCouplingRow ccr(var1, *coupling_matrix);
                  for (const auto & var2 : ccr)
                    {
                      unsigned int sub_m = context.get_elem_jacobian( var1, var2 ).m();
                      unsigned int sub_n = context.get_elem_jacobian( var1, var2 ).n();
                      DenseMatrix<Number> sub_jac(sub_m, sub_n);
                      for (unsigned int row=0; row<sub_m; row++)
                        for (unsigned int col=0; col<sub_n; col++)
                          {
                            sub_jac(row,col) = context.get_elem_jacobian( var1, var2 ).el(row,col);
                          }
                      input_matrix->add_matrix (sub_jac,
                                                context.get_dof_indices(var1),
                                                context.get_dof_indices(var2) );
                    }
                }
            }

        }

      if (assemble_vector)
        input_vector->add_vector (context.get_elem_residual(),
                                  context.get_dof_indices() );
    }

  if (assemble_matrix)
    input_matrix->close();
  if (assemble_vector)
    input_vector->close();
}

void RBConstruction::set_context_solution_vec(NumericVector<Number> & vec)
{
  // Set current_local_solution = vec so that we can access
  // vec from DGFEMContext during assembly
  vec.localize
    (*current_local_solution, this->get_dof_map().get_send_list());
}

void RBConstruction::truth_assembly()
{
  LOG_SCOPE("truth_assembly()", "RBConstruction");

  const RBParameters & mu = get_parameters();

  this->matrix->zero();
  this->rhs->zero();

  this->matrix->close();
  this->rhs->close();

  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
      {
        matrix->add(get_rb_theta_expansion().eval_A_theta(q_a, mu), *get_Aq(q_a));
      }

    std::unique_ptr<NumericVector<Number>> temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
      {
        *temp_vec = *get_Fq(q_f);
        temp_vec->scale( get_rb_theta_expansion().eval_F_theta(q_f, mu) );
        rhs->add(*temp_vec);
      }
  }

  this->matrix->close();
  this->rhs->close();
}

void RBConstruction::assemble_inner_product_matrix(SparseMatrix<Number> * input_matrix,
                                                   bool apply_dof_constraints)
{
  input_matrix->zero();

  if(!use_energy_inner_product)
  {
    add_scaled_matrix_and_vector(1.,
                                 inner_product_assembly,
                                 input_matrix,
                                 nullptr,
                                 false, /* symmetrize */
                                 apply_dof_constraints);
  }
  else
  {
    libmesh_error_msg_if(energy_inner_product_coeffs.size() != get_rb_theta_expansion().get_n_A_terms(),
                         "Error: invalid number of entries in energy_inner_product_coeffs.");

    // We symmetrize below so that we may use the energy inner-product even in cases
    // where the A_q are not symmetric.
    for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
      {
        add_scaled_matrix_and_vector(energy_inner_product_coeffs[q_a],
                                     &rb_assembly_expansion->get_A_assembly(q_a),
                                     input_matrix,
                                     nullptr,
                                     true, /* symmetrize */
                                     apply_dof_constraints);
      }
  }
}

void RBConstruction::assemble_Aq_matrix(unsigned int q,
                                        SparseMatrix<Number> * input_matrix,
                                        bool apply_dof_constraints)
{
  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_A_terms(),
                       "Error: We must have q < Q_a in assemble_Aq_matrix.");

  input_matrix->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_A_assembly(q),
                               input_matrix,
                               nullptr,
                               false, /* symmetrize */
                               apply_dof_constraints);
}

void RBConstruction::add_scaled_Aq(Number scalar,
                                   unsigned int q_a,
                                   SparseMatrix<Number> * input_matrix,
                                   bool symmetrize)
{
  LOG_SCOPE("add_scaled_Aq()", "RBConstruction");

  libmesh_error_msg_if(q_a >= get_rb_theta_expansion().get_n_A_terms(),
                       "Error: We must have q < Q_a in add_scaled_Aq.");

  if (!symmetrize)
    {
      input_matrix->add(scalar, *get_Aq(q_a));
      input_matrix->close();
    }
  else
    {
      add_scaled_matrix_and_vector(scalar,
                                   &rb_assembly_expansion->get_A_assembly(q_a),
                                   input_matrix,
                                   nullptr,
                                   symmetrize);
    }
}

void RBConstruction::assemble_misc_matrices()
{
  libMesh::out << "Assembling inner product matrix" << std::endl;
  assemble_inner_product_matrix(inner_product_matrix.get());

  if (store_non_dirichlet_operators)
    {
      libMesh::out << "Assembling non-Dirichlet inner product matrix" << std::endl;
      assemble_inner_product_matrix(non_dirichlet_inner_product_matrix.get(), false);
    }
}

void RBConstruction::assemble_all_affine_operators()
{
  for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      libMesh::out << "Assembling affine operator " << (q_a+1) << " of "
                   << get_rb_theta_expansion().get_n_A_terms() << std::endl;
      assemble_Aq_matrix(q_a, get_Aq(q_a));
    }

  if (store_non_dirichlet_operators)
    {
      for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
        {
          libMesh::out << "Assembling non-Dirichlet affine operator " << (q_a+1) << " of "
                       << get_rb_theta_expansion().get_n_A_terms() << std::endl;
          assemble_Aq_matrix(q_a, get_non_dirichlet_Aq(q_a), false);
        }
    }
}

void RBConstruction::assemble_all_affine_vectors()
{
  for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      libMesh::out << "Assembling affine vector " << (q_f+1) << " of "
                   << get_rb_theta_expansion().get_n_F_terms() << std::endl;
      assemble_Fq_vector(q_f, get_Fq(q_f));
    }

  if (store_non_dirichlet_operators)
    {
      for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          libMesh::out << "Assembling non-Dirichlet affine vector " << (q_f+1) << " of "
                       << get_rb_theta_expansion().get_n_F_terms() << std::endl;
          assemble_Fq_vector(q_f, get_non_dirichlet_Fq(q_f), false);
        }
    }

}

void RBConstruction::assemble_Fq_vector(unsigned int q,
                                        NumericVector<Number> * input_vector,
                                        bool apply_dof_constraints)
{
  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_F_terms(),
                       "Error: We must have q < Q_f in assemble_Fq_vector.");

  input_vector->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_F_assembly(q),
                               nullptr,
                               input_vector,
                               false,             /* symmetrize */
                               apply_dof_constraints /* apply_dof_constraints */);
}

void RBConstruction::assemble_all_output_vectors()
{
  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        libMesh::out << "Assembling output vector, (" << (n+1) << "," << (q_l+1)
                     << ") of (" << get_rb_theta_expansion().get_n_outputs()
                     << "," << get_rb_theta_expansion().get_n_output_terms(n) << ")"
                     << std::endl;
        get_output_vector(n, q_l)->zero();
        add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                     nullptr,
                                     get_output_vector(n,q_l),
                                     false, /* symmetrize */
                                     true   /* apply_dof_constraints */);
      }

  if (store_non_dirichlet_operators)
    {
      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
          {
            libMesh::out << "Assembling non-Dirichlet output vector, (" << (n+1) << "," << (q_l+1)
                         << ") of (" << get_rb_theta_expansion().get_n_outputs()
                         << "," << get_rb_theta_expansion().get_n_output_terms(n) << ")"
                         << std::endl;
            get_non_dirichlet_output_vector(n, q_l)->zero();
            add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                         nullptr,
                                         get_non_dirichlet_output_vector(n,q_l),
                                         false, /* symmetrize */
                                         false  /* apply_dof_constraints */);
          }
    }
}

Real RBConstruction::train_reduced_basis(const bool resize_rb_eval_data)
{
  if(get_RB_training_type() == "Greedy")
    {
      return train_reduced_basis_with_greedy(resize_rb_eval_data);
    }
  else if (get_RB_training_type() == "POD")
    {
      train_reduced_basis_with_POD();
      return 0.;
    }
  else
    {
      libmesh_error_msg("RB training type not recognized: " + get_RB_training_type());
    }

  return 0.;
}

Real RBConstruction::train_reduced_basis_with_greedy(const bool resize_rb_eval_data)
{
  LOG_SCOPE("train_reduced_basis_with_greedy()", "RBConstruction");

  int count = 0;

  RBEvaluation & rbe = get_rb_evaluation();

  // initialize rbe's parameters
  rbe.initialize_parameters(*this);

  // possibly resize data structures according to Nmax
  if (resize_rb_eval_data)
    rbe.resize_data_structures(get_Nmax());

  // Clear the Greedy param list
  for (auto & plist : rbe.greedy_param_list)
    plist.clear();

  rbe.greedy_param_list.clear();

  Real training_greedy_error = 0.;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if (rbe.get_n_basis_functions() >= get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << std::endl;
      return 0.;
    }

  // Optionally pre-evaluate the theta functions on the entire (local) training parameter set.
  if (get_preevaluate_thetas_flag())
    preevaluate_thetas();

  if(!skip_residual_in_train_reduced_basis)
    {
      // Compute the dual norms of the outputs if we haven't already done so.
      compute_output_dual_innerprods();

      // Compute the Fq Riesz representor dual norms if we haven't already done so.
      compute_Fq_representor_innerprods();
    }

  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
  Real initial_greedy_error = 0.;
  bool initial_greedy_error_initialized = false;
  while (true)
    {
      libMesh::out << std::endl << "---- Basis dimension: "
                   << rbe.get_n_basis_functions() << " ----" << std::endl;

      if (count > 0 || (count==0 && use_empty_rb_solve_in_greedy))
        {
          libMesh::out << "Performing RB solves on training set" << std::endl;
          training_greedy_error = compute_max_error_bound();

          libMesh::out << "Maximum error bound is " << training_greedy_error << std::endl << std::endl;

          // record the initial error
          if (!initial_greedy_error_initialized)
            {
              initial_greedy_error = training_greedy_error;
              initial_greedy_error_initialized = true;
            }

          // Break out of training phase if we have reached Nmax
          // or if the training_tolerance is satisfied.
          if (greedy_termination_test(training_greedy_error, initial_greedy_error, count))
            break;
        }

      libMesh::out << "Performing truth solve at parameter:" << std::endl;
      print_parameters();

      // Update the list of Greedily selected parameters
      this->update_greedy_param_list();

      // Perform an Offline truth solve for the current parameter
      truth_solve(-1);

      if (check_if_zero_truth_solve())
        {
          libMesh::out << "Zero basis function encountered hence ending basis enrichment" << std::endl;
          break;
        }

      // Add orthogonal part of the snapshot to the RB space
      libMesh::out << "Enriching the RB space" << std::endl;
      enrich_RB_space();

      update_system();

      // Check if we've reached Nmax now. We do this before calling
      // update_residual_terms() since we can skip that step if we've
      // already reached Nmax.
      if (rbe.get_n_basis_functions() >= this->get_Nmax())
      {
        libMesh::out << "Maximum number of basis functions reached: Nmax = "
                     << get_Nmax() << std::endl;
        break;
      }

      if(!skip_residual_in_train_reduced_basis)
        {
          update_residual_terms();
        }

      // Increment counter
      count++;
    }
  this->update_greedy_param_list();

  return training_greedy_error;
}

void RBConstruction::enrich_basis_from_rhs_terms(const bool resize_rb_eval_data)
{
  LOG_SCOPE("enrich_basis_from_rhs_terms()", "RBConstruction");

  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  // possibly resize data structures according to Nmax
  if (resize_rb_eval_data)
    {
      get_rb_evaluation().resize_data_structures(get_Nmax());
    }

  libMesh::out << std::endl << "---- Enriching basis from rhs terms ----" << std::endl;

  truth_assembly();

  for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      libMesh::out << std::endl << "Performing truth solve with rhs from rhs term " << q_f << std::endl;

      *rhs = *get_Fq(q_f);

      if (rhs->l2_norm() == 0)
      {
        // Skip enrichment if the rhs is zero
        continue;
      }

      // truth_assembly assembles into matrix and rhs, so use those for the solve
      if (extra_linear_solver)
        {
          // If extra_linear_solver has been initialized, then we use it for the
          // truth solves.
          solve_for_matrix_and_rhs(*extra_linear_solver, *matrix, *rhs);

          if (assert_convergence)
            check_convergence(*extra_linear_solver);
        }
      else
        {
          solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);

          if (assert_convergence)
            check_convergence(*get_linear_solver());
        }

      // Debugging: enable this code to print the rhs that was used in
      // the most recent truth solve to a uniquely-named file.
      if (false)
        {
          char temp_file[] = "truth_rhs_XXXXXX.dat";
          int fd = mkstemps(temp_file, 4);
          if (fd != -1)
            {
              libMesh::out << "Writing truth system rhs to file: " << temp_file << std::endl;
              rhs->print_matlab(std::string(temp_file));
            }
        }

      // Debugging: enable this code to print the most recent truth
      // solution to a uniquely-named file.
#ifdef LIBMESH_HAVE_EXODUS_API
      if (false)
        {
          // Note: mkstemps creates a file and returns an open file descriptor to it.
          // The filename is created from a template which must have 6 'X' characters followed
          // by a suffix having the specified length (in this case 4, for ".exo").
          char temp_file[] = "truth_XXXXXX.exo";
          int fd = mkstemps(temp_file, 4);
          if (fd != -1)
            {
              libMesh::out << "Writing truth solution to file: " << temp_file << std::endl;
              ExodusII_IO exo_io(this->get_mesh());
              std::set<std::string> system_names = {this->name()};
              exo_io.write_equation_systems(std::string(temp_file),
                                            this->get_equation_systems(), &system_names);
            }
        }
#endif

      // Call user-defined post-processing routines on the truth solution.
      post_process_truth_solution();

      // Add orthogonal part of the snapshot to the RB space
      libMesh::out << "Enriching the RB space" << std::endl;
      enrich_RB_space();

      update_system();
    }
}

void RBConstruction::train_reduced_basis_with_POD()
{
  // We need to use the same training set on all processes so that
  // the truth solves below work correctly in parallel.
  libmesh_error_msg_if(!serial_training_set, "We must use a serial training set with POD");
  libmesh_error_msg_if(get_rb_evaluation().get_n_basis_functions() > 0, "Basis should not already be initialized");

  get_rb_evaluation().initialize_parameters(*this);
  get_rb_evaluation().resize_data_structures(get_Nmax());

  // Storage for the POD snapshots
  unsigned int n_snapshots = get_n_training_samples();

  if (get_n_params() == 0)
    {
      // In this case we should have generated an empty training set
      // so assert this
      libmesh_assert(n_snapshots == 0);

      // If we have no parameters, then we should do exactly one "truth solve"
      n_snapshots = 1;
    }

  std::vector<std::unique_ptr<NumericVector<Number>>> POD_snapshots(n_snapshots);
  for (unsigned int i=0; i<n_snapshots; i++)
    {
      POD_snapshots[i] = NumericVector<Number>::build(this->comm());
      POD_snapshots[i]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    }

  // We use the same training set on all processes
  libMesh::out << std::endl;
  for (unsigned int i=0; i<n_snapshots; i++)
    {
      if (get_n_params() > 0)
        {
          set_params_from_training_set(i);
        }

      libMesh::out << "Truth solve " << (i+1) << " of " << n_snapshots << std::endl;

      truth_solve(-1);

      *POD_snapshots[i] = *solution;
    }
  libMesh::out << std::endl;

  // Set up the "correlation matrix"
  DenseMatrix<Number> correlation_matrix(n_snapshots,n_snapshots);
  for (unsigned int i=0; i<n_snapshots; i++)
    {
      get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(
        *inner_product_storage_vector, *POD_snapshots[i]);

      for (unsigned int j=0; j<=i; j++)
        {
          Number inner_prod = (POD_snapshots[j]->dot(*inner_product_storage_vector));

          correlation_matrix(i,j) = inner_prod;
          if(i != j)
            {
              correlation_matrix(j,i) = libmesh_conj(inner_prod);
            }
        }
    }

  // compute SVD of correlation matrix
  DenseVector<Real> sigma( n_snapshots );
  DenseMatrix<Number> U( n_snapshots, n_snapshots );
  DenseMatrix<Number> VT( n_snapshots, n_snapshots );
  correlation_matrix.svd(sigma, U, VT );

  libmesh_error_msg_if(sigma(0) == 0., "Zero singular value encountered in POD construction");

  // Add dominant vectors from the POD as basis functions.
  unsigned int j = 0;
  while (true)
    {
      if (j >= get_Nmax() || j >= n_snapshots)
        {
          libMesh::out << "Maximum number of basis functions (" << j << ") reached." << std::endl;
          break;
        }

      // The "energy" error in the POD approximation is determined by the first omitted
      // singular value, i.e. sigma(j). We normalize by sigma(0), which gives the total
      // "energy", in order to obtain a relative error.
      const Real rel_err = std::sqrt(sigma(j)) / std::sqrt(sigma(0));

      libMesh::out << "Number of basis functions: " << j
                   << ", POD error norm: " << rel_err << std::endl;

      if (rel_err < this->rel_training_tolerance)
        {
          libMesh::out << "Training tolerance reached." << std::endl;
          break;
        }

      std::unique_ptr< NumericVector<Number> > v = POD_snapshots[j]->zero_clone();
      for ( unsigned int i=0; i<n_snapshots; ++i )
        {
          v->add( U.el(i, j), *POD_snapshots[i] );
        }

      Real norm_v = std::sqrt(sigma(j));
      v->scale( 1./norm_v );

      get_rb_evaluation().basis_functions.emplace_back( std::move(v) );

      j++;
    }
  libMesh::out << std::endl;

  this->delta_N = get_rb_evaluation().get_n_basis_functions();
  update_system();
}

bool RBConstruction::greedy_termination_test(Real abs_greedy_error,
                                             Real initial_error,
                                             int)
{
  if (abs_greedy_error < this->abs_training_tolerance)
    {
      libMesh::out << "Absolute error tolerance reached." << std::endl;
      return true;
    }

  Real rel_greedy_error = abs_greedy_error/initial_error;
  if (rel_greedy_error < this->rel_training_tolerance)
    {
      libMesh::out << "Relative error tolerance reached." << std::endl;
      return true;
    }

  RBEvaluation & rbe = get_rb_evaluation();

  if (rbe.get_n_basis_functions() >= this->get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << std::endl;
      return true;
    }

  if (exit_on_repeated_greedy_parameters)
    for (auto & plist : rbe.greedy_param_list)
      if (plist == get_parameters())
        {
          libMesh::out << "Exiting greedy because the same parameters were selected twice" << std::endl;
          return true;
        }

  return false;
}

void RBConstruction::update_greedy_param_list()
{
  get_rb_evaluation().greedy_param_list.push_back( get_parameters() );
}

const RBParameters & RBConstruction::get_greedy_parameter(unsigned int i)
{
  libmesh_error_msg_if(i >= get_rb_evaluation().greedy_param_list.size(),
                       "Error: Argument in RBConstruction::get_greedy_parameter is too large.");

  return get_rb_evaluation().greedy_param_list[i];
}

Real RBConstruction::truth_solve(int plot_solution)
{
  LOG_SCOPE("truth_solve()", "RBConstruction");

  if(store_untransformed_basis && !_untransformed_solution)
  {
    _untransformed_solution = solution->zero_clone();
  }

  truth_assembly();

  // truth_assembly assembles into matrix and rhs, so use those for the solve
  if (extra_linear_solver)
    {
      // If extra_linear_solver has been initialized, then we use it for the
      // truth solves.
      solve_for_matrix_and_rhs(*extra_linear_solver, *matrix, *rhs);

      if (assert_convergence)
        check_convergence(*extra_linear_solver);
    }
  else
    {
      solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);

      if (assert_convergence)
        check_convergence(*get_linear_solver());
    }

  if (store_untransformed_basis)
    {
      *_untransformed_solution = *solution;
    }

  // Call user-defined post-processing routines on the truth solution.
  post_process_truth_solution();

  const RBParameters & mu = get_parameters();

  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      truth_outputs[n] = 0.;
      for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
        truth_outputs[n] += get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
          get_output_vector(n,q_l)->dot(*solution);
    }

#ifdef LIBMESH_HAVE_EXODUS_API
  if (plot_solution > 0)
    {
      ExodusII_IO exo_io(this->get_mesh());
      std::set<std::string> system_names = {this->name()};
      exo_io.write_equation_systems("truth.exo", this->get_equation_systems(), &system_names);
    }
#endif

  // Get the X norm of the truth solution
  // Useful for normalizing our true error data
  get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *solution);
  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));

  return libmesh_real(truth_X_norm);
}

void RBConstruction::set_RB_training_type(const std::string & RB_training_type_in)
{
  this->RB_training_type = RB_training_type_in;

  if(this->RB_training_type == "POD")
    {
      // We need to use a serial training set (so that the training
      // set is the same on all processes) if we're using POD
      this->serial_training_set = true;
    }
}

const std::string & RBConstruction::get_RB_training_type() const
{
  return RB_training_type;
}

void RBConstruction::set_Nmax(unsigned int Nmax_in)
{
  this->Nmax = Nmax_in;
}

void RBConstruction::load_basis_function(unsigned int i)
{
  LOG_SCOPE("load_basis_function()", "RBConstruction");

  libmesh_assert_less (i, get_rb_evaluation().get_n_basis_functions());

  *solution = get_rb_evaluation().get_basis_function(i);

  this->update();
}

void RBConstruction::enrich_RB_space()
{
  LOG_SCOPE("enrich_RB_space()", "RBConstruction");

  auto new_bf = solution->clone();

  std::unique_ptr<NumericVector<Number>> new_untransformed_bf;
  if (store_untransformed_basis)
    {
      new_untransformed_bf = _untransformed_solution->clone();
      libmesh_assert_equal_to(_untransformed_basis_functions.size(), get_rb_evaluation().get_n_basis_functions());
    }

  for (unsigned int index=0; index<get_rb_evaluation().get_n_basis_functions(); index++)
    {
      get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *new_bf);

      Number scalar =
        inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(index));
      new_bf->add(-scalar, get_rb_evaluation().get_basis_function(index));

      if (store_untransformed_basis)
        {
          new_untransformed_bf->add(-scalar, *_untransformed_basis_functions[index]);
        }
    }

  // Normalize new_bf
  get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *new_bf);
  Number new_bf_norm = std::sqrt( inner_product_storage_vector->dot(*new_bf) );

  if (new_bf_norm == 0.)
    {
      new_bf->zero(); // avoid potential nan's

      if (store_untransformed_basis)
        {
          new_untransformed_bf->zero();
        }
    }
  else
    {
      new_bf->scale(1./new_bf_norm);

      if (store_untransformed_basis)
        {
          new_untransformed_bf->scale(1./new_bf_norm);
        }
    }

  // load the new basis function into the basis_functions vector.
  get_rb_evaluation().basis_functions.emplace_back( std::move(new_bf) );

  if (store_untransformed_basis)
    {
      _untransformed_basis_functions.emplace_back( std::move(new_untransformed_bf) );
    }
}

void RBConstruction::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();
}

Real RBConstruction::get_RB_error_bound()
{
  get_rb_evaluation().set_parameters( get_parameters() );

  Real error_bound = 0.;
  if (get_preevaluate_thetas_flag())
    {
      // Obtain the pre-evaluated theta functions from the current training parameter index
      const auto & evaluated_thetas = get_evaluated_thetas(get_current_training_parameter_index());
      error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions(),
                                                 &evaluated_thetas);
    }
  else
    error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());


  if (normalize_rb_bound_in_greedy)
    {
      Real error_bound_normalization = get_rb_evaluation().get_error_bound_normalization();

      if ((error_bound < abs_training_tolerance) ||
          (error_bound_normalization < abs_training_tolerance))
        {
          // We don't want to normalize this error bound if the bound or the
          // normalization value are below the absolute tolerance. Hence do nothing
          // in this case.
        }
      else
        error_bound /= error_bound_normalization;
    }

  return error_bound;
}

void RBConstruction::recompute_all_residual_terms(bool compute_inner_products)
{
  // Compute the basis independent terms
  Fq_representor_innerprods_computed = false;
  compute_Fq_representor_innerprods(compute_inner_products);

  // and all the basis dependent terms
  unsigned int saved_delta_N = delta_N;
  delta_N = get_rb_evaluation().get_n_basis_functions();

  update_residual_terms(compute_inner_products);

  delta_N = saved_delta_N;
}

Real RBConstruction::compute_max_error_bound()
{
  LOG_SCOPE("compute_max_error_bound()", "RBConstruction");

  // Treat the case with no parameters in a special way
  if (get_n_params() == 0)
    {
      Real max_val;
      if (std::numeric_limits<Real>::has_infinity)
        {
          max_val = std::numeric_limits<Real>::infinity();
        }
      else
        {
          max_val = std::numeric_limits<Real>::max();
        }

      // Make sure we do at least one solve, but otherwise return a zero error bound
      // when we have no parameters
      return (get_rb_evaluation().get_n_basis_functions() == 0) ? max_val : 0.;
    }

  training_error_bounds.resize(this->get_local_n_training_samples());

  // keep track of the maximum error
  unsigned int max_err_index = 0;
  Real max_err = 0.;

  numeric_index_type first_index = get_first_local_training_index();
  for (unsigned int i=0; i<get_local_n_training_samples(); i++)
    {
      // Load training parameter i, this is only loaded
      // locally since the RB solves are local.
      set_params_from_training_set( first_index+i );

      // In case we pre-evaluate the theta functions,
      // also keep track of the current training parameter index.
      if (get_preevaluate_thetas_flag())
        set_current_training_parameter_index(first_index+i);


      training_error_bounds[i] = get_RB_error_bound();

      if (training_error_bounds[i] > max_err)
        {
          max_err_index = i;
          max_err = training_error_bounds[i];
        }
    }

  std::pair<numeric_index_type, Real> error_pair(first_index+max_err_index, max_err);
  get_global_max_error_pair(this->comm(),error_pair);

  // If we have a serial training set (i.e. a training set that is the same on all processors)
  // just set the parameters on all processors
  if (serial_training_set)
    {
      set_params_from_training_set( error_pair.first );
    }
  // otherwise, broadcast the parameter that produced the maximum error
  else
    {
      unsigned int root_id=0;
      if ((get_first_local_training_index() <= error_pair.first) &&
          (error_pair.first < get_last_local_training_index()))
        {
          set_params_from_training_set( error_pair.first );
          root_id = this->processor_id();
        }

      this->comm().sum(root_id); // root_id is only non-zero on one processor
      broadcast_parameters(root_id);
    }

  return error_pair.second;
}

void RBConstruction::update_RB_system_matrices()
{
  LOG_SCOPE("update_RB_system_matrices()", "RBConstruction");

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  std::unique_ptr<NumericVector<Number>> temp = NumericVector<Number>::build(this->comm());
  temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          get_rb_evaluation().RB_Fq_vector[q_f](i) = get_non_dirichlet_Fq_if_avail(q_f)->dot(get_rb_evaluation().get_basis_function(i));
        }
    }

  for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
          {
            get_rb_evaluation().RB_output_vectors[n][q_l](i) =
              get_output_vector(n,q_l)->dot(get_rb_evaluation().get_basis_function(i));
          }

      for (unsigned int j=0; j<RB_size; j++)
        {
          Number value = 0.;

          if (compute_RB_inner_product)
            {
              // Compute reduced inner_product_matrix
              temp->zero();
              get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

              value = temp->dot( get_rb_evaluation().get_basis_function(i) );
              get_rb_evaluation().RB_inner_product_matrix(i,j) = value;
              if (i!=j)
                {
                  // The inner product matrix is assumed
                  // to be hermitian
                  get_rb_evaluation().RB_inner_product_matrix(j,i) = libmesh_conj(value);
                }
            }

          for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
            {
              // Compute reduced Aq matrix
              temp->zero();
              get_non_dirichlet_Aq_if_avail(q_a)->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

              value = (*temp).dot(get_rb_evaluation().get_basis_function(i));
              get_rb_evaluation().RB_Aq_vector[q_a](i,j) = value;

              if (i!=j)
                {
                  temp->zero();
                  get_non_dirichlet_Aq_if_avail(q_a)->vector_mult(*temp, get_rb_evaluation().get_basis_function(i));

                  value = (*temp).dot(get_rb_evaluation().get_basis_function(j));
                  get_rb_evaluation().RB_Aq_vector[q_a](j,i) = value;
                }
            }
        }
    }
}


void RBConstruction::update_residual_terms(bool compute_inner_products)
{
  LOG_SCOPE("update_residual_terms()", "RBConstruction");

  libMesh::out << "Updating RB residual terms" << std::endl;

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  if (store_untransformed_basis)
    {
      libmesh_assert_equal_to(_untransformed_basis_functions.size(), get_rb_evaluation().get_n_basis_functions());
    }

  for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          // Initialize the vector in which we'll store the representor
          if (!get_rb_evaluation().Aq_representor[q_a][i])
            {
              get_rb_evaluation().Aq_representor[q_a][i] = NumericVector<Number>::build(this->comm());
              get_rb_evaluation().Aq_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }

          libmesh_assert(get_rb_evaluation().Aq_representor[q_a][i]->size()       == this->n_dofs()       &&
                         get_rb_evaluation().Aq_representor[q_a][i]->local_size() == this->n_local_dofs() );

          rhs->zero();
          if (!store_untransformed_basis)
            {
              get_Aq(q_a)->vector_mult(*rhs, get_rb_evaluation().get_basis_function(i));
            }
          else
            {
              get_Aq(q_a)->vector_mult(*rhs, *_untransformed_basis_functions[i]);
            }
          rhs->scale(-1.);

          if (!is_quiet())
            {
              libMesh::out << "Starting solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBConstruction::update_residual_terms() at "
                           << Utility::get_timestamp() << std::endl;
            }

          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);

          if (assert_convergence)
            check_convergence(*inner_product_solver);

          if (!is_quiet())
            {
              libMesh::out << "Finished solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBConstruction::update_residual_terms() at "
                           << Utility::get_timestamp() << std::endl;
              libMesh::out << this->n_linear_iterations() << " iterations, final residual "
                           << this->final_linear_residual() << std::endl;
            }

          // Store the representor
          *get_rb_evaluation().Aq_representor[q_a][i] = *solution;
        }
    }

  // Now compute and store the inner products (if requested)
  if (compute_inner_products)
    {

      for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector,*Fq_representor[q_f]);

          for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
            {
              for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  get_rb_evaluation().Fq_Aq_representor_innerprods[q_f][q_a][i] =
                    inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a][i]);
                }
            }
        }

      unsigned int q=0;
      for (unsigned int q_a1=0; q_a1<get_rb_theta_expansion().get_n_A_terms(); q_a1++)
        {
          for (unsigned int q_a2=q_a1; q_a2<get_rb_theta_expansion().get_n_A_terms(); q_a2++)
            {
              for (unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  for (unsigned int j=0; j<RB_size; j++)
                    {
                      get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][j]);
                      get_rb_evaluation().Aq_Aq_representor_innerprods[q][i][j] =
                        inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a1][i]);

                      if (i != j)
                        {
                          get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][i]);
                          get_rb_evaluation().Aq_Aq_representor_innerprods[q][j][i] =
                            inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a1][j]);
                        }
                    }
                }
              q++;
            }
        }
    } // end if (compute_inner_products)
}

SparseMatrix<Number> & RBConstruction::get_matrix_for_output_dual_solves()
{
  return *inner_product_matrix;
}

void RBConstruction::compute_output_dual_innerprods()
{
  // Skip calculations if we've already computed the output dual norms
  if (!output_dual_innerprods_computed)
    {
      // Short circuit if we don't have any outputs
      if (get_rb_theta_expansion().get_n_outputs() == 0)
        {
          output_dual_innerprods_computed = true;
          return;
        }

      // Only log if we get to here
      LOG_SCOPE("compute_output_dual_innerprods()", "RBConstruction");

      libMesh::out << "Compute output dual inner products" << std::endl;

      // Find out the largest value of Q_l
      unsigned int max_Q_l = 0;
      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        max_Q_l = (get_rb_theta_expansion().get_n_output_terms(n) > max_Q_l) ? get_rb_theta_expansion().get_n_output_terms(n) : max_Q_l;

      std::vector<std::unique_ptr<NumericVector<Number>>> L_q_representor(max_Q_l);
      for (unsigned int q=0; q<max_Q_l; q++)
        {
          L_q_representor[q] = NumericVector<Number>::build(this->comm());
          L_q_representor[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
        }

      for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        {
          for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
            {
              rhs->zero();
              rhs->add(1., *get_output_vector(n,q_l));

              if (!is_quiet())
                libMesh::out << "Starting solve n=" << n << ", q_l=" << q_l
                             << " in RBConstruction::compute_output_dual_innerprods() at "
                             << Utility::get_timestamp() << std::endl;

              // Use the main linear solver here instead of the inner_product solver, since
              // get_matrix_for_output_dual_solves() may not return the inner product matrix.
              solve_for_matrix_and_rhs(*get_linear_solver(), get_matrix_for_output_dual_solves(), *rhs);

              // We possibly perform multiple solves here with the same matrix, hence
              // set reuse_preconditioner(true) (and set it back to false again below
              // at the end of this function).
              linear_solver->reuse_preconditioner(true);

              if (assert_convergence)
                check_convergence(*get_linear_solver());

              if (!is_quiet())
                {
                  libMesh::out << "Finished solve n=" << n << ", q_l=" << q_l
                               << " in RBConstruction::compute_output_dual_innerprods() at "
                               << Utility::get_timestamp() << std::endl;

                  libMesh::out << this->n_linear_iterations()
                               << " iterations, final residual "
                               << this->final_linear_residual() << std::endl;
                }

              *L_q_representor[q_l] = *solution;
            }

          unsigned int q=0;
          for (unsigned int q_l1=0; q_l1<get_rb_theta_expansion().get_n_output_terms(n); q_l1++)
            {
              get_matrix_for_output_dual_solves().vector_mult(*inner_product_storage_vector, *L_q_representor[q_l1]);

              for (unsigned int q_l2=q_l1; q_l2<get_rb_theta_expansion().get_n_output_terms(n); q_l2++)
                {
                  output_dual_innerprods[n][q] = L_q_representor[q_l2]->dot(*inner_product_storage_vector);
                  libMesh::out << "output_dual_innerprods[" << n << "][" << q << "] = " << output_dual_innerprods[n][q] << std::endl;

                  q++;
                }
            }
        }

      // We may not need to use linear_solver again (e.g. this would happen if we use
      // extra_linear_solver for the truth_solves). As a result, let's clear linear_solver
      // to release any memory it may be taking up. If we do need it again, it will
      // be initialized when necessary.
      linear_solver->clear();
      linear_solver->reuse_preconditioner(false);

      output_dual_innerprods_computed = true;
    }

  get_rb_evaluation().output_dual_innerprods = output_dual_innerprods;
}

void RBConstruction::compute_Fq_representor_innerprods(bool compute_inner_products)
{

  // Skip calculations if we've already computed the Fq_representors
  if (!Fq_representor_innerprods_computed)
    {
      // Only log if we get to here
      LOG_SCOPE("compute_Fq_representor_innerprods()", "RBConstruction");

      for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          if (!Fq_representor[q_f])
            {
              Fq_representor[q_f] = NumericVector<Number>::build(this->comm());
              Fq_representor[q_f]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }

          libmesh_assert(Fq_representor[q_f]->size()       == this->n_dofs()       &&
                         Fq_representor[q_f]->local_size() == this->n_local_dofs() );

          rhs->zero();
          rhs->add(1., *get_Fq(q_f));

          if (!is_quiet())
            libMesh::out << "Starting solve q_f=" << q_f
                         << " in RBConstruction::update_residual_terms() at "
                         << Utility::get_timestamp() << std::endl;

          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);

          if (assert_convergence)
            check_convergence(*inner_product_solver);

          if (!is_quiet())
            {
              libMesh::out << "Finished solve q_f=" << q_f
                           << " in RBConstruction::update_residual_terms() at "
                           << Utility::get_timestamp() << std::endl;

              libMesh::out << this->n_linear_iterations()
                           << " iterations, final residual "
                           << this->final_linear_residual() << std::endl;
            }

          *Fq_representor[q_f] = *solution;
        }

      if (compute_inner_products)
        {
          unsigned int q=0;

          for (unsigned int q_f1=0; q_f1<get_rb_theta_expansion().get_n_F_terms(); q_f1++)
            {
              get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *Fq_representor[q_f1]);

              for (unsigned int q_f2=q_f1; q_f2<get_rb_theta_expansion().get_n_F_terms(); q_f2++)
                {
                  Fq_representor_innerprods[q] = inner_product_storage_vector->dot(*Fq_representor[q_f2]);

                  q++;
                }
            }
        } // end if (compute_inner_products)

      Fq_representor_innerprods_computed = true;
    }

  get_rb_evaluation().Fq_representor_innerprods = Fq_representor_innerprods;
}

void RBConstruction::load_rb_solution()
{
  LOG_SCOPE("load_rb_solution()", "RBConstruction");

  solution->zero();

  libmesh_error_msg_if(get_rb_evaluation().RB_solution.size() > get_rb_evaluation().get_n_basis_functions(),
                       "ERROR: System contains " << get_rb_evaluation().get_n_basis_functions() << " basis functions."
                       << " RB_solution vector constains " << get_rb_evaluation().RB_solution.size() << " entries."
                       << " RB_solution in RBConstruction::load_rb_solution is too long!");

  for (auto i : make_range(get_rb_evaluation().RB_solution.size()))
    solution->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));

  update();
}

// The slow (but simple, non-error prone) way to compute the residual dual norm
// Useful for error checking
Real RBConstruction::compute_residual_dual_norm_slow(const unsigned int N)
{
  LOG_SCOPE("compute_residual_dual_norm_slow()", "RBConstruction");

  // Put the residual in rhs in order to compute the norm of the Riesz representor
  // Note that this only works in serial since otherwise each processor will
  // have a different parameter value during the Greedy training.

  std::unique_ptr<NumericVector<Number>> RB_sol = NumericVector<Number>::build(comm());
  RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  std::unique_ptr<NumericVector<Number>> temp = NumericVector<Number>::build(comm());
  temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  if (store_untransformed_basis)
    {
      libmesh_assert_equal_to(_untransformed_basis_functions.size(), get_rb_evaluation().get_n_basis_functions());
    }

  for (unsigned int i=0; i<N; i++)
  {
    if (!store_untransformed_basis)
      {
        RB_sol->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
      }
    else
      {
        RB_sol->add(get_rb_evaluation().RB_solution(i), *_untransformed_basis_functions[i]);
      }
  }

  this->truth_assembly();
  matrix->vector_mult(*temp, *RB_sol);
  rhs->add(-1., *temp);

  // Then solve to get the Reisz representor
  matrix->zero();
  matrix->add(1., *inner_product_matrix);

  solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);
  get_non_dirichlet_inner_product_matrix_if_avail()->vector_mult(*inner_product_storage_vector, *solution);
  Number slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);

  return std::sqrt( libmesh_real(slow_residual_norm_sq) );
}

SparseMatrix<Number> * RBConstruction::get_inner_product_matrix()
{
  return inner_product_matrix.get();
}

SparseMatrix<Number> * RBConstruction::get_non_dirichlet_inner_product_matrix()
{
  libmesh_error_msg_if(!store_non_dirichlet_operators,
                       "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_inner_product_matrix.");

  return non_dirichlet_inner_product_matrix.get();
}

SparseMatrix<Number> * RBConstruction::get_non_dirichlet_inner_product_matrix_if_avail()
{
  if (store_non_dirichlet_operators)
    {
      return get_non_dirichlet_inner_product_matrix();
    }

  return get_inner_product_matrix();
}

SparseMatrix<Number> * RBConstruction::get_Aq(unsigned int q)
{
  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_A_terms(),
                       "Error: We must have q < Q_a in get_Aq.");

  return Aq_vector[q].get();
}

SparseMatrix<Number> * RBConstruction::get_non_dirichlet_Aq(unsigned int q)
{
  libmesh_error_msg_if(!store_non_dirichlet_operators,
                       "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Aq.");

  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_A_terms(),
                       "Error: We must have q < Q_a in get_Aq.");

  return non_dirichlet_Aq_vector[q].get();
}

SparseMatrix<Number> * RBConstruction::get_non_dirichlet_Aq_if_avail(unsigned int q)
{
  if (store_non_dirichlet_operators)
    {
      return get_non_dirichlet_Aq(q);
    }

  return get_Aq(q);
}

NumericVector<Number> * RBConstruction::get_Fq(unsigned int q)
{
  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_F_terms(),
                       "Error: We must have q < Q_f in get_Fq.");

  return Fq_vector[q].get();
}

NumericVector<Number> * RBConstruction::get_non_dirichlet_Fq(unsigned int q)
{
  libmesh_error_msg_if(!store_non_dirichlet_operators,
                       "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Fq.");

  libmesh_error_msg_if(q >= get_rb_theta_expansion().get_n_F_terms(),
                       "Error: We must have q < Q_f in get_Fq.");

  return non_dirichlet_Fq_vector[q].get();
}

NumericVector<Number> * RBConstruction::get_non_dirichlet_Fq_if_avail(unsigned int q)
{
  if (store_non_dirichlet_operators)
    {
      return get_non_dirichlet_Fq(q);
    }

  return get_Fq(q);
}

NumericVector<Number> * RBConstruction::get_output_vector(unsigned int n, unsigned int q_l)
{
  libmesh_error_msg_if((n >= get_rb_theta_expansion().get_n_outputs()) ||
                       (q_l >= get_rb_theta_expansion().get_n_output_terms(n)),
                       "Error: We must have n < n_outputs and "
                       "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_output_vector.");

  return outputs_vector[n][q_l].get();
}

NumericVector<Number> * RBConstruction::get_non_dirichlet_output_vector(unsigned int n, unsigned int q_l)
{
  libmesh_error_msg_if((n >= get_rb_theta_expansion().get_n_outputs()) ||
                       (q_l >= get_rb_theta_expansion().get_n_output_terms(n)),
                       "Error: We must have n < n_outputs and "
                       "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_non_dirichlet_output_vector.");

  return non_dirichlet_outputs_vector[n][q_l].get();
}

void RBConstruction::get_all_matrices(std::map<std::string, SparseMatrix<Number> *> & all_matrices)
{
  all_matrices.clear();

  all_matrices["inner_product"] = get_inner_product_matrix();

  if (store_non_dirichlet_operators)
    {
      all_matrices["non_dirichlet_inner_product"] = get_non_dirichlet_inner_product_matrix();
    }

  for (unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      std::stringstream matrix_name;
      matrix_name << "A" << q_a;
      all_matrices[matrix_name.str()] = get_Aq(q_a);

      if (store_non_dirichlet_operators)
        {
          matrix_name << "_non_dirichlet";
          all_matrices[matrix_name.str()] = get_non_dirichlet_Aq(q_a);
        }
    }
}

void RBConstruction::get_all_vectors(std::map<std::string, NumericVector<Number> *> & all_vectors)
{
  all_vectors.clear();

  get_output_vectors(all_vectors);

  for (unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      std::stringstream F_vector_name;
      F_vector_name << "F" << q_f;
      all_vectors[F_vector_name.str()] = get_Fq(q_f);

      if (store_non_dirichlet_operators)
        {
          F_vector_name << "_non_dirichlet";
          all_vectors[F_vector_name.str()] = get_non_dirichlet_Fq(q_f);
        }
    }
}

void RBConstruction::get_output_vectors(std::map<std::string, NumericVector<Number> *> & output_vectors)
{
  output_vectors.clear();

  for (unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for (unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        std::stringstream output_name;
        output_name << "output_" << n << "_"<< q_l;
        output_vectors[output_name.str()] = get_output_vector(n,q_l);

        if (store_non_dirichlet_operators)
          {
            output_name << "_non_dirichlet";
            output_vectors[output_name.str()] = get_non_dirichlet_output_vector(n,q_l);
          }
      }
}

#ifdef LIBMESH_ENABLE_DIRICHLET

std::unique_ptr<DirichletBoundary> RBConstruction::build_zero_dirichlet_boundary_object()
{
  ZeroFunction<> zf;

  std::set<boundary_id_type> dirichlet_ids;
  std::vector<unsigned int> variables;

  // The DirichletBoundary constructor clones zf, so it's OK that zf is only in local scope
  return libmesh_make_unique<DirichletBoundary>(dirichlet_ids, variables, &zf);
}

#endif

void RBConstruction::write_riesz_representors_to_files(const std::string & riesz_representors_dir,
                                                       const bool write_binary_residual_representors)
{
  LOG_SCOPE("write_riesz_representors_to_files()", "RBConstruction");

  // Write out Riesz representors. These are useful to have when restarting,
  // so you don't have to recompute them all over again.

  // First we write out the Fq representors, these are independent of an RBEvaluation object.
  libMesh::out << "Writing out the Fq_representors..." << std::endl;

  std::ostringstream file_name;
  const std::string riesz_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");
  struct stat stat_info;

  // Residual representors written out to their own separate directory
  if ( this->processor_id() == 0)
    if ( Utility::mkdir(riesz_representors_dir.c_str()) != 0)
      libMesh::out << "Skipping creating residual_representors directory: " << strerror(errno) << std::endl;

  for (auto i : index_range(Fq_representor))
    {
      if (Fq_representor[i])
        {
          file_name.str(""); // reset filename
          file_name << riesz_representors_dir << "/Fq_representor" << i << riesz_representor_suffix;

          // Check to see if file exists, if so, don't overwrite it, we assume it was
          // there from a previous call to this function.  Note: if stat returns zero
          // it means it successfully got the file attributes (and therefore the file
          // exists).  Because of the following factors:
          // 1.) write_serialized_data takes longer for proc 0 than it does for others,
          //     so processors can get out of sync.
          // 2.) The constructor for Xdr opens a (0 length) file on *all* processors,
          // there are typically hundreds of 0-length files created during this loop,
          // and that screws up checking for the existence of files.  One way to stay
          // in sync is to but a barrier at each iteration of the loop -- not sure how
          // bad this will affect performance, but it can't be much worse than serialized
          // I/O already is :)
          int stat_result = stat(file_name.str().c_str(), &stat_info);

          if ( (stat_result != 0) ||     // file definitely doesn't already exist
               (stat_info.st_size == 0)) // file exists, but has zero length (can happen if another proc already opened it!)
            {
              // No need to copy!
              // *solution = *(Fq_representor[i]);
              // std::swap doesn't work on pointers
              //std::swap(solution.get(), Fq_representor[i]);
              Fq_representor[i]->swap(*solution);

              Xdr fqr_data(file_name.str(),
                           write_binary_residual_representors ? ENCODE : WRITE);

              write_serialized_data(fqr_data, false);

              // Synchronize before moving on
              this->comm().barrier();

              // Swap back.
              Fq_representor[i]->swap(*solution);

              // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
              // for the system call, be sure to do it only on one processor, etc.
            }
        }
    }


  // Next, write out the Aq representors associated with rb_eval.
  libMesh::out << "Writing out the Aq_representors..." << std::endl;

  const unsigned int jstop  = get_rb_evaluation().get_n_basis_functions();
  const unsigned int jstart = jstop-get_delta_N();

  RBEvaluation & rbe = get_rb_evaluation();
  for (auto i : index_range(rbe.Aq_representor))
    for (unsigned int j=jstart; j<jstop; ++j)
      {
        libMesh::out << "Writing out Aq_representor[" << i << "][" << j << "]..." << std::endl;
        libmesh_assert(get_rb_evaluation().Aq_representor[i][j]);

        file_name.str(""); // reset filename
        file_name << riesz_representors_dir
                  << "/Aq_representor" << i << "_" << j << riesz_representor_suffix;

        {
          // No need to copy! Use swap instead.
          // *solution = *(Aq_representor[i][j]);
          rbe.Aq_representor[i][j]->swap(*solution);

          Xdr aqr_data(file_name.str(),
                       write_binary_residual_representors ? ENCODE : WRITE);

          write_serialized_data(aqr_data, false);

          // Synchronize before moving on
          this->comm().barrier();

          // Swap back.
          rbe.Aq_representor[i][j]->swap(*solution);

          // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
          // for the system call, be sure to do it only on one processor, etc.
        }
      }
}



void RBConstruction::read_riesz_representors_from_files(const std::string & riesz_representors_dir,
                                                        const bool read_binary_residual_representors)
{
  LOG_SCOPE("read_riesz_representors_from_files()", "RBConstruction");

  libMesh::out << "Reading in the Fq_representors..." << std::endl;

  const std::string riesz_representor_suffix = (read_binary_residual_representors ? ".xdr" : ".dat");
  std::ostringstream file_name;
  struct stat stat_info;

  // Read in the Fq_representors.  There should be Q_f of these.  FIXME:
  // should we be worried about leaks here?
  for (const auto & rep : Fq_representor)
    libmesh_error_msg_if(rep, "Error, must delete existing Fq_representor before reading in from file.");

  for (auto i : index_range(Fq_representor))
    {
      file_name.str(""); // reset filename
      file_name << riesz_representors_dir
                << "/Fq_representor" << i << riesz_representor_suffix;

      // On processor zero check to be sure the file exists
      if (this->processor_id() == 0)
        {
          int stat_result = stat(file_name.str().c_str(), &stat_info);

          libmesh_error_msg_if(stat_result != 0, "File does not exist: " << file_name.str());
        }

      Xdr fqr_data(file_name.str(),
                   read_binary_residual_representors ? DECODE : READ);

      read_serialized_data(fqr_data, false);

      Fq_representor[i] = NumericVector<Number>::build(this->comm());
      Fq_representor[i]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

      // No need to copy, just swap
      // *Fq_representor[i] = *solution;
      Fq_representor[i]->swap(*solution);
    }

  // Alert the update_residual_terms() function that we don't need to recompute
  // the Fq_representors as we have already read them in from file!
  Fq_representor_innerprods_computed = true;


  libMesh::out << "Reading in the Aq_representors..." << std::endl;

  // Read in the Aq representors.  The class makes room for [Q_a][Nmax] of these.  We are going to
  // read in [Q_a][get_rb_evaluation().get_n_basis_functions()].  FIXME:
  // should we be worried about leaks in the locations where we're about to fill entries?
  RBEvaluation & rbe = get_rb_evaluation();
  for (const auto & row : rbe.Aq_representor)
    for (const auto & rep : row)
      libmesh_error_msg_if(rep, "Error, must delete existing Aq_representor before reading in from file.");

  // Now ready to read them in from file!
  for (auto i : index_range(rbe.Aq_representor))
    for (unsigned int j=0; j<rbe.get_n_basis_functions(); ++j)
      {
        file_name.str(""); // reset filename
        file_name << riesz_representors_dir
                  << "/Aq_representor" << i << "_" << j << riesz_representor_suffix;

        // On processor zero check to be sure the file exists
        if (this->processor_id() == 0)
          {
            int stat_result = stat(file_name.str().c_str(), &stat_info);

            libmesh_error_msg_if(stat_result != 0, "File does not exist: " << file_name.str());
          }

        Xdr aqr_data(file_name.str(), read_binary_residual_representors ? DECODE : READ);

        read_serialized_data(aqr_data, false);

        rbe.Aq_representor[i][j] = NumericVector<Number>::build(this->comm());
        rbe.Aq_representor[i][j]->init (n_dofs(), n_local_dofs(),
                                            false, PARALLEL);

        // No need to copy, just swap
        rbe.Aq_representor[i][j]->swap(*solution);
      }
}

void RBConstruction::check_convergence(LinearSolver<Number> & input_solver)
{
  libMesh::LinearConvergenceReason conv_flag;

  conv_flag = input_solver.get_converged_reason();

  libmesh_error_msg_if(conv_flag < 0, "Convergence error. Error id: " << conv_flag);
}

bool RBConstruction::get_convergence_assertion_flag() const
{
  return assert_convergence;
}

void RBConstruction::set_convergence_assertion_flag(bool flag)
{
  assert_convergence = flag;
}

bool RBConstruction::get_preevaluate_thetas_flag() const
{
  return _preevaluate_thetas_flag;
}

void RBConstruction::set_preevaluate_thetas_flag(bool flag)
{
  _preevaluate_thetas_flag = flag;
}

unsigned int RBConstruction::get_current_training_parameter_index() const
{
  return _current_training_parameter_index;
}

void RBConstruction::set_current_training_parameter_index(unsigned int index)
{
  _current_training_parameter_index = index;
}

const std::vector<Number> &
RBConstruction::get_evaluated_thetas(unsigned int training_parameter_index) const
{
  const numeric_index_type first_index = get_first_local_training_index();
  libmesh_assert(training_parameter_index >= first_index);

  const numeric_index_type local_index = training_parameter_index - first_index;
  libmesh_assert(local_index < _evaluated_thetas.size());

  return _evaluated_thetas[local_index];
}

void RBConstruction::preevaluate_thetas()
{
  LOG_SCOPE("preevaluate_thetas()", "RBConstruction");

  _evaluated_thetas.resize(get_local_n_training_samples());

  if ( get_local_n_training_samples() == 0 )
    return;

  auto & rb_theta_expansion = get_rb_evaluation().get_rb_theta_expansion();
  const unsigned int n_A_terms = rb_theta_expansion.get_n_A_terms();
  const unsigned int n_F_terms = rb_theta_expansion.get_n_F_terms();

  // Collect all training parameters
  std::vector<RBParameters> mus(get_local_n_training_samples());
  const numeric_index_type first_index = get_first_local_training_index();
  for (unsigned int i=0; i<get_local_n_training_samples(); i++)
    {
      // Load training parameter i, this is only loaded
      // locally since the RB solves are local.
      set_params_from_training_set( first_index+i );
      mus[i] = get_parameters();
      _evaluated_thetas[i].resize(n_A_terms + n_F_terms);
    }

  // Evaluate thetas for all training parameters simultaneously
  for (unsigned int q_a=0; q_a<n_A_terms; q_a++)
    {
      const auto A_vals = rb_theta_expansion.eval_A_theta(q_a, mus);
      for (auto i : make_range(get_local_n_training_samples()))
        _evaluated_thetas[i][q_a] = A_vals[i];
    }

  for (unsigned int q_f=0; q_f<n_F_terms; q_f++)
    {
      const auto F_vals = rb_theta_expansion.eval_F_theta(q_f, mus);
      for (auto i : make_range(get_local_n_training_samples()))
        _evaluated_thetas[i][n_A_terms + q_f] = F_vals[i];
    }
}

} // namespace libMesh
