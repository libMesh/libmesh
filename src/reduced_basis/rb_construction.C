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

// C++ includes
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <limits>

namespace libMesh
{

RBConstruction::RBConstruction (EquationSystems & es,
                                const std::string & name_in,
                                const unsigned int number_in)
  : Parent(es, name_in, number_in),
    inner_product_solver(LinearSolver<Number>::build(es.comm())),
    extra_linear_solver(libmesh_nullptr),
    inner_product_matrix(SparseMatrix<Number>::build(es.comm())),
    exit_on_repeated_greedy_parameters(true),
    impose_internal_fluxes(false),
    compute_RB_inner_product(false),
    store_non_dirichlet_operators(false),
    use_empty_rb_solve_in_greedy(true),
    Fq_representor_innerprods_computed(false),
    Nmax(0),
    delta_N(1),
    quiet_mode(true),
    output_dual_innerprods_computed(false),
    assert_convergence(true),
    rb_eval(libmesh_nullptr),
    inner_product_assembly(libmesh_nullptr),
    rel_training_tolerance(1.e-4),
    abs_training_tolerance(1.e-12),
    normalize_rb_bound_in_greedy(false)
{
  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;
}


RBConstruction::~RBConstruction ()
{
  this->clear();
}


void RBConstruction::clear()
{
  LOG_SCOPE("clear()", "RBConstruction");

  Parent::clear();

  for (std::size_t q=0; q<Aq_vector.size(); q++)
    {
      delete Aq_vector[q];
      Aq_vector[q] = libmesh_nullptr;
    }

  for (std::size_t q=0; q<Fq_vector.size(); q++)
    {
      delete Fq_vector[q];
      Fq_vector[q] = libmesh_nullptr;
    }

  for (std::size_t i=0; i<outputs_vector.size(); i++)
    for (std::size_t q_l=0; q_l<outputs_vector[i].size(); q_l++)
      {
        delete outputs_vector[i][q_l];
        outputs_vector[i][q_l] = libmesh_nullptr;
      }

  if(store_non_dirichlet_operators)
    {
      for (std::size_t q=0; q<non_dirichlet_Aq_vector.size(); q++)
        {
          delete non_dirichlet_Aq_vector[q];
          non_dirichlet_Aq_vector[q] = libmesh_nullptr;
        }

      for (std::size_t q=0; q<non_dirichlet_Fq_vector.size(); q++)
        {
          delete non_dirichlet_Fq_vector[q];
          non_dirichlet_Fq_vector[q] = libmesh_nullptr;
        }

      for (std::size_t i=0; i<non_dirichlet_outputs_vector.size(); i++)
        for (std::size_t q_l=0; q_l<non_dirichlet_outputs_vector[i].size(); q_l++)
          {
            delete non_dirichlet_outputs_vector[i][q_l];
            non_dirichlet_outputs_vector[i][q_l] = libmesh_nullptr;
          }
    }

  // Also delete the Fq representors
  for (std::size_t q_f=0; q_f<Fq_representor.size(); q_f++)
    {
      delete Fq_representor[q_f];
      Fq_representor[q_f] = libmesh_nullptr;
    }
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
  const Real tol  =
    es.parameters.get<Real>("linear solver tolerance");

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

  // Update the system after the solve
  this->update();
}

void RBConstruction::set_rb_evaluation(RBEvaluation & rb_eval_in)
{
  rb_eval = &rb_eval_in;
}

RBEvaluation & RBConstruction::get_rb_evaluation()
{
  if(!rb_eval)
    libmesh_error_msg("Error: RBEvaluation object hasn't been initialized yet");

  return *rb_eval;
}

bool RBConstruction::is_rb_eval_initialized() const
{
  return (rb_eval != libmesh_nullptr);
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
  bool normalize_rb_bound_in_greedy = false;
  const bool normalize_rb_bound_in_greedy_in = infile("normalize_rb_bound_in_greedy",
                                                      normalize_rb_bound_in_greedy);

  // Read in the parameters from the input file too
  unsigned int n_continuous_parameters = infile.vector_variable_size("parameter_names");
  RBParameters mu_min_in;
  RBParameters mu_max_in;
  for(unsigned int i=0; i<n_continuous_parameters; i++)
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

  std::map< std::string, std::vector<Real> > discrete_parameter_values_in;

  unsigned int n_discrete_parameters = infile.vector_variable_size("discrete_parameter_names");
  for(unsigned int i=0; i<n_discrete_parameters; i++)
    {
      std::string param_name = infile("discrete_parameter_names", "NONE", i);

      unsigned int n_vals_for_param = infile.vector_variable_size(param_name);
      std::vector<Real> vals_for_param(n_vals_for_param);
      for (std::size_t j=0; j<vals_for_param.size(); j++)
        vals_for_param[j] = infile(param_name, 0., j);

      discrete_parameter_values_in[param_name] = vals_for_param;
    }

  std::map<std::string,bool> log_scaling_in;
  RBParameters::const_iterator it     = mu_min_in.begin();
  RBParameters::const_iterator it_end = mu_min_in.end();
  for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;

      // For now, just set all entries to false.
      // TODO: Implement a decent way to specify log-scaling true/false
      // in the input text file
      log_scaling_in[param_name] = false;
    }

  // Set the parameters that have been read in
  set_rb_construction_parameters(n_training_samples,
                                 deterministic_training,
                                 training_parameters_random_seed_in,
                                 quiet_mode_in,
                                 Nmax_in,
                                 rel_training_tolerance_in,
                                 abs_training_tolerance_in,
                                 normalize_rb_bound_in_greedy_in,
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
                                                    RBParameters mu_min_in,
                                                    RBParameters mu_max_in,
                                                    std::map< std::string, std::vector<Real> > discrete_parameter_values_in,
                                                    std::map<std::string,bool> log_scaling_in)
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

  // Initialize the parameter ranges and the parameters themselves
  initialize_parameters(mu_min_in, mu_max_in, discrete_parameter_values_in);

  initialize_training_parameters(this->get_parameters_min(),
                                 this->get_parameters_max(),
                                 n_training_samples_in,
                                 log_scaling_in,
                                 deterministic_training_in);   // use deterministic parameters
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
  if( is_rb_eval_initialized() )
    {
      libMesh::out << "Aq operators attached: " << get_rb_theta_expansion().get_n_A_terms() << std::endl;
      libMesh::out << "Fq functions attached: " << get_rb_theta_expansion().get_n_F_terms() << std::endl;
      libMesh::out << "n_outputs: " << get_rb_theta_expansion().get_n_outputs() << std::endl;
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        libMesh::out << "output " << n << ", Q_l = " << get_rb_theta_expansion().get_n_output_terms(n) << std::endl;
    }
  else
    {
      libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
    }
  libMesh::out << "Number of parameters: " << get_n_params() << std::endl;
  RBParameters::const_iterator it     = get_parameters().begin();
  RBParameters::const_iterator it_end = get_parameters().end();
  for( ; it != it_end; ++it)
    {
      std::string param_name = it->first;
      if(!is_discrete_parameter(param_name))
        {
          libMesh::out <<   "Parameter " << param_name
                       << ": Min = " << get_parameter_min(param_name)
                       << ", Max = " << get_parameter_max(param_name) << std::endl;
        }
    }
  print_discrete_parameter_values();
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << "quiet mode? " << is_quiet() << std::endl;
  libMesh::out << std::endl;
}

void RBConstruction::print_basis_function_orthogonality()
{
  UniquePtr< NumericVector<Number> > temp = solution->clone();

  for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
    {
      for(unsigned int j=0; j<get_rb_evaluation().get_n_basis_functions(); j++)
        {
          inner_product_matrix->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));
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
  if(!rb_assembly_expansion)
    libmesh_error_msg("Error: RBAssemblyExpansion object hasn't been initialized yet");

  return *rb_assembly_expansion;
}

void RBConstruction::set_inner_product_assembly(ElemAssembly & inner_product_assembly_in)
{
  inner_product_assembly = &inner_product_assembly_in;
}

ElemAssembly & RBConstruction::get_inner_product_assembly()
{
  if(!inner_product_assembly)
    libmesh_error_msg("Error: inner_product_assembly hasn't been initialized yet");

  return *inner_product_assembly;
}

void RBConstruction::zero_constrained_dofs_on_vector(NumericVector<Number> & vector)
{
  const DofMap & dof_map = get_dof_map();

  for(dof_id_type i=dof_map.first_dof(); i<dof_map.end_dof(); i++)
    {
      if(get_dof_map().is_constrained_dof(i))
        {
          vector.set(i, 0.);
        }
    }
  vector.close();
}

void RBConstruction::initialize_rb_construction(bool skip_matrix_assembly,
                                                bool skip_vector_assembly)
{
  if(!skip_matrix_assembly && !skip_vector_assembly)
    {
      // Check that the theta and assembly objects are consistently sized
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_A_terms(), get_rb_assembly_expansion().get_n_A_terms());
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_F_terms(), get_rb_assembly_expansion().get_n_F_terms());
      libmesh_assert_equal_to (get_rb_theta_expansion().get_n_outputs(), get_rb_assembly_expansion().get_n_outputs());
      for(unsigned int i=0; i<get_rb_theta_expansion().get_n_outputs(); i++)
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

  // Resize the Fq_representors and initialize each to NULL
  // These are basis independent and hence stored here, whereas
  // the Aq_representors are stored in RBEvaluation
  Fq_representor.resize(get_rb_theta_expansion().get_n_F_terms());

  // Initialize vectors for the inner products of the Fq representors
  // These are basis independent and therefore stored here.
  unsigned int Q_f_hat = get_rb_theta_expansion().get_n_F_terms()*(get_rb_theta_expansion().get_n_F_terms()+1)/2;
  Fq_representor_innerprods.resize(Q_f_hat);

  // Resize the output vectors
  outputs_vector.resize(get_rb_theta_expansion().get_n_outputs());
  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    outputs_vector[n].resize( get_rb_theta_expansion().get_n_output_terms(n) );

  // Resize the output dual norm vectors
  output_dual_innerprods.resize(get_rb_theta_expansion().get_n_outputs());
  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      unsigned int Q_l_hat = get_rb_theta_expansion().get_n_output_terms(n)*(get_rb_theta_expansion().get_n_output_terms(n)+1)/2;
      output_dual_innerprods[n].resize(Q_l_hat);
    }

  {
    DofMap & dof_map = this->get_dof_map();

    dof_map.attach_matrix(*inner_product_matrix);
    inner_product_matrix->init();
    inner_product_matrix->zero();

    for(unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
      {
        // Initialize the memory for the matrices
        Aq_vector[q] = SparseMatrix<Number>::build(this->comm()).release();
        dof_map.attach_matrix(*Aq_vector[q]);
        Aq_vector[q]->init();
        Aq_vector[q]->zero();
      }

    // We also need to initialize a second set of non-Dirichlet operators
    if(store_non_dirichlet_operators)
      {
        non_dirichlet_Aq_vector.resize(get_rb_theta_expansion().get_n_A_terms());
        for(unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
          {
            // Initialize the memory for the matrices
            non_dirichlet_Aq_vector[q] = SparseMatrix<Number>::build(this->comm()).release();
            dof_map.attach_matrix(*non_dirichlet_Aq_vector[q]);
            non_dirichlet_Aq_vector[q]->init();
            non_dirichlet_Aq_vector[q]->zero();
          }
      }
  }

  // Initialize the vectors
  for(unsigned int q=0; q<get_rb_theta_expansion().get_n_F_terms(); q++)
    {
      // Initialize the memory for the vectors
      Fq_vector[q] = NumericVector<Number>::build(this->comm()).release();
      Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    }

  // We also need to initialize a second set of non-Dirichlet operators
  if(store_non_dirichlet_operators)
    {
      non_dirichlet_Fq_vector.resize(get_rb_theta_expansion().get_n_F_terms());
      for(unsigned int q=0; q<get_rb_theta_expansion().get_n_F_terms(); q++)
        {
          // Initialize the memory for the vectors
          non_dirichlet_Fq_vector[q] = NumericVector<Number>::build(this->comm()).release();
          non_dirichlet_Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
        }
    }

  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        // Initialize the memory for the truth output vectors
        outputs_vector[n][q_l] = (NumericVector<Number>::build(this->comm()).release());
        outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      }

  if(store_non_dirichlet_operators)
    {
      non_dirichlet_outputs_vector.resize(get_rb_theta_expansion().get_n_outputs());
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        {
          non_dirichlet_outputs_vector[n].resize( get_rb_theta_expansion().get_n_output_terms(n) );
          for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
            {
              // Initialize the memory for the truth output vectors
              non_dirichlet_outputs_vector[n][q_l] = (NumericVector<Number>::build(this->comm()).release());
              non_dirichlet_outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }
        }
    }

  // Resize truth_outputs vector
  truth_outputs.resize(this->get_rb_theta_expansion().get_n_outputs());
}

UniquePtr<DGFEMContext> RBConstruction::build_context ()
{
  return UniquePtr<DGFEMContext>(new DGFEMContext(*this));
}

void RBConstruction::add_scaled_matrix_and_vector(Number scalar,
                                                  ElemAssembly * elem_assembly,
                                                  SparseMatrix<Number> * input_matrix,
                                                  NumericVector<Number> * input_vector,
                                                  bool symmetrize,
                                                  bool apply_dof_constraints)
{
  LOG_SCOPE("add_scaled_matrix_and_vector()", "RBConstruction");

  bool assemble_matrix = (input_matrix != libmesh_nullptr);
  bool assemble_vector = (input_vector != libmesh_nullptr);

  if(!assemble_matrix && !assemble_vector)
    return;

  const MeshBase & mesh = this->get_mesh();

  // First add any node-based terms (e.g. point loads)
  // We only enter this loop if we have at least one
  // nodeset, since we use nodesets to indicate
  // where to impose the node-based terms.
  if (mesh.get_boundary_info().n_nodeset_conds() > 0)
    {
      std::vector<numeric_index_type> node_id_list;
      std::vector<boundary_id_type> bc_id_list;

      // Get the list of nodes with boundary IDs
      mesh.get_boundary_info().build_node_list(node_id_list, bc_id_list);

      for (std::size_t i=0; i<node_id_list.size(); i++)
        {
          const Node & node = mesh.node_ref(node_id_list[i]);

          // If node is on this processor, then all dofs on node are too
          // so we can do the add below safely
          if (node.processor_id() == this->comm().rank())
            {
              // Get the values to add to the rhs vector
              std::map<numeric_index_type, Number> rhs_values;
              elem_assembly->get_nodal_rhs_values(rhs_values, *this, node);

              std::map<numeric_index_type, Number>::const_iterator it =
                rhs_values.begin();
              const std::map<numeric_index_type, Number>::const_iterator it_end =
                rhs_values.end();
              for ( ; it != it_end; ++it)
                {
                  numeric_index_type dof_index = it->first;
                  Number value = it->second;

                  input_vector->add( dof_index, value);
                }
            }
        }
    }

  UniquePtr<DGFEMContext> c = this->build_context();
  DGFEMContext & context  = cast_ref<DGFEMContext &>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Subdivision elements need special care:
      // - skip ghost elements
      // - init special quadrature rule
      const Elem * elem = *el;
      UniquePtr<QBase> qrule;
      if (elem->type() == TRI3SUBDIVISION)
        {
          const Tri3Subdivision * gh_elem = static_cast<const Tri3Subdivision *> (elem);
          if (gh_elem->is_ghost())
            continue ;
          // A Gauss quadrature rule for numerical integration.
          // For subdivision shell elements, a single Gauss point per
          // element is sufficient, hence we use extraorder = 0.
          const int extraorder = 0;
          FEBase * elem_fe = libmesh_nullptr;
          context.get_element_fe( 0, elem_fe );

          qrule = elem_fe->get_fe_type().default_quadrature_rule (2, extraorder);

          // Tell the finite element object to use our quadrature rule.
          elem_fe->attach_quadrature_rule (qrule.get());
        }

      context.pre_fe_reinit(*this, *el);
      context.elem_fe_reinit();
      elem_assembly->interior_assembly(context);

      for (context.side = 0;
           context.side != context.get_elem().n_sides();
           ++context.side )
        {
          // May not need to apply fluxes on non-boundary elements
          if( (context.get_elem().neighbor_ptr(context.get_side()) != libmesh_nullptr) && !impose_internal_fluxes )
            continue;

          // Impose boundary (e.g. Neumann) term
          context.side_fe_reinit();
          elem_assembly->boundary_assembly(context);

          if(context.dg_terms_are_active())
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

      // Need to symmetrize before imposing
      // periodic constraints
      if(assemble_matrix && symmetrize)
        {
          DenseMatrix<Number> Ke_transpose;
          context.get_elem_jacobian().get_transpose(Ke_transpose);
          context.get_elem_jacobian() += Ke_transpose;
          context.get_elem_jacobian() *= 0.5;
        }

      if(apply_dof_constraints)
        {
          // Apply constraints, e.g. Dirichlet and periodic constraints
          this->get_dof_map().constrain_element_matrix_and_vector
            (context.get_elem_jacobian(), context.get_elem_residual(), context.get_dof_indices() );
        }

      // Scale and add to global matrix and/or vector
      context.get_elem_jacobian() *= scalar;
      context.get_elem_residual() *= scalar;

      if(assemble_matrix)
        {

          CouplingMatrix * coupling_matrix = get_dof_map()._dof_coupling;
          if(!coupling_matrix)
            {
              // If we haven't defined a _dof_coupling matrix then just add
              // the whole matrix
              input_matrix->add_matrix (context.get_elem_jacobian(),
                                        context.get_dof_indices() );
            }
          else
            {
              // Otherwise we should only add the relevant submatrices
              for(unsigned int var1=0; var1<n_vars(); var1++)
                {
                  ConstCouplingRow ccr(var1, *coupling_matrix);
                  ConstCouplingRow::const_iterator end = ccr.end();
                  for (ConstCouplingRow::const_iterator it =
                         ccr.begin(); it != end; ++it)
                    {
                      unsigned int var2 = *it;

                      unsigned int sub_m = context.get_elem_jacobian( var1, var2 ).m();
                      unsigned int sub_n = context.get_elem_jacobian( var1, var2 ).n();
                      DenseMatrix<Number> sub_jac(sub_m, sub_n);
                      for(unsigned int row=0; row<sub_m; row++)
                        for(unsigned int col=0; col<sub_n; col++)
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

      if(assemble_vector)
        input_vector->add_vector (context.get_elem_residual(),
                                  context.get_dof_indices() );
    }

  if(assemble_matrix)
    input_matrix->close();
  if(assemble_vector)
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

    for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
      {
        matrix->add(get_rb_theta_expansion().eval_A_theta(q_a, mu), *get_Aq(q_a));
      }

    UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
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
  add_scaled_matrix_and_vector(1.,
                               inner_product_assembly,
                               input_matrix,
                               libmesh_nullptr,
                               false, /* symmetrize */
                               apply_dof_constraints);
}

void RBConstruction::assemble_Aq_matrix(unsigned int q,
                                        SparseMatrix<Number> * input_matrix,
                                        bool apply_dof_constraints)
{
  if(q >= get_rb_theta_expansion().get_n_A_terms())
    libmesh_error_msg("Error: We must have q < Q_a in assemble_Aq_matrix.");

  input_matrix->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_A_assembly(q),
                               input_matrix,
                               libmesh_nullptr,
                               false, /* symmetrize */
                               apply_dof_constraints);
}

void RBConstruction::add_scaled_Aq(Number scalar,
                                   unsigned int q_a,
                                   SparseMatrix<Number> * input_matrix,
                                   bool symmetrize)
{
  LOG_SCOPE("add_scaled_Aq()", "RBConstruction");

  if(q_a >= get_rb_theta_expansion().get_n_A_terms())
    libmesh_error_msg("Error: We must have q < Q_a in add_scaled_Aq.");

  if(!symmetrize)
    {
      input_matrix->add(scalar, *get_Aq(q_a));
      input_matrix->close();
    }
  else
    {
      add_scaled_matrix_and_vector(scalar,
                                   &rb_assembly_expansion->get_A_assembly(q_a),
                                   input_matrix,
                                   libmesh_nullptr,
                                   symmetrize);
    }
}

void RBConstruction::assemble_misc_matrices()
{
  libMesh::out << "Assembling inner product matrix" << std::endl;
  assemble_inner_product_matrix(inner_product_matrix.get());
}

void RBConstruction::assemble_all_affine_operators()
{
  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      libMesh::out << "Assembling affine operator " << (q_a+1) << " of "
                   << get_rb_theta_expansion().get_n_A_terms() << std::endl;
      assemble_Aq_matrix(q_a, get_Aq(q_a));
    }

  if(store_non_dirichlet_operators)
    {
      for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
        {
          libMesh::out << "Assembling non-Dirichlet affine operator " << (q_a+1) << " of "
                       << get_rb_theta_expansion().get_n_A_terms() << std::endl;
          assemble_Aq_matrix(q_a, get_non_dirichlet_Aq(q_a), false);
        }
    }
}

void RBConstruction::assemble_all_affine_vectors()
{
  for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      libMesh::out << "Assembling affine vector " << (q_f+1) << " of "
                   << get_rb_theta_expansion().get_n_F_terms() << std::endl;
      assemble_Fq_vector(q_f, get_Fq(q_f));
    }

  if(store_non_dirichlet_operators)
    {
      for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
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
  if(q >= get_rb_theta_expansion().get_n_F_terms())
    libmesh_error_msg("Error: We must have q < Q_f in assemble_Fq_vector.");

  input_vector->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_F_assembly(q),
                               libmesh_nullptr,
                               input_vector,
                               false,             /* symmetrize */
                               apply_dof_constraints /* apply_dof_constraints */);
}

void RBConstruction::assemble_all_output_vectors()
{
  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        libMesh::out << "Assembling output vector, (" << (n+1) << "," << (q_l+1)
                     << ") of (" << get_rb_theta_expansion().get_n_outputs()
                     << "," << get_rb_theta_expansion().get_n_output_terms(n) << ")"
                     << std::endl;
        get_output_vector(n, q_l)->zero();
        add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                     libmesh_nullptr,
                                     get_output_vector(n,q_l),
                                     false, /* symmetrize */
                                     true   /* apply_dof_constraints */);
      }

  if(store_non_dirichlet_operators)
    {
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
          {
            libMesh::out << "Assembling non-Dirichlet output vector, (" << (n+1) << "," << (q_l+1)
                         << ") of (" << get_rb_theta_expansion().get_n_outputs()
                         << "," << get_rb_theta_expansion().get_n_output_terms(n) << ")"
                         << std::endl;
            get_non_dirichlet_output_vector(n, q_l)->zero();
            add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                         libmesh_nullptr,
                                         get_non_dirichlet_output_vector(n,q_l),
                                         false, /* symmetrize */
                                         false  /* apply_dof_constraints */);
          }
    }
}

Real RBConstruction::train_reduced_basis(const bool resize_rb_eval_data)
{
  LOG_SCOPE("train_reduced_basis()", "RBConstruction");

  int count = 0;

  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  // possibly resize data structures according to Nmax
  if(resize_rb_eval_data)
    {
      get_rb_evaluation().resize_data_structures(get_Nmax());
    }

  // Clear the Greedy param list
  for (std::size_t i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
    get_rb_evaluation().greedy_param_list[i].clear();

  get_rb_evaluation().greedy_param_list.clear();

  Real training_greedy_error;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if(get_rb_evaluation().get_n_basis_functions() >= get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << std::endl;
      return 0.;
    }


  // Compute the dual norms of the outputs if we haven't already done so
  compute_output_dual_innerprods();

  // Compute the Fq Riesz representor dual norms if we haven't already done so
  compute_Fq_representor_innerprods();

  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
  Real initial_greedy_error = 0.;
  bool initial_greedy_error_initialized = false;
  while(true)
    {
      libMesh::out << std::endl << "---- Basis dimension: "
                   << get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;

      if( count > 0 || (count==0 && use_empty_rb_solve_in_greedy) )
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

      // Add orthogonal part of the snapshot to the RB space
      libMesh::out << "Enriching the RB space" << std::endl;
      enrich_RB_space();

      update_system();

      // Increment counter
      count++;
    }
  this->update_greedy_param_list();

  return training_greedy_error;
}

bool RBConstruction::greedy_termination_test(Real abs_greedy_error,
                                             Real initial_error,
                                             int)
{
  if(abs_greedy_error < this->abs_training_tolerance)
    {
      libMesh::out << "Absolute error tolerance reached." << std::endl;
      return true;
    }

  Real rel_greedy_error = abs_greedy_error/initial_error;
  if(rel_greedy_error < this->rel_training_tolerance)
    {
      libMesh::out << "Relative error tolerance reached." << std::endl;
      return true;
    }

  if(get_rb_evaluation().get_n_basis_functions() >= this->get_Nmax())
    {
      libMesh::out << "Maximum number of basis functions reached: Nmax = "
                   << get_Nmax() << std::endl;
      return true;
    }

  if(exit_on_repeated_greedy_parameters)
    {
      for (std::size_t i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
        {
          RBParameters & previous_parameters = get_rb_evaluation().greedy_param_list[i];
          if(previous_parameters == get_parameters())
            {
              libMesh::out << "Exiting greedy because the same parameters were selected twice"
                           << std::endl;
              return true;
            }
        }
    }

  return false;
}

void RBConstruction::update_greedy_param_list()
{
  get_rb_evaluation().greedy_param_list.push_back( get_parameters() );
}

const RBParameters & RBConstruction::get_greedy_parameter(unsigned int i)
{
  if( i >= get_rb_evaluation().greedy_param_list.size() )
    libmesh_error_msg("Error: Argument in RBConstruction::get_greedy_parameter is too large.");

  return get_rb_evaluation().greedy_param_list[i];
}

Real RBConstruction::truth_solve(int plot_solution)
{
  LOG_SCOPE("truth_solve()", "RBConstruction");

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



  const RBParameters & mu = get_parameters();

  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      truth_outputs[n] = 0.;
      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
        truth_outputs[n] += get_rb_theta_expansion().eval_output_theta(n, q_l, mu)*
          get_output_vector(n,q_l)->dot(*solution);
    }

  if(plot_solution > 0)
    {
#if defined(LIBMESH_USE_COMPLEX_NUMBERS)
      GMVIO(get_mesh()).write_equation_systems ("truth.gmv",
                                                this->get_equation_systems());
#else
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(get_mesh()).write_equation_systems ("truth.e",
                                                      this->get_equation_systems());
#endif
#endif
    }

  // Get the X norm of the truth solution
  // Useful for normalizing our true error data
  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));

  return libmesh_real(truth_X_norm);
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

  NumericVector<Number> * new_bf = NumericVector<Number>::build(this->comm()).release();
  new_bf->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  *new_bf = *solution;

  for(unsigned int index=0; index<get_rb_evaluation().get_n_basis_functions(); index++)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);

      Number scalar =
        inner_product_storage_vector->dot(get_rb_evaluation().get_basis_function(index));
      new_bf->add(-scalar, get_rb_evaluation().get_basis_function(index));
    }

  // Normalize new_bf
  inner_product_matrix->vector_mult(*inner_product_storage_vector, *new_bf);
  Number new_bf_norm = std::sqrt( inner_product_storage_vector->dot(*new_bf) );

  if(new_bf_norm == 0.)
    {
      new_bf->zero(); // avoid potential nan's
    }
  else
    {
      new_bf->scale(1./new_bf_norm);
    }

  // load the new basis function into the basis_functions vector.
  get_rb_evaluation().basis_functions.push_back( new_bf );
}

void RBConstruction::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();

  libMesh::out << "Updating RB residual terms" << std::endl;

  update_residual_terms();
}

Real RBConstruction::get_RB_error_bound()
{
  get_rb_evaluation().set_parameters( get_parameters() );

  Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());

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
  if(get_n_params() == 0)
    {
      Real max_val;
      if(std::numeric_limits<Real>::has_infinity)
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
  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
    {
      // Load training parameter i, this is only loaded
      // locally since the RB solves are local.
      set_params_from_training_set( first_index+i );

      training_error_bounds[i] = get_RB_error_bound();

      if(training_error_bounds[i] > max_err)
        {
          max_err_index = i;
          max_err = training_error_bounds[i];
        }
    }

  std::pair<numeric_index_type, Real> error_pair(first_index+max_err_index, max_err);
  get_global_max_error_pair(this->comm(),error_pair);

  // If we have a serial training set (i.e. a training set that is the same on all processors)
  // just set the parameters on all processors
  if(serial_training_set)
    {
      set_params_from_training_set( error_pair.first );
    }
  // otherwise, broadcast the parameter that produced the maximum error
  else
    {
      unsigned int root_id=0;
      if( (get_first_local_training_index() <= error_pair.first) &&
          (error_pair.first < get_last_local_training_index()) )
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

  UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
  temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          get_rb_evaluation().RB_Fq_vector[q_f](i) = get_Fq(q_f)->dot(get_rb_evaluation().get_basis_function(i));
        }
    }

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
          {
            get_rb_evaluation().RB_output_vectors[n][q_l](i) =
              get_output_vector(n,q_l)->dot(get_rb_evaluation().get_basis_function(i));
          }

      for(unsigned int j=0; j<RB_size; j++)
        {
          Number value = 0.;

          if(compute_RB_inner_product)
            {
              // Compute reduced inner_product_matrix
              temp->zero();
              inner_product_matrix->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

              value = temp->dot( get_rb_evaluation().get_basis_function(i) );
              get_rb_evaluation().RB_inner_product_matrix(i,j) = value;
              if(i!=j)
                {
                  // The inner product matrix is assumed
                  // to be hermitian
                  get_rb_evaluation().RB_inner_product_matrix(j,i) = libmesh_conj(value);
                }
            }

          for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
            {
              // Compute reduced Aq matrix
              temp->zero();
              get_Aq(q_a)->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

              value = (*temp).dot(get_rb_evaluation().get_basis_function(i));
              get_rb_evaluation().RB_Aq_vector[q_a](i,j) = value;

              if(i!=j)
                {
                  temp->zero();
                  get_Aq(q_a)->vector_mult(*temp, get_rb_evaluation().get_basis_function(i));

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

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    {
      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          // Initialize the vector in which we'll store the representor
          if(!get_rb_evaluation().Aq_representor[q_a][i])
            {
              get_rb_evaluation().Aq_representor[q_a][i] = (NumericVector<Number>::build(this->comm()).release());
              get_rb_evaluation().Aq_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }

          libmesh_assert(get_rb_evaluation().Aq_representor[q_a][i]->size()       == this->n_dofs()       &&
                         get_rb_evaluation().Aq_representor[q_a][i]->local_size() == this->n_local_dofs() );

          rhs->zero();
          get_Aq(q_a)->vector_mult(*rhs, get_rb_evaluation().get_basis_function(i));
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

      for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          inner_product_matrix->vector_mult(*inner_product_storage_vector,*Fq_representor[q_f]);

          for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
            {
              for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  get_rb_evaluation().Fq_Aq_representor_innerprods[q_f][q_a][i] =
                    inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a][i]);
                }
            }
        }

      unsigned int q=0;
      for(unsigned int q_a1=0; q_a1<get_rb_theta_expansion().get_n_A_terms(); q_a1++)
        {
          for(unsigned int q_a2=q_a1; q_a2<get_rb_theta_expansion().get_n_A_terms(); q_a2++)
            {
              for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  for(unsigned int j=0; j<RB_size; j++)
                    {
                      inner_product_matrix->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][j]);
                      get_rb_evaluation().Aq_Aq_representor_innerprods[q][i][j] =
                        inner_product_storage_vector->dot(*get_rb_evaluation().Aq_representor[q_a1][i]);

                      if(i != j)
                        {
                          inner_product_matrix->vector_mult(*inner_product_storage_vector, *get_rb_evaluation().Aq_representor[q_a2][i]);
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
  if(!output_dual_innerprods_computed)
    {
      // Short circuit if we don't have any outputs
      if( get_rb_theta_expansion().get_n_outputs() == 0 )
        {
          output_dual_innerprods_computed = true;
          return;
        }

      // Only log if we get to here
      LOG_SCOPE("compute_output_dual_innerprods()", "RBConstruction");

      libMesh::out << "Compute output dual inner products" << std::endl;

      // Find out the largest value of Q_l
      unsigned int max_Q_l = 0;
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        max_Q_l = (get_rb_theta_expansion().get_n_output_terms(n) > max_Q_l) ? get_rb_theta_expansion().get_n_output_terms(n) : max_Q_l;

      std::vector< NumericVector<Number> * > L_q_representor(max_Q_l);
      for(unsigned int q=0; q<max_Q_l; q++)
        {
          L_q_representor[q] = (NumericVector<Number>::build(this->comm()).release());
          L_q_representor[q]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
        }

      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        {
          for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
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
          for(unsigned int q_l1=0; q_l1<get_rb_theta_expansion().get_n_output_terms(n); q_l1++)
            {
              get_matrix_for_output_dual_solves().vector_mult(*inner_product_storage_vector, *L_q_representor[q_l1]);

              for(unsigned int q_l2=q_l1; q_l2<get_rb_theta_expansion().get_n_output_terms(n); q_l2++)
                {
                  output_dual_innerprods[n][q] = L_q_representor[q_l2]->dot(*inner_product_storage_vector);
                  libMesh::out << "output_dual_innerprods[" << n << "][" << q << "] = " << output_dual_innerprods[n][q] << std::endl;

                  q++;
                }
            }
        }

      // Finally clear the L_q_representor vectors
      for(unsigned int q=0; q<max_Q_l; q++)
        {
          if(L_q_representor[q])
            {
              delete L_q_representor[q];
              L_q_representor[q] = libmesh_nullptr;
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
  if(!Fq_representor_innerprods_computed)
    {
      // Only log if we get to here
      LOG_SCOPE("compute_Fq_representor_innerprods()", "RBConstruction");

      for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
        {
          if(!Fq_representor[q_f])
            {
              Fq_representor[q_f] = (NumericVector<Number>::build(this->comm()).release());
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

          for(unsigned int q_f1=0; q_f1<get_rb_theta_expansion().get_n_F_terms(); q_f1++)
            {
              inner_product_matrix->vector_mult(*inner_product_storage_vector, *Fq_representor[q_f1]);

              for(unsigned int q_f2=q_f1; q_f2<get_rb_theta_expansion().get_n_F_terms(); q_f2++)
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

  if(get_rb_evaluation().RB_solution.size() > get_rb_evaluation().get_n_basis_functions())
    libmesh_error_msg("ERROR: System contains " << get_rb_evaluation().get_n_basis_functions() << " basis functions." \
                      << " RB_solution vector constains " << get_rb_evaluation().RB_solution.size() << " entries." \
                      << " RB_solution in RBConstruction::load_rb_solution is too long!");

  for (std::size_t i=0; i<get_rb_evaluation().RB_solution.size(); i++)
    solution->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));

  update();
}

// The slow (but simple, non-error prone) way to compute the residual dual norm
// Useful for error checking
//Real RBConstruction::compute_residual_dual_norm(const unsigned int N)
//{
//   LOG_SCOPE("compute_residual_dual_norm()", "RBConstruction");
//
//   // Put the residual in rhs in order to compute the norm of the Riesz representor
//   // Note that this only works in serial since otherwise each processor will
//   // have a different parameter value during the Greedy training.
//
//   UniquePtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build();
//   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//
//   UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build();
//   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//
//   for(unsigned int i=0; i<N; i++)
//   {
//     RB_sol->add(RB_solution(i), get_rb_evaluation().get_basis_function(i));
//   }
//
//   this->truth_assembly();
//   matrix->vector_mult(*temp, *RB_sol);
//   rhs->add(-1., *temp);
//
//   // Then solve to get the Reisz representor
//   matrix->zero();
//   matrix->add(1., *inner_product_matrix);
//   if(constrained_problem)
//     matrix->add(1., *constraint_matrix);
//
//   solution->zero();
//   solve();
//   // Make sure we didn't max out the number of iterations
//   if( (this->n_linear_iterations() >=
//        this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
//       (this->final_linear_residual() >
//        this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
//   {
//     libmesh_error_msg("Warning: Linear solver may not have converged! Final linear residual = "
//                       << this->final_linear_residual() << ", number of iterations = "
//                       << this->n_linear_iterations());
//   }
//
//   inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//
//   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);
//
//   return std::sqrt( libmesh_real(slow_residual_norm_sq) );
//}

SparseMatrix<Number> * RBConstruction::get_inner_product_matrix()
{
  return inner_product_matrix.get();
}

SparseMatrix<Number> * RBConstruction::get_Aq(unsigned int q)
{
  if(q >= get_rb_theta_expansion().get_n_A_terms())
    libmesh_error_msg("Error: We must have q < Q_a in get_Aq.");

  return Aq_vector[q];
}

SparseMatrix<Number> * RBConstruction::get_non_dirichlet_Aq(unsigned int q)
{
  if(!store_non_dirichlet_operators)
    libmesh_error_msg("Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Aq.");

  if(q >= get_rb_theta_expansion().get_n_A_terms())
    libmesh_error_msg("Error: We must have q < Q_a in get_Aq.");

  return non_dirichlet_Aq_vector[q];
}

NumericVector<Number> * RBConstruction::get_Fq(unsigned int q)
{
  if(q >= get_rb_theta_expansion().get_n_F_terms())
    libmesh_error_msg("Error: We must have q < Q_f in get_Fq.");

  return Fq_vector[q];
}

NumericVector<Number> * RBConstruction::get_non_dirichlet_Fq(unsigned int q)
{
  if(!store_non_dirichlet_operators)
    libmesh_error_msg("Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Fq.");

  if(q >= get_rb_theta_expansion().get_n_F_terms())
    libmesh_error_msg("Error: We must have q < Q_f in get_Fq.");

  return non_dirichlet_Fq_vector[q];
}

NumericVector<Number> * RBConstruction::get_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_rb_theta_expansion().get_n_outputs()) || (q_l >= get_rb_theta_expansion().get_n_output_terms(n)) )
    libmesh_error_msg("Error: We must have n < n_outputs and "          \
                      << "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_output_vector.");

  return outputs_vector[n][q_l];
}

NumericVector<Number> * RBConstruction::get_non_dirichlet_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_rb_theta_expansion().get_n_outputs()) || (q_l >= get_rb_theta_expansion().get_n_output_terms(n)) )
    libmesh_error_msg("Error: We must have n < n_outputs and "          \
                      << "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_non_dirichlet_output_vector.");

  return non_dirichlet_outputs_vector[n][q_l];
}

void RBConstruction::get_all_matrices(std::map<std::string, SparseMatrix<Number> *> & all_matrices)
{
  all_matrices.clear();

  all_matrices["inner_product"] = get_inner_product_matrix();

  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
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

  for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
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

UniquePtr<DirichletBoundary> RBConstruction::build_zero_dirichlet_boundary_object()
{
  ZeroFunction<> zf;

  std::set<boundary_id_type> dirichlet_ids;
  std::vector<unsigned int> variables;

  // The DirichletBoundary constructor clones zf, so it's OK that zf is only in local scope
  return UniquePtr<DirichletBoundary> (new DirichletBoundary(dirichlet_ids, variables, &zf));
}

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
    if ( mkdir(riesz_representors_dir.c_str(), 0755) != 0)
      libMesh::out << "Skipping creating residual_representors directory: " << strerror(errno) << std::endl;

  for (std::size_t i=0; i<Fq_representor.size(); ++i)
    {
      if (Fq_representor[i] != libmesh_nullptr)
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
  for (std::size_t i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
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
          get_rb_evaluation().Aq_representor[i][j]->swap(*solution);

          Xdr aqr_data(file_name.str(),
                       write_binary_residual_representors ? ENCODE : WRITE);

          write_serialized_data(aqr_data, false);

          // Synchronize before moving on
          this->comm().barrier();

          // Swap back.
          get_rb_evaluation().Aq_representor[i][j]->swap(*solution);

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
  for (std::size_t i=0; i<Fq_representor.size(); ++i)
    {
      if (Fq_representor[i] != libmesh_nullptr)
        libmesh_error_msg("Error, must delete existing Fq_representor before reading in from file.");
    }

  for (std::size_t i=0; i<Fq_representor.size(); i++)
    {
      file_name.str(""); // reset filename
      file_name << riesz_representors_dir
                << "/Fq_representor" << i << riesz_representor_suffix;

      // On processor zero check to be sure the file exists
      if (this->processor_id() == 0)
        {
          int stat_result = stat(file_name.str().c_str(), &stat_info);

          if (stat_result != 0)
            libmesh_error_msg("File does not exist: " << file_name.str());
        }

      Xdr fqr_data(file_name.str(),
                   read_binary_residual_representors ? DECODE : READ);

      read_serialized_data(fqr_data, false);

      Fq_representor[i] = NumericVector<Number>::build(this->comm()).release();
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
  for (std::size_t i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
    for (std::size_t j=0; j<get_rb_evaluation().Aq_representor[i].size(); ++j)
      {
        if (get_rb_evaluation().Aq_representor[i][j] != libmesh_nullptr)
          libmesh_error_msg("Error, must delete existing Aq_representor before reading in from file.");
      }

  // Now ready to read them in from file!
  for (std::size_t i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
    for (unsigned int j=0; j<get_rb_evaluation().get_n_basis_functions(); ++j)
      {
        file_name.str(""); // reset filename
        file_name << riesz_representors_dir
                  << "/Aq_representor" << i << "_" << j << riesz_representor_suffix;

        // On processor zero check to be sure the file exists
        if (this->processor_id() == 0)
          {
            int stat_result = stat(file_name.str().c_str(), &stat_info);

            if (stat_result != 0)
              libmesh_error_msg("File does not exist: " << file_name.str());
          }

        Xdr aqr_data(file_name.str(), read_binary_residual_representors ? DECODE : READ);

        read_serialized_data(aqr_data, false);

        get_rb_evaluation().Aq_representor[i][j] = NumericVector<Number>::build(this->comm()).release();
        get_rb_evaluation().Aq_representor[i][j]->init (n_dofs(), n_local_dofs(),
                                                        false, PARALLEL);

        // No need to copy, just swap
        //*Aq_representor[i][j] = *solution;
        get_rb_evaluation().Aq_representor[i][j]->swap(*solution);
      }
}

void RBConstruction::check_convergence(LinearSolver<Number> & input_solver)
{
  libMesh::LinearConvergenceReason conv_flag;

  conv_flag = input_solver.get_converged_reason();

  if (conv_flag < 0)
    {
      std::stringstream err_msg;
      err_msg << "Convergence error. Error id: " << conv_flag;
      libmesh_error_msg(err_msg.str());
    }
}

bool RBConstruction::get_convergence_assertion_flag() const
{
  return assert_convergence;
}

void RBConstruction::set_convergence_assertion_flag(bool flag)
{
  assert_convergence = flag;
}

} // namespace libMesh
