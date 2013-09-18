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

// C++ includes
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <limits>

namespace libMesh
{

RBConstruction::RBConstruction (EquationSystems& es,
		                const std::string& name_in,
		                const unsigned int number_in)
  : Parent(es, name_in, number_in),
    inner_product_matrix(SparseMatrix<Number>::build(es.comm())),
    non_dirichlet_inner_product_matrix(SparseMatrix<Number>::build(es.comm())),
    constraint_matrix(SparseMatrix<Number>::build(es.comm())),
    constrained_problem(false),
    reuse_preconditioner(true),
    use_relative_bound_in_greedy(false),
    exit_on_repeated_greedy_parameters(true),
    write_data_during_training(false),
    impose_internal_dirichlet_BCs(false),
    impose_internal_fluxes(false),
    compute_RB_inner_product(false),
    store_non_dirichlet_operators(false),
    enforce_constraints_exactly(false),
    use_empty_rb_solve_in_greedy(true),
    Nmax(0),
    delta_N(1),
    quiet_mode(true),
    output_dual_innerprods_computed(false),
    Fq_representor_innerprods_computed(false),
    rb_eval(NULL),
    inner_product_assembly(NULL),
    constraint_assembly(NULL),
    training_tolerance(-1.)
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
  START_LOG("clear()", "RBConstruction");

  Parent::clear();

  for(unsigned int q=0; q<Aq_vector.size(); q++)
  {
    if(Aq_vector[q])
    {
      delete Aq_vector[q];
      Aq_vector[q] = NULL;
    }
  }

  for(unsigned int q=0; q<Fq_vector.size(); q++)
  {
    if(Fq_vector[q])
    {
      delete Fq_vector[q];
      Fq_vector[q] = NULL;
    }
  }

  for(unsigned int i=0; i<outputs_vector.size(); i++)
    for(unsigned int q_l=0; q_l<outputs_vector[i].size(); q_l++)
      if(outputs_vector[i][q_l])
      {
        delete outputs_vector[i][q_l];
        outputs_vector[i][q_l] = NULL;
      }

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q=0; q<non_dirichlet_Aq_vector.size(); q++)
    {
      if(non_dirichlet_Aq_vector[q])
      {
        delete non_dirichlet_Aq_vector[q];
        non_dirichlet_Aq_vector[q] = NULL;
      }
    }

    for(unsigned int q=0; q<non_dirichlet_Fq_vector.size(); q++)
    {
      if(non_dirichlet_Fq_vector[q])
      {
        delete non_dirichlet_Fq_vector[q];
        non_dirichlet_Fq_vector[q] = NULL;
      }
    }

    for(unsigned int i=0; i<non_dirichlet_outputs_vector.size(); i++)
      for(unsigned int q_l=0; q_l<non_dirichlet_outputs_vector[i].size(); q_l++)
        if(non_dirichlet_outputs_vector[i][q_l])
        {
          delete non_dirichlet_outputs_vector[i][q_l];
          non_dirichlet_outputs_vector[i][q_l] = NULL;
        }
  }

  // Also delete the Fq representors
  for(unsigned int q_f=0; q_f<Fq_representor.size(); q_f++)
  {
    if(Fq_representor[q_f])
    {
      delete Fq_representor[q_f];
      Fq_representor[q_f] = NULL;
    }
  }
  // Set Fq_representor_innerprods_computed flag to false now
  // that we've cleared the Fq representors
  Fq_representor_innerprods_computed = false;

  STOP_LOG("clear()", "RBConstruction");
}

std::string RBConstruction::system_type () const
{
  return "RBConstruction";
}

void RBConstruction::set_rb_evaluation(RBEvaluation& rb_eval_in)
{
  rb_eval = &rb_eval_in;
}

RBEvaluation& RBConstruction::get_rb_evaluation()
{
  if(!rb_eval)
  {
    libMesh::out << "Error: RBEvaluation object hasn't been initialized yet" << std::endl;
    libmesh_error();
  }

  return *rb_eval;
}

bool RBConstruction::is_rb_eval_initialized() const
{
  return (rb_eval != NULL);
}

RBThetaExpansion& RBConstruction::get_rb_theta_expansion()
{
  return get_rb_evaluation().get_rb_theta_expansion();
}

void RBConstruction::process_parameters_file (const std::string& parameters_filename)
{
  // First read in data from input_filename
  GetPot infile(parameters_filename);
  
  const unsigned int n_training_samples = infile("n_training_samples",0);
  const bool deterministic_training = infile("deterministic_training",false);
  std::string deterministic_training_parameter_name_in =
    infile("deterministic_training_parameter_name","NONE");
  const unsigned int deterministic_training_parameter_repeats_in =
    infile("deterministic_training_parameter_repeats",1);
  const std::string alternative_solver_in =
    infile("rb_alternative_solver",alternative_solver);
  const bool reuse_preconditioner_in = infile("reuse_preconditioner",
                                              reuse_preconditioner);
  const bool use_relative_bound_in_greedy_in = infile("use_relative_bound_in_greedy",
                                                      use_relative_bound_in_greedy);
  const bool write_data_during_training_in = infile("write_data_during_training",
                                                    write_data_during_training);
  unsigned int training_parameters_random_seed_in =
    static_cast<int>(-1);
  training_parameters_random_seed_in = infile("training_parameters_random_seed",
					   training_parameters_random_seed_in);
  const bool quiet_mode_in = infile("quiet_mode", quiet_mode);
  const unsigned int Nmax_in = infile("Nmax", Nmax);
  const Real training_tolerance_in = infile("training_tolerance",
                                            training_tolerance);

  // Read in the parameters from the input file too
  const unsigned int n_parameters = infile("n_parameters",1);
  RBParameters mu_min_in;
  RBParameters mu_max_in;
  RBParameters initial_mu_in;
  for(unsigned int i=0; i<n_parameters; i++)
  {
    // Read in the parameter names
    std::string param_name = infile("parameter_names", "NONE", i);

    for(unsigned int j=0; j<3; j++)
    {
      if(j==0)
      {
        Real min_val = infile(param_name, 0., j);
        mu_min_in.set_value(param_name, min_val);
      }
      else if(j==1)
      {
        Real max_val = infile(param_name, 0., j);
        mu_max_in.set_value(param_name, max_val);
      }
      else
      {
        Real init_val = infile(param_name, 0., j);
        initial_mu_in.set_value(param_name, init_val);
      }
    }
  }
  
  std::map<std::string,bool> log_scaling_in;
  RBParameters::const_iterator it     = mu_min_in.begin();
  RBParameters::const_iterator it_end = mu_min_in.end();
  unsigned int i=0;
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    log_scaling_in[param_name] = static_cast<bool>(infile("log_scaling", 0, i));
    i++;
  }

  // Set the parameters that have been read in
  set_rb_construction_parameters(n_training_samples,
                                 deterministic_training,
                                 deterministic_training_parameter_name_in,
                                 deterministic_training_parameter_repeats_in,
                                 alternative_solver_in,
                                 reuse_preconditioner_in,
                                 use_relative_bound_in_greedy_in,
                                 write_data_during_training_in,
                                 training_parameters_random_seed_in,
                                 quiet_mode_in,
                                 Nmax_in,
                                 training_tolerance_in,
                                 mu_min_in,
                                 mu_max_in,
                                 initial_mu_in,
                                 log_scaling_in);
}

void RBConstruction::set_rb_construction_parameters(
                                 unsigned int n_training_samples_in,
                                 bool deterministic_training_in,
                                 std::string deterministic_training_parameter_name_in,
                                 unsigned int deterministic_training_parameter_repeats_in,
                                 std::string alternative_solver_in,
                                 bool reuse_preconditioner_in,
                                 bool use_relative_bound_in_greedy_in,
                                 bool write_data_during_training_in,
                                 unsigned int training_parameters_random_seed_in,
                                 bool quiet_mode_in,
                                 unsigned int Nmax_in,
                                 Real training_tolerance_in,
                                 RBParameters mu_min_in,
                                 RBParameters mu_max_in,
                                 RBParameters initial_mu_in,
                                 std::map<std::string,bool> log_scaling_in)
{
  // Even if deterministic_training==false, we may specify one deterministic parameter
  set_deterministic_training_parameter_name(deterministic_training_parameter_name_in);

  // We also need to specify how many times each sample of the deterministic parameter is "repeated"
  set_deterministic_training_parameter_repeats(deterministic_training_parameter_repeats_in);

  // String which selects an alternate pc/solver combo for the update_residual_terms solves.
  // Possible values are:
  // "unchanged" -- use whatever we were using for truth solves
  // "amg" -- Use Boomeramg from Hypre.  DO NOT use on indefinite (Stokes) problems
  // "mumps" -- Use the sparse
  //update_residual_terms_solver = infile("update_residual_terms_solver",update_residual_terms_solver);
  alternative_solver = alternative_solver_in;

  // Tell the system that it is constrained (i.e. we want to use
  // the Stokes inner product matrix to compute Riesz representors)
  constrained_problem = (constraint_assembly != NULL);

  // Tell the system to reuse the preconditioner on consecutive
  // Offline solves to update residual data
  reuse_preconditioner = reuse_preconditioner_in;

  // Tell the system whether or not to use a relative error bound
  // in the Greedy algorithm
  use_relative_bound_in_greedy = use_relative_bound_in_greedy_in;

  // Tell the system whether or not to write out offline data during
  // train_reduced_basis. This allows us to continue from where the
  // training left off in case the computation stops for some reason.
  write_data_during_training = write_data_during_training_in;

  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  set_training_random_seed(training_parameters_random_seed_in);

  // Set quiet mode
  set_quiet_mode(quiet_mode_in);

  // Initialize RB parameters
  set_Nmax(Nmax_in);

  set_training_tolerance(training_tolerance_in);

  // Initialize the parameter ranges and the parameters themselves
  initialize_parameters(mu_min_in, mu_max_in, initial_mu_in);

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
  libMesh::out << "constrained_problem: " << constrained_problem << std::endl;
  libMesh::out << "Nmax: " << Nmax << std::endl;
  if(training_tolerance > 0.)
    libMesh::out << "Basis training error tolerance: " << get_training_tolerance() << std::endl;
  if( is_rb_eval_initialized() )
  {
    libMesh::out << "Aq operators attached: " << get_rb_theta_expansion().get_n_A_terms() << std::endl;
    libMesh::out << "Fq functions attached: " << get_rb_theta_expansion().get_n_F_terms() << std::endl;
    libMesh::out << "n_outputs: " << get_rb_theta_expansion().get_n_outputs() << std::endl;
    for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
      libMesh::out << "output " << n << ", Q_l = " << get_rb_theta_expansion().get_n_output_terms(n) << std::endl;
    libMesh::out << "Number of parameters: " << get_n_params() << std::endl;
  }
  else
  {
    libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
  }
  RBParameters::const_iterator it     = get_parameters().begin();
  RBParameters::const_iterator it_end = get_parameters().end();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    libMesh::out <<   "Parameter " << param_name
                 << ": Min = " << get_parameter_min(param_name)
                 << ", Max = " << get_parameter_max(param_name)
                 << ", value = " << get_parameters().get_value(param_name) << std::endl;
  }
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  if( get_deterministic_training_parameter_name() != "NONE" )
  {
    libMesh::out << "using partially random training set, deterministic parameter is: " << get_deterministic_training_parameter_name() << std::endl;
  }
  libMesh::out << "reuse preconditioner? " << reuse_preconditioner << std::endl;
  libMesh::out << "use a relative error bound in greedy? " << use_relative_bound_in_greedy << std::endl;
  libMesh::out << "write out data during basis training? " << write_data_during_training << std::endl;
  libMesh::out << "quiet mode? " << is_quiet() << std::endl;
  libMesh::out << std::endl;
}

void RBConstruction::print_basis_function_orthogonality()
{
  AutoPtr< NumericVector<Number> > temp = solution->clone();

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

void RBConstruction::set_rb_assembly_expansion(RBAssemblyExpansion& rb_assembly_expansion_in)
{
  rb_assembly_expansion = &rb_assembly_expansion_in;
}

RBAssemblyExpansion& RBConstruction::get_rb_assembly_expansion()
{
  if(!rb_assembly_expansion)
  {
    libMesh::out << "Error: RBAssemblyExpansion object hasn't been initialized yet" << std::endl;
    libmesh_error();
  }

  return *rb_assembly_expansion;
}

void RBConstruction::set_inner_product_assembly(ElemAssembly& inner_product_assembly_in)
{
  inner_product_assembly = &inner_product_assembly_in;
}

ElemAssembly& RBConstruction::get_inner_product_assembly()
{
  if(!inner_product_assembly)
  {
    libMesh::out << "Error: inner_product_assembly hasn't been initialized yet" << std::endl;
    libmesh_error();
  }

  return *inner_product_assembly;
}

void RBConstruction::set_constraint_assembly(ElemAssembly& constraint_assembly_in)
{
  constraint_assembly = &constraint_assembly_in;
}

ElemAssembly& RBConstruction::get_constraint_assembly()
{
  if(!constraint_assembly)
  {
    libMesh::out << "Error: constraint_assembly hasn't been initialized yet" << std::endl;
    libmesh_error();
  }

  return *constraint_assembly;
}

void RBConstruction::zero_constrained_dofs_on_vector(NumericVector<Number>& vector)
{
  const DofMap& dof_map = get_dof_map();

  for(dof_id_type i=dof_map.first_dof(); i<dof_map.end_dof(); i++)
  {
    if(get_dof_map().is_constrained_dof(i))
    {
      vector.set(i, 0.);
    }
  }
  vector.close();
}

void RBConstruction::initialize_rb_construction()
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


  // Perform the initialization
  allocate_data_structures();
  assemble_affine_expansion();
}

void RBConstruction::assemble_affine_expansion()
{
  // Assemble and store all of the matrices
  this->assemble_misc_matrices();
  this->assemble_all_affine_operators();
  
  // Assemble and store all of the vectors
  this->assemble_all_affine_vectors();
  this->assemble_all_output_vectors();
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
    DofMap& dof_map = this->get_dof_map();

    dof_map.attach_matrix(*inner_product_matrix);
    inner_product_matrix->init();
    inner_product_matrix->zero();

    if(this->constrained_problem)
    {
      dof_map.attach_matrix(*constraint_matrix);
      constraint_matrix->init();
      constraint_matrix->zero();
    }

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
      dof_map.attach_matrix(*non_dirichlet_inner_product_matrix);
      non_dirichlet_inner_product_matrix->init();
      non_dirichlet_inner_product_matrix->zero();

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
    Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  }

  // We also need to initialize a second set of non-Dirichlet operators
  if(store_non_dirichlet_operators)
  {
    non_dirichlet_Fq_vector.resize(get_rb_theta_expansion().get_n_F_terms());
    for(unsigned int q=0; q<get_rb_theta_expansion().get_n_F_terms(); q++)
    {
      // Initialize the memory for the vectors
      non_dirichlet_Fq_vector[q] = NumericVector<Number>::build(this->comm()).release();
      non_dirichlet_Fq_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }
  }

  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
    {
      // Initialize the memory for the truth output vectors
      outputs_vector[n][q_l] = (NumericVector<Number>::build(this->comm()).release());
      outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
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
        non_dirichlet_outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }
    }
  }

  // Resize truth_outputs vector
  truth_outputs.resize(this->get_rb_theta_expansion().get_n_outputs());
}

AutoPtr<DGFEMContext> RBConstruction::build_context ()
{
  return AutoPtr<DGFEMContext>(new DGFEMContext(*this));
}

void RBConstruction::add_scaled_matrix_and_vector(Number scalar,
                                                  ElemAssembly* elem_assembly,
                                                  SparseMatrix<Number>* input_matrix,
                                                  NumericVector<Number>* input_vector,
                                                  bool symmetrize,
                                                  bool apply_dof_constraints)
{
  START_LOG("add_scaled_matrix_and_vector()", "RBConstruction");

  bool assemble_matrix = (input_matrix != NULL);
  bool assemble_vector = (input_vector != NULL);

  if(!assemble_matrix && !assemble_vector)
    return;

  const MeshBase& mesh = this->get_mesh();

  AutoPtr<DGFEMContext> c = this->build_context();
  DGFEMContext &context  = libmesh_cast_ref<DGFEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
    elem_assembly->interior_assembly(context);

    for (context.side = 0;
         context.side != context.get_elem().n_sides();
         ++context.side )
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.get_elem().neighbor(context.get_side()) != NULL) && !impose_internal_fluxes )
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

      CouplingMatrix* coupling_matrix = get_dof_map()._dof_coupling;
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
          for(unsigned int var2=0; var2<n_vars(); var2++)
          {
            if( coupling_matrix->operator()(var1,var2) )
            {
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

  STOP_LOG("add_scaled_matrix_and_vector()", "RBConstruction");
}

void RBConstruction::set_context_solution_vec(NumericVector<Number>& vec)
{
  // Set current_local_solution = vec so that we can access
  // vec from DGFEMContext during assembly
  vec.localize
    (*current_local_solution, this->get_dof_map().get_send_list());
}

void RBConstruction::assemble_scaled_matvec(Number scalar,
                                            ElemAssembly* elem_assembly,
                                            NumericVector<Number>& dest,
                                            NumericVector<Number>& arg)
{
  START_LOG("assemble_scaled_matvec()", "RBConstruction");
  
  // This function isn't well tested lately, let's mark it as deprecated
  // In particular, it probably wouldn't work properly with DG terms
  libmesh_deprecated();

  dest.zero();

  // Set current_local_solution to be arg so that we
  // can access it from the DGFEMContext. Do this in a
  // function call so that it can be overloaded as
  // necessary in subclasses
  this->set_context_solution_vec(arg);

  const MeshBase& mesh = this->get_mesh();

  AutoPtr<DGFEMContext> c = this->build_context();
  DGFEMContext &context  = libmesh_cast_ref<DGFEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
    elem_assembly->interior_assembly(context);

    for (context.side = 0;
         context.side != context.get_elem().n_sides();
         ++context.side)
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.get_elem().neighbor(context.side) != NULL) && !impose_internal_fluxes )
        continue;

      context.side_fe_reinit();
      elem_assembly->boundary_assembly(context);
    }


    // Now perform the local matrix multiplcation
    context.get_elem_jacobian().vector_mult(context.get_elem_residual(), context.get_elem_solution());
    context.get_elem_residual() *= scalar;

    // Apply dof constraints, e.g. Dirichlet or periodic constraints
    this->get_dof_map().constrain_element_matrix_and_vector
      (context.get_elem_jacobian(), context.get_elem_residual(), context.get_dof_indices() );

    dest.add_vector (context.get_elem_residual(),
                     context.get_dof_indices() );
  }

  dest.close();

  STOP_LOG("assemble_scaled_matvec()", "RBConstruction");
}

void RBConstruction::truth_assembly()
{
  START_LOG("truth_assembly()", "RBConstruction");

  const RBParameters& mu = get_parameters();

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

    AutoPtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      *temp_vec = *get_Fq(q_f);
      temp_vec->scale( get_rb_theta_expansion().eval_F_theta(q_f, mu) );
      rhs->add(*temp_vec);
    }

    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }

  this->matrix->close();
  this->rhs->close();

  STOP_LOG("truth_assembly()", "RBConstruction");
}

void RBConstruction::assemble_inner_product_matrix(SparseMatrix<Number>* input_matrix, bool apply_dof_constraints)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1.,
                               inner_product_assembly,
                               input_matrix,
                               NULL,
                               false, /* symmetrize */
                               apply_dof_constraints);
}

void RBConstruction::assemble_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1., constraint_assembly, input_matrix, NULL);
}

void RBConstruction::assemble_and_add_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  add_scaled_matrix_and_vector(1., constraint_assembly, input_matrix, NULL);
}

void RBConstruction::assemble_Aq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix, bool apply_dof_constraints)
{
  if(q >= get_rb_theta_expansion().get_n_A_terms())
  {
    libMesh::err << "Error: We must have q < Q_a in assemble_Aq_matrix."
                 << std::endl;
    libmesh_error();
  }

  input_matrix->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_A_assembly(q),
                               input_matrix,
                               NULL,
                               false, /* symmetrize */
                               apply_dof_constraints);
}

void RBConstruction::add_scaled_Aq(Number scalar, unsigned int q_a, SparseMatrix<Number>* input_matrix, bool symmetrize)
{
  START_LOG("add_scaled_Aq()", "RBConstruction");

  if(q_a >= get_rb_theta_expansion().get_n_A_terms())
  {
    libMesh::err << "Error: We must have q < Q_a in add_scaled_Aq."
                 << std::endl;
    libmesh_error();
  }

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
                                 NULL,
                                 symmetrize);
  }

  STOP_LOG("add_scaled_Aq()", "RBConstruction");
}

void RBConstruction::assemble_misc_matrices()
{
  assemble_inner_product_matrix(inner_product_matrix.get());

  if(store_non_dirichlet_operators)
  {
    assemble_inner_product_matrix(non_dirichlet_inner_product_matrix.get(), /* apply_dof_constraints = */ false);
  }

  if( constrained_problem )
  {
    assemble_constraint_matrix(constraint_matrix.get());
  }
}

void RBConstruction::assemble_all_affine_operators()
{
  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
    assemble_Aq_matrix(q_a, get_Aq(q_a));

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
      assemble_Aq_matrix(q_a, get_non_dirichlet_Aq(q_a), false);
  }
}

void RBConstruction::assemble_all_affine_vectors()
{
  for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    assemble_Fq_vector(q_f, get_Fq(q_f));

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
      assemble_Fq_vector(q_f, get_non_dirichlet_Fq(q_f), false);
  }

}

void RBConstruction::assemble_Fq_vector(unsigned int q,
                                        NumericVector<Number>* input_vector,
                                        bool apply_dof_constraints)
{
  if(q >= get_rb_theta_expansion().get_n_F_terms())
  {
    libMesh::err << "Error: We must have q < Q_f in assemble_Fq_vector."
                 << std::endl;
    libmesh_error();
  }

  input_vector->zero();

  add_scaled_matrix_and_vector(1.,
                               &rb_assembly_expansion->get_F_assembly(q),
                               NULL,
                               input_vector,
                               false,             /* symmetrize */
                               apply_dof_constraints /* apply_dof_constraints */);
}

void RBConstruction::assemble_all_output_vectors()
{
  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
    {
      get_output_vector(n, q_l)->zero();
      add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                       NULL,
                                       get_output_vector(n,q_l),
                                       false, /* symmetrize */
                                       true   /* apply_dof_constraints */);
    }

  if(store_non_dirichlet_operators)
  {
    for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        get_non_dirichlet_output_vector(n, q_l)->zero();
        add_scaled_matrix_and_vector(1., &rb_assembly_expansion->get_output_assembly(n,q_l),
                                         NULL,
                                         get_non_dirichlet_output_vector(n,q_l),
                                         false, /* symmetrize */
                                         false  /* apply_dof_constraints */);
      }
  }
}

Real RBConstruction::train_reduced_basis(const std::string& directory_name,
                                         const bool resize_rb_eval_data)
{
  START_LOG("train_reduced_basis()", "RBConstruction");

  int count = 0;

  // initialize rb_eval's parameters
  get_rb_evaluation().initialize_parameters(*this);

  // possibly resize data structures according to Nmax
  if(resize_rb_eval_data)
  {
    get_rb_evaluation().resize_data_structures(get_Nmax());
  }

  // Clear the Greedy param list
  for(unsigned int i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
  {
    get_rb_evaluation().greedy_param_list[i].clear();
  }
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
  while(true)
  {
    libMesh::out << std::endl << "---- Basis dimension: "
                 << get_rb_evaluation().get_n_basis_functions() << " ----" << std::endl;

    if( count > 0 || (count==0 && use_empty_rb_solve_in_greedy) )
    {
      libMesh::out << "Performing RB solves on training set" << std::endl;
      training_greedy_error = compute_max_error_bound();

      libMesh::out << "Maximum " << (use_relative_bound_in_greedy ? "(relative)" : "(absolute)")
                   << " error bound is " << training_greedy_error << std::endl << std::endl;

      if(write_data_during_training)
      {
        std::stringstream new_dir_name;
        new_dir_name << directory_name << "_" << get_rb_evaluation().get_n_basis_functions();
        libMesh::out << "Writing out RB data to " << new_dir_name.str() << std::endl;
        get_rb_evaluation().write_offline_data_to_files(new_dir_name.str());
      }

      // Break out of training phase if we have reached Nmax
      // or if the training_tolerance is satisfied.
      if( greedy_termination_test(training_greedy_error, count) )
      {
        break;
      }
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
  STOP_LOG("train_reduced_basis()", "RBConstruction");

  return training_greedy_error;
}

bool RBConstruction::greedy_termination_test(Real training_greedy_error, int)
{
  if(training_greedy_error < this->training_tolerance)
  {
    libMesh::out << "Specified error tolerance reached." << std::endl;
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
    for(unsigned int i=0; i<get_rb_evaluation().greedy_param_list.size(); i++)
    {
      RBParameters& previous_parameters = get_rb_evaluation().greedy_param_list[i];
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

const RBParameters& RBConstruction::get_greedy_parameter(unsigned int i)
{
  if( i >= get_rb_evaluation().greedy_param_list.size() )
  {
    libMesh::out << "Error: Argument in RBConstruction::get_greedy_parameter is too large."
                 << std::endl;
    libmesh_error();
  }

  return get_rb_evaluation().greedy_param_list[i];
}

Real RBConstruction::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBConstruction");

  truth_assembly();

  // Safer to zero the solution first, especially when using iterative solvers
  solution->zero();
  solve();

  const RBParameters& mu = get_parameters();

  // Make sure we didn't max out the number of iterations
  if( (this->n_linear_iterations() >=
       this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
      (this->final_linear_residual() >
       this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
  {
      libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                   << this->final_linear_residual() << ", number of iterations = "
                   << this->n_linear_iterations() << std::endl << std::endl;
//     libmesh_error();
  }

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

  STOP_LOG("truth_solve()", "RBConstruction");

  return libmesh_real(truth_X_norm);
}

void RBConstruction::set_Nmax(unsigned int Nmax_in)
{
  this->Nmax = Nmax_in;
}

void RBConstruction::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "RBConstruction");

  libmesh_assert_less (i, get_rb_evaluation().get_n_basis_functions());

  *solution = get_rb_evaluation().get_basis_function(i);

  this->update();

  STOP_LOG("load_basis_function()", "RBConstruction");
}

void RBConstruction::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "RBConstruction");

  NumericVector<Number>* new_bf = NumericVector<Number>::build(this->comm()).release();
  new_bf->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
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

  STOP_LOG("enrich_RB_space()", "RBConstruction");
}

void RBConstruction::update_system()
{
  libMesh::out << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();

  libMesh::out << "Updating RB residual terms" << std::endl;

  // Note: the solves in this function employ a single system matrix and multiple
  // right-hand sides, so we may get better performance using a different
  // preconditioner, or even a direct solver.
  std::pair<std::string,std::string> orig_solver =
    this->set_alternative_solver(this->linear_solver);

  update_residual_terms();

  // Change the preconditioner, Krylov solver back to their original
  // value.  Note: does nothing if RBBase::alternative_solver ==
  // "unchanged".
  this->reset_alternative_solver(this->linear_solver, orig_solver);
}

Real RBConstruction::get_RB_error_bound()
{
  get_rb_evaluation().set_parameters( get_parameters() );

  Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());

  // Should we normalize the error bound to return a relative bound?
  if(use_relative_bound_in_greedy)
  {
    error_bound /= get_rb_evaluation().get_rb_solution_norm();
  }

  return error_bound;
}

void RBConstruction::recompute_all_residual_terms(bool compute_inner_products)
{
  // Use alternative solver for residual terms solves
  std::pair<std::string,std::string> orig_solver =
    this->set_alternative_solver(this->linear_solver);

  // Compute the basis independent terms
  Fq_representor_innerprods_computed = false;
  compute_Fq_representor_innerprods(compute_inner_products);

  // and all the basis dependent terms
  unsigned int saved_delta_N = delta_N;
  delta_N = get_rb_evaluation().get_n_basis_functions();

  update_residual_terms(compute_inner_products);

  delta_N = saved_delta_N;

  // Return to original solver
  this->reset_alternative_solver(this->linear_solver, orig_solver);
}

Real RBConstruction::compute_max_error_bound()
{
  START_LOG("compute_max_error_bound()", "RBConstruction");

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

    STOP_LOG("compute_max_error_bound()", "RBConstruction");

    // Make sure we do at least one solve, but otherwise return a zero error bound
    // when we have no parameters
    return (get_rb_evaluation().get_n_basis_functions() == 0) ? max_val : 0.;
  }

  training_error_bounds.resize(this->get_local_n_training_samples());

  // keep track of the maximum error
  unsigned int max_err_index = 0;
  Real max_err = 0.;

  unsigned int first_index = get_first_local_training_index();
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

  std::pair<unsigned int,Real> error_pair(first_index+max_err_index, max_err);
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

  STOP_LOG("compute_max_error_bound()", "RBConstruction");

  return error_pair.second;
}

void RBConstruction::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "RBConstruction");

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
  temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

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

  STOP_LOG("update_RB_system_matrices()", "RBConstruction");
}


void RBConstruction::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "RBConstruction");

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->reuse_preconditioner(false);
  }

  for(unsigned int q_a=0; q_a<get_rb_theta_expansion().get_n_A_terms(); q_a++)
  {
    for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      // Initialize the vector in which we'll store the representor
      if(!get_rb_evaluation().Aq_representor[q_a][i])
      {
        get_rb_evaluation().Aq_representor[q_a][i] = (NumericVector<Number>::build(this->comm()).release());
        get_rb_evaluation().Aq_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }

      libmesh_assert(get_rb_evaluation().Aq_representor[q_a][i]->size()       == this->n_dofs()       &&
                     get_rb_evaluation().Aq_representor[q_a][i]->local_size() == this->n_local_dofs() );

      rhs->zero();
      get_Aq(q_a)->vector_mult(*rhs, get_rb_evaluation().get_basis_function(i));
      rhs->scale(-1.);
//      zero_dirichlet_dofs_on_rhs();

      solution->zero();
      if (!is_quiet())
      {
        libMesh::out << "Starting solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBConstruction::update_residual_terms() at "
                     << Utility::get_timestamp() << std::endl;
      }

      solve();

      if (!is_quiet())
      {
        libMesh::out << "Finished solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBConstruction::update_residual_terms() at "
                     << Utility::get_timestamp() << std::endl;
        libMesh::out << this->n_linear_iterations() << " iterations, final residual "
                     << this->final_linear_residual() << std::endl;
      }

      // Make sure we didn't max out the number of iterations
      if( (this->n_linear_iterations() >=
          this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
          (this->final_linear_residual() >
          this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
      {
        libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                     << this->final_linear_residual() << ", number of iterations = "
                     << this->n_linear_iterations() << std::endl << std::endl;
//         libmesh_error();
      }

      // Store the representor
      *get_rb_evaluation().Aq_representor[q_a][i] = *solution;


      if(reuse_preconditioner)
      {
        // set this flag again in case we didn't do any F solves
        linear_solver->reuse_preconditioner(true);
      }
    }
  }

  if(reuse_preconditioner)
  {
    linear_solver->reuse_preconditioner(false);
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

  STOP_LOG("update_residual_terms()", "RBConstruction");
}

void RBConstruction::assemble_matrix_for_output_dual_solves()
{
  // By default we use the inner product matrix for steady problems

  {
    matrix->zero();
    matrix->close();
    matrix->add(1., *inner_product_matrix);
  }

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
    START_LOG("compute_output_dual_innerprods()", "RBConstruction");

    libMesh::out << "Compute output dual inner products" << std::endl;

    // Note: the solves in this function employ a single system matrix and multiple
    // right-hand sides, so we may get better performance using a different
    // preconditioner, or even a direct solver.
    std::pair<std::string,std::string> orig_solver =
      this->set_alternative_solver(this->linear_solver);

    // Find out the largest value of Q_l
    unsigned int max_Q_l = 0;
    for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
      max_Q_l = (get_rb_theta_expansion().get_n_output_terms(n) > max_Q_l) ? get_rb_theta_expansion().get_n_output_terms(n) : max_Q_l;

    std::vector< NumericVector<Number>* > L_q_representor(max_Q_l);
    for(unsigned int q=0; q<max_Q_l; q++)
    {
      L_q_representor[q] = (NumericVector<Number>::build(this->comm()).release());
      L_q_representor[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }

    if(reuse_preconditioner)
    {
      // For the first solve, make sure we generate a new preconditioner
      linear_solver->reuse_preconditioner(false);
    }

    for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      // If constrained_problem, we need to reassemble to add the constraint part back in
      if( (n==0) || constrained_problem)
      {
        assemble_matrix_for_output_dual_solves();

        if(constrained_problem)
        {
          matrix->add(1., *constraint_matrix);
        }
      }

      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
      {
        rhs->zero();
        rhs->add(1., *get_output_vector(n,q_l));
//        zero_dirichlet_dofs_on_rhs();

        solution->zero();

        if (!is_quiet())
          libMesh::out << "Starting solve n=" << n << ", q_l=" << q_l
                 << " in RBConstruction::compute_output_dual_innerprods() at "
                 << Utility::get_timestamp() << std::endl;

        solve();

        if (!is_quiet())
          {
            libMesh::out << "Finished solve n=" << n << ", q_l=" << q_l
                         << " in RBConstruction::compute_output_dual_innerprods() at "
                         << Utility::get_timestamp() << std::endl;

            libMesh::out << this->n_linear_iterations()
                         << " iterations, final residual "
                         << this->final_linear_residual() << std::endl;
          }

        // Make sure we didn't max out the number of iterations
        if( (this->n_linear_iterations() >=
             this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
            (this->final_linear_residual() >
             this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
           {
             libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                          << this->final_linear_residual() << ", number of iterations = "
                          << this->n_linear_iterations() << std::endl << std::endl;
             // libmesh_error();

           }
        *L_q_representor[q_l] = *solution;

        if(reuse_preconditioner)
          {
            // After we do a solve, tell PETSc we want to reuse the preconditioner
            // since the system matrix is not changing.
            linear_solver->reuse_preconditioner(true);
          }
      }

      // Get rid of the constraint part of the matrix before computing inner products
      if(constrained_problem)
        assemble_matrix_for_output_dual_solves();

      unsigned int q=0;
      for(unsigned int q_l1=0; q_l1<get_rb_theta_expansion().get_n_output_terms(n); q_l1++)
      {
        matrix->vector_mult(*inner_product_storage_vector, *L_q_representor[q_l1]);

        for(unsigned int q_l2=q_l1; q_l2<get_rb_theta_expansion().get_n_output_terms(n); q_l2++)
        {
          output_dual_innerprods[n][q] = L_q_representor[q_l2]->dot(*inner_product_storage_vector);
          libMesh::out << "output_dual_innerprods[" << n << "][" << q << "] = " << output_dual_innerprods[n][q] << std::endl;

          q++;
        }
      }
    }

    // reset same_preconditioner to false once all solves are finished
    if(reuse_preconditioner)
    {
      linear_solver->reuse_preconditioner(false);
    }

    // Finally clear the L_q_representor vectors
    for(unsigned int q=0; q<max_Q_l; q++)
    {
      if(L_q_representor[q])
      {
        delete L_q_representor[q];
        L_q_representor[q] = NULL;
      }
    }

    output_dual_innerprods_computed = true;

    // Change the preconditioner, Krylov solver back to their original
    // value.  Note: does nothing if RBBase::alternative_solver ==
    // "unchanged".
    this->reset_alternative_solver(this->linear_solver, orig_solver);

    STOP_LOG("compute_output_dual_innerprods()", "RBConstruction");

  }

  get_rb_evaluation().output_dual_innerprods = output_dual_innerprods;
}

void RBConstruction::compute_Fq_representor_innerprods(bool compute_inner_products)
{
  // Skip calculations if we've already computed the Fq_representors
  if(!Fq_representor_innerprods_computed)
  {
    // Only log if we get to here
    START_LOG("compute_Fq_representor_innerprods()", "RBConstruction");

    {
      matrix->zero();
      matrix->close();
      matrix->add(1., *inner_product_matrix);
      if(constrained_problem)
        matrix->add(1., *constraint_matrix);
    }

    if(reuse_preconditioner)
    {
      // For the first solve, make sure we generate a new preconditioner
      linear_solver->reuse_preconditioner(false);
    }

    for(unsigned int q_f=0; q_f<get_rb_theta_expansion().get_n_F_terms(); q_f++)
    {
      if(!Fq_representor[q_f])
      {
        Fq_representor[q_f] = (NumericVector<Number>::build(this->comm()).release());
        Fq_representor[q_f]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }

      libmesh_assert(Fq_representor[q_f]->size()       == this->n_dofs()       &&
                     Fq_representor[q_f]->local_size() == this->n_local_dofs() );

      rhs->zero();
      rhs->add(1., *get_Fq(q_f));
//      zero_dirichlet_dofs_on_rhs();

      solution->zero();

      if (!is_quiet())
        libMesh::out << "Starting solve q_f=" << q_f
		     << " in RBConstruction::update_residual_terms() at "
 		     << Utility::get_timestamp() << std::endl;

      solve();

      if (!is_quiet())
      {
        libMesh::out << "Finished solve q_f=" << q_f
		     << " in RBConstruction::update_residual_terms() at "
                     << Utility::get_timestamp() << std::endl;

        libMesh::out << this->n_linear_iterations()
                     << " iterations, final residual "
                     << this->final_linear_residual() << std::endl;
      }

      // Make sure we didn't max out the number of iterations
      if( (this->n_linear_iterations() >=
	   this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
          (this->final_linear_residual() >
           this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
      {
        libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
		     << this->final_linear_residual() << ", number of iterations = "
                     << this->n_linear_iterations() << std::endl << std::endl;
//      libmesh_error();

      }
      *Fq_representor[q_f] = *solution;

      if(reuse_preconditioner)
      {
        // After we do a solve, tell PETSc we want to reuse the preconditioner
        // since the system matrix is not changing.
        linear_solver->reuse_preconditioner(true);
      }
    }

    // Reset the same_preconditioner flag
    if(reuse_preconditioner)
    {
      linear_solver->reuse_preconditioner(false);
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

    STOP_LOG("compute_Fq_representor_innerprods()", "RBConstruction");
  }

  get_rb_evaluation().Fq_representor_innerprods = Fq_representor_innerprods;
}

void RBConstruction::load_rb_solution()
{
  START_LOG("load_rb_solution()", "RBConstruction");

  solution->zero();

  if(get_rb_evaluation().RB_solution.size() > get_rb_evaluation().get_n_basis_functions())
  {
    libMesh::err << "ERROR: System contains " << get_rb_evaluation().get_n_basis_functions() << " basis functions."
                 << " RB_solution vector constains " << get_rb_evaluation().RB_solution.size() << " entries."
                 << " RB_solution in RBConstruction::load_rb_solution is too long!" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<get_rb_evaluation().RB_solution.size(); i++)
  {
    solution->add(get_rb_evaluation().RB_solution(i), get_rb_evaluation().get_basis_function(i));
  }

  update();

  STOP_LOG("load_rb_solution()", "RBConstruction");
}

// The slow (but simple, non-error prone) way to compute the residual dual norm
// Useful for error checking
//Real RBConstruction::compute_residual_dual_norm(const unsigned int N)
//{
//  START_LOG("compute_residual_dual_norm()", "RBConstruction");
//
//   // Put the residual in rhs in order to compute the norm of the Riesz representor
//   // Note that this only works in serial since otherwise each processor will
//   // have a different parameter value during the Greedy training.
//
//   AutoPtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build();
//   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
//
//   AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
//   temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
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
//     libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
//                  << this->final_linear_residual() << ", number of iterations = "
//                  << this->n_linear_iterations() << std::endl << std::endl;
// //     libmesh_error();
//   }
//
//   inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//
//   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);
//
//  STOP_LOG("compute_residual_dual_norm()", "RBConstruction");
//
//  return std::sqrt( libmesh_real(slow_residual_norm_sq) );
//}

SparseMatrix<Number>* RBConstruction::get_inner_product_matrix()
{
  return inner_product_matrix.get();
}

SparseMatrix<Number>* RBConstruction::get_non_dirichlet_inner_product_matrix()
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true "
                 << "to access non_dirichlet_inner_product_matrix." << std::endl;
    libmesh_error();
  }

  return non_dirichlet_inner_product_matrix.get();
}

SparseMatrix<Number>* RBConstruction::get_Aq(unsigned int q)
{
  if(q >= get_rb_theta_expansion().get_n_A_terms())
  {
    libMesh::err << "Error: We must have q < Q_a in get_Aq."
                 << std::endl;
    libmesh_error();
  }

  return Aq_vector[q];
}

SparseMatrix<Number>* RBConstruction::get_non_dirichlet_Aq(unsigned int q)
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Aq." << std::endl;
    libmesh_error();
  }

  if(q >= get_rb_theta_expansion().get_n_A_terms())
  {
    libMesh::err << "Error: We must have q < Q_a in get_Aq."
                 << std::endl;
    libmesh_error();
  }

  return non_dirichlet_Aq_vector[q];
}

NumericVector<Number>* RBConstruction::get_Fq(unsigned int q)
{
  if(q >= get_rb_theta_expansion().get_n_F_terms())
  {
    libMesh::err << "Error: We must have q < Q_f in get_Fq."
                 << std::endl;
    libmesh_error();
  }

  return Fq_vector[q];
}

NumericVector<Number>* RBConstruction::get_non_dirichlet_Fq(unsigned int q)
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_Fq." << std::endl;
    libmesh_error();
  }

  if(q >= get_rb_theta_expansion().get_n_F_terms())
  {
    libMesh::err << "Error: We must have q < Q_f in get_Fq."
                 << std::endl;
    libmesh_error();
  }

  return non_dirichlet_Fq_vector[q];
}

NumericVector<Number>* RBConstruction::get_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_rb_theta_expansion().get_n_outputs()) || (q_l >= get_rb_theta_expansion().get_n_output_terms(n)) )
  {
    libMesh::err << "Error: We must have n < n_outputs and "
                 << "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_output_vector."
                 << std::endl;
    libmesh_error();
  }

  return outputs_vector[n][q_l];
}

NumericVector<Number>* RBConstruction::get_non_dirichlet_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_rb_theta_expansion().get_n_outputs()) || (q_l >= get_rb_theta_expansion().get_n_output_terms(n)) )
  {
    libMesh::err << "Error: We must have n < n_outputs and "
                 << "q_l < get_rb_theta_expansion().get_n_output_terms(n) in get_non_dirichlet_output_vector."
                 << std::endl;
    libmesh_error();
  }

  return non_dirichlet_outputs_vector[n][q_l];
}

AutoPtr<DirichletBoundary> RBConstruction::build_zero_dirichlet_boundary_object()
{
  ZeroFunction<> zf;

  std::set<boundary_id_type> dirichlet_ids;
  std::vector<unsigned int> variables;

  // The DirichletBoundary constructor clones zf, so it's OK that zf is only in local scope
  return AutoPtr<DirichletBoundary> (new DirichletBoundary(dirichlet_ids, variables, &zf));
}

void RBConstruction::write_riesz_representors_to_files(const std::string& riesz_representors_dir,
                                                       const bool write_binary_residual_representors)
{
  START_LOG("write_riesz_representors_to_files()", "RBConstruction");

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

  for (unsigned int i=0; i<Fq_representor.size(); ++i)
  {
    if (Fq_representor[i] != NULL)
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
  for (unsigned int i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
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

  STOP_LOG("write_riesz_representors_to_files()", "RBConstruction");
}



void RBConstruction::read_riesz_representors_from_files(const std::string& riesz_representors_dir,
                                                        const bool read_binary_residual_representors)
{
  START_LOG("read_riesz_representors_from_files()", "RBConstruction");

  libMesh::out << "Reading in the Fq_representors..." << std::endl;

  const std::string riesz_representor_suffix = (read_binary_residual_representors ? ".xdr" : ".dat");
  std::ostringstream file_name;
  struct stat stat_info;

  // Read in the Fq_representors.  There should be Q_f of these.  FIXME:
  // should we be worried about leaks here?
  for (unsigned int i=0; i<Fq_representor.size(); ++i)
  {
    if (Fq_representor[i] != NULL)
    {
      libMesh::out << "Error, must delete existing Fq_representor before reading in from file."
                   << std::endl;
      libmesh_error();
    }
  }

  for (unsigned int i=0; i<Fq_representor.size(); i++)
  {
    file_name.str(""); // reset filename
    file_name << riesz_representors_dir
              << "/Fq_representor" << i << riesz_representor_suffix;

    // On processor zero check to be sure the file exists
    if (this->processor_id() == 0)
    {
      int stat_result = stat(file_name.str().c_str(), &stat_info);

      if (stat_result != 0)
      {
        libMesh::out << "File does not exist: " << file_name.str() << std::endl;
        libmesh_error();
      }
    }

    Xdr fqr_data(file_name.str(),
                 read_binary_residual_representors ? DECODE : READ);

    read_serialized_data(fqr_data, false);

    Fq_representor[i] = NumericVector<Number>::build(this->comm()).release();
    Fq_representor[i]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

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
  for (unsigned int i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
    for (unsigned int j=0; j<get_rb_evaluation().Aq_representor[i].size(); ++j)
    {
      if (get_rb_evaluation().Aq_representor[i][j] != NULL)
      {
        libMesh::out << "Error, must delete existing Aq_representor before reading in from file."
                     << std::endl;
        libmesh_error();
      }
    }

  // Now ready to read them in from file!
  for (unsigned int i=0; i<get_rb_evaluation().Aq_representor.size(); ++i)
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
        {
          libMesh::out << "File does not exist: " << file_name.str() << std::endl;
          libmesh_error();
        }
      }

      Xdr aqr_data(file_name.str(), read_binary_residual_representors ? DECODE : READ);

      read_serialized_data(aqr_data, false);

      get_rb_evaluation().Aq_representor[i][j] = NumericVector<Number>::build(this->comm()).release();
      get_rb_evaluation().Aq_representor[i][j]->init (n_dofs(), n_local_dofs(),
                                            false, libMeshEnums::PARALLEL);

      // No need to copy, just swap
      //*Aq_representor[i][j] = *solution;
      get_rb_evaluation().Aq_representor[i][j]->swap(*solution);
    }

  STOP_LOG("read_riesz_representors_from_files()", "RBConstruction");
}


} // namespace libMesh
