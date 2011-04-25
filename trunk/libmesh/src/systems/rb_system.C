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

// LibMesh includes
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "dof_map.h"
#include "libmesh_logging.h"
#include "equation_systems.h"
#include "gmv_io.h"
#include "linear_solver.h"
#include "getpot.h"
#include "mesh_base.h"
#include "parallel.h"
#include "xdr_cxx.h"
#include "timestamp.h"
#include "petsc_linear_solver.h"

#include "fem_context.h"
#include "rb_system.h"
#include "rb_scm_system.h"
#include "rb_eim_system.h"
#include "rb_eim_evaluation.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "o_string_stream.h"
#include <fstream>
#include <sstream>

namespace libMesh
{

RBSystem::RBSystem (EquationSystems& es,
		    const std::string& name,
		    const unsigned int number)
  : Parent(es, name, number),
    rb_eval(NULL),
    inner_product_matrix(SparseMatrix<Number>::build()),
    non_dirichlet_inner_product_matrix(SparseMatrix<Number>::build()),
    constraint_matrix(SparseMatrix<Number>::build()),
    constrained_problem(false),
    store_basis_functions(false),
    store_representors(false),
    low_memory_mode(false),
    reuse_preconditioner(true),
    return_rel_error_bound(false),
    write_data_during_training(false),
    impose_internal_dirichlet_BCs(false),
    impose_internal_fluxes(false),
    compute_RB_inner_product(false),
    store_non_dirichlet_operators(false),
    parameters_filename(""),
    enforce_constraints_exactly(false),
    write_binary_basis_functions(true),
    read_binary_basis_functions(true),
    write_binary_residual_representors(true),
    read_binary_residual_representors(true),
    use_empty_RB_solve_in_greedy(true),
    Nmax(0),
    delta_N(1),
    quiet_mode(true),
    eigen_system_name(""),
    inner_prod_assembly(NULL),
    constraint_assembly(NULL),
    output_dual_norms_computed(false),
    Fq_representor_norms_computed(false),
    training_tolerance(-1.),
    _dirichlet_list_init(NULL),
    RB_system_initialized(false),
    current_EIM_system(NULL)
{
  rb_evaluation_objects.clear();

  // Clear the theta and assembly vectors so that we can push_back
  A_q_assembly_vector.clear();

  theta_q_f_vector.clear();
  F_q_assembly_vector.clear();

  // Make sure we clear EIM vectors so we can then push_back
  A_EIM_assembly_vector.clear();

  F_EIM_systems_vector.clear();
  F_EIM_assembly_vector.clear();

  theta_q_l_vector.clear();
  output_assembly_vector.clear();

  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;
}


RBSystem::~RBSystem ()
{
  this->clear();
}


void RBSystem::clear()
{
  START_LOG("clear()", "RBSystem");

  Parent::clear();

  for(unsigned int q=0; q<A_q_vector.size(); q++)
  {
    if(A_q_vector[q])
    {
      delete A_q_vector[q];
      A_q_vector[q] = NULL;
    }
  }

  for(unsigned int q=0; q<F_q_vector.size(); q++)
  {
    if(F_q_vector[q])
    {
      delete F_q_vector[q];
      F_q_vector[q] = NULL;
    }
  }

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q=0; q<non_dirichlet_A_q_vector.size(); q++)
    {
      if(non_dirichlet_A_q_vector[q])
      {
        delete non_dirichlet_A_q_vector[q];
        non_dirichlet_A_q_vector[q] = NULL;
      }
    }

    for(unsigned int q=0; q<non_dirichlet_F_q_vector.size(); q++)
    {
      if(non_dirichlet_F_q_vector[q])
      {
        delete non_dirichlet_F_q_vector[q];
        non_dirichlet_F_q_vector[q] = NULL;
      }
    }
  }

  for(unsigned int i=0; i<outputs_vector.size(); i++)
    for(unsigned int q_l=0; q_l<outputs_vector[i].size(); q_l++)
      if(outputs_vector[i][q_l])
      {
        delete outputs_vector[i][q_l];
        outputs_vector[i][q_l] = NULL;
      }

  // Also delete the Fq representors
  for(unsigned int q_f=0; q_f<F_q_representor.size(); q_f++)
  {
    if(F_q_representor[q_f])
    {
      delete F_q_representor[q_f];
      F_q_representor[q_f] = NULL;
    }
  }
  // Set Fq_representor_norms_computed flag to false now
  // that we've cleared the F_q representors
  Fq_representor_norms_computed = false;

  // Clear and delete all the RBEvaluation objects
  std::vector<RBEvaluation*>::iterator iter = rb_evaluation_objects.begin();
  for( ; iter != rb_evaluation_objects.end(); iter++)
  {
    RBEvaluation* eval = *iter;
    if(eval)
    {
      eval->clear(); 
      delete eval;
      eval = NULL;
    }
  }
  rb_evaluation_objects.clear();

  // Also, clear the rb_eval pointer. (The object should have been
  // deleted above, just need to set the pointer to NULL)
  if(rb_eval)
  {
    rb_eval = NULL;
  }

  STOP_LOG("clear()", "RBSystem");
}

std::string RBSystem::system_type () const
{
  return "RBSystem";
}

void RBSystem::process_parameters_file ()
{
  // First read in data from parameters_filename
  GetPot infile(parameters_filename);

  const unsigned int n_parameters = infile("n_parameters",1);
  const unsigned int n_training_samples_mu = infile("n_training_samples_mu",0);
  const bool deterministic_training = infile("deterministic_training",false);

  // String which selects an alternate pc/solver combo for the update_residual_terms solves.
  // Possible values are:
  // "unchanged" -- use whatever we were using for truth solves
  // "amg" -- Use Boomeramg from Hypre.  DO NOT use on indefinite (Stokes) problems
  // "mumps" -- Use the sparse
  //update_residual_terms_solver = infile("update_residual_terms_solver",update_residual_terms_solver);
  alternative_solver = infile("rb_alternative_solver",alternative_solver);

  // Set boolean which turns on/off storing the representor residuals
  store_representors = infile("store_representors", store_representors);

  // Booleans which control whether binary (instead of ASCII) residual
  // representors are written.  The store_representors boolean controls
  // whether these files are written out at all.
  write_binary_residual_representors = infile("write_binary_residual_representors",
					      write_binary_residual_representors);

  read_binary_residual_representors = infile("read_binary_residual_representors",
					      read_binary_residual_representors);


  // Set boolean which turns on/off writing the basis functions in
  // binary format.
  write_binary_basis_functions = infile("write_binary_basis_functions",
					write_binary_basis_functions);

  // Set boolean which turns on/off reading the basis functions in
  // binary format.
  read_binary_basis_functions = infile("read_binary_basis_functions",
					read_binary_basis_functions);

  // Tell the system that it is constrained (i.e. we want to use
  // the Stokes inner product matrix to compute Riesz representors)
  constrained_problem = (constraint_assembly != NULL);

  // Tell the system that we want to store or load the basis functions
  store_basis_functions = infile("store_basis_functions",
                                 store_basis_functions);

  // Tell the system if we're in low-memory mode
  low_memory_mode = infile("low_memory_mode",
                           low_memory_mode);

  // Tell the system to reuse the preconditioner on consecutive
  // Offline solves to update residual data
  reuse_preconditioner = infile("reuse_preconditioner",
                                reuse_preconditioner);

  // Tell the system whether or not to return a relative error bound
  // from each call to RB_solve
  return_rel_error_bound = infile("return_rel_error_bound",
                                  return_rel_error_bound);

  // Tell the system whether or not to write out offline data during
  // train_reduced_basis. This allows us to continue from where the
  // training left off in case the computation stops for some reason.
  write_data_during_training = infile("write_data_during_training",
                                      write_data_during_training);

  // Set boolean which turns on/off storing the representor residuals
  impose_internal_dirichlet_BCs = infile("impose_internal_dirichlet_BCs",
                                          impose_internal_dirichlet_BCs);

  // Set boolean which turns on/off storing the representor residuals
  impose_internal_fluxes = infile("impose_internal_fluxes",
                                   impose_internal_fluxes);

  // Set boolean flag to indicate whether or not we initialize
  // mesh dependent matrices and vectors when init_data
  // is called. Default value is true.
  initialize_mesh_dependent_data = infile("initialize_mesh_dependent_data",
                                          initialize_mesh_dependent_data);

  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  training_parameters_random_seed = infile("training_parameters_random_seed",
					   training_parameters_random_seed);
  
  // Set quiet mode
  const bool quiet_mode_in = infile("quiet_mode", quiet_mode);
  set_quiet_mode(quiet_mode_in);

  // Throw an error if we try to not initialize mesh-dependent data
  // when we also want to read in the basis functions.
  if(!initialize_mesh_dependent_data && store_basis_functions)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures if we want to read in basis "
                 << "functions."
                 << std::endl;
    libmesh_error();
  }

  // Initialize RB parameters
  const unsigned int Nmax_in = infile("Nmax", Nmax);
  set_Nmax(Nmax_in);
  
  const Real training_tolerance_in = infile("training_tolerance",
                                            training_tolerance);
  set_training_tolerance(training_tolerance_in);


  std::vector<Real> mu_min_vector(n_parameters);
  std::vector<Real> mu_max_vector(n_parameters);
  std::vector<bool> log_scaling(n_parameters);
  for(unsigned int i=0; i<n_parameters; i++)
  {
    // Read vector-based mu_min values.
    mu_min_vector[i] = infile("mu_min", mu_min_vector[i], i);

    // Read vector-based mu_max values.
    mu_max_vector[i] = infile("mu_max", mu_max_vector[i], i);

    // Read vector-based log scaling values.  Note the intermediate conversion to
    // int... this implies log_scaling = '1 1 1...' in the input file.
    log_scaling[i] = static_cast<bool>(infile("log_scaling", static_cast<int>(log_scaling[i]), i));
  }

  initialize_training_parameters(mu_min_vector,
                                 mu_max_vector,
                                 n_training_samples_mu,
                                 log_scaling,
                                 deterministic_training);   // use deterministic parameters



  // Set the initial parameter value to the minimum parameters
  set_current_parameters(mu_min_vector);

  libMesh::out << std::endl << "RBSystem parameters:" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "constrained_problem: " << constrained_problem << std::endl;
  libMesh::out << "Nmax: " << Nmax << std::endl;
  if(training_tolerance > 0.)
    libMesh::out << "Basis training error tolerance: " << get_training_tolerance() << std::endl;
  libMesh::out << "A_q operators attached: " << get_Q_a() << std::endl;
  libMesh::out << "F_q functions attached: " << get_Q_f() << std::endl;
  libMesh::out << "Number of A EIM systems: " << get_n_A_EIM_systems() << std::endl;
  libMesh::out << "Number of F EIM systems: " << get_n_F_EIM_systems() << std::endl;
  libMesh::out << "n_outputs: " << get_n_outputs() << std::endl;
  for(unsigned int n=0; n<get_n_outputs(); n++)
    libMesh::out << "output " << n << ", Q_l = " << get_Q_l(n) << std::endl;
  for(unsigned int i=0; i<n_parameters; i++)
  {
    libMesh::out <<   "Parameter " << i
                 << ": Min = " << get_parameter_min(i)
                 << ", Max = " << get_parameter_max(i)
                 << ", log scaling = " << log_scaling[i] << std::endl;
  }
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << "using deterministic training samples? " << deterministic_training << std::endl;
  libMesh::out << "store/load basis functions? " << store_basis_functions << std::endl;
  if(store_basis_functions)
  {
    libMesh::out << "  write out basis functions in binary format? "
                 << write_binary_basis_functions << std::endl;
    libMesh::out << "  read in basis functions in binary format? "
                 << read_binary_basis_functions << std::endl;
  }
  libMesh::out << "store/load residual representors? " << store_representors << std::endl;
  if(store_representors)
  {
    libMesh::out << "  write out residual representors in binary format? "
                 << write_binary_residual_representors << std::endl;
    libMesh::out << "  read in residual representors in binary format? "
                 << read_binary_residual_representors << std::endl;
  }
  libMesh::out << "low-memory mode? " << low_memory_mode << std::endl;
  libMesh::out << "reuse preconditioner? " << reuse_preconditioner << std::endl;
  libMesh::out << "return a relative error bound from RB_solve? " << return_rel_error_bound << std::endl;
  libMesh::out << "write out data during basis training? " << write_data_during_training << std::endl;
  libMesh::out << "initializing mesh-dependent data structures? "
               << initialize_mesh_dependent_data << std::endl;
  libMesh::out << "impose internal Dirichlet BCs? " << impose_internal_dirichlet_BCs << std::endl;
  libMesh::out << "impose internal fluxes? " << impose_internal_fluxes << std::endl;
  libMesh::out << "quiet mode? " << is_quiet() << std::endl;
  libMesh::out << "parameter initialized to mu_min: ";
  for(unsigned int i=0; i<n_parameters; i++)
  {
    libMesh::out << "mu[" << i << "] = " << get_current_parameters()[i];
    if(i < (n_parameters-1))
      libMesh::out << ", ";
    else
      libMesh::out << std::endl;
  }
  libMesh::out << std::endl;

  // We need Nmax to be initialized
  libmesh_assert(Nmax > 0);
}

void RBSystem::initialize_RB_system(bool do_not_assemble)
{
  process_parameters_file();
  allocate_data_structures();

  // Build a new RBEvaluation object
  libmesh_assert( rb_evaluation_objects.empty() );
  rb_eval = add_new_rb_evaluation_object();
  
  // And initialize rb_eval
  rb_eval->initialize();

  RB_system_initialized = true;

  if(!do_not_assemble)
  {
    // Initialize the non-Dirichlet and Dirichlet dofs lists
    this->initialize_dirichlet_dofs();

    // Assemble and store all of the matrices if we're
    // not in low-memory mode
    if(!low_memory_mode)
    {
      this->assemble_misc_matrices();
      this->assemble_all_affine_operators();
    }

    this->assemble_all_affine_vectors();
    this->assemble_all_output_vectors();
  }
}

void RBSystem::allocate_data_structures()
{
  // Resize vectors for storing mesh-dependent data but only
  // initialize if initialize_mesh_dependent_data == true
  A_q_vector.resize(get_Q_a());
  F_q_vector.resize(get_Q_f());
  
  // Resize the F_q_representors and initialize each to NULL
  // These are basis independent and hence stored here, whereas
  // the A_q_representors are stored in RBEvaluation
  F_q_representor.resize(get_Q_f());

  // Initialize vectors for the norms of the Fq representors
  // These are basis independent and therefore stored here.
  unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
  Fq_representor_norms.resize(Q_f_hat);

  // Resize the output vectors
  outputs_vector.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
    outputs_vector[n].resize( get_Q_l(n) );

  // Resize the output dual norm vectors
  output_dual_norms.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
    output_dual_norms[n].resize(Q_l_hat);
  }

  if(initialize_mesh_dependent_data)
  {
    // Only initialize matrices if we're not in low-memory mode
    if(!low_memory_mode)
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

      for(unsigned int q=0; q<get_Q_a(); q++)
      {
        // Initialize the memory for the matrices
        A_q_vector[q] = SparseMatrix<Number>::build().release();
        dof_map.attach_matrix(*A_q_vector[q]);
        A_q_vector[q]->init();
        A_q_vector[q]->zero();
      }
      
      // We also need to initialize a second set of non-Dirichlet operators
      if(store_non_dirichlet_operators)
      {
        dof_map.attach_matrix(*non_dirichlet_inner_product_matrix);
        non_dirichlet_inner_product_matrix->init();
        non_dirichlet_inner_product_matrix->zero();
        
        non_dirichlet_A_q_vector.resize(get_Q_a());
        for(unsigned int q=0; q<get_Q_a(); q++)
        {
          // Initialize the memory for the matrices
          non_dirichlet_A_q_vector[q] = SparseMatrix<Number>::build().release();
          dof_map.attach_matrix(*non_dirichlet_A_q_vector[q]);
          non_dirichlet_A_q_vector[q]->init();
          non_dirichlet_A_q_vector[q]->zero();
        }
      }
    }

    // Initialize the vectors even if we are in low-memory mode
    for(unsigned int q=0; q<get_Q_f(); q++)
    {
      // Initialize the memory for the vectors
      F_q_vector[q] = NumericVector<Number>::build().release();
      F_q_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }

    // We also need to initialize a second set of non-Dirichlet operators
    if(store_non_dirichlet_operators)
    {
      non_dirichlet_F_q_vector.resize(get_Q_f());
      for(unsigned int q=0; q<get_Q_f(); q++)
      {
        // Initialize the memory for the vectors
        non_dirichlet_F_q_vector[q] = NumericVector<Number>::build().release();
        non_dirichlet_F_q_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }
    }

    for(unsigned int n=0; n<get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        // Initialize the memory for the truth output vectors
        outputs_vector[n][q_l] = (NumericVector<Number>::build().release());
        outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }
  }

  // Resize truth_outputs vector
  truth_outputs.resize(this->get_n_outputs());
}

NumericVector<Number>& RBSystem::get_basis_function(unsigned int i)
{
  return rb_eval->get_basis_function(i);
}

RBEvaluation* RBSystem::add_new_rb_evaluation_object()
{
  RBEvaluation* e = new RBEvaluation(*this);
  rb_evaluation_objects.push_back(e);

  return e;
}

void RBSystem::attach_dirichlet_dof_initialization (DirichletDofAssembly* dirichlet_init)
{
  libmesh_assert (dirichlet_init != NULL);

  _dirichlet_list_init = dirichlet_init;
}

void RBSystem::initialize_dirichlet_dofs()
{
  // Short-circuit if _dirichlet_list_init is NULL
  if(!_dirichlet_list_init)
  {
    return;
  }

  START_LOG("initialize_dirichlet_dofs()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to initialize Dirichlet dofs."
                 << std::endl;
    libmesh_error();
  }

  // Initialize the lists of Dirichlet and non-Dirichlet degrees-of-freedom
  // Clear the set to store the Dirichlet dofs on this processor
  _dirichlet_list_init->dirichlet_dofs_set.clear();

  const MeshBase& mesh = this->get_equation_systems().get_mesh();

  AutoPtr<FEMContext> c = this->build_context();
  FEMContext &context  = libmesh_cast_ref<FEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();

    for (context.side = 0;
         context.side != context.elem->n_sides();
       ++context.side)
    {
      // Skip over non-boundary sides if we don't have internal Dirichlet BCs
      if ( (context.elem->neighbor(context.side) != NULL) && !impose_internal_dirichlet_BCs )
	continue;

      context.side_fe_reinit();
      _dirichlet_list_init->boundary_assembly(context);
    }
  }

  // Initialize the dirichlet dofs vector on each processor
  std::vector<unsigned int> dirichlet_dofs_vector;
  dirichlet_dofs_vector.clear();

  std::set<unsigned int>::iterator iter     = _dirichlet_list_init->dirichlet_dofs_set.begin();
  std::set<unsigned int>::iterator iter_end = _dirichlet_list_init->dirichlet_dofs_set.end();

  for ( ; iter != iter_end; ++iter)
  {
    unsigned int dirichlet_dof_index = *iter;
    dirichlet_dofs_vector.push_back(dirichlet_dof_index);
  }

  // Now take the union over all processors
  Parallel::allgather(dirichlet_dofs_vector);

  // Also, initialize the member data structure global_dirichlet_dofs_set
  global_dirichlet_dofs_set.clear();

  for (unsigned int ii=0; ii<dirichlet_dofs_vector.size(); ii++)
  {
    global_dirichlet_dofs_set.insert(dirichlet_dofs_vector[ii]);
  }

  STOP_LOG("initialize_dirichlet_dofs()", "RBSystem");
}

AutoPtr<FEMContext> RBSystem::build_context ()
{
  return AutoPtr<FEMContext>(new FEMContext(*this));
}

void RBSystem::add_scaled_matrix_and_vector(Number scalar,
                                            ElemAssembly* elem_assembly,
                                            SparseMatrix<Number>* input_matrix,
                                            NumericVector<Number>* input_vector,
                                            bool symmetrize,
                                            bool apply_dirichlet_bc)
{
  START_LOG("add_scaled_matrix_and_vector()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to perform add_scaled_matrix_and_vector."
                 << std::endl;
    libmesh_error();
  }

  bool assemble_matrix = (input_matrix != NULL);
  bool assemble_vector = (input_vector != NULL);

  if(!assemble_matrix && !assemble_vector)
    return;

  const MeshBase& mesh = this->get_mesh();

  AutoPtr<FEMContext> c = this->build_context();
  FEMContext &context  = libmesh_cast_ref<FEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
    elem_assembly->interior_assembly(context);

    for (context.side = 0;
          context.side != context.elem->n_sides();
          ++context.side)
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.elem->neighbor(context.side) != NULL) && !impose_internal_fluxes )
        continue;

      // Impose boundary (e.g. Neumann) term
      context.side_fe_reinit();
      elem_assembly->boundary_assembly(context);
    }

    // Need to symmetrize before imposing
    // periodic constraints
    if(assemble_matrix && symmetrize)
    {
      DenseMatrix<Number> Ke_transpose;
      context.elem_jacobian.get_transpose(Ke_transpose);
      context.elem_jacobian += Ke_transpose;
      context.elem_jacobian *= 0.5;
    }

    // Apply constraints, e.g. periodic constraints
    this->get_dof_map().constrain_element_matrix_and_vector
      (context.elem_jacobian, context.elem_residual, context.dof_indices);

    if(apply_dirichlet_bc)
    {
      // Apply Dirichlet boundary conditions, we assume zero Dirichlet BCs
      // Note that this cannot be inside the side-loop since non-boundary
      // elements may contain boundary dofs
      std::set<unsigned int>::const_iterator iter;
      for(unsigned int n=0; n<context.dof_indices.size(); n++)
      {
        iter = global_dirichlet_dofs_set.find( context.dof_indices[n] );
        if(iter != global_dirichlet_dofs_set.end())
        {
	  context.elem_jacobian.condense
	    (n,n,0.,context.elem_residual);
        }
      }
    }

    // Scale and add to global matrix and/or vector
    context.elem_jacobian *= scalar;
    context.elem_residual *= scalar;

    if(assemble_matrix)
      input_matrix->add_matrix (context.elem_jacobian,
                                context.dof_indices);
    if(assemble_vector)
      input_vector->add_vector (context.elem_residual,
                                context.dof_indices);
  }

  if(assemble_matrix)
    input_matrix->close();
  if(assemble_vector)
    input_vector->close();

  STOP_LOG("add_scaled_matrix_and_vector()", "RBSystem");
}

void RBSystem::set_context_solution_vec(NumericVector<Number>& vec)
{
  // Set current_local_solution = vec so that we can access
  // vec from FEMContext during assembly
  vec.localize
    (*current_local_solution, this->get_dof_map().get_send_list());
}

void RBSystem::assemble_scaled_matvec(Number scalar,
                                      ElemAssembly* elem_assembly,
                                      NumericVector<Number>& dest,
                                      NumericVector<Number>& arg)
{
  START_LOG("assemble_scaled_matvec()", "RBSystem");

  dest.zero();

  // Set current_local_solution to be arg so that we
  // can access it from the FEMContext. Do this in a
  // function call so that it can be overloaded as
  // necessary (e.g. for QNTransientRBSystem)
  this->set_context_solution_vec(arg);

  const MeshBase& mesh = this->get_mesh();

  AutoPtr<FEMContext> c = this->build_context();
  FEMContext &context  = libmesh_cast_ref<FEMContext&>(*c);

  this->init_context(context);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    context.pre_fe_reinit(*this, *el);
    context.elem_fe_reinit();
    elem_assembly->interior_assembly(context);

    for (context.side = 0;
         context.side != context.elem->n_sides();
         ++context.side)
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.elem->neighbor(context.side) != NULL) && !impose_internal_fluxes )
        continue;

      context.side_fe_reinit();
      elem_assembly->boundary_assembly(context);
    }


    // Now perform the local matrix multiplcation
    context.elem_jacobian.vector_mult(context.elem_residual, context.elem_solution);
    context.elem_residual *= scalar;

    // Apply constraints, e.g. periodic constraints
    this->get_dof_map().constrain_element_matrix_and_vector
      (context.elem_jacobian, context.elem_residual, context.dof_indices);
      
    // Apply Dirichlet boundary conditions, we assume zero Dirichlet BCs
    // This zeros the Dirichlet dofs in context.elem_residual
    std::set<unsigned int>::const_iterator iter;
    for(unsigned int n=0; n<context.dof_indices.size(); n++)
    {
      iter = global_dirichlet_dofs_set.find( context.dof_indices[n] );
      if(iter != global_dirichlet_dofs_set.end())
      {
        context.elem_jacobian.condense
          (n,n,0.,context.elem_residual);
      }
    }

    dest.add_vector (context.elem_residual,
                      context.dof_indices);
  }

  dest.close();

  STOP_LOG("assemble_scaled_matvec()", "RBSystem");
}

void RBSystem::truth_assembly()
{
  START_LOG("truth_assembly()", "RBSystem");

  this->matrix->zero();
  this->rhs->zero();

  if(!low_memory_mode)
  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      matrix->add(eval_theta_q_a(q_a), *get_A_q(q_a));
    }

    AutoPtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build();
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      *temp_vec = *get_F_q(q_f);
      temp_vec->scale( eval_theta_q_f(q_f) );
      rhs->add(*temp_vec);
    }

    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }
  else
  {
    // In low memory mode we do not store the matrices
    // from the affine expansion, so need to assemble

    // For efficiency (i.e. to avoid doing Q_a+Q_f loops
    // over the mesh) we do not use add_scaled_matrix_and_vector
    // here

    const MeshBase& mesh = this->get_mesh();

    std::vector<FEMContext*> Aq_context(get_Q_a());
    for(unsigned int q_a=0; q_a<Aq_context.size(); q_a++)
    {
      Aq_context[q_a] = this->build_context().release();
      this->init_context(*Aq_context[q_a]);
    }

    std::vector<FEMContext*> Fq_context(get_Q_f());
    for(unsigned int q_f=0; q_f<Fq_context.size(); q_f++)
    {
      Fq_context[q_f] = this->build_context().release();
      this->init_context(*Fq_context[q_f]);
    }

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        Aq_context[q_a]->pre_fe_reinit(*this, *el);
        Aq_context[q_a]->elem_fe_reinit();
        A_q_assembly_vector[q_a]->interior_assembly(*Aq_context[q_a]);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        Fq_context[q_f]->pre_fe_reinit(*this, *el);
        Fq_context[q_f]->elem_fe_reinit();
        F_q_assembly_vector[q_f]->interior_assembly(*Fq_context[q_f]);
      }

      for (Aq_context[0]->side = 0;
            Aq_context[0]->side != Aq_context[0]->elem->n_sides();
            ++Aq_context[0]->side)
      {
        // May not need to apply fluxes on non-boundary elements
        if( (Aq_context[0]->elem->neighbor(Aq_context[0]->side) != NULL) && !impose_internal_fluxes )
          continue;

        // Update the side information for all contexts
        for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
        {
          // Update the side information for all contexts
          Aq_context[q_a]->side = Aq_context[0]->side;

          Aq_context[q_a]->side_fe_reinit();
          A_q_assembly_vector[q_a]->boundary_assembly(*Aq_context[q_a]);
        }

        // Impose boundary terms, e.g. Neuman BCs
        for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
        {
          // Update the side information for all contexts
          Fq_context[q_f]->side = Aq_context[0]->side;

          Fq_context[q_f]->side_fe_reinit();
          F_q_assembly_vector[q_f]->boundary_assembly(*Fq_context[q_f]);
        }
      }

      // Constrain the dofs to impose hanging node or periodic constraints
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        this->get_dof_map().constrain_element_matrix
          (Aq_context[q_a]->elem_jacobian, Aq_context[q_a]->dof_indices);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        this->get_dof_map().constrain_element_vector
          (Fq_context[q_f]->elem_residual, Fq_context[q_f]->dof_indices);
      }
      
      // Apply Dirichlet boundary conditions, we assume zero Dirichlet BCs
      // Note that this cannot be inside the side-loop since non-boundary
      // elements may contain boundary dofs
      std::set<unsigned int>::const_iterator iter;
      for(unsigned int n=0; n<Aq_context[0]->dof_indices.size(); n++)
      {
	iter = global_dirichlet_dofs_set.find( Aq_context[0]->dof_indices[n] );
	if(iter != global_dirichlet_dofs_set.end())
	{
	  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
	  {
	    Aq_context[q_a]->elem_jacobian.condense
	      (n,n,0.,Aq_context[q_a]->elem_residual);
	  }

	  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	  {
	    Fq_context[q_f]->elem_jacobian.condense
	      (n,n,0.,Fq_context[q_f]->elem_residual);
	  }
	}
      }
      
      // Finally add local matrices/vectors to global system
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        // Scale by theta_q_a
        Aq_context[q_a]->elem_jacobian *= eval_theta_q_a(q_a);
        this->matrix->add_matrix (Aq_context[q_a]->elem_jacobian,
                                  Aq_context[q_a]->dof_indices);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        // Scale by theta_q_f
        Fq_context[q_f]->elem_residual *= eval_theta_q_f(q_f);
        this->rhs->add_vector (Fq_context[q_f]->elem_residual,
                              Fq_context[q_f]->dof_indices);
      }
    }

    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, matrix, NULL);

    // Delete all the ptrs to FEMContexts!
    for(unsigned int q_a=0; q_a<Aq_context.size(); q_a++)
    {
      delete Aq_context[q_a];
      Aq_context[q_a] = NULL;
    }
    Aq_context.clear();

    for(unsigned int q_f=0; q_f<Fq_context.size(); q_f++)
    {
      delete Fq_context[q_f];
      Fq_context[q_f];
    }
    Fq_context.clear();
  }

  this->matrix->close();
  this->rhs->close();

  STOP_LOG("truth_assembly()", "RBSystem");
}

void RBSystem::assemble_inner_product_matrix(SparseMatrix<Number>* input_matrix, bool apply_dirichlet_bc)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1.,
                               inner_prod_assembly,
                               input_matrix,
                               NULL,
                               false, /* symmetrize */
                               apply_dirichlet_bc);
}

void RBSystem::assemble_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1., constraint_assembly, input_matrix, NULL);
}

void RBSystem::assemble_and_add_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  add_scaled_matrix_and_vector(1., constraint_assembly, input_matrix, NULL);
}

void RBSystem::assemble_Aq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix, bool apply_dirichlet_bc)
{
  if(q >= get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in assemble_Aq_matrix."
                 << std::endl;
    libmesh_error();
  }

  input_matrix->zero();

  if(!is_A_EIM_operator(q))
  {
    add_scaled_matrix_and_vector(1.,
                                 A_q_assembly_vector[q],
                                 input_matrix,
                                 NULL,
                                 false, /* symmetrize */
                                 apply_dirichlet_bc);
  }
  else // We have an EIM function
  {
    // Find out the LHS EIM system index and the associated function index that corresponds to q
    std::pair<unsigned int, unsigned int> A_EIM_indices = get_A_EIM_indices(q);
    
    unsigned int system_id = A_EIM_indices.first;
    current_EIM_system = A_EIM_systems_vector[system_id];
    
    // Cache a GHOSTED version of the EIM basis function corresponding to the current EIM function index
    current_EIM_system->cache_ghosted_basis_function(A_EIM_indices.second);
    
    add_scaled_matrix_and_vector(1.,
                                 A_EIM_assembly_vector[system_id],
                                 input_matrix,
                                 NULL,
                                 false, /* symmetrize */
                                 apply_dirichlet_bc);
  }
}

void RBSystem::add_scaled_Aq(Number scalar, unsigned int q_a, SparseMatrix<Number>* input_matrix, bool symmetrize)
{
  START_LOG("add_scaled_Aq()", "RBSystem");

  if(q_a >= get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in add_scaled_Aq."
                 << std::endl;
    libmesh_error();
  }

  if(!low_memory_mode && !symmetrize)
  {
    input_matrix->add(scalar, *get_A_q(q_a));
    input_matrix->close();
  }
  else
  {
    if(!is_A_EIM_operator(q_a))
    {
      add_scaled_matrix_and_vector(scalar,
                                   A_q_assembly_vector[q_a],
                                   input_matrix,
                                   NULL,
                                   symmetrize);
    }
    else // We have an EIM function
    {
      // Find out the LHS EIM system index and the associated function index that corresponds to q_a
      std::pair<unsigned int, unsigned int> A_EIM_indices = get_A_EIM_indices(q_a);

      unsigned int system_id = A_EIM_indices.first;
      current_EIM_system = A_EIM_systems_vector[system_id];
      
      // Cache a GHOSTED version of the EIM basis function corresponding to the current EIM function index
      current_EIM_system->cache_ghosted_basis_function(A_EIM_indices.second);

      add_scaled_matrix_and_vector(scalar,
                                   A_EIM_assembly_vector[system_id],
                                   input_matrix,
                                   NULL,
                                   symmetrize);
    }
  }

  STOP_LOG("add_scaled_Aq()", "RBSystem");
}

void RBSystem::assemble_misc_matrices()
{
  if(low_memory_mode)
  {
    libMesh::out << "Error: Cannot store misc matrices in low-memory mode." << std::endl;
    libmesh_error();
  }

  assemble_inner_product_matrix(inner_product_matrix.get());

  if(store_non_dirichlet_operators)
  {
    assemble_inner_product_matrix(non_dirichlet_inner_product_matrix.get(), /* apply_dirichlet_bc = */ false);
  }

  if( constrained_problem )
    assemble_constraint_matrix(constraint_matrix.get());
}

void RBSystem::assemble_all_affine_operators()
{
  if(low_memory_mode)
  {
    libMesh::out << "Error: Cannot store affine matrices in low-memory mode." << std::endl;
    libmesh_error();
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    assemble_Aq_matrix(q_a, get_A_q(q_a));

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      assemble_Aq_matrix(q_a, get_non_dirichlet_A_q(q_a), false);
  }
}

void RBSystem::assemble_all_affine_vectors()
{
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    assemble_Fq_vector(q_f, get_F_q(q_f));

  if(store_non_dirichlet_operators)
  {
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      assemble_Fq_vector(q_f, get_non_dirichlet_F_q(q_f), false);
  }

}

void RBSystem::assemble_Fq_vector(unsigned int q,
                                  NumericVector<Number>* input_vector,
                                  bool apply_dirichlet_bc)
{
  if(q >= get_Q_f())
  {
    libMesh::err << "Error: We must have q < Q_f in assemble_Fq_vector."
                 << std::endl;
    libmesh_error();
  }

  input_vector->zero();
    
  if(!is_F_EIM_function(q))
  {
    add_scaled_matrix_and_vector(1.,
                                 F_q_assembly_vector[q],
                                 NULL,
                                 input_vector,
                                 false,             /* symmetrize */
                                 apply_dirichlet_bc /* apply_dirichlet_bc */);
  }
  else // We have an EIM function
  {
    // Find out the LHS EIM system index and the associated function index that corresponds to q
    std::pair<unsigned int, unsigned int> F_EIM_indices = get_F_EIM_indices(q);
    
    unsigned int system_id = F_EIM_indices.first;
    current_EIM_system = F_EIM_systems_vector[system_id];
      
    // Cache a GHOSTED version of the EIM basis function corresponding to the current EIM function index
    current_EIM_system->cache_ghosted_basis_function(F_EIM_indices.second);
    
    add_scaled_matrix_and_vector(1.,
                                 F_EIM_assembly_vector[system_id],
                                 NULL,
                                 input_vector,
                                 false,             /* symmetrize */
                                 apply_dirichlet_bc /* apply_dirichlet_bc */);
  }
}

void RBSystem::assemble_all_output_vectors()
{
  for(unsigned int n=0; n<get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      get_output_vector(n, q_l)->zero();
      add_scaled_matrix_and_vector(1., output_assembly_vector[n][q_l],
                                       NULL,
                                       get_output_vector(n,q_l));
    }
}

Real RBSystem::train_reduced_basis(const std::string& directory_name)
{
  START_LOG("train_reduced_basis()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to train reduced basis."
                 << std::endl;
    libmesh_error();
  }

  int count = 0;

  // Clear the Greedy param list
  for(unsigned int i=0; i<rb_eval->greedy_param_list.size(); i++)
    rb_eval->greedy_param_list[i].clear();
  rb_eval->greedy_param_list.clear();

  Real training_greedy_error;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if (this->get_n_basis_functions() >= Nmax) return 0.;
  

  // Compute the dual norms of the outputs if we haven't already done so
  compute_output_dual_norms();

  // Compute the Fq Riesz representor dual norms if we haven't already done so
  compute_Fq_representor_norms();

  libMesh::out << std::endl << "---- Performing Greedy basis enrichment ----" << std::endl;
  while(true)
  {
    libMesh::out << std::endl << "---- Basis dimension: "
                 << get_n_basis_functions() << " ----" << std::endl;

    if( count > 0 || (count==0 && use_empty_RB_solve_in_greedy) )
    {
      libMesh::out << "Performing RB solves on training set" << std::endl;
      training_greedy_error = compute_max_error_bound();

      libMesh::out << "Maximum " << (return_rel_error_bound ? "(relative)" : "(absolute)")
                   << " error bound is " << training_greedy_error << std::endl << std::endl;

      if(write_data_during_training)
      {
        OStringStream new_dir_name;
        new_dir_name << directory_name << "_" << get_n_basis_functions();
        libMesh::out << "Writing out RB data to " << new_dir_name.str() << std::endl;
        write_offline_data_to_files(new_dir_name.str());
      }

      // Break out of training phase if we have reached Nmax
      // or if the training_tolerance is satisfied.
      if( greedy_termination_test(training_greedy_error, count) )
      {
        break;
      }
    }

    libMesh::out << "Performing truth solve at parameter:" << std::endl;
    print_current_parameters();

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
  STOP_LOG("train_reduced_basis()", "RBSystem");

  return training_greedy_error;
}

Real RBSystem::eval_output_dual_norm(unsigned int n)
{
  Number output_bound_sq = 0.;    
  unsigned int q=0;
  for(unsigned int q_l1=0; q_l1<get_Q_l(n); q_l1++)
  {
    for(unsigned int q_l2=q_l1; q_l2<get_Q_l(n); q_l2++)
    {
      Real delta = (q_l1==q_l2) ? 1. : 2.;
      output_bound_sq += delta * libmesh_real(
        libmesh_conj(eval_theta_q_l(n,q_l1))*eval_theta_q_l(n,q_l2) * output_dual_norms[n][q] );
      q++;
    }
  }
    
  return libmesh_real(std::sqrt( output_bound_sq ));
}

bool RBSystem::greedy_termination_test(Real training_greedy_error, int)
{
  if(training_greedy_error < this->training_tolerance)
  {
    libMesh::out << "Specified error tolerance reached." << std::endl;
    return true;
  }

  if(get_n_basis_functions() >= this->get_Nmax())
  {
    libMesh::out << "Maximum number of basis functions reached: Nmax = "
                 << get_Nmax() << std::endl;
    return true;
  }

  return false;
}

void RBSystem::update_greedy_param_list()
{
  rb_eval->greedy_param_list.push_back( get_current_parameters() );
}

std::vector<Real> RBSystem::get_greedy_parameter(unsigned int i)
{
  if( i >= rb_eval->greedy_param_list.size() )
  {
    libMesh::out << "Error: Argument in RBSystem::get_greedy_parameter is too large."
                 << std::endl;
    libmesh_error();
  }

  return rb_eval->greedy_param_list[i];
}

Real RBSystem::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to do a truth solve."
                 << std::endl;
    libmesh_error();
  }

  truth_assembly();
  
  // Safer to zero the solution first, especially when using iterative solvers
  solution->zero();
  solve();

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

  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    truth_outputs[n] = 0.;
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      truth_outputs[n] += libmesh_conj(eval_theta_q_l(n, q_l))*get_output_vector(n,q_l)->dot(*solution);
  }

  if(plot_solution > 0)
  {
    const MeshBase& mesh = get_mesh();
    GMVIO(mesh).write_equation_systems ("truth.gmv",
                                        this->get_equation_systems());
  }

  // Get the X norm of the truth solution
  // Useful for normalizing our true error data
  if(!low_memory_mode)
  {
    inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
    matrix->vector_mult(*inner_product_storage_vector, *solution);
  }

  Number truth_X_norm = std::sqrt(inner_product_storage_vector->dot(*solution));

  STOP_LOG("truth_solve()", "RBSystem");

  return libmesh_real(truth_X_norm);
}

void RBSystem::set_eigen_system_name(const std::string& name)
{
  eigen_system_name = name;
}

void RBSystem::set_Nmax(unsigned int Nmax_in)
{
  this->Nmax = Nmax_in;
}

unsigned int RBSystem::get_n_F_EIM_systems() const
{
  return F_EIM_systems_vector.size();
}

unsigned int RBSystem::get_n_F_EIM_functions() const
{
  unsigned int count = 0;
  
  for(unsigned int i=0; i<F_EIM_systems_vector.size(); i++)
    count += F_EIM_systems_vector[i]->get_n_affine_functions();

  return count;
}

void RBSystem::attach_A_q(RBTheta* theta_q_a,
                          ElemAssembly* A_q_assembly)
{
  theta_q_a_vector.push_back(theta_q_a);
  A_q_assembly_vector.push_back(A_q_assembly);
}

void RBSystem::attach_F_q(RBTheta* theta_q_f,
                          ElemAssembly* F_q_assembly)
{
  theta_q_f_vector.push_back(theta_q_f);
  F_q_assembly_vector.push_back(F_q_assembly);
}

void RBSystem::attach_A_EIM_operators(RBEIMSystem* eim_system,
                                      ElemAssembly* EIM_assembly)
{
  A_EIM_systems_vector.push_back( eim_system );
  A_EIM_assembly_vector.push_back(EIM_assembly);
}

void RBSystem::attach_F_EIM_vectors(RBEIMSystem* eim_system,
                                    ElemAssembly* EIM_assembly)
{
  F_EIM_systems_vector.push_back( eim_system );
  F_EIM_assembly_vector.push_back(EIM_assembly);
}

void RBSystem::attach_output(std::vector<RBTheta*> theta_q_l,
                             std::vector<ElemAssembly*> output_assembly)
{
  // Make sure the input vectors are all the same size!
  if( theta_q_l.size() == output_assembly.size() )
  {
    theta_q_l_vector.push_back(theta_q_l);
    output_assembly_vector.push_back(output_assembly);
  }
  else
  {
    libMesh::out << "Error: The input vectors in attach_output must be the same size in attach_output"
                 << std::endl;
    libmesh_error();
  }
}

void RBSystem::attach_output(RBTheta* theta_q_l,
                             ElemAssembly* output_assembly)
{
  std::vector<RBTheta*> theta_l_vector(1); theta_l_vector[0] = theta_q_l;
  std::vector<ElemAssembly*> L_vector(1); L_vector[0] = output_assembly;

  attach_output(theta_l_vector, L_vector);
}

bool RBSystem::is_F_EIM_function(unsigned int q)
{
  libmesh_assert(q < get_Q_f());

  return (q >= theta_q_f_vector.size());
}

void RBSystem::attach_inner_prod_assembly(ElemAssembly* IP_assembly_in)
{
  inner_prod_assembly = IP_assembly_in;
}

void RBSystem::attach_constraint_assembly(ElemAssembly* constraint_assembly_in)
{
  constraint_assembly = constraint_assembly_in;
}

unsigned int RBSystem::get_Q_l(unsigned int index) const
{
  if(index >= get_n_outputs())
  {
    libMesh::err << "Error: We must have index < n_outputs in get_Q_l."
                 << std::endl;
    libmesh_error();
  }
  return theta_q_l_vector[index].size();
}

Number RBSystem::eval_theta_q_f(unsigned int q)
{
  if(q >= get_Q_f())
  {
    libMesh::err << "Error: We must have q < Q_f in eval_theta_q_f."
                 << std::endl;
    libmesh_error();
  }

  if( q < theta_q_f_vector.size() )
  {
    libmesh_assert(theta_q_f_vector[q] != NULL);
    return theta_q_f_vector[q]->evaluate( get_current_parameters() );
  }
  else
  {
    // Find out the RHS EIM system index and the associated function index that corresponds to q
    std::pair<unsigned int, unsigned int> F_EIM_indices = get_F_EIM_indices(q);

    RBEIMSystem& eim_system = *F_EIM_systems_vector[F_EIM_indices.first];
    eim_system.set_current_parameters(get_current_parameters());
    
    RBEIMEvaluation& eim_eval = libmesh_cast_ref<RBEIMEvaluation&>(*eim_system.rb_eval);
    eim_eval.RB_solve(eim_system.get_n_basis_functions());

    return eim_system.rb_eval->RB_solution(F_EIM_indices.second);
  }
}

std::pair<unsigned int, unsigned int> RBSystem::get_F_EIM_indices(unsigned int q)
{
  // Find out which EIM system q refers to
  int function_index = q - theta_q_f_vector.size();
  unsigned int system_index = 0;

  if( function_index != 0)
  {
    for(system_index = 0; system_index<F_EIM_systems_vector.size(); system_index++)
    {
      unsigned int increment = F_EIM_systems_vector[system_index]->get_n_affine_functions();
      if( (static_cast<int>(function_index - increment)) < 0)
        break;
    
      function_index -= increment;
    }
  }

  std::pair<unsigned int, unsigned int> EIM_indices(system_index, function_index);
  return EIM_indices;
}

Number RBSystem::eval_theta_q_l(unsigned int output_index, unsigned int q_l)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_Q_l(output_index)) )
  {
    libMesh::err << "Error: We must have output_index < n_outputs and "
                 << "q_l < get_Q_l(output_index) in eval_theta_q_l."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_l_vector[output_index][q_l] != NULL);

  return theta_q_l_vector[output_index][q_l]->evaluate( get_current_parameters() );
}

void RBSystem::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to load basis function."
                 << std::endl;
    libmesh_error();
  }

  libmesh_assert(i < get_n_basis_functions());

  *solution = get_basis_function(i);

  this->update();

  STOP_LOG("load_basis_function()", "RBSystem");
}

void RBSystem::enrich_RB_space()
{
  START_LOG("enrich_RB_space()", "RBSystem");

  NumericVector<Number>* new_bf = NumericVector<Number>::build().release();
  new_bf->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  *new_bf = *solution;

  // compute orthogonalization
  AutoPtr< NumericVector<Number> > proj_index = NumericVector<Number>::build();
  AutoPtr< NumericVector<Number> > proj_sum   = NumericVector<Number>::build();
  proj_index->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  proj_sum->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  if(low_memory_mode)
    assemble_inner_product_matrix(matrix);

  for(unsigned int index=0; index<get_n_basis_functions(); index++)
  {
    // invoke copy constructor for NumericVector
    *proj_index = get_basis_function(index);
    if(!low_memory_mode)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector,*proj_index);
    }
    else
    {
      matrix->vector_mult(*inner_product_storage_vector,*proj_index);
    }
    Number scalar = inner_product_storage_vector->dot(*new_bf);
    proj_sum->add( scalar, *proj_index);
  }
  new_bf->add(-1.,*proj_sum);

  // Normalize new_bf
  if(!low_memory_mode)
  {
    inner_product_matrix->vector_mult(*inner_product_storage_vector,*new_bf);
  }
  else
  {
    matrix->vector_mult(*inner_product_storage_vector,*new_bf);
  }
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
  rb_eval->basis_functions.push_back( new_bf );

  STOP_LOG("enrich_RB_space()", "RBSystem");
}

void RBSystem::update_system()
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

void RBSystem::recompute_all_residual_terms(bool compute_inner_products)
{
  // Use alternative solver for residual terms solves
  std::pair<std::string,std::string> orig_solver =
    this->set_alternative_solver(this->linear_solver);

  // Compute the basis independent terms
  Fq_representor_norms_computed = false;
  compute_Fq_representor_norms(compute_inner_products);

  // and all the basis dependent terms
  unsigned int saved_delta_N = delta_N;
  delta_N = get_n_basis_functions();
  
  update_residual_terms(compute_inner_products);

  delta_N = saved_delta_N;

  // Return to original solver
  this->reset_alternative_solver(this->linear_solver, orig_solver);
}

Real RBSystem::compute_max_error_bound()
{
  START_LOG("compute_max_error_bound()", "RBSystem");

  training_error_bounds.resize(this->get_local_n_training_samples());

  // keep track of the maximum error
  unsigned int max_err_index = 0;
  Real max_err = 0.;

  unsigned int first_index = get_first_local_training_index();
  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
  {
    // Load training parameter i, this is only loaded
    // locally since the RB solves are local.
    set_current_parameters( get_training_parameter(first_index+i) );

    training_error_bounds[i] = get_RB_error_bound();
//     libMesh::out << "Error bound at training index " << first_index+i << " is "
//                  << training_error_bounds[i] << std::endl;

    if(training_error_bounds[i] > max_err)
    {
      max_err_index = i;
      max_err = training_error_bounds[i];
    }
  }

  std::pair<unsigned int,Real> error_pair(first_index+max_err_index, max_err);
  get_global_max_error_pair(error_pair);

  // Now broadcast the parameter that produced the maximum error
  unsigned int root_id=0;
  if( (get_first_local_training_index() <= error_pair.first) &&
      (error_pair.first < get_last_local_training_index()) )
  {
    set_current_parameters( get_training_parameter(error_pair.first) );
    root_id = libMesh::processor_id();
  }

  Parallel::sum(root_id); // root_id is only non-zero on one processor
  broadcast_current_parameters(root_id);

  STOP_LOG("compute_max_error_bound()", "RBSystem");

  return error_pair.second;
}

Real RBSystem::get_SCM_lower_bound()
{
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
  // Get the SCM lower bound from eigen_system
  EquationSystems& es = this->get_equation_systems();
  RBSCMSystem& eigen_system = es.get_system<RBSCMSystem>(eigen_system_name);

  eigen_system.set_current_parameters( this->get_current_parameters() );
  return eigen_system.get_SCM_LB();
#else
  libMesh::out << "SLEPc and GLPK must be installed for SCM functions to work." << std::endl;
  libmesh_error();
#endif // defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
}

Real RBSystem::get_SCM_upper_bound()
{
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
  // Get the SCM lower bound from eigen_system
  EquationSystems& es = this->get_equation_systems();
  RBSCMSystem& eigen_system = es.get_system<RBSCMSystem>(eigen_system_name);

  eigen_system.set_current_parameters( this->get_current_parameters() );
  return eigen_system.get_SCM_UB();
#else
  libMesh::out << "SLEPc and GLPK must be installed for SCM functions to work." << std::endl;
  libmesh_error();
#endif // defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
}

Real RBSystem::residual_scaling_denom(Real alpha_LB)
{
  // Here we implement the residual scaling for a coercive
  // problem.
  return alpha_LB;
}

void RBSystem::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "RBSystem");

  unsigned int RB_size = get_n_basis_functions();

  AutoPtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  temp->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      rb_eval->RB_F_q_vector[q_f](i) = get_F_q(q_f)->dot(get_basis_function(i));
    }
  }

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
  {
    for(unsigned int n=0; n<get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        rb_eval->RB_output_vectors[n][q_l](i) = get_output_vector(n,q_l)->dot(get_basis_function(i));
      }

    for(unsigned int j=0; j<RB_size; j++)
    {
      Number value = 0.;
      
      if(compute_RB_inner_product)
      {
        // Compute reduced inner_product_matrix
        temp->zero();
        if(!low_memory_mode)
        {
          inner_product_matrix->vector_mult(*temp, get_basis_function(j));
        }
        else
        {
          assemble_inner_product_matrix(matrix);
          matrix->vector_mult(*temp, get_basis_function(j));
        }

        value = get_basis_function(i).dot(*temp);
        rb_eval->RB_inner_product_matrix(i,j) = value;
        if(i!=j)
        {
          // The inner product matrix is assumed
          // to be symmetric
          rb_eval->RB_inner_product_matrix(j,i) = value;
        }
      }

      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        // Compute reduced A_q matrix
        temp->zero();
        if(!low_memory_mode)
        {
          get_A_q(q_a)->vector_mult(*temp, get_basis_function(j));
        }
        else
        {
          assemble_Aq_matrix(q_a,matrix);
          matrix->vector_mult(*temp, get_basis_function(j));
        }

        value = (*temp).dot(get_basis_function(i));
        rb_eval->RB_A_q_vector[q_a](i,j) = value;

        if(i!=j)
        {
          temp->zero();
          if(!low_memory_mode)
          {
            get_A_q(q_a)->vector_mult(*temp, get_basis_function(i));
          }
          else
          {
            // matrix should still hold affine matrix q_a
            matrix->vector_mult(*temp, get_basis_function(i));
          }

          value = (*temp).dot(get_basis_function(j));
          rb_eval->RB_A_q_vector[q_a](j,i) = value;
        }
      }
    }
  }

  STOP_LOG("update_RB_system_matrices()", "RBSystem");
}


void RBSystem::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "RBSystem");

  unsigned int RB_size = get_n_basis_functions();

  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }

  if(low_memory_mode)
  {
    assemble_inner_product_matrix(matrix);
    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, matrix, NULL);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->same_preconditioner = false;
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      // Initialize the vector in which we'll store the representor
      if(!rb_eval->A_q_representor[q_a][i])
      {
        rb_eval->A_q_representor[q_a][i] = (NumericVector<Number>::build().release());
        rb_eval->A_q_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }

      libmesh_assert(rb_eval->A_q_representor[q_a][i]->size()       == this->n_dofs()       && 
                     rb_eval->A_q_representor[q_a][i]->local_size() == this->n_local_dofs() );

      rhs->zero();
      if(!low_memory_mode)
      {
        get_A_q(q_a)->vector_mult(*rhs, get_basis_function(i));
      }
      else
      {
        assemble_scaled_matvec(1.,
                               A_q_assembly_vector[q_a],
                               *rhs,
                               get_basis_function(i));
      }
      rhs->scale(-1.);
      zero_dirichlet_dofs_on_rhs();

      solution->zero();
      if (!is_quiet())
	    {
        libMesh::out << "Starting solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBSystem::update_residual_terms() at "
                     << Utility::get_timestamp() << std::endl;
	    }

      solve();

      if (!is_quiet())
	    {
        libMesh::out << "Finished solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBSystem::update_residual_terms() at "
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
      *rb_eval->A_q_representor[q_a][i] = *solution;


      if(reuse_preconditioner)
      {
        // set this flag again in case we didn't do any F solves
        linear_solver->same_preconditioner = true;
      }
    }
  }

  if(reuse_preconditioner)
  {
    linear_solver->same_preconditioner = false;
  }

  // Now compute and store the inner products (if requested)
  if (compute_inner_products)
    {
      if(low_memory_mode && constrained_problem)
	assemble_inner_product_matrix(matrix);

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	{
	  if(!low_memory_mode)
	    {
	      inner_product_matrix->vector_mult(*inner_product_storage_vector,*F_q_representor[q_f]);
	    }
	  else
	    {
	      matrix->vector_mult(*inner_product_storage_vector,*F_q_representor[q_f]);
	    }

	  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
	    {
	      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
		{
		  rb_eval->Fq_Aq_representor_norms[q_f][q_a][i] =
		    inner_product_storage_vector->dot(*rb_eval->A_q_representor[q_a][i]);
		}
	    }
	}

      unsigned int q=0;
      for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
	{
	  for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
	    {
	      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
		{
		  for(unsigned int j=0; j<RB_size; j++)
		    {
		      if(!low_memory_mode)
			{
			  inner_product_matrix->vector_mult(*inner_product_storage_vector, *rb_eval->A_q_representor[q_a2][j]);
			}
		      else
			{
			  matrix->vector_mult(*inner_product_storage_vector, *rb_eval->A_q_representor[q_a2][j]);
			}
		      rb_eval->Aq_Aq_representor_norms[q][i][j] = inner_product_storage_vector->dot(*rb_eval->A_q_representor[q_a1][i]);

		      if(i != j)
			{
			  if(!low_memory_mode)
			    {
			      inner_product_matrix->vector_mult(*inner_product_storage_vector, *rb_eval->A_q_representor[q_a2][i]);
			    }
			  else
			    {
			      matrix->vector_mult(*inner_product_storage_vector, *rb_eval->A_q_representor[q_a2][i]);
			    }
			  rb_eval->Aq_Aq_representor_norms[q][j][i] = inner_product_storage_vector->dot(*rb_eval->A_q_representor[q_a1][j]);
			}
		    }
		}
	      q++;
	    }
	}
    } // end if (compute_inner_products)

  STOP_LOG("update_residual_terms()", "RBSystem");
}

void RBSystem::assemble_matrix_for_output_dual_solves()
{
  // By default we use the inner product matrix for steady problems
  
  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
  }
  else
  {
    assemble_inner_product_matrix(matrix);
  }
}

void RBSystem::compute_output_dual_norms()
{
  // short-circuit if we've already computed the output dual norms
  if(output_dual_norms_computed)
  {
    return;
  }

  // Short circuit if we don't have any outputs
  if( get_n_outputs() == 0 )
  {
    output_dual_norms_computed = true;
    return;
  }

  START_LOG("compute_output_dual_norms()", "RBSystem");
  
  libMesh::out << "Compute output dual norms" << std::endl;
  
  // Note: the solves in this function employ a single system matrix and multiple
  // right-hand sides, so we may get better performance using a different
  // preconditioner, or even a direct solver.
  std::pair<std::string,std::string> orig_solver =
    this->set_alternative_solver(this->linear_solver);

  // Find out the largest value of Q_l
  unsigned int max_Q_l = 0;
  for(unsigned int n=0; n<get_n_outputs(); n++)
    max_Q_l = (get_Q_l(n) > max_Q_l) ? get_Q_l(n) : max_Q_l;

  std::vector< NumericVector<Number>* > L_q_representor(max_Q_l);
  for(unsigned int q=0; q<max_Q_l; q++)
  {
    L_q_representor[q] = (NumericVector<Number>::build().release());
    L_q_representor[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->same_preconditioner = false;
  }
  
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    // If constrained_problem, we need to reassemble to add the constraint part back in
    if( (n==0) || constrained_problem)
    {
      assemble_matrix_for_output_dual_solves();
      
      if(constrained_problem)
      {
        if(!low_memory_mode)
        {
          matrix->add(1., *constraint_matrix);
        }
        else
        {
          add_scaled_matrix_and_vector(1., constraint_assembly, matrix, NULL);
        }
      }
    }

    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      rhs->zero();
            rhs->add(1., *get_output_vector(n,q_l));
      zero_dirichlet_dofs_on_rhs();

      solution->zero();

      if (!is_quiet())
        libMesh::out << "Starting solve n=" << n << ", q_l=" << q_l
               << " in RBSystem::compute_output_dual_norms() at "
               << Utility::get_timestamp() << std::endl;

      solve();

      if (!is_quiet())
        {
          libMesh::out << "Finished solve n=" << n << ", q_l=" << q_l
                       << " in RBSystem::compute_output_dual_norms() at "
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
          linear_solver->same_preconditioner = true;
        }
    }

    // Get rid of the constraint part of the matrix before computing inner products
    if(constrained_problem)
      assemble_matrix_for_output_dual_solves();

    unsigned int q=0;
    for(unsigned int q_l1=0; q_l1<get_Q_l(n); q_l1++)
    {
      matrix->vector_mult(*inner_product_storage_vector, *L_q_representor[q_l1]);

      for(unsigned int q_l2=q_l1; q_l2<get_Q_l(n); q_l2++)
      {
        output_dual_norms[n][q] = L_q_representor[q_l2]->dot(*inner_product_storage_vector);
        libMesh::out << "output_dual_norms[" << n << "][" << q << "] = " << output_dual_norms[n][q] << std::endl;
        
        q++;
      }
    }
  }

  // reset same_preconditioner to false once all solves are finished
  if(reuse_preconditioner)
  {
    linear_solver->same_preconditioner = false;
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
  
  output_dual_norms_computed = true;

  // Change the preconditioner, Krylov solver back to their original
  // value.  Note: does nothing if RBBase::alternative_solver ==
  // "unchanged".
  this->reset_alternative_solver(this->linear_solver, orig_solver);

  STOP_LOG("compute_output_dual_norms()", "RBSystem");
}

void RBSystem::compute_Fq_representor_norms(bool compute_inner_products)
{
  // Short-circuit if we've already computed the Fq_representors
  if(Fq_representor_norms_computed)
  {
    return;
  }

  START_LOG("compute_Fq_representor_norms()", "RBSystem");

  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }

  if(low_memory_mode)
  {
    assemble_inner_product_matrix(matrix);
    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, matrix, NULL);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->same_preconditioner = false;
  }

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    if(!F_q_representor[q_f])
    {
      F_q_representor[q_f] = (NumericVector<Number>::build().release());
      F_q_representor[q_f]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }

    libmesh_assert(F_q_representor[q_f]->size()       == this->n_dofs()       && 
                   F_q_representor[q_f]->local_size() == this->n_local_dofs() );

    rhs->zero();
    rhs->add(1., *get_F_q(q_f));
    zero_dirichlet_dofs_on_rhs();

    solution->zero();

    if (!is_quiet())
      libMesh::out << "Starting solve q_f=" << q_f
		   << " in RBSystem::update_residual_terms() at "
		   << Utility::get_timestamp() << std::endl;

    solve();

    if (!is_quiet())
    {
      libMesh::out << "Finished solve q_f=" << q_f
		   << " in RBSystem::update_residual_terms() at "
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
    *F_q_representor[q_f] = *solution;

    if(reuse_preconditioner)
    {
      // After we do a solve, tell PETSc we want to reuse the preconditioner
      // since the system matrix is not changing.
      linear_solver->same_preconditioner = true;
    }
  }

  // Reset the same_preconditioner flag
  if(reuse_preconditioner)
  {
    linear_solver->same_preconditioner = false;
  }

  if (compute_inner_products)
  {
    unsigned int q=0;
    if(low_memory_mode && constrained_problem)
      assemble_inner_product_matrix(matrix);

    for(unsigned int q_f1=0; q_f1<get_Q_f(); q_f1++)
    {
      if(!low_memory_mode)
      {
        inner_product_matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f1]);
      }
      else
      {
        matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f1]);
      }

      for(unsigned int q_f2=q_f1; q_f2<get_Q_f(); q_f2++)
      {
        Fq_representor_norms[q] = inner_product_storage_vector->dot(*F_q_representor[q_f2]);

        q++;
      }
    }
  } // end if (compute_inner_products)
  
  Fq_representor_norms_computed = true;

  STOP_LOG("compute_Fq_representor_norms()", "RBSystem");
}

void RBSystem::load_RB_solution()
{
  START_LOG("load_RB_solution()", "RBSystem");

  if(!initialize_mesh_dependent_data)
  {
    libMesh::err << "Error: We must initialize the mesh dependent "
                 << "data structures in order to load RB solution."
                 << std::endl;
    libmesh_error();
  }

  solution->zero();

  if(rb_eval->RB_solution.size() > get_n_basis_functions())
  {
    libMesh::err << "ERROR: System contains " << get_n_basis_functions() << " basis functions."
                 << " RB_solution vector constains " << rb_eval->RB_solution.size() << " entries."
                 << " RB_solution in RBSystem::load_RB_solution is too long!" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<rb_eval->RB_solution.size(); i++)
  {
    solution->add(rb_eval->RB_solution(i), get_basis_function(i));
  }

  update();

  STOP_LOG("load_RB_solution()", "RBSystem");
}

// The slow (but simple, non-error prone) way to compute the residual dual norm
// Useful for error checking
//Real RBSystem::compute_residual_dual_norm(const unsigned int N)
//{
//  START_LOG("compute_residual_dual_norm()", "RBSystem");
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
//     RB_sol->add(RB_solution(i), get_basis_function(i));
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
//   if(!low_memory_mode)
//   {
//     inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//   }
//   else
//   {
//     assemble_inner_product_matrix(matrix);
//     matrix->vector_mult(*inner_product_storage_vector, *solution);
//   }
//
//   Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);
//
//  STOP_LOG("compute_residual_dual_norm()", "RBSystem");
//
//  return std::sqrt( libmesh_real(slow_residual_norm_sq) );
//}

SparseMatrix<Number>* RBSystem::get_inner_product_matrix()
{
  if(low_memory_mode)
  {
    libMesh::err << "Error: The inner-product matrix is not stored in low-memory mode." << std::endl;
    libmesh_error();
  }

  return inner_product_matrix.get();
}

SparseMatrix<Number>* RBSystem::get_non_dirichlet_inner_product_matrix()
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true "
                 << "to access non_dirichlet_inner_product_matrix." << std::endl;
    libmesh_error();
  }
  if(low_memory_mode)
  {
    libMesh::err << "Error: The non-Dirichlet inner-product matrix is not stored in low-memory mode." << std::endl;
    libmesh_error();
  }

  return non_dirichlet_inner_product_matrix.get();
}

SparseMatrix<Number>* RBSystem::get_A_q(unsigned int q)
{
  if(low_memory_mode)
  {
    libMesh::err << "Error: The affine matrices are not stored in low-memory mode." << std::endl;
    libmesh_error();
  }

  if(q >= get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in get_A_q."
                 << std::endl;
    libmesh_error();
  }

  return A_q_vector[q];
}

SparseMatrix<Number>* RBSystem::get_non_dirichlet_A_q(unsigned int q)
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_A_q." << std::endl;
    libmesh_error();
  }

  if(low_memory_mode)
  {
    libMesh::err << "Error: The affine matrices are not stored in low-memory mode." << std::endl;
    libmesh_error();
  }

  if(q >= get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in get_A_q."
                 << std::endl;
    libmesh_error();
  }

  return non_dirichlet_A_q_vector[q];
}

RBEIMSystem& RBSystem::get_A_EIM_system(unsigned int index)
{
  if(index >= A_EIM_systems_vector.size())
  {
    libMesh::err << "Error: We must have index < get_n_A_EIM_systems() in get_A_EIM_system."
                 << std::endl;
    libmesh_error();
  }
  
  return *A_EIM_systems_vector[index];
}

RBEIMSystem& RBSystem::get_F_EIM_system(unsigned int index)
{
  if(index >= F_EIM_systems_vector.size())
  {
    libMesh::err << "Error: We must have index < get_n_F_EIM_systems() in get_F_EIM_system."
                 << std::endl;
    libmesh_error();
  }
  
  return *F_EIM_systems_vector[index];
}

std::vector<Number> RBSystem::evaluate_current_EIM_function(Elem& element, const std::vector<Point>& qpoints)
{
  return current_EIM_system->evaluate_current_affine_function(element, qpoints);
}

NumericVector<Number>* RBSystem::get_F_q(unsigned int q)
{
  if(q >= get_Q_f())
  {
    libMesh::err << "Error: We must have q < Q_f in get_F_q."
                 << std::endl;
    libmesh_error();
  }

  return F_q_vector[q];
}

NumericVector<Number>* RBSystem::get_non_dirichlet_F_q(unsigned int q)
{
  if(!store_non_dirichlet_operators)
  {
    libMesh::err << "Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_F_q." << std::endl;
    libmesh_error();
  }

  if(q >= get_Q_f())
  {
    libMesh::err << "Error: We must have q < Q_f in get_F_q."
                 << std::endl;
    libmesh_error();
  }

  return non_dirichlet_F_q_vector[q];
}

NumericVector<Number>* RBSystem::get_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_n_outputs()) || (q_l >= get_Q_l(n)) )
  {
    libMesh::err << "Error: We must have n < n_outputs and "
                 << "q_l < get_Q_l(n) in get_output_vector."
                 << std::endl;
    libmesh_error();
  }

  return outputs_vector[n][q_l];
}

void RBSystem::zero_dirichlet_dofs_on_rhs()
{
  this->zero_dirichlet_dofs_on_vector(*rhs);
}

void RBSystem::zero_dirichlet_dofs_on_vector(NumericVector<Number>& temp)
{
  START_LOG("zero_dirichlet_dofs_on_vector()", "RBSystem");

  std::set<unsigned int>::iterator iter     = global_dirichlet_dofs_set.begin();
  std::set<unsigned int>::iterator iter_end = global_dirichlet_dofs_set.end();

  DofMap& dof_map = this->get_dof_map();
  for ( ; iter != iter_end; ++iter)
  {
    unsigned int index = *iter;
    if( (dof_map.first_dof() <= index) && (index < dof_map.end_dof()) )
    {
      temp.set(index, 0.);
    }
  }
  temp.close();

  STOP_LOG("zero_dirichlet_dofs_on_vector()", "RBSystem");
}




void RBSystem::write_offline_data_to_files(const std::string& directory_name,
                                           const RBDataIO io_flag)
{
  START_LOG("write_offline_data_to_files()", "RBSystem");

  if( (io_flag == ALL_DATA) || (io_flag == BASIS_DEPENDENT) )
  {
    rb_eval->write_offline_data_to_files(directory_name);
  }
  
  // return here if we only want basis dependent data
  if( io_flag == BASIS_DEPENDENT )
  {
    STOP_LOG("write_offline_data_to_files()", "RBSystem");
    return;
  }

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Make a directory to store all the data files
    mkdir(directory_name.c_str(), 0777);

    // Write out F_q representor norm data
    std::ofstream RB_Fq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_norms.dat";
      RB_Fq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_norms_out.good() )
    {
      libMesh::err << "Error opening Fq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_norms_out.precision(precision_level);
    unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
    for(unsigned int i=0; i<Q_f_hat; i++)
    {
      RB_Fq_norms_out << std::scientific << Fq_representor_norms[i] << " ";
    }
    RB_Fq_norms_out.close();

    // Write out output data
    for(unsigned int n=0; n<get_n_outputs(); n++)
    {
      std::ofstream output_dual_norms_out;
      {
        OStringStream file_name;
        file_name << directory_name << "/output_";
        OSSRealzeroright(file_name,3,0,n);
        file_name << "_dual_norms.dat";
        output_dual_norms_out.open(file_name.str().c_str());
      }
      if ( !output_dual_norms_out.good() )
      {
        libMesh::err << "Error opening output " << n << " dual norms file" << std::endl;
        libmesh_error();
      }
      output_dual_norms_out.precision(precision_level);
      
      unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
      for(unsigned int q=0; q<Q_l_hat; q++)
      {
        output_dual_norms_out << std::scientific << output_dual_norms[n][q] << " ";
      }
      output_dual_norms_out.close();
    }
  }

  // Write out residual representors if requested
  if (store_representors)
    {
      // Write out F_q_representors.  These are useful to have when restarting,
      // so you don't have to recompute them all over again.  There should be
      // Q_f of these.
      if (!is_quiet())
	libMesh::out << "Writing out the F_q_representors..." << std::endl;

      std::ostringstream file_name;
      const std::string residual_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");
      struct stat stat_info;

      // Residual representors written out to their own separate directory
      std::string residual_representors_dir = "residual_representors";
      if ( libMesh::processor_id() == 0)
	if ( mkdir(residual_representors_dir.c_str(), 0755) != 0)
	  libMesh::out << "Skipping creating residual_representors directory: " << strerror(errno) << std::endl;

      for (unsigned int i=0; i<F_q_representor.size(); ++i)
	{
	  if (F_q_representor[i] != NULL)
	    {

	      file_name.str(""); // reset filename
	      file_name << residual_representors_dir << "/F_q_representor" << i << residual_representor_suffix;

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
		  // *solution = *(F_q_representor[i]);
		  // std::swap doesn't work on pointers
		  //std::swap(solution.get(), F_q_representor[i]);
		  F_q_representor[i]->swap(*solution);

		  Xdr fqr_data(file_name.str(),
			       write_binary_residual_representors ? ENCODE : WRITE);

		  write_serialized_data(fqr_data, false);

		  // Synchronize before moving on
		  Parallel::barrier();

		  // Swap back.
		  F_q_representor[i]->swap(*solution);

		  // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
		  // for the system call, be sure to do it only on one processor, etc.
		}
	    }
	}
    } // end if (store_representors)

  STOP_LOG("write_offline_data_to_files()", "RBSystem");
}



void RBSystem::read_offline_data_from_files(const std::string& directory_name,
                                            const RBDataIO io_flag)
{
  START_LOG("read_offline_data_from_files()", "RBSystem");

  if( (io_flag == ALL_DATA) || (io_flag == BASIS_DEPENDENT) )
  {
    rb_eval->read_offline_data_from_files(directory_name);
  }

  // Return here if we only want basis dependent data
  if( io_flag == BASIS_DEPENDENT )
  {
    STOP_LOG("read_offline_data_from_files()", "RBSystem");
    return;
  }

  // Next read in F_q representor norm data
  std::ifstream RB_Fq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Fq_norms.dat";
    RB_Fq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Fq_norms_in.good() )
  {
    libMesh::err << "Error opening Fq_norms.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
  for(unsigned int i=0; i<Q_f_hat; i++)
  {
    RB_Fq_norms_in >> Fq_representor_norms[i];
  }
  RB_Fq_norms_in.close();


  // Read in output data
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    std::ifstream output_dual_norms_in;
    {
      OStringStream file_name;
      file_name << directory_name << "/output_";
      OSSRealzeroright(file_name,3,0,n);
      file_name << "_dual_norms.dat";
      output_dual_norms_in.open(file_name.str().c_str());
    }
    if ( !output_dual_norms_in.good() )
    {
      libMesh::err << "Error opening input " << n << " dual norms file" << std::endl;
      libmesh_error();
    }
    
    unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
    for(unsigned int q=0; q<Q_l_hat; q++)
    {
      output_dual_norms_in >> output_dual_norms[n][q];
    }
    output_dual_norms_in.close();
    output_dual_norms_computed = true;
  }

  // Read in the representor vectors if requested
  if (store_representors)
  {
    libMesh::out << "Reading in the F_q_representors..." << std::endl;

    const std::string residual_representors_dir = "residual_representors";
    const std::string residual_representor_suffix = (read_binary_residual_representors ? ".xdr" : ".dat");
    std::ostringstream file_name;
    struct stat stat_info;

    // Read in the F_q_representors.  There should be Q_f of these.  FIXME:
    // should we be worried about leaks here?
    for (unsigned int i=0; i<F_q_representor.size(); ++i)
    {
      if (F_q_representor[i] != NULL)
      {
        libMesh::out << "Error, must delete existing F_q_representor before reading in from file."
                     << std::endl;
        libmesh_error();
      }
    }

    for (unsigned int i=0; i<F_q_representor.size(); i++)
    {
      file_name.str(""); // reset filename
      file_name << residual_representors_dir
		<< "/F_q_representor" << i << residual_representor_suffix;

      // On processor zero check to be sure the file exists
      if (libMesh::processor_id() == 0)
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

      F_q_representor[i] = NumericVector<Number>::build().release();
      F_q_representor[i]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

      // No need to copy, just swap
      // *F_q_representor[i] = *solution;
      F_q_representor[i]->swap(*solution);
    }

    // Alert the update_residual_terms() function that we don't need to recompute
    // the F_q_representors as we have already read them in from file!
    Fq_representor_norms_computed = true;

  } // end if (store_representors)

  STOP_LOG("read_offline_data_from_files()", "RBSystem");
}


} // namespace libMesh



