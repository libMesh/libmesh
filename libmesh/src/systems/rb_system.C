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
#include "parallel.h"

#include "fem_context.h"
#include "rb_system.h"
#include "rb_scm_system.h"

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
    inner_product_matrix(SparseMatrix<Number>::build()),
    constraint_matrix(SparseMatrix<Number>::build()),
    constrained_problem(false),
    store_basis_functions(false),
    store_representors(false),
    low_memory_mode(false),
    reuse_preconditioner(true),
    return_rel_error_bound(true),
    write_data_during_training(false),
    impose_internal_dirichlet_BCs(false),
    impose_internal_fluxes(false),
    parameters_filename(""),
    extra_quadrature_order(0),
    enforce_constraints_exactly(false),
    Nmax(0),
    delta_N(1),
    write_binary_basis_functions(true),
    read_binary_basis_functions(true),
    write_binary_residual_representors(true),
    read_binary_residual_representors(true),
    quiet(true),
    eigen_system_name(""),
    inner_prod_assembly(NULL),
    constraint_assembly(NULL),
    training_tolerance(-1.),
    update_residual_terms_called(false),
    _dirichlet_list_init(NULL),
    initial_Nmax(0),
    RB_system_initialized(false)
{
  RB_solution.resize(0);

  Fq_representor_norms.clear();
  Fq_Aq_representor_norms.clear();
  Aq_Aq_representor_norms.clear();

  // Clear the theta and assembly vectors so that we can push_back
  theta_q_a_vector.clear();
  A_q_intrr_assembly_vector.clear();
  A_q_bndry_assembly_vector.clear();

  theta_q_f_vector.clear();
  F_q_intrr_assembly_vector.clear();
  F_q_bndry_assembly_vector.clear();

  theta_q_l_vector.clear();
  output_intrr_assembly_vector.clear();
  output_bndry_assembly_vector.clear();

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
  Parent::clear();

  for(unsigned int q=0; q<get_Q_a(); q++)
  {
    if(A_q_vector[q])
    {
      delete A_q_vector[q];
      A_q_vector[q] = NULL;
    }
  }

  for(unsigned int q=0; q<get_Q_f(); q++)
  {
    if(F_q_vector[q])
    {
      delete F_q_vector[q];
      F_q_vector[q] = NULL;
    }
  }

  for(unsigned int i=0; i<get_n_outputs(); i++)
    for(unsigned int q_l=0; q_l<get_Q_l(i); q_l++)
      if(outputs_vector[i][q_l])
      {
        delete outputs_vector[i][q_l];
        outputs_vector[i][q_l] = NULL;
      }

  // Clear the basis functions and the
  // basis-function-dependent data using
  // the non-virtual helper function
  this->clear_basis_helper();
}

void RBSystem::clear_basis_function_dependent_data()
{
  update_residual_terms_called = false;

  // Clear the Greedy param list
  for(unsigned int i=0; i<greedy_param_list.size(); i++)
    greedy_param_list[i].clear();
  greedy_param_list.clear();

  // Call non-virtual helper class to clear the
  // basis related data
  clear_basis_helper();
}

void RBSystem::clear_basis_helper()
{
  // Clear the basis functions
  for(unsigned int i=0; i<basis_functions.size(); i++)
    {
      if (basis_functions[i])
	{
	  basis_functions[i]->clear();
	  delete basis_functions[i];
	  basis_functions[i] = NULL;
	}
    }
  basis_functions.resize(0);

  // Also delete the representors
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    if(F_q_representor[q_f])
    {
      delete F_q_representor[q_f];
      F_q_representor[q_f] = NULL;
    }
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int i=0; i<A_q_representor[q_a].size(); i++)
    {
      if(A_q_representor[q_a][i])
      {
        delete A_q_representor[q_a][i];
        A_q_representor[q_a][i] = NULL;
      }
    }
  }

}

std::string RBSystem::system_type () const
{
  return "RBSystem";
}

void RBSystem::init_data ()
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

  // Tell the system whether or not to initialize \calN dependent data
  // structures.
  initialize_calN_dependent_data = infile("initialize_calN_dependent_data",
                                          initialize_calN_dependent_data);

  // Set boolean which turns on/off storing the representor residuals
  impose_internal_dirichlet_BCs = infile("impose_internal_dirichlet_BCs",
                                          impose_internal_dirichlet_BCs);

  // Set boolean which turns on/off storing the representor residuals
  impose_internal_fluxes = infile("impose_internal_fluxes",
                                   impose_internal_fluxes);

  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  training_parameters_random_seed = infile("training_parameters_random_seed",
					   training_parameters_random_seed);
  
  // Set quiet mode
  const bool quiet_in = infile("quiet", quiet);
  set_quiet(quiet_in);

  // Throw an error if we try to not initialize calN-dependent data
  // when we also want to read in the basis functions.
  if(!initialize_calN_dependent_data && store_basis_functions)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures if we want to read in basis "
              << "functions."
              << std::endl;
    libmesh_error();
  }

  // Initialize RB parameters
  const unsigned int Nmax_in = infile("Nmax", Nmax);
  set_Nmax(Nmax_in);
  initial_Nmax = get_Nmax();
  
  const Real training_tolerance_in = infile("training_tolerance",
                                            training_tolerance);
  set_training_tolerance(training_tolerance_in);


  std::vector<Real> mu_min_vector(n_parameters);
  std::vector<Real> mu_max_vector(n_parameters);
  std::vector<bool> log_scaling(n_parameters);
  std::vector<Real> init_mu_vector(n_parameters);
  for(unsigned int i=0; i<n_parameters; i++)
  {
    // Read vector-based mu_min values.
    mu_min_vector[i] = infile("mu_min", mu_min_vector[i], i);

    // Read vector-based mu_max values.
    mu_max_vector[i] = infile("mu_max", mu_max_vector[i], i);

    // Read vector-based log scaling values.  Note the intermediate conversion to
    // int... this implies log_scaling = '1 1 1...' in the input file.
    log_scaling[i] = static_cast<bool>(infile("log_scaling", static_cast<int>(log_scaling[i]), i));

    // Read vector-based init_mu values.
    init_mu_vector[i] = infile("init_mu", init_mu_vector[i], i);
  }

  initialize_training_parameters(mu_min_vector,
                                 mu_max_vector,
                                 n_training_samples_mu,
                                 log_scaling,
                                 deterministic_training);   // use deterministic parameters



  // Set the initial parameter value
  set_current_parameters(init_mu_vector);

  std::cout << std::endl << "RBSystem parameters:" << std::endl;
  std::cout << "system name: " << this->name() << std::endl;
  std::cout << "constrained_problem: " << constrained_problem << std::endl;
  std::cout << "Nmax: " << Nmax << std::endl;
  if(training_tolerance > 0.)
    std::cout << "Basis training error tolerance: " << get_training_tolerance() << std::endl;
  std::cout << "Q_a: " << get_Q_a() << std::endl;
  std::cout << "Q_f: " << get_Q_f() << std::endl;
  std::cout << "n_outputs: " << get_n_outputs() << std::endl;
  for(unsigned int n=0; n<get_n_outputs(); n++)
    std::cout << "output " << n << ", Q_l = " << get_Q_l(n) << std::endl;
  for(unsigned int i=0; i<n_parameters; i++)
  {
    std::cout <<   "Parameter " << i
              << ": Min = " << get_parameter_min(i)
              << ", Max = " << get_parameter_max(i)
              << ", log scaling = " << log_scaling[i] << std::endl;
  }
  std::cout << "n_training_samples: " << get_n_training_samples() << std::endl;
  std::cout << "using deterministic training samples? " << deterministic_training << std::endl;
  std::cout << "store/load basis functions? " << store_basis_functions << std::endl;
  if(store_basis_functions)
  {
    std::cout << "  write out basis functions in binary format? "
              << write_binary_basis_functions << std::endl;
    std::cout << "  read in basis functions in binary format? "
              << read_binary_basis_functions << std::endl;
  }
  std::cout << "store/load residual representors? " << store_representors << std::endl;
  if(store_representors)
  {
    std::cout << "  write out residual representors in binary format? "
              << write_binary_residual_representors << std::endl;
    std::cout << "  read in residual representors in binary format? "
              << read_binary_residual_representors << std::endl;
  }
  std::cout << "low-memory mode? " << low_memory_mode << std::endl;
  std::cout << "reuse preconditioner? " << reuse_preconditioner << std::endl;
  std::cout << "return a relative error bound from RB_solve? " << return_rel_error_bound << std::endl;
  std::cout << "write out data during basis training? " << write_data_during_training << std::endl;
  std::cout << "initializing calN-dependent data structures? "
            << initialize_calN_dependent_data << std::endl;
  std::cout << "impose internal Dirichlet BCs? " << impose_internal_dirichlet_BCs << std::endl;
  std::cout << "impose internal fluxes? " << impose_internal_fluxes << std::endl;
  std::cout << "quiet mode? " << quiet << std::endl;
  std::cout << "initial parameter: ";
  for(unsigned int i=0; i<n_parameters; i++)
  {
    std::cout << "mu[" << i << "] = " << get_current_parameters()[i];
    if(i < (n_parameters-1))
      std::cout << ", ";
    else
      std::cout << std::endl;
  }
  std::cout << std::endl;

  // We need Nmax to be initialized
  libmesh_assert(Nmax > 0);

  // Now that input parameters (e.g. Q_a) have been read in,
  // call the Parent's initialization routine.
  Parent::init_data();

  // Resize vectors for storing calN-dependent data but only
  // initialize if initialize_calN_dependent_data == true
  A_q_vector.resize(get_Q_a());
  F_q_vector.resize(get_Q_f());
  F_q_representor.resize(get_Q_f());
  A_q_representor.resize(get_Q_a());
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    A_q_representor[q_a].resize(Nmax);
  }
  outputs_vector.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
    outputs_vector[n].resize( get_Q_l(n) );

  if(initialize_calN_dependent_data)
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
    }

    // Initialize the vectors even if we are in low-memory mode
    for(unsigned int q=0; q<get_Q_f(); q++)
    {
      // Initialize the memory for the vectors
      F_q_vector[q] = NumericVector<Number>::build().release();
      F_q_vector[q]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
    }

    for(unsigned int n=0; n<get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        // Initialize the memory for the truth output vectors
        outputs_vector[n][q_l] = (NumericVector<Number>::build().release());
        outputs_vector[n][q_l]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
      }
  }

  // Now allocate the N (i.e. RB) dependent data structures

  // Allocate dense matrices for RB solves
  RB_A_q_vector.resize(get_Q_a());

  for(unsigned int q=0; q<get_Q_a(); q++)
  {
    // Initialize the memory for the RB matrices
    RB_A_q_vector[q].resize(Nmax,Nmax);
  }

  RB_F_q_vector.resize(get_Q_f());

  for(unsigned int q=0; q<get_Q_f(); q++)
  {
    // Initialize the memory for the RB vectors
    RB_F_q_vector[q].resize(Nmax);
  }

  // Initialize vectors for the norms of the representors
  unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
  Fq_representor_norms.resize(Q_f_hat);

  Fq_Aq_representor_norms.resize(get_Q_f());
  for(unsigned int i=0; i<get_Q_f(); i++)
  {
    Fq_Aq_representor_norms[i].resize(get_Q_a());
    for(unsigned int j=0; j<get_Q_a(); j++)
    {
      Fq_Aq_representor_norms[i][j].resize(Nmax);
    }
  }

  unsigned int Q_a_hat = get_Q_a()*(get_Q_a()+1)/2;
  Aq_Aq_representor_norms.resize(Q_a_hat);
  for(unsigned int i=0; i<Q_a_hat; i++)
  {
    Aq_Aq_representor_norms[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      Aq_Aq_representor_norms[i][j].resize(Nmax);
    }
  }

  // Initialize the RB output vectors
  RB_output_vectors.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    RB_output_vectors[n].resize(get_Q_l(n));
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].resize(Nmax);
    }
  }

  // Initialize vectors storing output data
  output_dual_norms.resize(get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
    output_dual_norms[n].resize(Q_l_hat);
  }
  RB_outputs.resize(get_n_outputs());
  RB_output_error_bounds.resize(get_n_outputs());

  RB_system_initialized = true;
}

void RBSystem::attach_dirichlet_dof_initialization (dirichlet_list_fptr dirichlet_init)
{
  libmesh_assert (dirichlet_init != NULL);

  _dirichlet_list_init = dirichlet_init;
}

void RBSystem::initialize_dirichlet_dofs()
{
  START_LOG("initialize_dirichlet_dofs()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to initialize Dirichlet dofs."
              << std::endl;
    libmesh_error();
  }

  // Create a set to store the Dirichlet dofs on this processor
  std::set<unsigned int> dirichlet_dofs_set;
  dirichlet_dofs_set.clear();

  // Initialize the lists of Dirichlet and non-Dirichlet degrees-of-freedom
  if (_dirichlet_list_init != NULL)
    {
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
	      _dirichlet_list_init(context, *this, dirichlet_dofs_set);
	    }
	}
    }

  // Initialize the dirichlet dofs vector on each processor
  std::vector<unsigned int> dirichlet_dofs_vector;
  dirichlet_dofs_vector.clear();

  std::set<unsigned int>::iterator iter     = dirichlet_dofs_set.begin();
  std::set<unsigned int>::iterator iter_end = dirichlet_dofs_set.end();

  for ( ; iter != iter_end; ++iter)
  {
    unsigned int dirichlet_dof_index = *iter;
    dirichlet_dofs_vector.push_back(dirichlet_dof_index);
  }

  // Now take the union over all processors
  Parallel::allgather(dirichlet_dofs_vector);

  // Put all local dofs into non_dirichlet_dofs_set and
  // then erase the Dirichlet dofs
  // Note that this approach automatically ignores non-local Dirichlet dofs
  std::set<unsigned int> non_dirichlet_dofs_set;
  for(unsigned int i=this->get_dof_map().first_dof(); i<this->get_dof_map().end_dof(); i++)
    non_dirichlet_dofs_set.insert(i);

  // Also, initialize the member data structure global_dirichlet_dofs_set
  global_dirichlet_dofs_set.clear();

  for (unsigned int ii=0; ii<dirichlet_dofs_vector.size(); ii++)
  {
    non_dirichlet_dofs_set.erase(dirichlet_dofs_vector[ii]);
    global_dirichlet_dofs_set.insert(dirichlet_dofs_vector[ii]);
  }

  // Finally, load the non-Dirichlet dofs into the system
  iter     = non_dirichlet_dofs_set.begin();
  iter_end = non_dirichlet_dofs_set.end();

  this->non_dirichlet_dofs_vector.clear();

  for ( ; iter != iter_end; ++iter)
    {
      unsigned int non_dirichlet_dof_index = *iter;

      this->non_dirichlet_dofs_vector.push_back(non_dirichlet_dof_index);
    }

  STOP_LOG("initialize_dirichlet_dofs()", "RBSystem");
}

void RBSystem::perform_initial_assembly()
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

  // Compute the dual norms of the outputs
  this->compute_output_dual_norms();
}

AutoPtr<FEMContext> RBSystem::build_context ()
{
  return AutoPtr<FEMContext>(new FEMContext(*this));
}

void RBSystem::add_scaled_matrix_and_vector(Number scalar,
                                 affine_assembly_fptr intrr_assembly,
                                 affine_assembly_fptr bndry_assembly,
                                 SparseMatrix<Number>* input_matrix,
                                 NumericVector<Number>* input_vector,
                                 bool symmetrize)
{
  START_LOG("add_scaled_matrix_and_vector()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
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
    if(intrr_assembly != NULL)
      intrr_assembly(context, *this);

    for (context.side = 0;
          context.side != context.elem->n_sides();
          ++context.side)
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.elem->neighbor(context.side) != NULL) && !impose_internal_fluxes )
        continue;

      // Impose boundary (e.g. Neumann) term
      if( bndry_assembly != NULL )
      {
        context.side_fe_reinit();
        bndry_assembly(context, *this);
      }
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
                                      affine_assembly_fptr intrr_assembly,
                                      affine_assembly_fptr bndry_assembly,
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
    if(intrr_assembly != NULL)
      intrr_assembly(context, *this);

    for (context.side = 0;
         context.side != context.elem->n_sides();
         ++context.side)
    {
      // May not need to apply fluxes on non-boundary elements
      if( (context.elem->neighbor(context.side) != NULL) && !impose_internal_fluxes )
        continue;

      if( bndry_assembly != NULL )
        {
          context.side_fe_reinit();
          bndry_assembly(context, *this);
        }
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
        if(A_q_intrr_assembly_vector[q_a] != NULL)
          this->A_q_intrr_assembly_vector[q_a](*Aq_context[q_a], *this);
      }

      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
        Fq_context[q_f]->pre_fe_reinit(*this, *el);
        Fq_context[q_f]->elem_fe_reinit();
        if(F_q_intrr_assembly_vector[q_f] != NULL)
          this->F_q_intrr_assembly_vector[q_f](*Fq_context[q_f], *this);
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

          if( A_q_bndry_assembly_vector[q_a] != NULL )
          {
            Aq_context[q_a]->side_fe_reinit();
            this->A_q_bndry_assembly_vector[q_a](*Aq_context[q_a], *this);
          }
        }

        // Impose boundary terms, e.g. Neuman BCs
        for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
        {
          // Update the side information for all contexts
          Fq_context[q_f]->side = Aq_context[0]->side;

          if( F_q_bndry_assembly_vector[q_f] != NULL )
          {
            Fq_context[q_f]->side_fe_reinit();
            this->F_q_bndry_assembly_vector[q_f](*Fq_context[q_f], *this);
          }
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
      add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);

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

void RBSystem::assemble_inner_product_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1., inner_prod_assembly, NULL, input_matrix, NULL);
}

void RBSystem::assemble_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1., constraint_assembly, NULL, input_matrix, NULL);
}

void RBSystem::assemble_and_add_constraint_matrix(SparseMatrix<Number>* input_matrix)
{
  add_scaled_matrix_and_vector(1., constraint_assembly, NULL, input_matrix, NULL);
}

void RBSystem::assemble_Aq_matrix(unsigned int q, SparseMatrix<Number>* input_matrix)
{
  if(q >= get_Q_a())
  {
    std::cerr << "Error: We must have q < Q_a in assemble_Aq_matrix."
              << std::endl;
    libmesh_error();
  }

  input_matrix->zero();
  add_scaled_matrix_and_vector(1., A_q_intrr_assembly_vector[q], A_q_bndry_assembly_vector[q], input_matrix, NULL);
}

void RBSystem::add_scaled_Aq(Number scalar, unsigned int q_a, SparseMatrix<Number>* input_matrix, bool symmetrize)
{
  START_LOG("add_scaled_Aq()", "RBSystem");

  if(q_a >= get_Q_a())
  {
    std::cerr << "Error: We must have q < Q_a in add_scaled_Aq."
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
    add_scaled_matrix_and_vector(scalar,
                                 A_q_intrr_assembly_vector[q_a],
                                 A_q_bndry_assembly_vector[q_a],
                                 input_matrix,
                                 NULL,
                                 symmetrize);
  }

  STOP_LOG("add_scaled_Aq()", "RBSystem");
}

void RBSystem::assemble_misc_matrices()
{
  if(low_memory_mode)
  {
    std::cout << "Error: Cannot store misc matrices in low-memory mode." << std::endl;
    libmesh_error();
  }

  assemble_inner_product_matrix(inner_product_matrix.get());

  if( constrained_problem )
    assemble_constraint_matrix(constraint_matrix.get());
}

void RBSystem::assemble_all_affine_operators()
{
  if(low_memory_mode)
  {
    std::cout << "Error: Cannot store affine matrices in low-memory mode." << std::endl;
    libmesh_error();
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    assemble_Aq_matrix(q_a, get_A_q(q_a));
}

void RBSystem::assemble_all_affine_vectors()
{
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    get_F_q(q_f)->zero();
    add_scaled_matrix_and_vector(1., F_q_intrr_assembly_vector[q_f],
                                 F_q_bndry_assembly_vector[q_f], NULL, get_F_q(q_f));
  }
}

void RBSystem::assemble_all_output_vectors()
{
  for(unsigned int n=0; n<get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      get_output_vector(n, q_l)->zero();
      add_scaled_matrix_and_vector(1., output_intrr_assembly_vector[n][q_l],
                                       output_bndry_assembly_vector[n][q_l],
                                       NULL,
                                       get_output_vector(n,q_l));
    }
}

Real RBSystem::train_reduced_basis(const std::string& directory_name)
{
  START_LOG("train_reduced_basis()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to train reduced basis."
              << std::endl;
    libmesh_error();
  }

  int count = 1;

  // Clear the Greedy param list
  for(unsigned int i=0; i<greedy_param_list.size(); i++)
    greedy_param_list[i].clear();
  greedy_param_list.clear();

  Real training_greedy_error;


  // If we are continuing from a previous training run,
  // we might already be at the max number of basis functions.
  // If so, we can just return.
  if (this->get_n_basis_functions() >= Nmax) return 0.;

  while(true)
  {
    std::cout << std::endl << "---- Training solve " << count << " ----" << std::endl;
    print_current_parameters();

    // Update the list of Greedily selected parameters
    this->update_greedy_param_list();

    // Perform an Offline truth solve for the current parameter
    truth_solve(-1);

    // Add orthogonal part of the snapshot to the RB space
    std::cout << std::endl << "Enriching the RB space" << std::endl;
    enrich_RB_space();

    unsigned int RB_size = get_n_basis_functions();
    std::cout << "Reduced basis dimension = " << RB_size << std::endl;

    update_system();

    std::cout << "Performing RB solves on training set" << std::endl;
    training_greedy_error = compute_a_posteriori_bounds();


    std::cout << "Maximum a posteriori error is "
              << training_greedy_error << std::endl << std::endl;

    if(write_data_during_training)
    {
      OStringStream new_dir_name;
      new_dir_name << directory_name << "_" << get_n_basis_functions();
      std::cout << "Writing out RB data to " << new_dir_name.str() << std::endl;
      write_offline_data_to_files(new_dir_name.str());
    }

    // Break out of training phase if we have reached Nmax
    // or if the training_tolerance is satisfied.
    if( greedy_termination_test(training_greedy_error, count) )
    {
      break;
    }

    // Increment counter
    count++;
  }
  this->update_greedy_param_list();
  STOP_LOG("train_reduced_basis()", "RBSystem");

  return training_greedy_error;
}

bool RBSystem::greedy_termination_test(Real training_greedy_error, int)
{
  if(training_greedy_error < this->training_tolerance)
  {
    std::cout << "Specified error tolerance reached." << std::endl;
    return true;
  }

  if(get_n_basis_functions() >= this->get_Nmax())
  {
    std::cout << "Maximum number of basis functions reached: Nmax = "
              << get_Nmax() << std::endl;
    return true;
  }

  return false;
}


void RBSystem::update_greedy_param_list()
{
  greedy_param_list.push_back( get_current_parameters() );
}

std::vector<Real> RBSystem::get_greedy_parameter(unsigned int i)
{
  if( i >= greedy_param_list.size() )
  {
    std::cout << "Error: Argument in RBSystem::get_greedy_parameter is too large."
              << std::endl;
    libmesh_error();
  }

  return greedy_param_list[i];
}

Real RBSystem::truth_solve(int plot_solution)
{
  START_LOG("truth_solve()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to do a truth solve."
              << std::endl;
    libmesh_error();
  }

  truth_assembly();
  solve();

  // Make sure we didn't max out the number of iterations
  if( (this->n_linear_iterations() >=
       this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
      (this->final_linear_residual() >
       this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
  {
      std::cout << "Warning: Linear solver may not have converged! Final linear residual = "
                << this->final_linear_residual() << ", number of iterations = "
                << this->n_linear_iterations() << std::endl << std::endl;
//     libmesh_error();
  }

  truth_outputs.resize(this->get_n_outputs());
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    truth_outputs[n] = 0.;
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      truth_outputs[n] += eval_theta_q_l(n, q_l)*get_output_vector(n,q_l)->dot(*solution);
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
  // We assume that Nmax is initialized in
  // init_data to initial_Nmax, and that it
  // cannot be increased above initial_Nmax
  if(RB_system_initialized && (Nmax_in > initial_Nmax))
  {
    std::cerr << "Error: System was initialized with Nmax = " << initial_Nmax
              << ", cannot set Nmax higher than this value."  << std::endl;
    libmesh_error();
  }

  if(RB_system_initialized && (Nmax_in < this->get_n_basis_functions()))
  {
    std::cerr << "Error: Cannot set Nmax to be less than the "
              << "current number of basis functions."  << std::endl;
    libmesh_error();
  }

  this->Nmax = Nmax_in;
}

void RBSystem::attach_A_q(theta_q_fptr theta_q_a,
                                   affine_assembly_fptr A_q_intrr_assembly,
                                   affine_assembly_fptr A_q_bndry_assembly)
{
  theta_q_a_vector.push_back(theta_q_a);
  A_q_intrr_assembly_vector.push_back(A_q_intrr_assembly);
  A_q_bndry_assembly_vector.push_back(A_q_bndry_assembly);
}

void RBSystem::attach_F_q(theta_q_fptr theta_q_f,
                          affine_assembly_fptr F_q_intrr_assembly,
                          affine_assembly_fptr F_q_bndry_assembly)
{
  theta_q_f_vector.push_back(theta_q_f);
  F_q_intrr_assembly_vector.push_back(F_q_intrr_assembly);
  F_q_bndry_assembly_vector.push_back(F_q_bndry_assembly);
}

void RBSystem::attach_output(std::vector<theta_q_fptr> theta_q_l,
                             std::vector<affine_assembly_fptr> output_intrr_assembly,
                             std::vector<affine_assembly_fptr> output_bndry_assembly)
{
  // Make sure the input vectors are all the same size!
  if( (theta_q_l.size() == output_intrr_assembly.size()) &&
      (output_intrr_assembly.size() == output_bndry_assembly.size()) )
  {
    theta_q_l_vector.push_back(theta_q_l);
    output_intrr_assembly_vector.push_back(output_intrr_assembly);
    output_bndry_assembly_vector.push_back(output_bndry_assembly);
  }
  else
  {
    std::cout << "Error: The input vectors in attach_output must all be the same size in attach_output"
              << std::endl;
    libmesh_error();
  }
}

void RBSystem::attach_inner_prod_assembly(affine_assembly_fptr IP_assembly)
{
  inner_prod_assembly = IP_assembly;
}

void RBSystem::attach_constraint_assembly(affine_assembly_fptr constraint_assembly_in)
{
  constraint_assembly = constraint_assembly_in;
}

unsigned int RBSystem::get_Q_l(unsigned int index) const
{
  if(index >= get_n_outputs())
  {
    std::cerr << "Error: We must have index < n_outputs in get_Q_l."
              << std::endl;
    libmesh_error();
  }
  return theta_q_l_vector[index].size();
}

Number RBSystem::eval_theta_q_f(unsigned int q)
{
  if(q >= get_Q_f())
  {
    std::cerr << "Error: We must have q < Q_f in eval_theta_q_f."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_f_vector[q] != NULL);

  return theta_q_f_vector[q](current_parameters);
}

Number RBSystem::eval_theta_q_l(unsigned int output_index, unsigned int q_l)
{
  if( (output_index >= get_n_outputs()) || (q_l >= get_Q_l(output_index)) )
  {
    std::cerr << "Error: We must have output_index < n_outputs and "
              << "q_l < get_Q_l(output_index) in eval_theta_q_l."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(theta_q_l_vector[output_index][q_l] != NULL);

  return theta_q_l_vector[output_index][q_l](current_parameters);
}

void RBSystem::load_basis_function(unsigned int i)
{
  START_LOG("load_basis_function()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to load basis function."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(i < basis_functions.size());

  *solution = *basis_functions[i];

  // synchronise solution and current_local_solution
  this->update();

  STOP_LOG("load_basis_function()", "RBSystem");
}

void RBSystem::load_basis_function(unsigned int i, NumericVector<Number>& vec)
{
  START_LOG("load_basis_function()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to load basis function."
              << std::endl;
    libmesh_error();
  }

  libmesh_assert(i < basis_functions.size());

  vec = *basis_functions[i];

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

  for(unsigned int index=0; index<basis_functions.size(); index++)
  {
    // invoke copy constructor for NumericVector
    *proj_index = *basis_functions[index];
    if(!low_memory_mode)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector,*proj_index);
    }
    else
    {
      matrix->vector_mult(*inner_product_storage_vector,*proj_index);
    }
    Number scalar = new_bf->dot(*inner_product_storage_vector);
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
  Number new_bf_norm = std::sqrt( new_bf->dot(*inner_product_storage_vector) );
  new_bf->scale(1./new_bf_norm);

  // load the new basis function into the basis_functions vector.
  basis_functions.push_back( new_bf );

  STOP_LOG("enrich_RB_space()", "RBSystem");
}

void RBSystem::update_system()
{
  std::cout << "Updating RB matrices" << std::endl;
  update_RB_system_matrices();

  std::cout << "Updating RB residual terms" << std::endl;

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

void RBSystem::recompute_all_residual_terms()
{
  unsigned int saved_delta_N = delta_N;
  delta_N = get_n_basis_functions();

  // Use alternative solver for residual terms solves
  std::pair<std::string,std::string> orig_solver =
    this->set_alternative_solver(this->linear_solver);
  
  update_residual_terms(/*compute_inner_products=*/false);

  // Return to original solver
  this->reset_alternative_solver(this->linear_solver, orig_solver);
  
  delta_N = saved_delta_N;
}

Real RBSystem::compute_a_posteriori_bounds()
{
  START_LOG("compute_a_posteriori_bounds()", "RBSystem");

  unsigned int RB_size = basis_functions.size();


  training_error_bounds.resize(this->get_local_n_training_samples());

  // keep track of the maximum error
  unsigned int max_err_index = 0.;
  Real max_err = 0.;

  unsigned int first_index = get_first_local_training_index();
  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
  {
    // Load training parameter i, this is only loaded
    // locally since the RB solves are local.
    set_current_parameters( get_training_parameter(first_index+i) );

    training_error_bounds[i] = RB_solve(RB_size);
//     std::cout << "Error bound at training index " << first_index+i << " is "
//               << training_error_bounds[i] << std::endl;

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

  STOP_LOG("compute_a_posteriori_bounds()", "RBSystem");

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
  std::cout << "SLEPc and GLPK must be installed for SCM functions to work." << std::endl;
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
  std::cout << "SLEPc and GLPK must be installed for SCM functions to work." << std::endl;
  libmesh_error();
#endif // defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)
}


Real RBSystem::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "RBSystem");

  if(N > get_n_basis_functions())
  {
    std::cerr << "ERROR: N cannot be larger than the number "
              << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(N==0)
  {
    std::cerr << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }

  // Resize (and clear) the solution vector
  RB_solution.resize(N);

  // Assemble the RB system
  DenseMatrix<Number> RB_system_matrix(N,N);
  RB_system_matrix.zero();

  DenseMatrix<Number> RB_A_q_a;
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    RB_A_q_vector[q_a].get_principal_submatrix(N, RB_A_q_a);

    RB_system_matrix.add(eval_theta_q_a(q_a), RB_A_q_a);
  }

  // Assemble the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  DenseVector<Number> RB_F_q_f;
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    RB_F_q_vector[q_f].get_principal_subvector(N, RB_F_q_f);

    RB_rhs.add(eval_theta_q_f(q_f), RB_F_q_f);
  }
  
  // Solve the linear system
  RB_system_matrix.lu_solve(RB_rhs, RB_solution);

  // Evaluate the dual norm of the residual for RB_solution_vector
  Real epsilon_N = compute_residual_dual_norm(N);

  // Get lower bound for coercivity constant
  const Real alpha_LB = get_SCM_lower_bound();
  // alpha_LB needs to be positive to get a valid error bound
  libmesh_assert( alpha_LB > 0. );

  // Store (absolute) error bound
  Real abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

  // Compute the norm of RB_solution
  Real RB_solution_norm = RB_solution.l2_norm();

  // Now compute the outputs and associated errors
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<get_n_outputs(); n++)
  {
    RB_outputs[n] = 0.;
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs[n] += eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_solution);
    }
    
    Number output_bound_sq = 0.;    
    unsigned int q=0;
    for(unsigned int q_l1=0; q_l1<get_Q_l(n); q_l1++)
    {
      for(unsigned int q_l2=q_l1; q_l2<get_Q_l(n); q_l2++)
      {
        Real delta = (q_l1==q_l2) ? 1. : 2.;
        output_bound_sq += delta*eval_theta_q_l(n,q_l1)*eval_theta_q_l(n,q_l2) * output_dual_norms[n][q];
        q++;
      }
    }
    
    RB_output_error_bounds[n] = abs_error_bound * libmesh_real(std::sqrt( output_bound_sq ));
  }

  STOP_LOG("RB_solve()", "RBSystem");

  return ( return_rel_error_bound ? abs_error_bound/RB_solution_norm : abs_error_bound );
}

Real RBSystem::residual_scaling_denom(Real alpha_LB)
{
  // Here we implement the residual scaling for a coercive
  // problem.
  return std::sqrt(alpha_LB);
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
      RB_F_q_vector[q_f](i) = get_F_q(q_f)->dot(*basis_functions[i]);
    }
  }

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
  {
    for(unsigned int n=0; n<get_n_outputs(); n++)
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l](i) = get_output_vector(n,q_l)->dot(*basis_functions[i]);
      }

    for(unsigned int j=0; j<RB_size; j++)
    {
      Number value = 0.;

      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        // Compute reduced A_q matrix
        temp->zero();
        if(!low_memory_mode)
        {
          get_A_q(q_a)->vector_mult(*temp, *basis_functions[j]);
        }
        else
        {
          assemble_Aq_matrix(q_a,matrix);
          matrix->vector_mult(*temp, *basis_functions[j]);
        }

        value = (*basis_functions[i]).dot(*temp);
        RB_A_q_vector[q_a](i,j) = value;

        if(i!=j)
        {
          temp->zero();
          if(!low_memory_mode)
          {
            get_A_q(q_a)->vector_mult(*temp, *basis_functions[i]);
          }
          else
          {
            // matrix should still hold affine matrix q_a
            matrix->vector_mult(*temp, *basis_functions[i]);
          }

          value = (*basis_functions[j]).dot(*temp);
          RB_A_q_vector[q_a](j,i) = value;
        }
      }
    }
  }

  STOP_LOG("update_RB_system_matrices()", "RBSystem");
}


void RBSystem::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "RBSystem");

  // First we need to compute the representors of
  // each basis function using inner_product_matrix

  unsigned int RB_size = get_n_basis_functions();

  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    linear_solver->same_preconditioner = false;
  }

  // We only need to compute the representors for
  // the right-hand side once.  Note: this will be called at least
  // once even after a restart, in which case F_q_representor[0]
  // will already be set.
  if(!update_residual_terms_called)
  {
    if(low_memory_mode)
    {
      assemble_inner_product_matrix(matrix);
      if(constrained_problem)
        add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);
    }

    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
      {
	F_q_representor[q_f] = (NumericVector<Number>::build().release());
	F_q_representor[q_f]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

	rhs->zero();
	rhs->add(1., *get_F_q(q_f));
	zero_dirichlet_dofs_on_rhs();

	solution->zero();

	if (!quiet)
	  std::cout << "Starting solve q_f=" << q_f
		    << " in RBSystem::update_residual_terms() at "
		    << Utility::get_timestamp() << std::endl;

	solve();

	if (!quiet)
	  {
	    std::cout << "Finished solve q_f=" << q_f
		      << " in RBSystem::update_residual_terms() at "
		      << Utility::get_timestamp() << std::endl;

	    std::cout << this->n_linear_iterations()
		      << " iterations, final residual "
		      << this->final_linear_residual() << std::endl;
	  }

	// Make sure we didn't max out the number of iterations
	if( (this->n_linear_iterations() >=
	     this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
	    (this->final_linear_residual() >
	     this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
	  {
	    std::cout << "Warning: Linear solver may not have converged! Final linear residual = "
		      << this->final_linear_residual() << ", number of iterations = "
		      << this->n_linear_iterations() << std::endl << std::endl;
	    //         libmesh_error();

	  }
	*F_q_representor[q_f] = *solution;

	if(reuse_preconditioner)
	  {
	    // After we do a solve, tell PETSc we want to reuse the preconditioner
	    // since the system matrix is not changing.
	    linear_solver->same_preconditioner = true;
	  }
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
		Fq_representor_norms[q] = F_q_representor[q_f2]->dot(*inner_product_storage_vector);

		q++;
	      }
	  }
      } // end if (compute_inner_products)

    
  } // end if(!update_residual_terms_called)


  if(low_memory_mode)
  {
    assemble_inner_product_matrix(matrix);
    if(constrained_problem)
      add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);
  }

  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      // Initialize the vector in which we'll store the representor
      A_q_representor[q_a][i] = (NumericVector<Number>::build().release());
      A_q_representor[q_a][i]->init(this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

      rhs->zero();
      if(!low_memory_mode)
      {
        get_A_q(q_a)->vector_mult(*rhs, *basis_functions[i]);
      }
      else
      {
        assemble_scaled_matvec(1.,
                               A_q_intrr_assembly_vector[q_a],
                               A_q_bndry_assembly_vector[q_a],
                               *rhs,
                               *basis_functions[i]);
      }
      rhs->scale(-1.);
      zero_dirichlet_dofs_on_rhs();

      solution->zero();
      if (!quiet)
	    {
        std::cout << "Starting solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBSystem::update_residual_terms() at "
                  << Utility::get_timestamp() << std::endl;
	    }

      solve();

      if (!quiet)
	    {
        std::cout << "Finished solve [q_a][i]=[" << q_a <<"]["<< i << "] in RBSystem::update_residual_terms() at "
                  << Utility::get_timestamp() << std::endl;
        std::cout << this->n_linear_iterations() << " iterations, final residual "
                  << this->final_linear_residual() << std::endl;
	    }

      // Make sure we didn't max out the number of iterations
      if( (this->n_linear_iterations() >=
          this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
          (this->final_linear_residual() >
          this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
      {
        std::cout << "Warning: Linear solver may not have converged! Final linear residual = "
                  << this->final_linear_residual() << ", number of iterations = "
                  << this->n_linear_iterations() << std::endl << std::endl;
//         libmesh_error();
      }

      // Store the representor
      *A_q_representor[q_a][i] = *solution;


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
		  Fq_Aq_representor_norms[q_f][q_a][i] =
		    A_q_representor[q_a][i]->dot(*inner_product_storage_vector);
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
			  inner_product_matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a2][j]);
			}
		      else
			{
			  matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a2][j]);
			}
		      Aq_Aq_representor_norms[q][i][j] = A_q_representor[q_a1][i]->dot(*inner_product_storage_vector);

		      if(i != j)
			{
			  if(!low_memory_mode)
			    {
			      inner_product_matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a2][i]);
			    }
			  else
			    {
			      matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a2][i]);
			    }
			  Aq_Aq_representor_norms[q][j][i] = A_q_representor[q_a1][j]->dot(*inner_product_storage_vector);
			}
		    }
		}
	      q++;
	    }
	}
    } // end if (compute_inner_products)
  
  update_residual_terms_called = true;

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
  START_LOG("compute_output_dual_norms()", "RBSystem");
  
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
          add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);
        }
      }
    }

    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      rhs->zero();
            rhs->add(1., *get_output_vector(n,q_l));
      zero_dirichlet_dofs_on_rhs();

      solution->zero();

      if (!quiet)
        std::cout << "Starting solve n=" << n << ", q_l=" << q_l
            << " in RBSystem::compute_output_dual_norms() at "
            << Utility::get_timestamp() << std::endl;

      solve();

      if (!quiet)
        {
          std::cout << "Finished solve n=" << n << ", q_l=" << q_l
                    << " in RBSystem::compute_output_dual_norms() at "
                    << Utility::get_timestamp() << std::endl;

          std::cout << this->n_linear_iterations()
                    << " iterations, final residual "
                    << this->final_linear_residual() << std::endl;
        }

      // Make sure we didn't max out the number of iterations
      if( (this->n_linear_iterations() >=
           this->get_equation_systems().parameters.get<unsigned int>("linear solver maximum iterations")) &&
          (this->final_linear_residual() >
           this->get_equation_systems().parameters.get<Real>("linear solver tolerance")) )
         {
           std::cout << "Warning: Linear solver may not have converged! Final linear residual = "
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
        std::cout << "output_dual_norms[" << n << "][" << q << "] = " << output_dual_norms[n][q] << std::endl;
        
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

  // Change the preconditioner, Krylov solver back to their original
  // value.  Note: does nothing if RBBase::alternative_solver ==
  // "unchanged".
  this->reset_alternative_solver(this->linear_solver, orig_solver);

  STOP_LOG("compute_output_dual_norms()", "RBSystem");
}

void RBSystem::load_RB_solution()
{
  START_LOG("load_RB_solution()", "RBSystem");

  if(!initialize_calN_dependent_data)
  {
    std::cerr << "Error: We must initialize the calN dependent "
              << "data structures in order to load RB solution."
              << std::endl;
    libmesh_error();
  }

  solution->zero();

  if(RB_solution.size() > basis_functions.size())
  {
    std::cerr << "ERROR: System contains " << basis_functions.size() << " basis functions."
              << " RB_solution vector constains " << RB_solution.size() << " entries."
              << " RB_solution in RBSystem::load_RB_solution is too long!" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<RB_solution.size(); i++)
  {
    solution->add(RB_solution(i), *basis_functions[i]);
  }

  update();

  STOP_LOG("load_RB_solution()", "RBSystem");
}

Real RBSystem::compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("compute_residual_dual_norm()", "RBSystem");

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
//     RB_sol->add(RB_solution(i), *basis_functions[i]);
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
//     std::cout << "Warning: Linear solver may not have converged! Final linear residual = "
//               << this->final_linear_residual() << ", number of iterations = "
//               << this->n_linear_iterations() << std::endl << std::endl;
// //     libmesh_error();
//   }
//
//   Real slow_residual_norm_sq = inner_product( *solution, *solution );


  // Use the stored representor inner product values
  // to evaluate the residual norm
  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<get_Q_f(); q_f1++)
  {
    for(unsigned int q_f2=q_f1; q_f2<get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      residual_norm_sq += delta*eval_theta_q_f(q_f1)*eval_theta_q_f(q_f2) * Fq_representor_norms[q];

      q++;
    }
  }

  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      for(unsigned int i=0; i<N; i++)
      {
        Real delta = 2.;
        residual_norm_sq +=
          RB_solution(i)*delta*eval_theta_q_f(q_f)*eval_theta_q_a(q_a) * Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  q=0;
  for(unsigned int q_a1=0; q_a1<get_Q_a(); q_a1++)
  {
    for(unsigned int q_a2=q_a1; q_a2<get_Q_a(); q_a2++)
    {
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq +=
            RB_solution(i)*RB_solution(j)* delta * eval_theta_q_a(q_a1) * eval_theta_q_a(q_a2) *
            Aq_Aq_representor_norms[q][i][j];
        }
      }

      q++;
    }
  }

//  if(libmesh_real(residual_norm_sq) < 0.)
//  {
//    std::cout << "Warning: Square of residual norm is negative "
//              << "in RBSystem::compute_residual_dual_norm()" << std::endl;

    // Sometimes this is negative due to rounding error,
    // but error is on the order of 1.e-10, so shouldn't
    // affect error bound much...
//     libmesh_error();
//    residual_norm_sq = std::abs(residual_norm_sq);
//  }

//   std::cout << "Slow residual norm squared = " << slow_residual_norm_sq
//             << ", fast residual norm squared = " << residual_norm_sq << std::endl;

  STOP_LOG("compute_residual_dual_norm()", "RBSystem");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}

SparseMatrix<Number>* RBSystem::get_A_q(unsigned int q)
{
  if(low_memory_mode)
  {
    std::cerr << "Error: The affine matrices are not store in low-memory mode." << std::endl;
    libmesh_error();
  }

  if(q >= get_Q_a())
  {
    std::cerr << "Error: We must have q < Q_a in get_A_q."
              << std::endl;
    libmesh_error();
  }

  return A_q_vector[q];
}

NumericVector<Number>* RBSystem::get_F_q(unsigned int q)
{
  if(q >= get_Q_f())
  {
    std::cerr << "Error: We must have q < Q_f in get_F_q."
              << std::endl;
    libmesh_error();
  }

  return F_q_vector[q];
}

NumericVector<Number>* RBSystem::get_output_vector(unsigned int n, unsigned int q_l)
{
  if( (n >= get_n_outputs()) || (q_l >= get_Q_l(n)) )
  {
    std::cerr << "Error: We must have n < n_outputs and "
              << "q_l < get_Q_l(n) in get_output_vector."
              << std::endl;
    libmesh_error();
  }

  return outputs_vector[n][q_l];
}

void RBSystem::zero_dirichlet_dofs_on_rhs()
{
  START_LOG("zero_dirichlet_dofs_on_rhs()", "RBSystem");

  this->zero_dirichlet_dofs_on_vector(*rhs);

  STOP_LOG("zero_dirichlet_dofs_on_rhs()", "RBSystem");
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
      temp.set(index, 0.);
  }
  temp.close();

  STOP_LOG("zero_dirichlet_dofs_on_vector()", "RBSystem");
}




void RBSystem::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBSystem");

  const unsigned int precision_level = 14;

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= Nmax );

  if(libMesh::processor_id() == 0)
  {

    // Make a directory to store all the data files
    if( mkdir(directory_name.c_str(), 0777) == -1)
    {
      std::cout << "In RBSystem::write_offline_data_to_files, directory "
                << directory_name << " already exists, overwriting contents." << std::endl;
    }

    // First, write out how many basis functions we have generated
    {
      std::ofstream n_bfs_out;
      {
        OStringStream file_name;
        file_name << directory_name << "/n_bfs.dat";
        n_bfs_out.open(file_name.str().c_str());
      }
      if ( !n_bfs_out.good() )
      {
        std::cerr << "Error opening n_bfs.dat" << std::endl;
        libmesh_error();
      }
      n_bfs_out << n_bfs;
      n_bfs_out.close();
    }

    // Also, write out the greedily selected parameters
    {
      std::ofstream greedy_params_out;
      {
        OStringStream file_name;
        file_name << directory_name << "/greedy_params.dat";
        greedy_params_out.open(file_name.str().c_str());
      }
      if ( !greedy_params_out.good() )
      {
        std::cerr << "Error opening greedy_params.dat" << std::endl;
        libmesh_error();
      }
      for(unsigned int i=0; i<greedy_param_list.size(); i++)
      {
        for(unsigned int j=0; j<get_n_params(); j++)
        {
          greedy_params_out << greedy_param_list[i][j] << " ";
        }
        greedy_params_out << std::endl;
      }
      greedy_params_out.close();
    }

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
        std::cerr << "Error opening output " << n << " dual norms file" << std::endl;
        libmesh_error();
      }
      output_dual_norms_out.precision(precision_level);
      
      unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
      for(unsigned int q=0; q<Q_l_hat; q++)
      {
        output_dual_norms_out << std::scientific << output_dual_norms[n][q] << " ";
      }
      output_dual_norms_out.close();
      
      for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
      {
        std::ofstream output_n_out;
        {
          OStringStream file_name;
          file_name << directory_name << "/output_";
          OSSRealzeroright(file_name,3,0,n);
          file_name << "_";
          OSSRealzeroright(file_name,3,0,q_l);
          file_name << ".dat";
          output_n_out.open(file_name.str().c_str());
        }
        if( !output_n_out.good() )
        {
          std::cerr << "Error opening output file for output " << n << std::endl;
          libmesh_error();
        }
        output_n_out.precision(precision_level);

        for(unsigned int j=0; j<n_bfs; j++)
        {
          output_n_out << std::scientific << RB_output_vectors[n][q_l](j) << " ";
        }
        output_n_out.close();
      }
    }

    // Next write out the F_q vectors
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_F_";
      OSSRealzeroright(file_name,3,0,q_f);
      file_name << ".dat";
      std::ofstream RB_F_q_f_out(file_name.str().c_str());

      if ( !RB_F_q_f_out.good() )
      {
        std::cerr << "Error opening RB_F_" << q_f << ".dat" << std::endl;
        libmesh_error();
      }

      RB_F_q_f_out.precision(precision_level);
      for(unsigned int i=0; i<n_bfs; i++)
      {
        RB_F_q_f_out << std::scientific << RB_F_q_vector[q_f](i) << " ";
      }
      RB_F_q_f_out.close();
    }

    // Next write out the A_q matrices
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_A_";
      OSSRealzeroright(file_name,3,0,q_a);
      file_name << ".dat";
      std::ofstream RB_A_q_a_out(file_name.str().c_str());

      if ( !RB_A_q_a_out.good() )
      {
        std::cerr << "Error opening RB_A_" << q_a << ".dat" << std::endl;
        libmesh_error();
      }

      RB_A_q_a_out.precision(precision_level);
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_A_q_a_out << std::scientific << RB_A_q_vector[q_a](i,j) << " ";
        }
      }
      RB_A_q_a_out.close();
    }


    // Next write out F_q representor norm data
    std::ofstream RB_Fq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_norms.dat";
      RB_Fq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_norms_out.good() )
    {
      std::cerr << "Error opening Fq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_norms_out.precision(precision_level);
    unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
    for(unsigned int i=0; i<Q_f_hat; i++)
    {
      RB_Fq_norms_out << std::scientific << Fq_representor_norms[i] << " ";
    }
    RB_Fq_norms_out.close();

    // Next write out Fq_Aq representor norm data
    std::ofstream RB_Fq_Aq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_Aq_norms.dat";
      RB_Fq_Aq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_Aq_norms_out.good() )
    {
      std::cerr << "Error opening Fq_Aq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_Aq_norms_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
    {
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
      {
        for(unsigned int i=0; i<n_bfs; i++)
        {
          RB_Fq_Aq_norms_out << std::scientific << Fq_Aq_representor_norms[q_f][q_a][i] << " ";
        }
      }
    }
    RB_Fq_Aq_norms_out.close();

    // Next write out Aq_Aq representor norm data
    std::ofstream RB_Aq_Aq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_Aq_norms.dat";
      RB_Aq_Aq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Aq_Aq_norms_out.good() )
    {
      std::cerr << "Error opening Aq_Aq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Aq_Aq_norms_out.precision(precision_level);
    unsigned int Q_a_hat = get_Q_a()*(get_Q_a()+1)/2;
    for(unsigned int i=0; i<Q_a_hat; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        for(unsigned int l=0; l<n_bfs; l++)
        {
          RB_Aq_Aq_norms_out << std::scientific << Aq_Aq_representor_norms[i][j][l] << " ";
        }
      }
    }
    RB_Aq_Aq_norms_out.close();

  }

  // Now write out the basis functions if requested
  if(store_basis_functions)
  {
    std::cout << "Writing out the basis functions..." << std::endl;

    std::ostringstream file_name;
    const std::string basis_function_suffix = (write_binary_basis_functions ? ".xdr" : ".dat");

    // Use System::write_serialized_data to write out the basis functions
    // by copying them into this->solution one at a time.
    for(unsigned int i=0; i<n_bfs; i++)
    {
      // No need to copy, just swap
      // *solution = *basis_functions[i];
      basis_functions[i]->swap(*solution);

      file_name.str(""); // reset the string
      file_name << directory_name << "/bf" << i << basis_function_suffix;

      Xdr bf_data(file_name.str(),
		  write_binary_basis_functions ? ENCODE : WRITE);

      write_serialized_data(bf_data, false);

      // Synchronize before moving on
      Parallel::barrier();

      // Swap back
      basis_functions[i]->swap(*solution);
    }
  }

  // Write out residual representors if requested
  if (store_representors)
    {
      // Write out F_q_representors.  These are useful to have when restarting,
      // so you don't have to recompute them all over again.  There should be
      // Q_f of these.
      if (!quiet)
	std::cout << "Writing out the F_q_representors..." << std::endl;

      std::ostringstream file_name;
      const std::string residual_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");
      struct stat stat_info;

      // Residual representors written out to their own separate directory
      std::string residual_representors_dir = "residual_representors";
      if ( libMesh::processor_id() == 0)
	if ( mkdir(residual_representors_dir.c_str(), 0755) != 0)
	  std::cout << "Skipping creating residual_representors directory: " << strerror(errno) << std::endl;

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

      // Write out A_q_representors.  These are useful to have when restarting,
      // so you don't have to recompute them all over again.  There should be
      // Q_a * this->get_n_basis_functions() of these.
      if (!quiet)
	std::cout << "Writing out the A_q_representors..." << std::endl;

      const unsigned int jstop  = this->get_n_basis_functions();
      const unsigned int jstart = jstop-delta_N;
      for (unsigned int i=0; i<A_q_representor.size(); ++i)
	for (unsigned int j=jstart; j<jstop; ++j)
	  {
	    std::cout << "Writing out A_q_representor[" << i << "][" << j << "]..." << std::endl;
	    libmesh_assert(A_q_representor[i][j] != NULL);

	    file_name.str(""); // reset filename
	    file_name << residual_representors_dir
		      << "/A_q_representor" << i << "_" << j << residual_representor_suffix;

	    {
	      // No need to copy!
	      // *solution = *(A_q_representor[i][j]);
	      A_q_representor[i][j]->swap(*solution);

	      Xdr aqr_data(file_name.str(),
			   write_binary_residual_representors ? ENCODE : WRITE);

	      write_serialized_data(aqr_data, false);

	      // Synchronize before moving on
	      Parallel::barrier();

	      // Swap back.
	      A_q_representor[i][j]->swap(*solution);

	      // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
	      // for the system call, be sure to do it only on one processor, etc.
	    }
	  }
    } // end if (store_representors)

  STOP_LOG("write_offline_data_to_files()", "RBSystem");
}



void RBSystem::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBSystem");

  // First, find out how many basis functions we had when Greedy terminated
  unsigned int n_bfs;
  {
    OStringStream file_name;
    file_name << directory_name << "/n_bfs.dat";
    std::ifstream n_bfs_in(file_name.str().c_str());

    if ( !n_bfs_in.good() )
    {
      std::cerr << "Error opening n_bfs.dat" << std::endl;
      libmesh_error();
    }

    n_bfs_in >> n_bfs;
    n_bfs_in.close();
  }
  libmesh_assert( n_bfs <= Nmax );

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
      std::cerr << "Error opening input " << n << " dual norms file" << std::endl;
      libmesh_error();
    }
    
    unsigned int Q_l_hat = get_Q_l(n)*(get_Q_l(n)+1)/2;
    for(unsigned int q=0; q<Q_l_hat; q++)
    {
      output_dual_norms_in >> output_dual_norms[n][q];
    }
    output_dual_norms_in.close();
    
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      std::ifstream output_n_in;
      {
        OStringStream file_name;
        file_name << directory_name << "/output_";
        OSSRealzeroright(file_name,3,0,n);
        file_name << "_";
        OSSRealzeroright(file_name,3,0,q_l);
        file_name << ".dat";
        output_n_in.open(file_name.str().c_str());
      }
      if( !output_n_in.good() )
      {
        std::cerr << "Error opening input file for output " << n << std::endl;
        libmesh_error();
      }

      for(unsigned int j=0; j<n_bfs; j++)
      {
        Number  value;
        output_n_in >> value;
        RB_output_vectors[n][q_l](j) = value;
      }
      output_n_in.close();
    }
  }

  // Next read in the F_q vectors
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_F_";
    OSSRealzeroright(file_name,3,0,q_f);
    file_name << ".dat";
    std::ifstream RB_F_q_f_in(file_name.str().c_str());

    if ( !RB_F_q_f_in.good() )
    {
      std::cerr << "Error opening RB_F_" << q_f << ".dat" << std::endl;
      libmesh_error();
    }

    for(unsigned int i=0; i<n_bfs; i++)
    {
      Number  value;
      RB_F_q_f_in >> value;
      RB_F_q_vector[q_f](i) = value;
    }
    RB_F_q_f_in.close();
  }

  // Next read in the A_q matrices
  for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_A_";
    OSSRealzeroright(file_name,3,0,q_a);
    file_name << ".dat";
    std::ifstream RB_A_q_a_in(file_name.str().c_str());

    if ( !RB_A_q_a_in.good() )
    {
      std::cerr << "Error opening RB_A_" << q_a << ".dat" << std::endl;
      libmesh_error();
    }

    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        Number  value;
        RB_A_q_a_in >> value;
        RB_A_q_vector[q_a](i,j) = value;
      }
    }
    RB_A_q_a_in.close();
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
    std::cerr << "Error opening Fq_norms.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_f_hat = get_Q_f()*(get_Q_f()+1)/2;
  for(unsigned int i=0; i<Q_f_hat; i++)
  {
    RB_Fq_norms_in >> Fq_representor_norms[i];
  }
  RB_Fq_norms_in.close();

  // Next read in Fq_Aq representor norm data
  std::ifstream RB_Fq_Aq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Fq_Aq_norms.dat";
    RB_Fq_Aq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Fq_Aq_norms_in.good() )
  {
    std::cerr << "Error opening Fq_Aq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
  {
    for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        RB_Fq_Aq_norms_in >> Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }
  RB_Fq_Aq_norms_in.close();

  // Next read in Aq_Aq representor norm data
  std::ifstream RB_Aq_Aq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_Aq_norms.dat";
    RB_Aq_Aq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Aq_Aq_norms_in.good() )
  {
    std::cerr << "Error opening Aq_Aq_norms.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_a_hat = get_Q_a()*(get_Q_a()+1)/2;
  for(unsigned int i=0; i<Q_a_hat; i++)
  {
    for(unsigned int j=0; j<n_bfs; j++)
    {
      for(unsigned int l=0; l<n_bfs; l++)
      {
        RB_Aq_Aq_norms_in >> Aq_Aq_representor_norms[i][j][l];
      }
    }
  }
  RB_Aq_Aq_norms_in.close();

  // Resize basis_functions even if we don't read them in so that
  // get_n_bfs() returns the correct value. Initialize the pointers
  // to NULL
  basis_functions.resize(n_bfs);
  for(unsigned int i=0; i<n_bfs; i++)
    {
      basis_functions[i] = NULL;
    }

  if(store_basis_functions)
  {
    std::cout << "Reading in the basis functions..." << std::endl;

    std::ostringstream file_name;
    const std::string basis_function_suffix = (read_binary_basis_functions ? ".xdr" : ".dat");
    struct stat stat_info;

    // Use System::read_serialized_data to read in the basis functions
    // into this->solution and then swap with the appropriate
    // of basis function.
    for(unsigned int i=0; i<n_bfs; i++)
    {
      file_name.str(""); // reset the string
      file_name << directory_name << "/bf" << i << basis_function_suffix;

      // On processor zero check to be sure the file exists
      if (libMesh::processor_id() == 0)
	{
	  int stat_result = stat(file_name.str().c_str(), &stat_info);

	  if (stat_result != 0)
	    {
	      std::cout << "File does not exist: " << file_name.str() << std::endl;
	      libmesh_error();
	    }
	}

      Xdr bf_data(file_name.str(),
		  read_binary_basis_functions ? DECODE : READ);

      read_serialized_data(bf_data, false);

      basis_functions[i] = NumericVector<Number>::build().release();
      basis_functions[i]->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

      // No need to copy, just swap
      // *basis_functions[i] = *solution;
      basis_functions[i]->swap(*solution);
    }
  }

  // Read in the representor vectors if requested
  if (store_representors)
    {
      std::cout << "Reading in the F_q_representors..." << std::endl;

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
	      std::cout << "Error, must delete existing F_q_representor before reading in from file."
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
		  std::cout << "File does not exist: " << file_name.str() << std::endl;
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
      update_residual_terms_called=true;

      std::cout << "Reading in the A_q_representors..." << std::endl;

      // Read in the A_q representors.  The class makes room for [Q_a][Nmax] of these.  We are going to
      // read in [Q_a][this->get_n_basis_functions()].  FIXME:
      // should we be worried about leaks in the locations where we're about to fill entries?
      for (unsigned int i=0; i<A_q_representor.size(); ++i)
	for (unsigned int j=0; j<A_q_representor[i].size(); ++j)
	  {
	    if (A_q_representor[i][j] != NULL)
	      {
		std::cout << "Error, must delete existing A_q_representor before reading in from file."
			  << std::endl;
		libmesh_error();
	      }
	  }

      // Now ready to read them in from file!
      for (unsigned int i=0; i<A_q_representor.size(); ++i)
	for (unsigned int j=0; j<this->get_n_basis_functions(); ++j)
	  {
	    file_name.str(""); // reset filename
	    file_name << residual_representors_dir
		      << "/A_q_representor" << i << "_" << j << residual_representor_suffix;

	    // On processor zero check to be sure the file exists
	    if (libMesh::processor_id() == 0)
	    {
	      int stat_result = stat(file_name.str().c_str(), &stat_info);

	      if (stat_result != 0)
		{
		  std::cout << "File does not exist: " << file_name.str() << std::endl;
		  libmesh_error();
		}
	    }

	    Xdr aqr_data(file_name.str(),
			 read_binary_residual_representors ? DECODE : READ);

	    read_serialized_data(aqr_data, false);

	    A_q_representor[i][j] = NumericVector<Number>::build().release();
	    A_q_representor[i][j]->init (this->n_dofs(), this->n_local_dofs(),
					 false, libMeshEnums::PARALLEL);

	    // No need to copy, just swap
	    //*A_q_representor[i][j] = *solution;
	    A_q_representor[i][j]->swap(*solution);
	  }
    } // end if (store_representors)
  STOP_LOG("read_offline_data_from_files()", "RBSystem");
}

} // namespace libMesh



