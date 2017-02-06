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
#include "libmesh/transient_rb_construction.h"
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/transient_rb_theta_expansion.h"
#include "libmesh/transient_rb_assembly_expansion.h"

// LibMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/linear_solver.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/getpot.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/timestamp.h"
#include "libmesh/parallel.h"

// For checking for the existence of files
#include <sys/stat.h>

#include <fstream>
#include <sstream>

// Need SLEPc to get the POD eigenvalues
#if defined(LIBMESH_HAVE_SLEPC)
// LAPACK include (via SLEPc)
#include <petscsys.h>
#include <slepcblaslapack.h>
#endif // LIBMESH_HAVE_SLEPC

namespace libMesh
{

TransientRBConstruction::TransientRBConstruction (EquationSystems & es,
                                                  const std::string & name_in,
                                                  const unsigned int number_in)
  : Parent(es, name_in, number_in),
    L2_matrix(SparseMatrix<Number>::build(es.comm())),
    non_dirichlet_L2_matrix(SparseMatrix<Number>::build(es.comm())),
    nonzero_initialization(false),
    compute_truth_projection_error(false),
    init_filename(""),
    POD_tol(-1.),
    max_truth_solves(-1),
    L2_assembly(libmesh_nullptr)
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;

  temporal_data.resize(0);

  // We should not necessarily exit the greedy due to repeated parameters in
  // the transient case
  exit_on_repeated_greedy_parameters = false;
}



TransientRBConstruction::~TransientRBConstruction ()
{
  this->clear();
}


void TransientRBConstruction::clear()
{
  Parent::clear();

  // clear the mass matrices
  for (std::size_t q=0; q<M_q_vector.size(); q++)
    {
      delete M_q_vector[q];
      M_q_vector[q] = libmesh_nullptr;
    }

  if(store_non_dirichlet_operators)
    {
      for (std::size_t q=0; q<non_dirichlet_M_q_vector.size(); q++)
        {
          delete non_dirichlet_M_q_vector[q];
          non_dirichlet_M_q_vector[q] = libmesh_nullptr;
        }
    }

  // clear the temporal_data
  for (std::size_t i=0; i<temporal_data.size(); i++)
    {
      if(temporal_data[i])
        {
          temporal_data[i]->clear();
          delete temporal_data[i];
          temporal_data[i] = libmesh_nullptr;
        }
    }
  temporal_data.resize(0);
}

void TransientRBConstruction::initialize_rb_construction(bool skip_matrix_assembly,
                                                         bool skip_vector_assembly)
{
  // Check that the theta and assembly objects are consistently sized
#ifndef NDEBUG
  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  TransientRBAssemblyExpansion & trans_assembly_expansion =
    cast_ref<TransientRBAssemblyExpansion &>(get_rb_assembly_expansion());
#endif
  // This assert only gets called if DEBUG is on
  libmesh_assert_equal_to (trans_theta_expansion.get_n_M_terms(), trans_assembly_expansion.get_n_M_terms());

  Parent::initialize_rb_construction(skip_matrix_assembly, skip_vector_assembly);
}

void TransientRBConstruction::process_parameters_file (const std::string & parameters_filename)
{
  Parent::process_parameters_file(parameters_filename);

  // Read in data from parameters_filename
  GetPot infile(parameters_filename);

  // Read in the generic temporal discretization data
  process_temporal_parameters_file(parameters_filename);

  // Read in the data specific to Construction
  nonzero_initialization = infile("nonzero_initialization",nonzero_initialization);
  init_filename = infile("init_filename",init_filename);

  const Real POD_tol_in         = infile("POD_tol", POD_tol);
  const int max_truth_solves_in = infile("max_truth_solves", max_truth_solves);
  const unsigned int delta_N_in = infile("delta_N", delta_N);

  set_POD_tol(POD_tol_in);
  set_max_truth_solves(max_truth_solves_in);
  set_delta_N(delta_N_in);

  // Pass the temporal discretization data to the RBEvaluation
  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());
  trans_rb_eval.pull_temporal_discretization_data( *this );
}

void TransientRBConstruction::print_info()
{
  Parent::print_info();

  libMesh::out << std::endl << "TransientRBConstruction parameters:" << std::endl;

  if( is_rb_eval_initialized() )
    {
      // Print out info that describes the current setup
      TransientRBThetaExpansion & trans_theta_expansion =
        cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
      libMesh::out << "Q_m: " << trans_theta_expansion.get_n_M_terms() << std::endl;
    }
  else
    {
      libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
    }
  libMesh::out << "Number of time-steps: " << get_n_time_steps() << std::endl;
  libMesh::out << "dt: " << get_delta_t() << std::endl;
  libMesh::out << "euler_theta (time discretization parameter): " << get_euler_theta() << std::endl;
  if(get_POD_tol() > 0.)
    libMesh::out << "POD_tol: " << get_POD_tol() << std::endl;
  if(max_truth_solves > 0)
    libMesh::out << "Maximum number of truth solves: " << max_truth_solves << std::endl;
  libMesh::out << "delta_N (number of basis functions to add each POD-Greedy step): " << get_delta_N() << std::endl;
  if(nonzero_initialization)
    {
      libMesh::out << "Reading initial condition from " << init_filename << std::endl;
    }
  else
    {
      libMesh::out << "Using zero initial condition" << std::endl;
    }
  libMesh::out << std::endl;
}

void TransientRBConstruction::allocate_data_structures()
{
  Parent::allocate_data_structures();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m       = trans_theta_expansion.get_n_M_terms();
  const unsigned int n_outputs = trans_theta_expansion.get_n_outputs();

  // Resize and allocate vectors for storing mesh-dependent data
  const unsigned int n_time_levels = get_n_time_steps()+1;
  temporal_data.resize(n_time_levels);

  // Resize vectors for storing mesh-dependent data but only
  // initialize if initialize_mesh_dependent_data == true
  M_q_vector.resize(Q_m);

  // Only initialize the mass matrices if we
  // are not in single-matrix mode
  {
    DofMap & dof_map = this->get_dof_map();

    dof_map.attach_matrix(*L2_matrix);
    L2_matrix->init();
    L2_matrix->zero();

    for(unsigned int q=0; q<Q_m; q++)
      {
        // Initialize the memory for the matrices
        M_q_vector[q] = SparseMatrix<Number>::build(this->comm()).release();
        dof_map.attach_matrix(*M_q_vector[q]);
        M_q_vector[q]->init();
        M_q_vector[q]->zero();
      }

    // We also need to initialize a second set of non-Dirichlet operators
    if(store_non_dirichlet_operators)
      {
        dof_map.attach_matrix(*non_dirichlet_L2_matrix);
        non_dirichlet_L2_matrix->init();
        non_dirichlet_L2_matrix->zero();

        non_dirichlet_M_q_vector.resize(Q_m);
        for(unsigned int q=0; q<Q_m; q++)
          {
            // Initialize the memory for the matrices
            non_dirichlet_M_q_vector[q] = SparseMatrix<Number>::build(this->comm()).release();
            dof_map.attach_matrix(*non_dirichlet_M_q_vector[q]);
            non_dirichlet_M_q_vector[q]->init();
            non_dirichlet_M_q_vector[q]->zero();
          }
      }
  }

  for(unsigned int i=0; i<n_time_levels; i++)
    {
      temporal_data[i] = (NumericVector<Number>::build(this->comm()).release());
      temporal_data[i]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
    }

  // and the truth output vectors
  truth_outputs_all_k.resize(n_outputs);
  for(unsigned int n=0; n<n_outputs; n++)
    {
      truth_outputs_all_k[n].resize(n_time_levels);
    }

  // This vector is for storing rhs entries for
  // computing the projection of the initial condition
  // into the RB space
  RB_ic_proj_rhs_all_N.resize(Nmax);
}

void TransientRBConstruction::assemble_affine_expansion(bool skip_matrix_assembly,
                                                        bool skip_vector_assembly)
{
  // Call parent's assembly functions
  Parent::assemble_affine_expansion(skip_matrix_assembly, skip_vector_assembly);

  // Now update RB_ic_proj_rhs_all_N if necessary.
  // This allows us to compute the L2 projection
  // of the initial condition into the RB space
  // so that we can continue to enrich a given RB
  // space.
  if(get_rb_evaluation().get_n_basis_functions() > 0)
    {
      // Load the initial condition into the solution vector
      initialize_truth();

      UniquePtr< NumericVector<Number> > temp1 = NumericVector<Number>::build(this->comm());
      temp1->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

      // First compute the right-hand side vector for the L2 projection
      L2_matrix->vector_mult(*temp1, *solution);

      for(unsigned int i=0; i<get_rb_evaluation().get_n_basis_functions(); i++)
        {
          RB_ic_proj_rhs_all_N(i) = temp1->dot(get_rb_evaluation().get_basis_function(i));
        }
    }
}

Real TransientRBConstruction::train_reduced_basis(const bool resize_rb_eval_data)
{
  compute_truth_projection_error = true;
  Real value = Parent::train_reduced_basis(resize_rb_eval_data);
  compute_truth_projection_error = false;

  return value;
}

SparseMatrix<Number> * TransientRBConstruction::get_M_q(unsigned int q)
{
  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  if(q >= trans_theta_expansion.get_n_M_terms())
    libmesh_error_msg("Error: We must have q < Q_m in get_M_q.");

  return M_q_vector[q];
}

SparseMatrix<Number> * TransientRBConstruction::get_non_dirichlet_M_q(unsigned int q)
{
  if(!store_non_dirichlet_operators)
    libmesh_error_msg("Error: Must have store_non_dirichlet_operators==true to access non_dirichlet_M_q.");

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  if(q >= trans_theta_expansion.get_n_M_terms())
    libmesh_error_msg("Error: We must have q < Q_m in get_M_q.");

  return non_dirichlet_M_q_vector[q];
}

void TransientRBConstruction::get_all_matrices(std::map<std::string, SparseMatrix<Number> *> & all_matrices)
{
  Parent::get_all_matrices(all_matrices);

  all_matrices["L2_matrix"] = L2_matrix.get();

  if (store_non_dirichlet_operators)
    all_matrices["L2_matrix_non_dirichlet"] = non_dirichlet_L2_matrix.get();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();

  for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      std::stringstream matrix_name;
      matrix_name << "M" << q_m;
      all_matrices[matrix_name.str()] = get_M_q(q_m);

      if (store_non_dirichlet_operators)
        {
          matrix_name << "_non_dirichlet";
          all_matrices[matrix_name.str()] = get_non_dirichlet_M_q(q_m);
        }
    }
}

void TransientRBConstruction::assemble_L2_matrix(SparseMatrix<Number> * input_matrix, bool apply_dirichlet_bc)
{
  input_matrix->zero();
  add_scaled_matrix_and_vector(1.,
                               L2_assembly,
                               input_matrix,
                               libmesh_nullptr,
                               false, /* symmetrize */
                               apply_dirichlet_bc);
}

void TransientRBConstruction::assemble_mass_matrix(SparseMatrix<Number> * input_matrix)
{
  input_matrix->zero();
  add_scaled_mass_matrix(1., input_matrix);
}

void TransientRBConstruction::add_scaled_mass_matrix(Number scalar, SparseMatrix<Number> * input_matrix)
{
  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();

  for(unsigned int q=0; q<Q_m; q++)
    input_matrix->add(scalar * trans_theta_expansion.eval_M_theta(q,mu), *get_M_q(q));
}

void TransientRBConstruction::mass_matrix_scaled_matvec(Number scalar,
                                                        NumericVector<Number> & dest,
                                                        NumericVector<Number> & arg)
{
  LOG_SCOPE("mass_matrix_scaled_matvec()", "TransientRBConstruction");

  dest.zero();

  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();

  UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
  temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  for(unsigned int q=0; q<Q_m; q++)
    {
      get_M_q(q)->vector_mult(*temp_vec, arg);
      dest.add(scalar * trans_theta_expansion.eval_M_theta(q,mu), *temp_vec);
    }
}

void TransientRBConstruction::truth_assembly()
{
  LOG_SCOPE("truth_assembly()", "TransientRBConstruction");

  this->matrix->close();

  this->matrix->zero();
  this->rhs->zero();

  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  const Real dt          = get_delta_t();
  const Real euler_theta = get_euler_theta();

  {
    // We should have already assembled the matrices
    // and vectors in the affine expansion, so
    // just use them

    add_scaled_mass_matrix(1./dt, matrix);
    mass_matrix_scaled_matvec(1./dt, *rhs, *current_local_solution);

    UniquePtr< NumericVector<Number> > temp_vec = NumericVector<Number>::build(this->comm());
    temp_vec->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

    for(unsigned int q_a=0; q_a<Q_a; q_a++)
      {
        matrix->add(euler_theta*trans_theta_expansion.eval_A_theta(q_a,mu), *get_Aq(q_a));

        get_Aq(q_a)->vector_mult(*temp_vec, *current_local_solution);
        temp_vec->scale( -(1.-euler_theta)*trans_theta_expansion.eval_A_theta(q_a,mu) );
        rhs->add(*temp_vec);
      }

    for(unsigned int q_f=0; q_f<Q_f; q_f++)
      {
        *temp_vec = *get_Fq(q_f);
        temp_vec->scale( get_control(get_time_step())*trans_theta_expansion.eval_F_theta(q_f,mu) );
        rhs->add(*temp_vec);
      }

  }

  this->matrix->close();
  this->rhs->close();
}

void TransientRBConstruction::set_L2_assembly(ElemAssembly & L2_assembly_in)
{
  L2_assembly = &L2_assembly_in;
}

ElemAssembly & TransientRBConstruction::get_L2_assembly()
{
  if(!L2_assembly)
    libmesh_error_msg("Error: L2_assembly hasn't been initialized yet");

  return *L2_assembly;
}

void TransientRBConstruction::assemble_Mq_matrix(unsigned int q, SparseMatrix<Number> * input_matrix, bool apply_dirichlet_bc)
{
  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  TransientRBAssemblyExpansion & trans_assembly_expansion =
    cast_ref<TransientRBAssemblyExpansion &>(get_rb_assembly_expansion());

  if(q >= trans_theta_expansion.get_n_M_terms())
    libmesh_error_msg("Error: We must have q < Q_m in assemble_Mq_matrix.");

  input_matrix->zero();
  add_scaled_matrix_and_vector(1.,
                               &trans_assembly_expansion.get_M_assembly(q),
                               input_matrix,
                               libmesh_nullptr,
                               false, /* symmetrize */
                               apply_dirichlet_bc);
}

void TransientRBConstruction::assemble_all_affine_operators()
{
  Parent::assemble_all_affine_operators();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  for(unsigned int q=0; q<trans_theta_expansion.get_n_M_terms(); q++)
    assemble_Mq_matrix(q, get_M_q(q));

  if(store_non_dirichlet_operators)
    {
      for(unsigned int q=0; q<trans_theta_expansion.get_n_M_terms(); q++)
        assemble_Mq_matrix(q, get_non_dirichlet_M_q(q), false);
    }
}

void TransientRBConstruction::assemble_misc_matrices()
{
  libMesh::out << "Assembling L2 matrix" << std::endl;
  assemble_L2_matrix(L2_matrix.get());

  if(store_non_dirichlet_operators)
    {
      libMesh::out << "Assembling non-Dirichlet L2 matrix" << std::endl;
      assemble_L2_matrix(non_dirichlet_L2_matrix.get(), /* apply_dirichlet_bc = */ false);
    }

  Parent::assemble_misc_matrices();
}

Real TransientRBConstruction::truth_solve(int write_interval)
{
  LOG_SCOPE("truth_solve()", "TransientRBConstruction");

  const RBParameters & mu = get_parameters();
  const unsigned int n_time_steps = get_n_time_steps();

  //   // NumericVector for computing true L2 error
  //   UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build();
  //   temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  // Apply initial condition again.
  initialize_truth();
  set_time_step(0);

  // Now compute the truth outputs
  for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
    {
      truth_outputs_all_k[n][0] = 0.;
      for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
        {
          truth_outputs_all_k[n][0] += get_rb_theta_expansion().eval_output_theta(n,q_l,mu)*
            get_output_vector(n,q_l)->dot(*solution);
        }
    }

  // Load initial projection error into temporal_data dense matrix
  if(compute_truth_projection_error)
    set_error_temporal_data();

  for(unsigned int time_level=1; time_level<=n_time_steps; time_level++)
    {
      set_time_step(time_level);

      *old_local_solution = *current_local_solution;

      // We assume that the truth assembly has been attached to the system
      truth_assembly();

      // truth_assembly assembles into matrix and rhs, so use those for the solve
      solve_for_matrix_and_rhs(*get_linear_solver(), *matrix, *rhs);

      // The matrix doesn't change at each timestep, so we
      // can set reuse_preconditioner == true
      linear_solver->reuse_preconditioner(true);

      if (assert_convergence)
        {
          check_convergence(*get_linear_solver());
        }

      // Now compute the truth outputs
      for(unsigned int n=0; n<get_rb_theta_expansion().get_n_outputs(); n++)
        {
          truth_outputs_all_k[n][time_level] = 0.;
          for(unsigned int q_l=0; q_l<get_rb_theta_expansion().get_n_output_terms(n); q_l++)
            {
              truth_outputs_all_k[n][time_level] +=
                get_rb_theta_expansion().eval_output_theta(n,q_l,mu)*get_output_vector(n,q_l)->dot(*solution);
            }
        }

      // load projection error into column _k of temporal_data matrix
      if(compute_truth_projection_error)
        set_error_temporal_data();

      if ( (write_interval > 0) && (time_level%write_interval == 0) )
        {
          libMesh::out << std::endl << "Truth solve, plotting time step " << time_level << std::endl;

          std::ostringstream file_name;

          file_name << "truth.e.";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << time_level;

#ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(get_mesh()).write_equation_systems (file_name.str(),
                                                          this->get_equation_systems());
#endif
        }
    }

  // Set reuse_preconditioner back to false for subsequent solves.
  linear_solver->reuse_preconditioner(false);

  // Get the L2 norm of the truth solution at time-level _K
  // Useful for normalizing our true error data
  L2_matrix->vector_mult(*inner_product_storage_vector, *solution);
  Real final_truth_L2_norm = libmesh_real(std::sqrt(inner_product_storage_vector->dot(*solution)));


  return final_truth_L2_norm;
}

bool
TransientRBConstruction::greedy_termination_test(Real abs_greedy_error,
                                                 Real initial_greedy_error,
                                                 int count)
{
  if ( (get_max_truth_solves()>0) && (count >= get_max_truth_solves()) )
    {
      libMesh::out << "Maximum number of truth solves reached: max = "
                   << count << std::endl;
      return true;
    }

  return Parent::greedy_termination_test(abs_greedy_error, initial_greedy_error, count);
}

Number TransientRBConstruction::set_error_temporal_data()
{
  LOG_SCOPE("set_error_temporal_data()", "TransientRBConstruction");

  // first compute the projection of solution onto the current
  // RB space

  const unsigned int time_step = get_time_step();

  if(get_rb_evaluation().get_n_basis_functions() == 0)
    {
      // If the basis is empty, then the error is the solution itself
      temporal_data[time_step]->zero();
      temporal_data[time_step]->add(1., *solution);
    }
  else
    {
      unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

      UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
      temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

      // First compute the right-hand side vector for the projection
      inner_product_matrix->vector_mult(*temp, *solution);

      //    zero_dirichlet_dofs_on_vector(*temp);

      // Do not assume that RB_stiffness matrix is diagonal,
      // diagonality degrades as N increases

      // Get an appropriately sized copy of RB_inner_product_matrix
      DenseMatrix<Number> RB_inner_product_matrix_N(RB_size,RB_size);
      for(unsigned int i=0; i<RB_size; i++)
        for(unsigned int j=0; j<RB_size; j++)
          {
            RB_inner_product_matrix_N(i,j) = get_rb_evaluation().RB_inner_product_matrix(i,j);
          }

      // Compute the projection RHS
      DenseVector<Number> RB_proj_rhs(RB_size);
      for(unsigned int i=0; i<RB_size; i++)
        {
          RB_proj_rhs(i) = temp->dot(get_rb_evaluation().get_basis_function(i));
        }

      DenseVector<Number> RB_proj(RB_size);

      // Now solve the linear system
      RB_inner_product_matrix_N.lu_solve(RB_proj_rhs, RB_proj);

      // Load the RB projection into temp
      temp->zero();
      for(unsigned int i=0; i<RB_size; i++)
        {
          temp->add(RB_proj(i), get_rb_evaluation().get_basis_function(i));
        }

      temp->add(-1., *solution);

      // Now temp holds the projection error, store in temporal_data
      *(temporal_data[time_step]) = *temp;
    }

  // return the square of the X norm of the truth solution
  inner_product_matrix->vector_mult(*inner_product_storage_vector,*solution);

  return solution->dot(*inner_product_storage_vector);
}

const NumericVector<Number> & TransientRBConstruction::get_error_temporal_data()
{
  LOG_SCOPE("get_error_temporal_data()", "TransientRBConstruction");

  const unsigned int time_step = get_time_step();

  return *temporal_data[time_step];
}

void TransientRBConstruction::initialize_truth ()
{
  if (nonzero_initialization)
    {
      // Use System::read_serialized_data to read the initial condition
      // into this->solution
      Xdr IC_data(init_filename, READ);
      read_serialized_data(IC_data, false);
    }
  else
    {
      // Otherwise zero out the solution as a default
      this->solution->zero();
    }
  this->solution->close();
  this->update();
}

void TransientRBConstruction::add_IC_to_RB_space()
{
  LOG_SCOPE("add_IC_to_RB_space()", "TransientRBConstruction");

  if (get_rb_evaluation().get_n_basis_functions() > 0)
    libmesh_error_msg("Error: Should not call TransientRBConstruction::add_IC_to_RB_space() " \
                      << "on a system that already contains basis functions.");

  if (!nonzero_initialization)
    libmesh_error_msg("Error: Should not call TransientRBConstruction::add_IC_to_RB_space() " \
                      << "when nonzero_initialization==false.");

  initialize_truth();

  // load the new basis function into the basis_functions vector.
  get_rb_evaluation().basis_functions.push_back( NumericVector<Number>::build(this->comm()).release() );
  NumericVector<Number> & current_bf = get_rb_evaluation().get_basis_function(get_rb_evaluation().get_n_basis_functions()-1);
  current_bf.init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
  current_bf = *solution;

  // We can just set the norm to 1.
  inner_product_matrix->vector_mult(*inner_product_storage_vector,*solution);

  Real current_bf_norm = libmesh_real(std::sqrt( current_bf.dot(*inner_product_storage_vector) ));
  current_bf.scale(1./current_bf_norm);

  unsigned int saved_delta_N = get_delta_N();
  set_delta_N(1);
  update_system();
  set_delta_N(saved_delta_N);
}

void TransientRBConstruction::enrich_RB_space()
{
  // Need SLEPc to get the POD eigenvalues
#if defined(LIBMESH_HAVE_SLEPC)
  LOG_SCOPE("enrich_RB_space()", "TransientRBConstruction");

  // With the "method of snapshots", the size of
  // the eigenproblem is determined by the number
  // of time-steps (rather than the number of spatial dofs).
  PetscBLASInt eigen_size = temporal_data.size();
  PetscBLASInt LDA = eigen_size; // The leading order of correlation_matrix
  std::vector<Number> correlation_matrix(LDA*eigen_size);

  for(int i=0; i<eigen_size; i++)
    {
      inner_product_matrix->vector_mult(*inner_product_storage_vector, *temporal_data[i]);

      for(int j=i; j<eigen_size; j++)
        {
          // Scale the inner products by the number of time-steps to normalize the
          // POD energy norm appropriately
          Number inner_prod = (temporal_data[j]->dot(*inner_product_storage_vector)) /
            (Real)(get_n_time_steps()+1);

          // Fill upper triangular part of correlation_matrix
          correlation_matrix[j*eigen_size+i] = inner_prod;
        }
    }

  // Call the LAPACK eigensolver
  char JOBZ = 'V'; // Compute eigenvectors and eigenvalues
  char RANGE = 'A'; // Compute all eigenvalues
  char UPLO = 'U'; // Upper triangular symmetric matrix

  Real VL = 0.; // Not used when RANGE = A
  Real VU = 0.; // Not used when RANGE = A

  PetscBLASInt IL = 0; // Not used when RANGE = A
  PetscBLASInt IU = 0; // Not used when RANGE = A

  Real ABSTOL = 1.e-14; // Absolute tolerance for eigensolver

  PetscBLASInt M = 0; // (output) The total number of evals found

  std::vector<Real> W(eigen_size); // (output) the eigenvalues

  PetscBLASInt LDZ = eigen_size; // The leading order of Z
  std::vector<Number> Z(LDZ*eigen_size); // (output) the eigenvectors

  std::vector<PetscBLASInt> ISUPPZ(2*eigen_size); // Indicates which evecs in Z are nonzero

  // Work array, sized according to lapack documentation
  PetscBLASInt LWORK = 26*eigen_size;
  std::vector<Number> WORK(LWORK);

  // Work array, sized according to lapack documentation
  PetscBLASInt LIWORK = 10*eigen_size;
  std::vector<PetscBLASInt> IWORK(LIWORK);

  PetscBLASInt INFO = 0;

#ifdef LIBMESH_USE_REAL_NUMBERS // Use real numbers
  // Call the eigensolver for symmetric eigenvalue problems.
  // NOTE: evals in W are in ascending order
  LAPACKsyevr_(&JOBZ, &RANGE, &UPLO, &eigen_size, &correlation_matrix[0],
               &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ,
               &ISUPPZ[0], &WORK[0], &LWORK, &IWORK[0], &LIWORK, &INFO );
#elif LIBMESH_USE_COMPLEX_NUMBERS // Use complex numbers
  // Need some extra data in the complex case

  // Work array, sized according to lapack documentation
  PetscBLASInt LRWORK = 24*eigen_size;
  std::vector<Real> RWORK(LRWORK);

  // Call the eigensolver for symmetric eigenvalue problems.
  // NOTE: evals in W are in ascending order
  LAPACKsyevr_(&JOBZ, &RANGE, &UPLO, &eigen_size, &correlation_matrix[0],
               &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ,
               &ISUPPZ[0], &WORK[0], &LWORK, &RWORK[0], &LRWORK, &IWORK[0],
               &LIWORK, &INFO );
#else
#error libMesh does not yet support quaternions!
#endif

  if (INFO != 0)
    libmesh_error_msg("Error in LAPACK syev eigensolver routine, INFO = " << INFO);

  // eval and evec now hold the sorted eigenvalues/eigenvectors
  libMesh::out << std::endl << "POD Eigenvalues:" << std::endl;
  for(unsigned int i=0; i<=1; i++)
    {
      libMesh::out << "eigenvalue " << i << " = " << W[eigen_size-1-i] << std::endl;
    }
  libMesh::out << "..." << std::endl;// << "." << std::endl << "." << std::endl;
  libMesh::out << "last eigenvalue = " << W[0] << std::endl;
  libMesh::out << std::endl;

  // Now load the new basis functions
  unsigned int count = 0;
  while (true)
    {
      // load the new basis function into the basis_functions vector.
      get_rb_evaluation().basis_functions.push_back( NumericVector<Number>::build(this->comm()).release() );
      NumericVector<Number> & current_bf = get_rb_evaluation().get_basis_function(get_rb_evaluation().get_n_basis_functions()-1);
      current_bf.init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      current_bf.zero();

      // Perform the matrix multiplication of temporal data with
      // the next POD eigenvector
      for(int j=0; j<eigen_size; j++)
        {
          int Z_row = j;
          int Z_col = (eigen_size-1-count);
          int Z_index = Z_col*eigen_size + Z_row;
          current_bf.add(Z[Z_index], *temporal_data[j]);
        }

      // We just set the norm to 1.
      inner_product_matrix->vector_mult(*inner_product_storage_vector,current_bf);

      Real current_bf_norm = std::abs( std::sqrt( current_bf.dot(*inner_product_storage_vector) ) );
      current_bf.scale(1./current_bf_norm);

      // Increment count here since we use the incremented counter
      // in the if clauses below
      count++;

      // If positive POD_tol, we use it to determine the number of basis functions
      // to add, and then break the loop when POD_tol is satisfied, or after Nmax
      // basis functions have been added. Else we break the loop after delta_N
      // (or Nmax) new basis functions.
      if (POD_tol > 0.)
        {
          set_delta_N(1);

          // We need to define the updated RB system matrices before the RB solve
          update_system();
          Real error_bound = get_rb_evaluation().rb_solve(get_rb_evaluation().get_n_basis_functions());

          if ( (error_bound <= POD_tol) || (get_rb_evaluation().get_n_basis_functions()==get_Nmax()) )
            {
              set_delta_N(0);
              break;
            }
        }
      else
        {
          if (count == get_delta_N())
            {
              break;
            }
          else
            if (get_rb_evaluation().get_n_basis_functions()==get_Nmax())
              {
                set_delta_N(count);
                break;
              }
        }
    }

#else
  libmesh_not_implemented();
#endif
}


void TransientRBConstruction::update_system()
{
  // If delta_N is set to zero, there is nothing to update
  if(get_delta_N() == 0)
    return;

  Parent::update_system();

  libMesh::out << "Updating RB initial conditions" << std::endl;
  update_RB_initial_condition_all_N();
}

SparseMatrix<Number> & TransientRBConstruction::get_matrix_for_output_dual_solves()
{
  return *L2_matrix;
}

void TransientRBConstruction::load_rb_solution()
{
  LOG_SCOPE("load_rb_solution()", "TransientRBConstruction");

  solution->zero();

  const unsigned int time_step = get_time_step();

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());
  DenseVector<Number> RB_solution_vector_k = trans_rb_eval.RB_temporal_solution_data[time_step];

  if(RB_solution_vector_k.size() > get_rb_evaluation().get_n_basis_functions())
    libmesh_error_msg("ERROR: rb_eval object contains " \
                      << get_rb_evaluation().get_n_basis_functions() \
                      << " basis functions. RB_solution vector constains " \
                      << RB_solution_vector_k.size() \
                      << " entries. RB_solution in TransientRBConstruction::load_rb_solution is too long!");

  for (unsigned int i=0; i<RB_solution_vector_k.size(); i++)
    solution->add(RB_solution_vector_k(i), get_rb_evaluation().get_basis_function(i));

  update();
}

void TransientRBConstruction::update_RB_system_matrices()
{
  LOG_SCOPE("update_RB_system_matrices()", "TransientRBConstruction");

  Parent::update_RB_system_matrices();

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  UniquePtr< NumericVector<Number> > temp = NumericVector<Number>::build(this->comm());
  temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      for(unsigned int j=0; j<RB_size; j++)
        {
          // Compute reduced L2 matrix
          temp->zero();
          L2_matrix->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

          Number value = get_rb_evaluation().get_basis_function(i).dot(*temp);
          trans_rb_eval.RB_L2_matrix(i,j) = value;
          if(i!=j)
            {
              // The L2 matrix is assumed
              // to be symmetric
              trans_rb_eval.RB_L2_matrix(j,i) = value;
            }

          for(unsigned int q_m=0; q_m<Q_m; q_m++)
            {
              // Compute reduced M_q matrix
              temp->zero();
              get_M_q(q_m)->vector_mult(*temp, get_rb_evaluation().get_basis_function(j));

              value = (get_rb_evaluation().get_basis_function(i)).dot(*temp);
              trans_rb_eval.RB_M_q_vector[q_m](i,j) = value;

              if(i!=j)
                {
                  // Each mass matrix term is assumed
                  // to be symmetric
                  trans_rb_eval.RB_M_q_vector[q_m](j,i) = value;
                }
            }

        }
    }
}




void TransientRBConstruction::update_residual_terms(bool compute_inner_products)
{
  LOG_SCOPE("update_residual_terms()", "TransientRBConstruction");

  Parent::update_residual_terms(compute_inner_products);

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());

  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          // Initialize the vectors when we need them
          if(!trans_rb_eval.M_q_representor[q_m][i])
            {
              trans_rb_eval.M_q_representor[q_m][i] = (NumericVector<Number>::build(this->comm()).release());
              trans_rb_eval.M_q_representor[q_m][i]->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
            }

          libmesh_assert(trans_rb_eval.M_q_representor[q_m][i]->size()       == this->n_dofs()       &&
                         trans_rb_eval.M_q_representor[q_m][i]->local_size() == this->n_local_dofs() );

          rhs->zero();
          M_q_vector[q_m]->vector_mult(*rhs, get_rb_evaluation().get_basis_function(i));

          if (!is_quiet())
            libMesh::out << "Starting solve i="
                         << i << " in TransientRBConstruction::update_residual_terms() at "
                         << Utility::get_timestamp() << std::endl;

          solve_for_matrix_and_rhs(*inner_product_solver, *inner_product_matrix, *rhs);

          if (assert_convergence)
            check_convergence(*inner_product_solver);

          if (!is_quiet())
            {
              libMesh::out << "Finished solve i="
                           << i << " in TransientRBConstruction::update_residual_terms() at "
                           << Utility::get_timestamp() << std::endl;

              libMesh::out << this->n_linear_iterations()
                           << " iterations, final residual "
                           << this->final_linear_residual() << std::endl;
            }

          *trans_rb_eval.M_q_representor[q_m][i] = *solution;
        }
    }

  // Now compute and store the inner products if requested
  if (compute_inner_products)
    {
      for(unsigned int q_f=0; q_f<Q_f; q_f++)
        {
          inner_product_matrix->vector_mult(*inner_product_storage_vector, *Fq_representor[q_f]);

          for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
            {
              for(unsigned int q_m=0; q_m<Q_m; q_m++)
                {
                  trans_rb_eval.Fq_Mq_representor_innerprods[q_f][q_m][i] =
                    trans_rb_eval.M_q_representor[q_m][i]->dot(*inner_product_storage_vector);
                } // end for q_m
            } // end for i
        } // end for q_f

      unsigned int q=0;
      for(unsigned int q_m1=0; q_m1<Q_m; q_m1++)
        {
          for(unsigned int q_m2=q_m1; q_m2<Q_m; q_m2++)
            {
              for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
                {
                  for(unsigned int j=0; j<RB_size; j++)
                    {
                      inner_product_matrix->vector_mult(*inner_product_storage_vector, *trans_rb_eval.M_q_representor[q_m2][j]);

                      trans_rb_eval.Mq_Mq_representor_innerprods[q][i][j] =
                        trans_rb_eval.M_q_representor[q_m1][i]->dot(*inner_product_storage_vector);

                      if(i != j)
                        {
                          inner_product_matrix->vector_mult(*inner_product_storage_vector,
                                                            *trans_rb_eval.M_q_representor[q_m2][i]);

                          trans_rb_eval.Mq_Mq_representor_innerprods[q][j][i] =
                            trans_rb_eval.M_q_representor[q_m1][j]->dot(*inner_product_storage_vector);
                        }
                    } // end for j
                } // end for i
              q++;
            } // end for q_m2
        } // end for q_m1


      for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
        {
          for(unsigned int j=0; j<RB_size; j++)
            {
              for(unsigned int q_a=0; q_a<Q_a; q_a++)
                {
                  for(unsigned int q_m=0; q_m<Q_m; q_m++)
                    {
                      inner_product_matrix->vector_mult(*inner_product_storage_vector,
                                                        *trans_rb_eval.M_q_representor[q_m][j]);

                      trans_rb_eval.Aq_Mq_representor_innerprods[q_a][q_m][i][j] =
                        trans_rb_eval.Aq_representor[q_a][i]->dot(*inner_product_storage_vector);

                      if(i != j)
                        {
                          inner_product_matrix->vector_mult(*inner_product_storage_vector,
                                                            *trans_rb_eval.M_q_representor[q_m][i]);

                          trans_rb_eval.Aq_Mq_representor_innerprods[q_a][q_m][j][i] =
                            trans_rb_eval.Aq_representor[q_a][j]->dot(*inner_product_storage_vector);
                        }
                    } // end for q_m
                } // end for q_a
            } // end for j
        } // end for i
    } // end if (compute_inner_products)
}


void TransientRBConstruction::update_RB_initial_condition_all_N()
{
  LOG_SCOPE("update_RB_initial_condition_all_N()", "TransientRBConstruction");

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());

  // Load the initial condition into the solution vector
  initialize_truth();

  UniquePtr< NumericVector<Number> > temp1 = NumericVector<Number>::build(this->comm());
  temp1->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  UniquePtr< NumericVector<Number> > temp2 = NumericVector<Number>::build(this->comm());
  temp2->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);


  unsigned int RB_size = get_rb_evaluation().get_n_basis_functions();

  // First compute the right-hand side vector for the L2 projection
  L2_matrix->vector_mult(*temp1, *solution);

  for(unsigned int i=(RB_size-delta_N); i<RB_size; i++)
    {
      RB_ic_proj_rhs_all_N(i) = temp1->dot(get_rb_evaluation().get_basis_function(i));
    }


  // Now compute the projection for each N
  DenseMatrix<Number> RB_L2_matrix_N;
  DenseVector<Number> RB_rhs_N;
  for(unsigned int N=(RB_size-delta_N); N<RB_size; N++)
    {
      // We have to index here by N+1 since the loop index is zero-based.
      trans_rb_eval.RB_L2_matrix.get_principal_submatrix(N+1, RB_L2_matrix_N);

      RB_ic_proj_rhs_all_N.get_principal_subvector(N+1, RB_rhs_N);

      DenseVector<Number> RB_ic_N(N+1);

      // Now solve the linear system
      RB_L2_matrix_N.lu_solve(RB_rhs_N, RB_ic_N);

      // Load RB_ic_N into RB_initial_condition_all_N
      trans_rb_eval.RB_initial_condition_all_N[N] = RB_ic_N;

      // Compute the L2 error for the RB initial condition
      // This part is dependent on the truth space.

      // load the RB solution into temp1
      temp1->zero();
      for(unsigned int i=0; i<N+1; i++)
        {
          temp1->add(RB_ic_N(i), get_rb_evaluation().get_basis_function(i));
        }

      // subtract truth initial condition from RB_ic_N
      temp1->add(-1., *solution);

      // Compute L2 norm error, i.e. sqrt(M(solution,solution))
      temp2->zero();
      L2_matrix->vector_mult(*temp2, *temp1);

      trans_rb_eval.initial_L2_error_all_N[N] = libmesh_real(std::sqrt(temp2->dot(*temp1)));
    }
}

//Real TransientRBConstruction::uncached_compute_residual_dual_norm(const unsigned int N)
//{
//   LOG_SCOPE("uncached_compute_residual_dual_norm()", "TransientRBConstruction");
//
//   // This is the "slow" way of computing the residual, but it is useful
//   // for validating the "fast" way.
//   // Note that this only works in serial since otherwise each processor will
//   // have a different parameter value during the Greedy training.
//
//   // Assemble the right-hand side to find the Reisz representor
//   // of the residual in the X norm
//   UniquePtr< NumericVector<Number> > RB_sol = NumericVector<Number>::build();
//   RB_sol->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//   RB_sol->zero();
//
//   UniquePtr< NumericVector<Number> > ghosted_temp = NumericVector<Number>::build();
//   ghosted_temp->init (this->n_dofs(), this->n_local_dofs(),
//                       this->get_dof_map().get_send_list(), false,
//                       GHOSTED);
//
//   UniquePtr< NumericVector<Number> > parallel_temp = NumericVector<Number>::build();
//   parallel_temp->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
//
//   // Store current_local_solution, since we don't want to corrupt it
//   *ghosted_temp = *current_local_solution;
//
//   for(unsigned int i=0; i<N; i++)
//   {
//     RB_sol->add(RB_solution(i), rb_eval->get_basis_function(i));
//     parallel_temp->add(old_RB_solution(i), rb_eval->get_basis_function(i));
//   }
//
//   // Load parallel_temp into current_local_solution in order to do assembly
//   const std::vector<unsigned int> & send_list = this->get_dof_map().get_send_list();
//   parallel_temp->localize (*current_local_solution, send_list);
//
//   // Load the system_matrix
//   this->truth_assembly();
//
//   // Restore current_local_solution
//   *current_local_solution = *ghosted_temp;
//
//   matrix->vector_mult(*parallel_temp, *RB_sol);
//   rhs->add(-1., *parallel_temp);
//   rhs->close();
//
//   zero_dirichlet_dofs_on_rhs();
//
//   // Then solve the system to get the Reisz representor
//   matrix->zero();
//   {
//     matrix->add(1., *inner_product_matrix);
//     if(constrained_problem)
//       matrix->add(1., *constraint_matrix);
//   }
//
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
//                  << this->final_linear_residual() << ", number of iterations = "
//                  << this->n_linear_iterations());
//   }
//
//  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//
//  Real slow_residual_norm_sq = solution->dot(*inner_product_storage_vector);
//
//  return libmesh_real(std::sqrt( slow_residual_norm_sq ));
//}

void TransientRBConstruction::write_riesz_representors_to_files(const std::string & riesz_representors_dir,
                                                                const bool write_binary_residual_representors)
{
  LOG_SCOPE("write_riesz_representors_to_files()", "TransientRBConstruction");

  // Write out the M_q_representors.  These are useful to have when restarting,
  // so you don't have to recompute them all over again.  There should be
  // this->rb_eval->get_n_basis_functions() of these.
  libMesh::out << "Writing out the M_q_representors..." << std::endl;

  std::ostringstream file_name;
  const std::string riesz_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());

  const unsigned int istop  = trans_rb_eval.get_n_basis_functions();
  const unsigned int istart = istop-get_delta_N();

  for (std::size_t q=0; q<trans_rb_eval.M_q_representor.size(); ++q)
    for (unsigned int i=istart; i<istop; ++i)
      {
        libMesh::out << "Writing out M_q_representor[" << q << "][" << i << "]..." << std::endl;
        libmesh_assert(trans_rb_eval.M_q_representor[q][i]);

        file_name.str(""); // reset filename
        file_name << riesz_representors_dir << "/M_q_representor" << i << riesz_representor_suffix;

        {
          // No need to copy!
          //*solution = *(M_q_representor[q][i]);
          trans_rb_eval.M_q_representor[q][i]->swap(*solution);

          Xdr mr_data(file_name.str(),
                      write_binary_residual_representors ? ENCODE : WRITE);

          write_serialized_data(mr_data, false);

          // Synchronize before moving on
          this->comm().barrier();

          // Swap back.
          trans_rb_eval.M_q_representor[q][i]->swap(*solution);

          // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
          // for the system call, be sure to do it only on one processor, etc.
        }
      }
}

void TransientRBConstruction::read_riesz_representors_from_files(const std::string & riesz_representors_dir,
                                                                 const bool read_binary_residual_representors)
{
  LOG_SCOPE("read_riesz_representors_from_files()", "TransientRBConstruction");

  const std::string riesz_representor_suffix =
    (read_binary_residual_representors ? ".xdr" : ".dat");

  std::ostringstream file_name;
  struct stat stat_info;

  TransientRBEvaluation & trans_rb_eval = cast_ref<TransientRBEvaluation &>(get_rb_evaluation());

  libMesh::out << "Reading in the M_q_representors..." << std::endl;

  // Read in the Aq representors.  The class makes room for [Q_m][Nmax] of these.  We are going to
  // read in [Q_m][this->rb_eval->get_n_basis_functions()].  FIXME:
  // should we be worried about leaks in the locations where we're about to fill entries?
  for (std::size_t i=0; i<trans_rb_eval.M_q_representor.size(); ++i)
    for (std::size_t j=0; j<trans_rb_eval.M_q_representor[i].size(); ++j)
      {
        if (trans_rb_eval.M_q_representor[i][j] != libmesh_nullptr)
          libmesh_error_msg("Error, must delete existing M_q_representor before reading in from file.");
      }

  // Now ready to read them in from file!
  for (std::size_t i=0; i<trans_rb_eval.M_q_representor.size(); ++i)
    for (std::size_t j=0; j<trans_rb_eval.get_n_basis_functions(); ++j)
      {
        file_name.str(""); // reset filename
        file_name << riesz_representors_dir
                  << "/M_q_representor" << i << "_" << j << riesz_representor_suffix;

        // On processor zero check to be sure the file exists
        if (this->processor_id() == 0)
          {
            int stat_result = stat(file_name.str().c_str(), &stat_info);

            if (stat_result != 0)
              libmesh_error_msg("File does not exist: " << file_name.str());
          }

        Xdr aqr_data(file_name.str(),
                     read_binary_residual_representors ? DECODE : READ);

        read_serialized_data(aqr_data, false);

        trans_rb_eval.M_q_representor[i][j] = NumericVector<Number>::build(this->comm()).release();
        trans_rb_eval.M_q_representor[i][j]->init (n_dofs(), n_local_dofs(),
                                                   false, PARALLEL);

        // No need to copy, just swap
        //*M_q_representor[i][j] = *solution;
        trans_rb_eval.M_q_representor[i][j]->swap(*solution);
      }
}

} // namespace libMesh
