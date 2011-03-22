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


// Configuration data
#include "libmesh_config.h"

// This class requires QNTransientSCMSystem, which is not
// defined if SLEPc is not present.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "qn_transient_rb_system.h"
#include "qn_transient_rb_context.h"
#include "qn_transient_scm_system.h"
#include "qn_transient_rb_evaluation.h"

#include "equation_systems.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "gmv_io.h"
#include "linear_solver.h"
#include "getpot.h"
#include "timestamp.h"
#include "xdr_cxx.h"
#include "parallel.h"

// For checking for the existence of files
#include <sys/stat.h>

namespace libMesh
{

QNTransientRBSystem::QNTransientRBSystem (EquationSystems& es,
                                          const std::string& name,
                                          const unsigned int number)
  : Parent(es, name, number),
    current_newton_iterate(NumericVector<Number>::build()),
    use_nominal_rho_LB(true),
    nonlinear_tolerance(1.e-8),
    n_newton_steps(15),
    use_eisenstat(false),
    eisenstat_eta0(0.01),
    eisenstat_gamma(0.9),
    eisenstat_alpha(2.),
    C_assembly(NULL),
    truth_intrr_assembly(NULL),
    truth_bndry_assembly(NULL)
{
}

QNTransientRBSystem::~QNTransientRBSystem ()
{
  this->clear();
}

void QNTransientRBSystem::clear()
{
  Parent::clear();

  // Don't need to call clear_basis_function_dependent_data
  // since that just clears matrices/vectors, but we need
  // to free the memory.
  for(unsigned int n=0; n<C_n_vector.size(); n++)
  {
    if(C_n_vector[n])
    {
      delete C_n_vector[n];
      C_n_vector[n] = NULL;
    }
  }

  // Delete the Nmax^2 representors
  for(unsigned int i=0; i<C_representor.size(); i++)
  {
    for(unsigned int j=0; j<C_representor[i].size(); j++)
    {
      if(C_representor[i][j])
      {
        delete C_representor[i][j];
        C_representor[i][j] = NULL;
      }
    }
  }

}

void QNTransientRBSystem::clear_basis_function_dependent_data()
{
  // Clear the vector storing the C_n matrices
  for(unsigned int n=0; n<C_n_vector.size(); n++)
  {
    if(C_n_vector[n])
    {
      C_n_vector[n]->clear();
      C_n_vector[n]->init();
    }
  }

  // Delete the Nmax^2 representors
  for(unsigned int i=0; i<C_representor.size(); i++)
  {
    for(unsigned int j=0; j<C_representor[i].size(); j++)
    {
      if(C_representor[i][j])
      {
        // Need to NULL this out so that it gets reinitialized properly
        delete C_representor[i][j];
        C_representor[i][j] = NULL;
      }
    }
  }

  Parent::clear_basis_function_dependent_data();
}

void QNTransientRBSystem::process_parameters_file ()
{
  // Call the Parent's initialization routine.
  Parent::process_parameters_file();

  // Read in data from parameters_filename
  GetPot infile(parameters_filename);

  const Real nonlinear_tolerance_in    = infile("nonlinear_tolerance",
                                                nonlinear_tolerance);
  const unsigned int n_newton_steps_in = infile("n_newton_steps",
                                                n_newton_steps);

  set_nonlinear_tolerance(nonlinear_tolerance_in);
  set_n_newton_steps(n_newton_steps_in);

  // Read in and set Eisenstat formula parameters
  use_eisenstat = infile("use_eisenstat", false);
  eisenstat_eta0 = infile("eisenstat_eta0", 0.01);
  eisenstat_gamma = infile("eisenstat_gamma", 0.9);
  eisenstat_alpha = infile("eisenstat_alpha", 2.);
  
  // By default we enforce constraints exactly in the nonlinear case
  enforce_constraints_exactly = true;

  libMesh::out << std::endl << "QNTransientRBSystem parameters:" << std::endl;
  libMesh::out << "Nonlinear tolerance: " << nonlinear_tolerance << std::endl;
  libMesh::out << "Maximum number of Newton steps: " << n_newton_steps << std::endl;
  libMesh::out << std::endl;
}

void QNTransientRBSystem::add_new_rb_evaluation_object()
{
  QNTransientRBEvaluation* e = new QNTransientRBEvaluation(*this);
  rb_evaluation_objects.push_back(e);
}

void QNTransientRBSystem::allocate_data_structures()
{
  Parent::allocate_data_structures();

  // Resize C_n_vector
  C_n_vector.resize(Nmax);

  C_representor.resize(Nmax);
  for(unsigned int i=0; i<Nmax; i++)
  {
    C_representor[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      C_representor[i][j] = NULL;
    }
  }

  // Optionally initialize the mesh dependent data structures
  if(initialize_mesh_dependent_data)
  {
    DofMap& dof_map = this->get_dof_map();

    // Initialize current_newton_iterate
    current_newton_iterate->init (this->n_dofs(), this->n_local_dofs(),
                                  dof_map.get_send_list(), false,
                                  GHOSTED);

    // If we're not in low memory mode, initialize the memory for the C_n matrices
    if(!low_memory_mode)
      for(unsigned int n=0; n<Nmax; n++)
      {
        C_n_vector[n] = SparseMatrix<Number>::build().release();
        dof_map.attach_matrix(*C_n_vector[n]);

        // Initialize and zero the C_n matrices
        C_n_vector[n]->init();
      }
  }

}

AutoPtr<FEMContext> QNTransientRBSystem::build_context ()
{
  return AutoPtr<FEMContext>(new QNTransientRBContext(*this));
}

void QNTransientRBSystem::truth_assembly()
{
  START_LOG("truth_assembly()", "QNTransientRBSystem");

  // Here we directly assemble the Jacobian and residual for
  // the nonlinear system using user-specified assembly functions
  // Just call RBSystem::add_scaled_matrix_and_vector

  if( (truth_intrr_assembly == NULL) && (truth_bndry_assembly == NULL) )
  {
    libMesh::out << "Error: Truth assembly not attached to the system." << std::endl;
    libmesh_error();
  }

  this->matrix->zero();
  this->rhs->zero();

  RBSystem::add_scaled_matrix_and_vector(1., truth_intrr_assembly, truth_bndry_assembly,
                                         this->matrix, this->rhs);

  STOP_LOG("truth_assembly()", "QNTransientRBSystem");
}

Number QNTransientRBSystem::eval_theta_c()
{
  libmesh_assert(theta_c != NULL);

  return theta_c( *theta_data );
}

void QNTransientRBSystem::initialize_truth ()
{
  START_LOG("initialize_truth()", "QNTransientRBSystem");

  if (nonzero_initialization)
  {
    // Put the initial condition in solution and current_local_solution
    Parent::initialize_truth();

    // We have loaded the initialization into the solution
    // vector, now copy to current_newton_iterate
    *current_newton_iterate = *current_local_solution;
  }
  else
  {
    // Otherwise zero out current_newton_iterate as a default
    current_newton_iterate->zero();

    // Also, set the solution to zero
    // since we use the solution vector in
    // TransientRBSystem::update_RB_initial_condition_all_N
    solution->zero();
  }

  STOP_LOG("initialize_truth()", "QNTransientRBSystem");
}

Real QNTransientRBSystem::truth_solve(int write_interval)
{
  START_LOG("truth_solve()", "QNTransientRBSystem");

  const MeshBase& mesh = get_mesh();

  initialize_truth();

  set_time_level(0);

  // Now compute the truth outputs
  for(unsigned int n=0; n<get_n_outputs(); n++)
    for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
    {
      truth_outputs_all_k[n][_k] = eval_theta_q_l(n,q_l)*get_output_vector(n,q_l)->dot(*current_newton_iterate);
    }

  // Load initial projection error into temporal_data dense matrix
  solution->zero();
  solution->add(1., *current_newton_iterate);
  if(compute_truth_projection_error)
    set_error_temporal_data();

  for(unsigned int time_level=1; time_level<=_K; time_level++)
    {
      set_time_level(time_level); // update the member variable _k

      if(!quiet)
      {
        libMesh::out << std::endl << "Truth solve, time step " << _k << std::endl;
      }

      // Set the old solution to be the result of the previous Newton loop
      *old_local_solution = *current_newton_iterate;

      // Needed for eisenstat tolerance calculation
      Real old_rhs_norm = 0.;

      // std::ostringstream event_name;
      // event_name << "Newton iterations for timestep " << time_level;
      // PerfLog perf_log ("Newton Loop");
      // perf_log.start_event(event_name.str());

      // Now we begin the nonlinear loop
      for (unsigned int l=0; l<n_newton_steps; ++l)
        {
          // Assemble & solve the linear system.
          truth_assembly();

	  // Compute norm of rhs, this is the current nonlinear residual, and
	  // the linear solver tolerance should (in most cases) be set based on this value.
	  rhs->close(); // necessary?
	  const Real rhs_norm = rhs->l2_norm();

	  // If we print this out, we can compare how much the previous linear iterations
	  // have "oversolved" the system, if at all.
	  if(!quiet)
	    libMesh::out << "Current nonlinear residual: " << rhs_norm << std::endl;

	  // Assume fixed tolerance == eisenstat_eta0
	  Real current_linear_tolerance = eisenstat_eta0;

	  if (use_eisenstat)
	    {
	      // 1.) Tolerance selection method of Dembo and Steihaug.  This method
	      // is simple, but it's not robust in practice because it depends on the
	      // scale of the residual norm...
	      // tol = min{1/(l+2), rhs_norm}
	      //	  const Real current_linear_tolerance =
	      //	    std::min( 1./static_cast<Real>(l+2), rhs_norm );

	      // 2.) Eisenstat tolerance calculation
	      if (l!=0)
		{
		  current_linear_tolerance = eisenstat_gamma * std::pow(rhs_norm/old_rhs_norm, eisenstat_alpha);

		  // But don't let it get larger than the initial value!
		  current_linear_tolerance = std::min(eisenstat_eta0, current_linear_tolerance);

		  // Eisenstat doesn't take into account the desired *nonlinear* tolerance
		  // so you will almost always oversolve the final Newton step unless you set:
		  current_linear_tolerance = std::max(current_linear_tolerance, 0.1*nonlinear_tolerance);
		}

	      // Save rhs_norm, for next iteration!
	      old_rhs_norm = rhs_norm;
	    }

	  this->get_equation_systems().parameters.set<Real>("linear solver tolerance") = current_linear_tolerance;

	  if(!quiet)
	    libMesh::out << "Using relative linear tolerance: " << current_linear_tolerance << std::endl;

          // We're solving for delta_u, so set our intial guess to zero
          solution->zero(); // Do this after assembly to be sure
          solve();

#if defined(LIBMESH_ENABLE_AMR) || defined(LIBMESH_ENABLE_PERIODIC)
          if(enforce_constraints_exactly)
          {
            get_dof_map().enforce_constraints_exactly(*this); 
            update();
          }
#endif

          // solution now holds delta_u, so update the Newton iterate
          current_newton_iterate->add(1., *current_local_solution);

          // Compute the l2 norm of the difference
          const Real norm_delta = solution->l2_norm();

          // How many iterations were required to solve the linear system?
          const unsigned int n_linear_iterations = this->n_linear_iterations();

	  // What was the initial residual of the linear system?  Unfortunately,
	  // this is a little less straightforward to obtain: There is not a generic
	  // LinearSolver routine to do it, we have to cast to PetscLinearSolver...
	  // I like to have this value since the difference between this and
	  // the final linear residual from the previous step is a measure of the amount
	  // of oversolving.

          // What was the final residual of the linear system?
          const Real final_linear_residual = this->final_linear_residual();

          // Print out convergence information for the linear and
          // nonlinear iterations.
          if(!quiet)
          {
	    // Save cout flags
	    std::ios_base::fmtflags flags = libMesh::out.flags();

	    // Use scientific formatting for these numbers
	    libMesh::out.setf(std::ios::scientific);

            libMesh::out << "Linear solver converged at step: "
                      << n_linear_iterations
                      << ", final residual: "
                      << final_linear_residual
                      << "  Nonlinear convergence: ||delta_u|| = "
                      << norm_delta
                      << std::endl;

	    // reset the original format flags
	    libMesh::out.flags(flags);
          }

          if ((norm_delta < nonlinear_tolerance) &&
              (final_linear_residual < nonlinear_tolerance))
            {
              if(!quiet)
              {
                libMesh::out << " Nonlinear solver converged at step "
                          << l
                          << std::endl << std::endl << std::endl;
              }
              break;
            }

          // Make sure the solver converges correctly
          if( (l==(n_newton_steps-1)) &&
              ((norm_delta > nonlinear_tolerance) || (final_linear_residual > nonlinear_tolerance)) )
          {
            libMesh::out << "Warning: Linear solver may not have converged! Final linear residual = "
                      << this->final_linear_residual() << ", number of iterations = "
                      << this->n_linear_iterations() << std::endl << std::endl;
//             libmesh_error();
          }
        } // end nonlinear loop

      // perf_log.stop_event(event_name.str());

//       // Perform a truth SCM solve based on current_newton_iterate
//       EquationSystems& es = this->get_equation_systems();
//       QNTransientSCMSystem& eigen_system = es.get_system<QNTransientSCMSystem>(eigen_system_name);
//       eigen_system.compute_truth_stability_constant();

      // Now compute the truth outputs
      for(unsigned int n=0; n<get_n_outputs(); n++)
        for(unsigned int q_l=0; q_l<get_Q_l(n); q_l++)
        {
          truth_outputs_all_k[n][_k] = eval_theta_q_l(n,q_l)*get_output_vector(n,q_l)->dot(*current_newton_iterate);
        }

      if ( (write_interval > 0) && (_k%write_interval == 0) )
        {
          solution->zero();
          solution->add(1., *current_newton_iterate);

          OStringStream file_name;

          // We write the file name in the gmv auto-read format.
          file_name << "truth.gmv.";
          OSSRealzeroright(file_name,3,0, _k);

          GMVIO(mesh).write_equation_systems (file_name.str(),
                                              this->get_equation_systems());
        }

      // Load projection error into temporal_data dense matrix
      solution->zero();
      solution->add(1., *current_newton_iterate);
      if(compute_truth_projection_error)
        set_error_temporal_data();

    }

  // Get the L2 norm of the truth solution at time-level _K
  // Useful for normalizing our true error data
  if(!low_memory_mode)
  {
    L2_matrix->vector_mult(*inner_product_storage_vector, *current_newton_iterate);
  }
  else
  {
    assemble_L2_matrix(matrix);
    matrix->vector_mult(*inner_product_storage_vector, *current_newton_iterate);
  }

  Real final_truth_L2_norm = libmesh_real(std::sqrt( inner_product_storage_vector->dot(*current_newton_iterate)));

  STOP_LOG("truth_solve()", "QNTransientRBSystem");

  return final_truth_L2_norm;
}


void QNTransientRBSystem::update_system()
{
  if(get_delta_N() == 0)
    return;

  if(!low_memory_mode)
  {
    libMesh::out << "Updating trilinear form operators" << std::endl;
    update_trilinear_form_operators();
  }

  Parent::update_system();
}

// The slow (but simple, non-error prone) way to compute the residual dual norm
// Useful for error checking
//Real QNTransientRBSystem::compute_residual_dual_norm(const unsigned int N)
//{
//  START_LOG("compute_residual_dual_norm()", "QNTransientRBSystem");
//
//  // This is the "slow" way of computing the residual, but it is useful
//  // for validating the "fast" way.
//  // Note that this only works in serial since otherwise each processor will
//  // have a different parameter value during the Greedy training.
//
//  // Create temporary vectors to store the system data
//  AutoPtr< NumericVector<Number> > temp1 = NumericVector<Number>::build();
//  temp1->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
//  temp1->zero();
//
//  AutoPtr< NumericVector<Number> > temp2 = NumericVector<Number>::build();
//  temp2->init (this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);
//  temp2->zero();
//
//  AutoPtr< NumericVector<Number> > ghosted_temp1 = NumericVector<Number>::build();
//  ghosted_temp1->init (this->n_dofs(), this->n_local_dofs(),
//                       this->get_dof_map().get_send_list(), false,
//                       GHOSTED);
//
//  AutoPtr< NumericVector<Number> > ghosted_temp2 = NumericVector<Number>::build();
//  ghosted_temp2->init (this->n_dofs(), this->n_local_dofs(),
//                       this->get_dof_map().get_send_list(), false,
//                       GHOSTED);
//
//  // Store old_solution, since we don't want to corrupt it in case call
//  // this function in between truth solves
//  *ghosted_temp1 = *old_local_solution;
//  *ghosted_temp2 = *current_newton_iterate;
//
//  current_newton_iterate->zero();
//  old_local_solution->zero();
//  for(unsigned int i=0; i<N; i++)
//  {
//    temp1->add(RB_solution(i),     *basis_functions[i]);
//    temp2->add(old_RB_solution(i), *basis_functions[i]);
//  }
//  // Load temp1 and temp2 into current_newton_iterate and old_local_solution in order to do assembly
//  temp1->localize (*current_newton_iterate, this->get_dof_map().get_send_list());
//  temp2->localize (*old_local_solution, this->get_dof_map().get_send_list());
//
//  // Call truth assembly to put the residual in the rhs vector
//  this->truth_assembly();
//
//  // Restore current_newton_iterate and old_local_solution
//  *current_newton_iterate = *ghosted_temp1;
//  *old_local_solution     = *ghosted_temp2;
//
//  zero_dirichlet_dofs_on_rhs();
//
//  // Then solve the system to get the Reisz representor
//  matrix->zero();
//  matrix->add(1., *inner_product_matrix);
//
//  if(constrained_problem)
//    matrix->add(1., *constraint_matrix);
//
//  solution->zero();
//  solve();
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
//  inner_product_matrix->vector_mult(*inner_product_storage_vector, *solution);
//  Real slow_residual_norm_sq = inner_product_storage_vector->dot(*solution);
//
//  STOP_LOG("compute_residual_dual_norm()", "QNTransientRBSystem");
//
//  return libmesh_real(std::sqrt( slow_residual_norm_sq ));
//}

void QNTransientRBSystem::update_RB_system_matrices()
{
  START_LOG("update_RB_system_matrices()", "QNTransientRBSystem");

  Parent::update_RB_system_matrices();

  QNTransientRBEvaluation* qn_rb_eval = libmesh_cast_ptr<QNTransientRBEvaluation*>(rb_eval);

  // Now we just need to assemble the trilinear form
  unsigned int RB_size = get_n_basis_functions();

  for(unsigned int n=0; n<RB_size; n++)
  {
    if(low_memory_mode)
    {
      assemble_C_n_matrix(n, matrix);
    }

    for(unsigned int j=0; j<RB_size; j++)
    {
      if(!low_memory_mode)
      {
        get_C_n(n)->vector_mult(*inner_product_storage_vector, get_basis_function(j));
      }
      else
      {
        matrix->vector_mult(*inner_product_storage_vector, get_basis_function(j));
      }

      for(unsigned int i=0; i<RB_size; i++)
      {
        // Short-circuit all entries that have already been initialized
        if(qn_rb_eval->RB_trilinear_form[n][i][j] == 0.)
        {
          qn_rb_eval->RB_trilinear_form[n][i][j] = get_basis_function(i).dot(*inner_product_storage_vector);
        }
      }
    }
  }

  STOP_LOG("update_RB_system_matrices()", "QNTransientRBSystem");
}

void QNTransientRBSystem::set_context_solution_vec(NumericVector<Number>& vec)
{
  // Set current_local_solution = vec so that we can access
  // vec from QNTransientRBContext during assembly
  vec.localize
    (*current_newton_iterate, this->get_dof_map().get_send_list());
}

void QNTransientRBSystem::update_residual_terms(bool compute_inner_products)
{
  START_LOG("update_residual_terms()", "QNTransientRBSystem");

  Parent::update_residual_terms(compute_inner_products);

  QNTransientRBEvaluation* qn_rb_eval = libmesh_cast_ptr<QNTransientRBEvaluation*>(rb_eval);

  unsigned int RB_size = get_n_basis_functions();

  if(!low_memory_mode)
  {
    matrix->zero();
    matrix->add(1., *inner_product_matrix);
    if(constrained_problem)
      matrix->add(1., *constraint_matrix);
  }
  // Don't need to assemble here in low-memory mode
  // since we will assemble during the loop below

  if(reuse_preconditioner)
  {
    // For the first solve, make sure we generate a new preconditioner
    // Actually, this is not necessary as we can use the preconditioner
    // generated by RBSystem::update_residual_terms()
    linear_solver->same_preconditioner = true;
  }


  for(unsigned int n=0; n<RB_size; n++)
  {
    for(unsigned int j=0; j<RB_size; j++)
    {
      if(C_representor[n][j] == NULL) // Short-circuit entries that are already initialized
      {
        C_representor[n][j] = NumericVector<Number>::build().release();
        C_representor[n][j]->init(this->n_dofs(), this->n_local_dofs(), false, libMeshEnums::PARALLEL);

        if(!low_memory_mode)
        {
          get_C_n(n)->vector_mult(*rhs, get_basis_function(j));
        }
        else
        {
          assemble_C_n_matrix(n, matrix);
          matrix->vector_mult(*rhs, get_basis_function(j));

          // Now need to restore the inner product matrix
          assemble_inner_product_matrix(matrix);
          if(constrained_problem)
            add_scaled_matrix_and_vector(1., constraint_assembly, NULL, matrix, NULL);

// Can't use a matvec here since we need to access two vectors
// for the trilinear form and matvec just accesses one vector
//           assemble_scaled_matvec(1.,
//                                  C_assembly,
//                                  *rhs,
//                                  *basis_functions[j]);
        }

        rhs->scale(-1.); // Negate: more convenient for residual calculations
        zero_dirichlet_dofs_on_rhs();

        solution->zero();

	if (!quiet)
	  libMesh::out << "Starting solve [n][j]=["
		       << n <<"]["<< j << "] in QNTransientRBSystem::update_residual_terms() at "
		       << Utility::get_timestamp() << std::endl;

        solve();

	if (!quiet)
	  {
	    libMesh::out << "Finished solve [n][j]=["
		         << n <<"]["<< j << "] in QNTransientRBSystem::update_residual_terms() at "
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
          libMesh::out << "Warning: Linear solver may not have converged. Final linear residual = "
                       << this->final_linear_residual() << ", number of iterations = "
                       << this->n_linear_iterations() << std::endl << std::endl;
//           libmesh_error();
        }

        *C_representor[n][j] = *solution;

        if(reuse_preconditioner)
        {
          // After we do a solve, tell PETSc we want to reuse the preconditioner
          // since the system matrix is not changing.
          linear_solver->same_preconditioner = true;
        }
      }
    }
  }

  if(reuse_preconditioner)
  {
    // We no longer want to reuse the preconditioner
    linear_solver->same_preconditioner = false;
  }

  // Now compute the necessary inner products if requested
  if (compute_inner_products)
    {
      if(low_memory_mode && constrained_problem)
	assemble_inner_product_matrix(matrix);

      // We short-circuit all entries that are non-zero
      // to avoid re-computing inner products unnecessarily
      for(unsigned int q_f=0; q_f<get_Q_f(); q_f++)
	{
	  if(!low_memory_mode)
	    {
	      inner_product_matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f]);
	    }
	  else
	    {
	      matrix->vector_mult(*inner_product_storage_vector, *F_q_representor[q_f]);
	    }

	  for(unsigned int n1=0; n1<RB_size; n1++)
	    {
	      for(unsigned int j1=0; j1<RB_size; j1++)
		{
		  if( qn_rb_eval->Fq_C_representor_norms[q_f][n1][j1] == 0. )
		    {
		      qn_rb_eval->Fq_C_representor_norms[q_f][n1][j1] =
			C_representor[n1][j1]->dot(*inner_product_storage_vector);
		    }
		} // end for j1
	    } // end for n1
	} // end for q_f

      for(unsigned int q_m=0; q_m<get_Q_m(); q_m++)
	{
	  for(unsigned int i=0; i<RB_size; i++)
	    {
	      if(!low_memory_mode)
		{
		  inner_product_matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][i]);
		}
	      else
		{
		  matrix->vector_mult(*inner_product_storage_vector, *M_q_representor[q_m][i]);
		}

	      for(unsigned int n1=0; n1<RB_size; n1++)
		{
		  for(unsigned int j1=0; j1<RB_size; j1++)
		    {
		      if( qn_rb_eval->Mq_C_representor_norms[q_m][i][n1][j1] == 0. )
			{
			  qn_rb_eval->Mq_C_representor_norms[q_m][i][n1][j1] =
			    C_representor[n1][j1]->dot(*inner_product_storage_vector);
			}
		    } // end for j1
		} // end for n1
	    } // end for i
	} // end for q_m

      
      for(unsigned int q_a=0; q_a<get_Q_a(); q_a++)
	{
	  for(unsigned int i=0; i<RB_size; i++)
	    {
	      if(!low_memory_mode)
		{
		  inner_product_matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a][i]);
		}
	      else
		{
		  matrix->vector_mult(*inner_product_storage_vector, *A_q_representor[q_a][i]);
		}

	      for(unsigned int n1=0; n1<RB_size; n1++)
		{
		  for(unsigned int j1=0; j1<RB_size; j1++)
		    {
		      if( qn_rb_eval->Aq_C_representor_norms[q_a][i][n1][j1] == 0. )
			{
			  qn_rb_eval->Aq_C_representor_norms[q_a][i][n1][j1] =
			    C_representor[n1][j1]->dot(*inner_product_storage_vector);
			}
		    } // end for j1
		} // end for n1
	    } // end for i
	} // end for q_a

      
      for(unsigned int n1=0; n1<RB_size; n1++)
	{
	  for(unsigned int j1=0; j1<RB_size; j1++)
	    {
	      if(!low_memory_mode)
		{
		  inner_product_matrix->vector_mult(*inner_product_storage_vector, *C_representor[n1][j1]);
		}
	      else
		{
		  matrix->vector_mult(*inner_product_storage_vector, *C_representor[n1][j1]);
		}

	      for(unsigned int n2=0; n2<RB_size; n2++)
		{
		  for(unsigned int j2=0; j2<RB_size; j2++)
		    {
		      if( qn_rb_eval->C_C_representor_norms[n1][j1][n2][j2] == 0.)
			{
			  qn_rb_eval->C_C_representor_norms[n1][j1][n2][j2] =
			    C_representor[n2][j2]->dot(*inner_product_storage_vector);
			}
		    } // end for j2
		} // end for n2
	    } // end for j1
	} // end for n1
    } // end if (compute_inner_products)


  STOP_LOG("update_residual_terms()", "QNTransientRBSystem");
}

void QNTransientRBSystem::attach_C(theta_q_fptr theta_c_in, affine_assembly_fptr C_assembly_in)
{
  C_assembly = C_assembly_in;
  theta_c = theta_c_in;
}

void QNTransientRBSystem::attach_truth_assembly(affine_assembly_fptr truth_intrr_assembly_in,
                                               affine_assembly_fptr truth_bndry_assembly_in)
{
  truth_intrr_assembly = truth_intrr_assembly_in;
  truth_bndry_assembly = truth_bndry_assembly_in;
}

void QNTransientRBSystem::update_trilinear_form_operators()
{
  START_LOG("update_trilinear_form_operators()", "QNTransientRBSystem");

  unsigned int RB_size = get_n_basis_functions();
  for(unsigned int n=(RB_size-delta_N); n<RB_size; n++)
  {
    assemble_C_n_matrix(n, get_C_n(n));
  }
  STOP_LOG("update_trilinear_form_operators()", "QNTransientRBSystem");
}

void QNTransientRBSystem::update_all_trilinear_form_operators()
{
  if(!low_memory_mode && initialize_mesh_dependent_data)
  {
    unsigned int saved_delta_N = delta_N;
    delta_N = get_n_basis_functions();
    update_trilinear_form_operators();
    delta_N = saved_delta_N;
  }
}

SparseMatrix<Number>* QNTransientRBSystem::get_C_n(unsigned int n)
{
  if(low_memory_mode)
  {
    libMesh::err << "Error: The C matrices are not stored in low-memory mode." << std::endl;
    libmesh_error();
  }

  if(n >= get_n_basis_functions())
  {
    libMesh::err << "Error: We must have n < get_n_basis_functions() in get_C_n."
                 << std::endl;
    libmesh_error();
  }

  return C_n_vector[n];
}

void QNTransientRBSystem::assemble_C_n_matrix(unsigned int n, SparseMatrix<Number>* input_matrix)
{
  if(n >= get_n_basis_functions())
  {
    libMesh::err << "Error: We must have n < get_n_basis_functions() in assemble_C_n_matrix."
                 << std::endl;
    libmesh_error();
  }

  // Store basis_functions[n] in current_newton_iterate for assembly
  get_basis_function(n).localize
    (*current_newton_iterate, this->get_dof_map().get_send_list());

  input_matrix->zero();
  add_scaled_matrix_and_vector(1.,
                               C_assembly,
                               NULL,
                               input_matrix,
                               NULL);
}

void QNTransientRBSystem::add_scaled_Cn(Number scalar, unsigned int n, SparseMatrix<Number>* input_matrix, bool symmetrize)
{
  START_LOG("add_scaled_Cn()", "QNTransientRBSystem");

  if(n >= get_n_basis_functions())
  {
    libMesh::err << "Error: We must have n < get_n_basis_functions() in add_scaled_Cn."
                 << std::endl;
    libmesh_error();
  }

  if(!low_memory_mode && !symmetrize)
  {
    input_matrix->add(scalar, *get_C_n(n));
    input_matrix->close();
  }
  else
  {
    // Store basis_functions[n] in current_newton_iterate for assembly
    get_basis_function(n).localize
      (*current_newton_iterate, this->get_dof_map().get_send_list());

    add_scaled_matrix_and_vector(scalar,
                                C_assembly,
                                NULL,
                                input_matrix,
                                NULL,
                                symmetrize);
  }

  STOP_LOG("add_scaled_Cn()", "QNTransientRBSystem");
}

void QNTransientRBSystem::add_scaled_current_C(Number scalar, SparseMatrix<Number>* input_matrix, bool symmetrize)
{
  START_LOG("add_scaled_current_C()", "QNTransientRBSystem");

  // Assemble based on the current value of current_newton_iterate
  add_scaled_matrix_and_vector(scalar,
                               C_assembly,
                               NULL,
                               input_matrix,
                               NULL,
                               symmetrize);

  STOP_LOG("add_scaled_current_C()", "QNTransientRBSystem");
}

Real QNTransientRBSystem::get_SCM_lower_bound()
{
  // Get the SCM lower bound from eigen_system
  EquationSystems& es = this->get_equation_systems();
  QNTransientSCMSystem& eigen_system = es.get_system<QNTransientSCMSystem>(eigen_system_name);

  // Create a parameter vector in which the current time-level
  // is appended to current_parameters.
  std::vector<Real> params(get_n_params()+1);
  for(unsigned int i=0; i<get_n_params(); i++)
    params[i] = get_current_parameters()[i];

  params[params.size()-1] = _k;

  // Also, construct a vector storing the RB coefficients
  QNTransientRBEvaluation* qn_rb_eval = libmesh_cast_ptr<QNTransientRBEvaluation*>(rb_eval);
  std::vector<Number> RB_u_euler_theta(qn_rb_eval->RB_solution.size());
  for(unsigned int i=0; i<qn_rb_eval->RB_solution.size(); i++)
    RB_u_euler_theta[i] = euler_theta*qn_rb_eval->RB_solution(i) + (1.-euler_theta)*qn_rb_eval->old_RB_solution(i);

  // Pass params and RB_u_euler_theta to the associated SCM system
  eigen_system.set_current_parameters( params );
  eigen_system.set_current_RB_coeffs( RB_u_euler_theta );
  return eigen_system.get_SCM_LB();
}

Real QNTransientRBSystem::get_SCM_upper_bound()
{
  // Get the SCM lower bound from eigen_system
  EquationSystems& es = this->get_equation_systems();
  QNTransientSCMSystem& eigen_system = es.get_system<QNTransientSCMSystem>(eigen_system_name);

  // Create a parameter vector in which the current time-level
  // is appended to current_parameters.
  std::vector<Real> params(get_n_params()+1);
  for(unsigned int i=0; i<get_n_params(); i++)
    params[i] = get_current_parameters()[i];

  params[params.size()-1] = _k;

  QNTransientRBEvaluation* qn_rb_eval = libmesh_cast_ptr<QNTransientRBEvaluation*>(rb_eval);
  std::vector<Number> RB_u_euler_theta(qn_rb_eval->RB_solution.size());
  for(unsigned int i=0; i<qn_rb_eval->RB_solution.size(); i++)
    RB_u_euler_theta[i] = euler_theta*qn_rb_eval->RB_solution(i) + (1.-euler_theta)*qn_rb_eval->old_RB_solution(i);

  eigen_system.set_current_parameters( params );
  eigen_system.set_current_RB_coeffs( RB_u_euler_theta );
  return eigen_system.get_SCM_UB();
}



void QNTransientRBSystem::write_offline_data_to_files(const std::string& directory_name,
                                                      RBDataIO io_flag)
{
  START_LOG("write_offline_data_to_files()", "QNTransientRBSystem");

  Parent::write_offline_data_to_files(directory_name, io_flag);
  
  // return here if we only want basis dependent data
  if(io_flag == BASIS_DEPENDENT)
    return;

  if (store_representors)
    {
      // Write out C_representors.  These are useful to have when restarting,
      // so you don't have to recompute them all over again.  There should be
      // Nmax^2 of these, but we only write out up to N, where N is the current RB
      // dimension.
      libMesh::out << "Writing out the C_representors..." << std::endl;

      std::ostringstream file_name;
      const std::string residual_representor_suffix = (write_binary_residual_representors ? ".xdr" : ".dat");
      //  struct stat stat_info;

      // Residual representors written out to their own separate directory
      std::string residual_representors_dir = "residual_representors";

      unsigned int new_start_index = this->get_n_basis_functions() - delta_N;

      for (unsigned int i=0; i<this->get_n_basis_functions(); ++i)
	for (unsigned int j=0; j<this->get_n_basis_functions(); ++j)
	  {
	    // if (C_representor[i][j] != NULL)
	    if ( (i >= new_start_index) || // do all columns in the new rows
		 (j >= new_start_index) )  // in old rows, just do new columns
	      {
		libMesh::out << "Writing out C_representor[" << i << "][" << j << "]..." << std::endl;
		libmesh_assert(C_representor[i][j] != NULL);

		file_name.str(""); // reset filename
		file_name << residual_representors_dir
			  << "/C_representor" << i << "_" << j << residual_representor_suffix;

		{
		  // No need to copy!
		  // *solution = *(C_representor[i][j]);
		  C_representor[i][j]->swap(*solution);

		  Xdr c_data(file_name.str(),
			     write_binary_residual_representors ? ENCODE : WRITE);

		  write_serialized_data(c_data, false);

		  // Synchronize before moving on
		  Parallel::barrier();

		  // Swap back
		  C_representor[i][j]->swap(*solution);

		  // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
		  // for the system call, be sure to do it only on one processor, etc.
		}
	      }
	  }
    } // end if (store_representors)

  STOP_LOG("write_offline_data_to_files()", "QNTransientRBSystem");
}




void QNTransientRBSystem::read_offline_data_from_files(const std::string& directory_name,
                                                       const RBDataIO io_flag)
{
  START_LOG("read_offline_data_from_files()", "QNTransientRBSystem");

  Parent::read_offline_data_from_files(directory_name, io_flag);
  
  // return if we only want basis dependent data
  if(io_flag == BASIS_DEPENDENT)
    return;

  // Read in the residual representors from file if requested
  if (store_representors)
    {
      if (!quiet)
	libMesh::out << "Reading in the C_representors..." << std::endl;

      const std::string residual_representors_dir = "residual_representors";
      const std::string residual_representor_suffix =
	(read_binary_residual_representors ? ".xdr" : ".dat");

      std::ostringstream file_name;
      struct stat stat_info;

      // Read in the C_representors from file.  There should be
      // this->get_n_basis_functions()^2 of these.  FIXME: Should we worry
      // about memory leaks here?
      for (unsigned int i=0; i<this->get_n_basis_functions(); ++i)
	for (unsigned int j=0; j<this->get_n_basis_functions(); ++j)
	  {
            if (C_representor[i][j] != NULL)
	    {
	      libMesh::out << "Error, must delete existing C_representor before reading in from file."
			   << std::endl;
	      libmesh_error();
	    }

            file_name.str(""); // reset the filename
            file_name << residual_representors_dir
              << "/C_representor" << i << "_" << j << residual_representor_suffix;

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

            Xdr c_data(file_name.str(),
               read_binary_residual_representors ? DECODE : READ);

            read_serialized_data(c_data, false);

            C_representor[i][j] = NumericVector<Number>::build().release();
            C_representor[i][j]->init (this->n_dofs(),
                                       this->n_local_dofs(),
                                       false, libMeshEnums::PARALLEL);

            // No need to copy, just swap
            // *C_representor[i][j] = *solution;
            C_representor[i][j]->swap(*solution);
	  }

      // Finally, if we read in the representors then we must be planning to continue the basis training
      // therefore if we are in high-memory mode we need to recompute all the trilinear form matrices
      if(store_basis_functions && !low_memory_mode)
        {
          libMesh::out << "Updating trilinear form operators" << std::endl;
          update_all_trilinear_form_operators();
        }

    } // end if (store_representors)



  STOP_LOG("read_offline_data_from_files()", "QNTransientRBSystem");
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
