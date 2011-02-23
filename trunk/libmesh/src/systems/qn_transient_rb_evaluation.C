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

#include "qn_transient_rb_evaluation.h"
#include "qn_transient_rb_system.h"
#include "libmesh_logging.h"

#include "o_string_stream.h"
#include <fstream>
#include <sstream>

namespace libMesh
{
        
QNTransientRBEvaluation::QNTransientRBEvaluation (RBSystem& rb_sys_in)
  :
  Parent(rb_sys_in)
{
}

void QNTransientRBEvaluation::clear()
{
  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  // Zero RB_trilinear_form and the representor norms.
  // Do this _before_ calling Parent::clear()
  // since we rely on get_n_basis_functions.
  // It is important to be able to clear the data below
  // since we check whether these values are zero before
  // updating representor data in QNTransientRBSystem
  unsigned int RB_size = get_n_basis_functions();
  for(unsigned int n=0; n<RB_size; n++)
    for(unsigned int j=0; j<RB_size; j++)
      for(unsigned int i=0; i<RB_size; i++)
        RB_trilinear_form[n][i][j] = 0.;

  for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
    for(unsigned int n1=0; n1<RB_size; n1++)
      for(unsigned int j1=0; j1<RB_size; j1++)
        Fq_C_representor_norms[q_f][n1][j1] = 0.;

  for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
    for(unsigned int i=0; i<RB_size; i++)
      for(unsigned int n1=0; n1<RB_size; n1++)
        for(unsigned int j1=0; j1<RB_size; j1++)
          Mq_C_representor_norms[q_m][i][n1][j1] = 0.;

  for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
    for(unsigned int i=0; i<RB_size; i++)
      for(unsigned int n1=0; n1<RB_size; n1++)
        for(unsigned int j1=0; j1<RB_size; j1++)
          Aq_C_representor_norms[q_a][i][n1][j1] = 0.;

  for(unsigned int n1=0; n1<RB_size; n1++)
    for(unsigned int j1=0; j1<RB_size; j1++)
      for(unsigned int n2=0; n2<RB_size; n2++)
        for(unsigned int j2=0; j2<RB_size; j2++)
          C_C_representor_norms[n1][j1][n2][j2] = 0.;

  Parent::clear();
}

void QNTransientRBEvaluation::initialize()
{
  Parent::initialize();

  // Initialize the N (i.e. RB) dependent data structures
  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  // Resize the RB trilinear form
  const unsigned int Nmax = qn_rb_sys.get_Nmax();
  RB_trilinear_form.resize(Nmax);
  for(unsigned int i=0; i<Nmax; i++)
  {
    RB_trilinear_form[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      RB_trilinear_form[i][j].resize(Nmax);
      for(unsigned int l=0; l<Nmax; l++)
      {
        RB_trilinear_form[i][j][l] = 0.;
      }
    }
  }


  // Initialize vectors for the norms of the representors
  Fq_C_representor_norms.resize(qn_rb_sys.get_Q_f());
  for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
  {
    Fq_C_representor_norms[q_f].resize(Nmax);
    for(unsigned int n1=0; n1<Nmax; n1++)
    {
      Fq_C_representor_norms[q_f][n1].resize(Nmax);
      for(unsigned int j1=0; j1<Nmax; j1++)
      {
        Fq_C_representor_norms[q_f][n1][j1] = 0.;
      }
    }
  }

  Mq_C_representor_norms.resize(qn_rb_sys.get_Q_m());
  for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
  {
    Mq_C_representor_norms[q_m].resize(Nmax);
    for(unsigned int i=0; i<Nmax; i++)
    {
      Mq_C_representor_norms[q_m][i].resize(Nmax);
      for(unsigned int n1=0; n1<Nmax; n1++)
      {
        Mq_C_representor_norms[q_m][i][n1].resize(Nmax);
        for(unsigned int j1=0; j1<Nmax; j1++)
        {
          Mq_C_representor_norms[q_m][i][n1][j1] = 0.;
        }
      }
    }
  }

  Aq_C_representor_norms.resize(qn_rb_sys.get_Q_a());
  for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
  {
    Aq_C_representor_norms[q_a].resize(Nmax);
    for(unsigned int i=0; i<Nmax; i++)
    {
      Aq_C_representor_norms[q_a][i].resize(Nmax);
      for(unsigned int n1=0; n1<Nmax; n1++)
      {
        Aq_C_representor_norms[q_a][i][n1].resize(Nmax);
        for(unsigned int j1=0; j1<Nmax; j1++)
        {
          Aq_C_representor_norms[q_a][i][n1][j1] = 0.;
        }
      }
    }
  }

  C_C_representor_norms.resize(Nmax);
  for(unsigned int n1=0; n1<Nmax; n1++)
  {
    C_C_representor_norms[n1].resize(Nmax);
    for(unsigned int j1=0; j1<Nmax; j1++)
    {
      C_C_representor_norms[n1][j1].resize(Nmax);
      for(unsigned int n2=0; n2<Nmax; n2++)
      {
        C_C_representor_norms[n1][j1][n2].resize(Nmax);
        for(unsigned int j2=0; j2<Nmax; j2++)
        {
          C_C_representor_norms[n1][j1][n2][j2] = 0.;
        }
      }
    }
  }
}

Real QNTransientRBEvaluation::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "QNTransientRBEvaluation");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }
  if(N==0)
  {
    libMesh::err << "ERROR: N must be greater than 0 in RB_solve" << std::endl;
    libmesh_error();
  }
  
  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  DenseMatrix<Number> Base_RB_LHS_matrix(N,N);
  Base_RB_LHS_matrix.zero();

  DenseMatrix<Number> RB_RHS_Aq_matrix(N,N);
  RB_RHS_Aq_matrix.zero();

  DenseMatrix<Number> RB_mass_matrix_N(N,N);
  RB_mass_matrix_N.zero();
  DenseMatrix<Number> RB_M_q_m;
  for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
  {
    RB_M_q_vector[q_m].get_principal_submatrix(N, RB_M_q_m);
    RB_mass_matrix_N.add(qn_rb_sys.eval_theta_q_m(q_m), RB_M_q_m);
  }

  const Real dt = qn_rb_sys.get_dt();

  Base_RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);

  for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = qn_rb_sys.eval_theta_q_a(q_a);
    for(unsigned int i=0; i<N; i++)
    {
      for(unsigned int j=0; j<N; j++)
      {
        Base_RB_LHS_matrix(i,j) += qn_rb_sys.get_euler_theta()*cached_theta_q_a*RB_A_q_vector[q_a](i,j);
        RB_RHS_Aq_matrix(i,j) += -cached_theta_q_a*RB_A_q_vector[q_a](i,j);
      }
    }
  }

  // Set system time level to 0
  error_bound_all_k.resize(qn_rb_sys.get_K()+1);
  qn_rb_sys.set_time_level(0);

  // This is the actual LHS matrix
  DenseMatrix<Number> RB_LHS_matrix(N,N);
  RB_LHS_matrix.zero();

  // Initialize a vector to store our current Newton iterate
  DenseVector<Number> RB_u_bar(N);
  RB_u_bar.zero();

  // Load the initial condition into RB_u_bar
  RB_u_bar = RB_initial_condition_all_N[N-1];

  // Initialize solution storage vectors
  RB_solution.resize(N);
  old_RB_solution.resize(N);
  RB_temporal_solution_data.resize(qn_rb_sys.get_K()+1);
  for(unsigned int time_level=0; time_level<=qn_rb_sys.get_K(); time_level++)
  {
    RB_temporal_solution_data[time_level].resize(N);
  }
  // and load the initial data
  RB_solution = RB_u_bar;
  RB_temporal_solution_data[0] = RB_u_bar;

  Real error_bound_sum = pow( initial_L2_error_all_N[N-1], 2.);

  // Set error bound at _k=0
  error_bound_all_k[0] = std::sqrt(error_bound_sum);

  // Compute the outputs and associated error bounds at the initial time
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<qn_rb_sys.get_n_outputs(); n++)
  {
    RB_outputs_all_k[n][0] = 0.;
    for(unsigned int q_l=0; q_l<qn_rb_sys.get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs_all_k[n][0] += qn_rb_sys.eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_u_bar);
    }

    RB_output_error_bounds_all_k[n][0] = error_bound_all_k[0] * qn_rb_sys.eval_output_dual_norm(n);
  }

  // Initialize a vector to store the solution from the old time-step
  DenseVector<Number> RB_u_old(N);
  RB_u_old.zero();

  // Initialize a vector to store the Newton increment, RB_delta_u
  DenseVector<Number> RB_delta_u(N);
  RB_delta_u.zero();

  // Initialize the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  // Initialize the RB_u_euler_theta
  DenseVector<Number> RB_u_euler_theta(N);
  RB_u_euler_theta.zero();

  // Pre-compute eval_theta_c()
  Number cached_theta_c = qn_rb_sys.eval_theta_c();

  // These vectors allow us to plot the rho upper and lower bounds
  qn_rb_sys.rho_LB_vector.resize(qn_rb_sys.get_K()+1);
  qn_rb_sys.rho_UB_vector.resize(qn_rb_sys.get_K()+1);

  for(unsigned int time_level=1; time_level<=qn_rb_sys.get_K(); time_level++)
  {
    qn_rb_sys.set_time_level(time_level); // update the member variable _k

    // Set RB_u_old to be the result of the previous Newton loop
    RB_u_old = RB_u_bar;

    // Now we begin the nonlinear loop
    for (unsigned int l=0; l<qn_rb_sys.get_n_newton_steps(); ++l)
    {
      // Get u_euler_theta = euler_theta*RB_u_bar + (1-euler_theta)*RB_u_old
      RB_u_euler_theta.zero();
      for(unsigned int n=0; n<N; n++)
      {
        RB_u_euler_theta(n) += qn_rb_sys.get_euler_theta()*RB_u_bar(n) +
                              (1.-qn_rb_sys.get_euler_theta())*RB_u_old(n);
      }

      // Assemble the left-hand side for the RB linear system
      RB_LHS_matrix = Base_RB_LHS_matrix;

      // Add the trilinear term
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          for(unsigned int n=0; n<N; n++)
          {
            RB_LHS_matrix(i,j) += cached_theta_c*qn_rb_sys.get_euler_theta()*RB_u_euler_theta(n)*
                                  (RB_trilinear_form[n][i][j] + RB_trilinear_form[j][i][n]);
          }
        }
      }

      // Assemble the right-hand side for the RB linear system (the residual)
      // First add forcing terms
      RB_rhs.zero();
      for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
      {
        Number cached_theta_f = qn_rb_sys.eval_theta_q_f(q_f);
        for(unsigned int i=0; i<N; i++)
        {
          RB_rhs(i) += cached_theta_f*RB_F_q_vector[q_f](i);
        }
      }


      // Now add -1./dt * M * (RB_u_bar - RB_u_old)
      RB_mass_matrix_N.vector_mult_add(RB_rhs, -1./dt, RB_u_bar);
      RB_mass_matrix_N.vector_mult_add(RB_rhs,  1./dt, RB_u_old);

      // Now add -mu*A1*u_euler_theta
      RB_RHS_Aq_matrix.vector_mult_add(RB_rhs, 1., RB_u_euler_theta);

      // Finally add the trilinear term
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          Number RB_u_euler_theta_j = RB_u_euler_theta(j);
          for(unsigned int n=0; n<N; n++)
          {
            RB_rhs(i) += -cached_theta_c*RB_u_euler_theta(n)*RB_u_euler_theta_j*RB_trilinear_form[n][i][j];
          }
        }
      }

      RB_LHS_matrix.lu_solve(RB_rhs, RB_delta_u);

      // update the Newton iterate
      RB_u_bar += RB_delta_u;

      // Compute the l2 norm of RB_delta_u
      Real RB_delta_u_norm = RB_delta_u.l2_norm();

//      if(!quiet)
//      {
//        libMesh::out << "||RB_delta_u|| = " << RB_delta_u_norm << std::endl;
//      }

      if( RB_delta_u_norm < qn_rb_sys.get_nonlinear_tolerance())
      {
//        if(!quiet)
//        {
//          libMesh::out << "RB Newton solve converged at step " << l << std::endl << std::endl;
//        }
        break;
      }

      if( (l==(qn_rb_sys.get_n_newton_steps()-1)) &&
          (RB_delta_u_norm > qn_rb_sys.get_nonlinear_tolerance()) )
      {
        libMesh::out << "ERROR: RB Newton loop did not converge" << std::endl << std::endl;
        libmesh_error();
      }
    }

    // Load RB_solution into RB_solution_vector for residual computation
    RB_solution = RB_u_bar;
    old_RB_solution = RB_u_old;
    RB_temporal_solution_data[time_level] = RB_u_bar;

    Real rho_LB = qn_rb_sys.use_nominal_rho_LB ? qn_rb_sys.get_nominal_rho_LB() :
                                       qn_rb_sys.get_SCM_lower_bound();

    qn_rb_sys.rho_LB_vector[time_level] = rho_LB;
    qn_rb_sys.rho_UB_vector[time_level] = qn_rb_sys.use_nominal_rho_LB ? 0. : qn_rb_sys.get_SCM_upper_bound();

    // Evaluate the dual norm of the residual for RB_solution_vector
    Real epsilon_N = compute_residual_dual_norm(N);

    error_bound_sum += qn_rb_sys.residual_scaling_numer(rho_LB) * pow(epsilon_N, 2.);

    // store error bound at time-level
    error_bound_all_k[time_level] = std::sqrt(error_bound_sum/qn_rb_sys.residual_scaling_denom(rho_LB));

    // Now compute the outputs and associated errors
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<qn_rb_sys.get_n_outputs(); n++)
    {
      RB_outputs_all_k[n][time_level] = 0.;
      for(unsigned int q_l=0; q_l<qn_rb_sys.get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
        RB_outputs_all_k[n][time_level] += qn_rb_sys.eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_u_bar);
      }

      RB_output_error_bounds_all_k[n][time_level] = error_bound_all_k[time_level] * qn_rb_sys.eval_output_dual_norm(n);
    }
  }

  // Now compute the L2 norm of the RB solution at time-level _K
  // to normalize the error bound
  // We reuse RB_rhs here
  RB_mass_matrix_N.vector_mult(RB_rhs, RB_u_bar);
  Real final_RB_L2_norm = libmesh_real(std::sqrt( RB_u_bar.dot(RB_rhs) ));

  STOP_LOG("RB_solve()", "QNTransientRBEvaluation");

  return ( qn_rb_sys.return_rel_error_bound ? error_bound_all_k[qn_rb_sys.get_K()]/final_RB_L2_norm :
                                              error_bound_all_k[qn_rb_sys.get_K()] );
}

Real QNTransientRBEvaluation::compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("compute_residual_dual_norm()", "QNTransientRBEvaluation");

  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  // Use the stored representor inner product values
  // to evaluate the residual dual norm
  const Real dt = qn_rb_sys.get_dt();

  Number residual_norm_sq = 0.;

  // Use TransientRBSystem to compute all the linear terms
  residual_norm_sq += pow( Parent::uncached_compute_residual_dual_norm(N), 2.);

  // Now just need to add the terms involving the nonlinearity
  std::vector<Number> RB_u_euler_theta(N);
  std::vector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta[i]  = qn_rb_sys.get_euler_theta()*RB_solution(i) + (1.-qn_rb_sys.get_euler_theta())*old_RB_solution(i);
    mass_coeffs[i] = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  // Pre-compute eval_theta_c()
  Number cached_theta_c = qn_rb_sys.eval_theta_c();

  // All residual terms can be treated as positive quantities...
  for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = qn_rb_sys.eval_theta_q_f(q_f);
    for(unsigned int n1=0; n1<N; n1++)
    {
      for(unsigned int j1=0; j1<N; j1++)
      {
        residual_norm_sq += 2.*cached_theta_q_f*cached_theta_c*RB_u_euler_theta[n1]*RB_u_euler_theta[j1]
                              *Fq_C_representor_norms[q_f][n1][j1];
      }
    }
  }

  for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
  {
    Number cached_theta_q_m = qn_rb_sys.eval_theta_q_m(q_m);
    for(unsigned int i=0; i<N; i++)
    {
      for(unsigned int n1=0; n1<N; n1++)
      {
        for(unsigned int j1=0; j1<N; j1++)
        {
          residual_norm_sq += 2.*cached_theta_q_m*cached_theta_c*
                             mass_coeffs[i]*RB_u_euler_theta[n1]*RB_u_euler_theta[j1]*
                             Mq_C_representor_norms[q_m][i][n1][j1];
        }
      }
    }
  }

  for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = qn_rb_sys.eval_theta_q_a(q_a);
    for(unsigned int i=0; i<N; i++)
    {
      for(unsigned int n1=0; n1<N; n1++)
      {
        for(unsigned int j1=0; j1<N; j1++)
        {
          residual_norm_sq += 2.*cached_theta_q_a*cached_theta_c*
                             RB_u_euler_theta[i]*RB_u_euler_theta[n1]*RB_u_euler_theta[j1]*
                             Aq_C_representor_norms[q_a][i][n1][j1];
        }
      }
    }
  }

  for(unsigned int n1=0; n1<N; n1++)
  {
    for(unsigned int j1=0; j1<N; j1++)
    {
      Number RB_u_euler_theta_1 = RB_u_euler_theta[n1]*RB_u_euler_theta[j1];

      for(unsigned int n2=n1; n2<N; n2++)
      {
        unsigned int init_j2_index = (n2 == n1) ? j1 : 0;
        for(unsigned int j2=init_j2_index; j2<N; j2++)
        {
          Number RB_u_euler_theta_2 = RB_u_euler_theta[n2]*RB_u_euler_theta[j2];

          Real delta = ( (n2 == n1) && (j2 == j1) ) ? 1. : 2.;

          residual_norm_sq += delta*
                              cached_theta_c*cached_theta_c*
                              RB_u_euler_theta_1*RB_u_euler_theta_2*
                              C_C_representor_norms[n1][j1][n2][j2];
        }
      }

    }
  }

  if(libmesh_real(residual_norm_sq) < 0.)
  {
    libMesh::err << "Warning: Square of residual norm is negative "
                 << "in QNTransientRBEvaluation::compute_residual_dual_norm " << std::endl;

    // Sometimes this is negative due to rounding error,
    // but this error shouldn't affect the error bound
    // too much...
//     libmesh_error();
    residual_norm_sq = std::abs(residual_norm_sq);
  }

//  libMesh::out << "slow residual_sq = " << slow_residual_norm_sq
//               << ", fast residual_sq = " << residual_norm_sq << std::endl;

  STOP_LOG("compute_residual_dual_norm()", "QNTransientRBEvaluation");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}


void QNTransientRBEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "QNTransientRBEvaluation");

  Parent::write_offline_data_to_files(directory_name);
  
  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= qn_rb_sys.get_Nmax() );

  // We use a lower precision here (only 10 digits) since
  // the trilinear form terms can take up a _lot_ of disk
  // space
  const unsigned int precision_level = 10;

  if(libMesh::processor_id() == 0)
  {
    // Write out the RB trilinear form
    std::ofstream RB_trilinear_form_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_trilinear_form.dat";
      RB_trilinear_form_out.open(file_name.str().c_str());
    }
    if ( !RB_trilinear_form_out.good() )
    {
      libMesh::err << "Error opening RB_trilinear_form.dat" << std::endl;
      libmesh_error();
    }
    RB_trilinear_form_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        for(unsigned int l=0; l<n_bfs; l++)
        {
          RB_trilinear_form_out << std::scientific << RB_trilinear_form[i][j][l] << " ";
        }
      }
    }
    RB_trilinear_form_out << std::endl;
    RB_trilinear_form_out.close();

    // Next write out the Fq_C representor norm data
    std::ofstream RB_Fq_C_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_C_norms.dat";
      RB_Fq_C_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_C_norms_out.good() )
    {
      libMesh::err << "Error opening Fq_C_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_C_norms_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
    {
      for(unsigned int n1=0; n1<n_bfs; n1++)
      {
        for(unsigned int j1=0; j1<n_bfs; j1++)
        {
          RB_Fq_C_norms_out << std::scientific << Fq_C_representor_norms[q_f][n1][j1] << " ";
        }
      }
    }
    RB_Fq_C_norms_out.close();

    // Next write out the Mq_C representor norm data
    std::ofstream RB_Mq_C_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Mq_C_norms.dat";
      RB_Mq_C_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Mq_C_norms_out.good() )
    {
      libMesh::err << "Error opening Mq_C_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Mq_C_norms_out.precision(precision_level);
    for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int n1=0; n1<n_bfs; n1++)
        {
          for(unsigned int j1=0; j1<n_bfs; j1++)
          {
            RB_Mq_C_norms_out << std::scientific << Mq_C_representor_norms[q_m][i][n1][j1] << " ";
          }
        }
      }
    }
    RB_Mq_C_norms_out.close();

    // Next write out the Aq_C representor norm data
    std::ofstream RB_Aq_C_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_C_norms.dat";
      RB_Aq_C_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Aq_C_norms_out.good() )
    {
      libMesh::err << "Error opening Aq_C_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Aq_C_norms_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int n1=0; n1<n_bfs; n1++)
        {
          for(unsigned int j1=0; j1<n_bfs; j1++)
          {
            RB_Aq_C_norms_out << std::scientific << Aq_C_representor_norms[q_a][i][n1][j1] << " ";
          }
        }
      }
    }
    RB_Aq_C_norms_out.close();

    // Next write out the C_C representor norm data
    std::ofstream RB_C_C_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/C_C_norms.dat";
      RB_C_C_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_C_C_norms_out.good() )
    {
      libMesh::err << "Error opening C_C_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_C_C_norms_out.precision(precision_level);
    for(unsigned int n1=0; n1<n_bfs; n1++)
    {
      for(unsigned int j1=0; j1<n_bfs; j1++)
      {
        for(unsigned int n2=0; n2<n_bfs; n2++)
        {
          for(unsigned int j2=0; j2<n_bfs; j2++)
          {
            RB_C_C_norms_out << std::scientific << C_C_representor_norms[n1][j1][n2][j2] << " ";
          }
        }
      }
    }
    RB_C_C_norms_out << std::endl;
    RB_C_C_norms_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "QNTransientRBEvaluation");
}

void QNTransientRBEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "QNTransientRBEvaluation");

  Parent::read_offline_data_from_files(directory_name);
  
  QNTransientRBSystem& qn_rb_sys = libmesh_cast_ref<QNTransientRBSystem&>(rb_sys);

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();

  // Read in the RB trilinear form
  std::ifstream RB_trilinear_form_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_trilinear_form.dat";
    RB_trilinear_form_in.open(file_name.str().c_str());
  }
  if ( !RB_trilinear_form_in.good() )
  {
    libMesh::err << "Error opening RB_trilinear_form.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<n_bfs; j++)
    {
      for(unsigned int l=0; l<n_bfs; l++)
      {
        RB_trilinear_form_in >> RB_trilinear_form[i][j][l];
      }
    }
  }
  RB_trilinear_form_in.close();

  // Next read in the Fq_C representor norm data
  std::ifstream RB_Fq_C_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Fq_C_norms.dat";
    RB_Fq_C_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Fq_C_norms_in.good() )
  {
    libMesh::err << "Error opening Fq_C_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<qn_rb_sys.get_Q_f(); q_f++)
  {
    for(unsigned int n1=0; n1<n_bfs; n1++)
    {
      for(unsigned int j1=0; j1<n_bfs; j1++)
      {
        RB_Fq_C_norms_in >> Fq_C_representor_norms[q_f][n1][j1];
      }
    }
  }
  RB_Fq_C_norms_in.close();

  // Next read in the Mq_C representor norm data
  std::ifstream RB_Mq_C_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Mq_C_norms.dat";
    RB_Mq_C_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Mq_C_norms_in.good() )
  {
    libMesh::err << "Error opening Mq_C_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_m=0; q_m<qn_rb_sys.get_Q_m(); q_m++)
  {
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int n1=0; n1<n_bfs; n1++)
      {
        for(unsigned int j1=0; j1<n_bfs; j1++)
        {
          RB_Mq_C_norms_in >> Mq_C_representor_norms[q_m][i][n1][j1];
        }
      }
    }
  }
  RB_Mq_C_norms_in.close();

  // Next read in the Aq_C representor norm data
  std::ifstream RB_Aq_C_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_C_norms.dat";
    RB_Aq_C_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Aq_C_norms_in.good() )
  {
    libMesh::err << "Error opening Aq_C_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<qn_rb_sys.get_Q_a(); q_a++)
  {
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int n1=0; n1<n_bfs; n1++)
      {
        for(unsigned int j1=0; j1<n_bfs; j1++)
        {
          RB_Aq_C_norms_in >> Aq_C_representor_norms[q_a][i][n1][j1];
        }
      }
    }
  }
  RB_Aq_C_norms_in.close();

  // Next read in the C_C representor norm data
  std::ifstream RB_C_C_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/C_C_norms.dat";
    RB_C_C_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_C_C_norms_in.good() )
  {
    libMesh::err << "Error opening C_C_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int n1=0; n1<n_bfs; n1++)
  {
    for(unsigned int j1=0; j1<n_bfs; j1++)
    {
      for(unsigned int n2=0; n2<n_bfs; n2++)
      {
        for(unsigned int j2=0; j2<n_bfs; j2++)
        {
          RB_C_C_norms_in >> C_C_representor_norms[n1][j1][n2][j2];
        }
      }
    }
  }
  RB_C_C_norms_in.close();

  STOP_LOG("read_offline_data_from_files()", "QNTransientRBEvaluation");
}

} // namespace libMesh
