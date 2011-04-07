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

#include "transient_rb_evaluation.h"
#include "transient_rb_system.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"
#include "xdr_cxx.h"
#include "parallel.h"

// For checking for the existence of files
#include <sys/stat.h>

#include "o_string_stream.h"
#include <fstream>
#include <sstream>

namespace libMesh
{

TransientRBEvaluation::TransientRBEvaluation (TransientRBSystem& rb_sys_in)
  :
  Parent(rb_sys_in)
{
}

void TransientRBEvaluation::clear()
{
  Parent::clear();
  
  clear_riesz_representors();
}

void TransientRBEvaluation::clear_riesz_representors()
{
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  // Delete the M_q representors
  for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
  {
    for(unsigned int i=0; i<M_q_representor[q_m].size(); i++)
    {
      if(M_q_representor[q_m][i])
      {
        delete M_q_representor[q_m][i];
        M_q_representor[q_m][i] = NULL;
      }
    }
  }
}

void TransientRBEvaluation::initialize()
{
  // Now allocate the N (i.e. RB) dependent data structures
  Parent::initialize();
  
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);
  const unsigned int Nmax = trans_rb_sys.get_Nmax();

  RB_L2_matrix.resize(Nmax,Nmax);

  // Allocate dense matrices for RB solves
  RB_M_q_vector.resize(trans_rb_sys.get_Q_m());
  for(unsigned int q=0; q<trans_rb_sys.get_Q_m(); q++)
  {
    // Initialize the memory for the RB matrices
    RB_M_q_vector[q].resize(Nmax,Nmax);
  }

  // Initialize vectors for the norms of the representors
  Fq_Mq_representor_norms.resize(trans_rb_sys.get_Q_f());
  for(unsigned int i=0; i<trans_rb_sys.get_Q_f(); i++)
  {
    Fq_Mq_representor_norms[i].resize(trans_rb_sys.get_Q_m());
    for(unsigned int j=0; j<trans_rb_sys.get_Q_m(); j++)
    {
      Fq_Mq_representor_norms[i][j].resize(Nmax);
    }
  }

  unsigned int Q_m_hat = trans_rb_sys.get_Q_m()*(trans_rb_sys.get_Q_m()+1)/2;
  Mq_Mq_representor_norms.resize(Q_m_hat);
  for(unsigned int i=0; i<Q_m_hat; i++)
  {
    Mq_Mq_representor_norms[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      Mq_Mq_representor_norms[i][j].resize(Nmax);
    }
  }

  Aq_Mq_representor_norms.resize(trans_rb_sys.get_Q_a());
  for(unsigned int i=0; i<trans_rb_sys.get_Q_a(); i++)
  {
    Aq_Mq_representor_norms[i].resize(trans_rb_sys.get_Q_m());
    for(unsigned int j=0; j<trans_rb_sys.get_Q_m(); j++)
    {
      Aq_Mq_representor_norms[i][j].resize(Nmax);
      for(unsigned int k=0; k<Nmax; k++)
      {
        Aq_Mq_representor_norms[i][j][k].resize(Nmax);
      }
    }
  }

  // Initialize the initial condition storage
  RB_initial_condition_all_N.resize(Nmax);
  for(unsigned int i=0; i<RB_initial_condition_all_N.size(); i++)
  {
    // The i^th row holds a vector of lenght i+1
    RB_initial_condition_all_N[i].resize(i+1);
  }

  initial_L2_error_all_N.resize(Nmax);

  // Resize the RB output vectors
  RB_outputs_all_k.resize(trans_rb_sys.get_n_outputs());
  RB_output_error_bounds_all_k.resize(trans_rb_sys.get_n_outputs());
  for(unsigned int n=0; n<trans_rb_sys.get_n_outputs(); n++)
  {
    RB_outputs_all_k[n].resize(trans_rb_sys.get_K()+1);
    RB_output_error_bounds_all_k[n].resize(trans_rb_sys.get_K()+1);
  }

  // Resize M_q_representor
  M_q_representor.resize(trans_rb_sys.get_Q_m());
  for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
  {
    M_q_representor[q_m].resize(trans_rb_sys.get_Nmax());
  }
}

Real TransientRBEvaluation::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "TransientRBEvaluation");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in RB_solve" << std::endl;
    libmesh_error();
  }

  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);
  const Real dt = trans_rb_sys.get_dt();

  // First assemble the mass matrix
  DenseMatrix<Number> RB_mass_matrix_N(N,N);
  RB_mass_matrix_N.zero();
  DenseMatrix<Number> RB_M_q_m;
  for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
  {
    RB_M_q_vector[q_m].get_principal_submatrix(N, RB_M_q_m);
    RB_mass_matrix_N.add(trans_rb_sys.eval_theta_q_m(q_m), RB_M_q_m);
  }

  DenseMatrix<Number> RB_LHS_matrix(N,N);
  RB_LHS_matrix.zero();

  DenseMatrix<Number> RB_RHS_matrix(N,N);
  RB_RHS_matrix.zero();

  RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);
  RB_RHS_matrix.add(1./dt, RB_mass_matrix_N);

  DenseMatrix<Number> RB_A_q_a;
  for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
  {
    RB_A_q_vector[q_a].get_principal_submatrix(N, RB_A_q_a);

    RB_LHS_matrix.add(       trans_rb_sys.get_euler_theta()*trans_rb_sys.eval_theta_q_a(q_a), RB_A_q_a);
    RB_RHS_matrix.add( -(1.-trans_rb_sys.get_euler_theta())*trans_rb_sys.eval_theta_q_a(q_a), RB_A_q_a);
  }

  // Set system time level to 0
  error_bound_all_k.resize(trans_rb_sys.get_K()+1);
  trans_rb_sys.set_time_level(0); // Sets the member variable _k to zero

  // Resize/clear the solution vector
  RB_solution.resize(N);

  // Load the initial condition into RB_solution
  if(N > 0)
  {
    RB_solution = RB_initial_condition_all_N[N-1];
  }

  // Resize/clear the old solution vector
  old_RB_solution.resize(N);

  // Initialize the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  // Initialize the vectors storing solution data
  RB_temporal_solution_data.resize(trans_rb_sys.get_K()+1);
  for(unsigned int time_level=0; time_level<=trans_rb_sys.get_K(); time_level++)
  {
    RB_temporal_solution_data[time_level].resize(N);
  }
  // and load the initial data
  RB_temporal_solution_data[0] = RB_solution;

  Real error_bound_sum = 0.;
  if(N > 0)
  {
    error_bound_sum += pow( initial_L2_error_all_N[N-1], 2.);
  }

  // Set error bound at the initial time
  error_bound_all_k[trans_rb_sys.get_time_level()] = std::sqrt(error_bound_sum);

  // Compute the outputs and associated error bounds at the initial time
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<trans_rb_sys.get_n_outputs(); n++)
  {
    RB_outputs_all_k[n][0] = 0.;
    for(unsigned int q_l=0; q_l<trans_rb_sys.get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs_all_k[n][0] += trans_rb_sys.eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_solution);
    }

    RB_output_error_bounds_all_k[n][0] = error_bound_all_k[0] * trans_rb_sys.eval_output_dual_norm(n);
  }

  Real alpha_LB = trans_rb_sys.get_SCM_lower_bound();
  
  // Precompute time-invariant parts of the dual norm of the residual.
  cache_online_residual_terms(N);

  for(unsigned int time_level=1; time_level<=trans_rb_sys.get_K(); time_level++)
  {
    trans_rb_sys.set_time_level(time_level); // This updates the member variable _k
    old_RB_solution = RB_solution;

    // Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
    RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

    // Add forcing terms
    DenseVector<Number> RB_F_q_f;
    for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
    {
      RB_F_q_vector[q_f].get_principal_subvector(N, RB_F_q_f);
      RB_rhs.add(trans_rb_sys.eval_theta_q_f(q_f), RB_F_q_f);
    }

    if(N > 0)
    {
      RB_LHS_matrix.lu_solve(RB_rhs, RB_solution);
    }

    // Save RB_solution for current time level
    RB_temporal_solution_data[time_level] = RB_solution;

    // Evaluate the dual norm of the residual for RB_solution_vector
//    Real epsilon_N = uncached_compute_residual_dual_norm(N);
    Real epsilon_N = compute_residual_dual_norm(N);

    error_bound_sum += trans_rb_sys.residual_scaling_numer(alpha_LB) * pow(epsilon_N, 2.);

    // store error bound at time-level _k
    error_bound_all_k[time_level] = std::sqrt(error_bound_sum/trans_rb_sys.residual_scaling_denom(alpha_LB));

    // Now compute the outputs and associated errors
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<trans_rb_sys.get_n_outputs(); n++)
    {
      RB_outputs_all_k[n][time_level] = 0.;
      for(unsigned int q_l=0; q_l<trans_rb_sys.get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
        RB_outputs_all_k[n][time_level] += trans_rb_sys.eval_theta_q_l(n,q_l)*RB_output_vector_N.dot(RB_solution);
      }

      RB_output_error_bounds_all_k[n][time_level] = error_bound_all_k[time_level] * trans_rb_sys.eval_output_dual_norm(n);
    }
  }

  // Now compute the L2 norm of the RB solution at time-level _K
  // to normalize the error bound
  // We reuse RB_rhs here
  DenseMatrix<Number> RB_L2_matrix_N;
  RB_L2_matrix.get_principal_submatrix(N,RB_L2_matrix_N);
  RB_L2_matrix_N.vector_mult(RB_rhs, RB_solution);
  Real final_RB_L2_norm = libmesh_real(std::sqrt(RB_solution.dot(RB_rhs)));

  STOP_LOG("RB_solve()", "TransientRBEvaluation");

   return ( trans_rb_sys.return_rel_error_bound ? error_bound_all_k[trans_rb_sys.get_K()]/final_RB_L2_norm :
                                                  error_bound_all_k[trans_rb_sys.get_K()] );
}

void TransientRBEvaluation::cache_online_residual_terms(const unsigned int N)
{
  START_LOG("cache_online_residual_terms()", "TransientRBEvaluation");

  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  cached_Fq_term = 0.;
  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<trans_rb_sys.get_Q_f(); q_f1++)
  {
    Number cached_theta_q_f1 = trans_rb_sys.eval_theta_q_f(q_f1);
    for(unsigned int q_f2=q_f1; q_f2<trans_rb_sys.get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      cached_Fq_term += delta*cached_theta_q_f1*trans_rb_sys.eval_theta_q_f(q_f2) * trans_rb_sys.Fq_representor_norms[q];

      q++;
    }
  }

  cached_Fq_Aq_vector.resize(N);
  for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = trans_rb_sys.eval_theta_q_f(q_f);
    for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
    {
      Number cached_theta_q_a = trans_rb_sys.eval_theta_q_a(q_a);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Aq_vector(i) += 2.*cached_theta_q_f*cached_theta_q_a*
                                  Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  cached_Aq_Aq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_a1=0; q_a1<trans_rb_sys.get_Q_a(); q_a1++)
  {
    Number cached_theta_q_a1 = trans_rb_sys.eval_theta_q_a(q_a1);
    for(unsigned int q_a2=q_a1; q_a2<trans_rb_sys.get_Q_a(); q_a2++)
    {
      Number cached_theta_q_a2 = trans_rb_sys.eval_theta_q_a(q_a2);
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Aq_Aq_matrix(i,j) += delta*
                                      cached_theta_q_a1*cached_theta_q_a2*
                                      Aq_Aq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  cached_Fq_Mq_vector.resize(N);
  for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = trans_rb_sys.eval_theta_q_f(q_f);
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = trans_rb_sys.eval_theta_q_m(q_m);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Mq_vector(i) += 2.*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  cached_Aq_Mq_matrix.resize(N,N);
  for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = trans_rb_sys.eval_theta_q_a(q_a);
    
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = trans_rb_sys.eval_theta_q_m(q_m);
      
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Aq_Mq_matrix(i,j) += 2.*cached_theta_q_a*cached_theta_q_m*Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  
  cached_Mq_Mq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_m1=0; q_m1<trans_rb_sys.get_Q_m(); q_m1++)
  {
    Number cached_theta_q_m1 = trans_rb_sys.eval_theta_q_m(q_m1);
    for(unsigned int q_m2=q_m1; q_m2<trans_rb_sys.get_Q_m(); q_m2++)
    {
      Number cached_theta_q_m2 = trans_rb_sys.eval_theta_q_m(q_m2);
      Real delta = (q_m1==q_m2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          cached_Mq_Mq_matrix(i,j) += delta*
                                      cached_theta_q_m1*cached_theta_q_m2*
                                      Mq_Mq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }
  
  STOP_LOG("cache_online_residual_terms()", "TransientRBEvaluation");
}

Real TransientRBEvaluation::compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("compute_residual_dual_norm()", "TransientRBEvaluation");

  // This assembly assumes we have already called cache_online_residual_terms
  // and that the RB_solve parameter is constant in time
  
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  const Real dt = trans_rb_sys.get_dt();

  DenseVector<Number> RB_u_euler_theta(N);
  DenseVector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta(i)  = trans_rb_sys.get_euler_theta()*RB_solution(i) +
                          (1.-trans_rb_sys.get_euler_theta())*old_RB_solution(i);
    mass_coeffs(i) = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  Number residual_norm_sq = cached_Fq_term;
  
  residual_norm_sq += RB_u_euler_theta.dot(cached_Fq_Aq_vector);
  residual_norm_sq += mass_coeffs.dot(cached_Fq_Mq_vector);

  for(unsigned int i=0; i<N; i++)
    for(unsigned int j=0; j<N; j++)
    {
      residual_norm_sq += RB_u_euler_theta(i)*RB_u_euler_theta(j)*cached_Aq_Aq_matrix(i,j);
      residual_norm_sq += mass_coeffs(i)*mass_coeffs(j)*cached_Mq_Mq_matrix(i,j);
      residual_norm_sq += RB_u_euler_theta(i)*mass_coeffs(j)*cached_Aq_Mq_matrix(i,j);
    }


  if(libmesh_real(residual_norm_sq) < 0)
  {
    libMesh::out << "Warning: Square of residual norm is negative "
                 << "in TransientRBEvaluation::compute_residual_dual_norm()" << std::endl;

    // Sometimes this is negative due to rounding error,
    // but error is on the order of 1.e-10, so shouldn't
    // affect result
//    libmesh_error();
     residual_norm_sq = std::abs(residual_norm_sq);
  }

  STOP_LOG("compute_residual_dual_norm()", "TransientRBEvaluation");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}

Real TransientRBEvaluation::uncached_compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("uncached_compute_residual_dual_norm()", "TransientRBEvaluation");

  // Use the stored representor inner product values
  // to evaluate the residual norm
  
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  const Real dt = trans_rb_sys.get_dt();

  std::vector<Number> RB_u_euler_theta(N);
  std::vector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta[i]  = trans_rb_sys.get_euler_theta()*RB_solution(i) +
                          (1.-trans_rb_sys.get_euler_theta())*old_RB_solution(i);
    mass_coeffs[i] = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<trans_rb_sys.get_Q_f(); q_f1++)
  {
    Number cached_theta_q_f1 = trans_rb_sys.eval_theta_q_f(q_f1);
    for(unsigned int q_f2=q_f1; q_f2<trans_rb_sys.get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      residual_norm_sq += delta*cached_theta_q_f1*trans_rb_sys.eval_theta_q_f(q_f2) * trans_rb_sys.Fq_representor_norms[q];

      q++;
    }
  }

  for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = trans_rb_sys.eval_theta_q_f(q_f);
    for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
    {
      Number cached_theta_q_a = trans_rb_sys.eval_theta_q_a(q_a);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*RB_u_euler_theta[i]*cached_theta_q_f*cached_theta_q_a*
                               Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  q=0;
  for(unsigned int q_a1=0; q_a1<trans_rb_sys.get_Q_a(); q_a1++)
  {
    Number cached_theta_q_a1 = trans_rb_sys.eval_theta_q_a(q_a1);
    for(unsigned int q_a2=q_a1; q_a2<trans_rb_sys.get_Q_a(); q_a2++)
    {
      Number cached_theta_q_a2 = trans_rb_sys.eval_theta_q_a(q_a2);
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += delta*RB_u_euler_theta[i]*RB_u_euler_theta[j]*
                              cached_theta_q_a1*cached_theta_q_a2*
                              Aq_Aq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  // Now add the terms due to the time-derivative
  q=0;
  for(unsigned int q_m1=0; q_m1<trans_rb_sys.get_Q_m(); q_m1++)
  {
    Number cached_theta_q_m1 = trans_rb_sys.eval_theta_q_m(q_m1);
    for(unsigned int q_m2=q_m1; q_m2<trans_rb_sys.get_Q_m(); q_m2++)
    {
      Number cached_theta_q_m2 = trans_rb_sys.eval_theta_q_m(q_m2);
      Real delta = (q_m1==q_m2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += delta*mass_coeffs[i]*mass_coeffs[j]*
                              cached_theta_q_m1*cached_theta_q_m2*
                              Mq_Mq_representor_norms[q][i][j];
        }
      }
      q++;
    }
  }

  for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
  {
    Number cached_theta_q_f = trans_rb_sys.eval_theta_q_f(q_f);
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = trans_rb_sys.eval_theta_q_m(q_m);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*mass_coeffs[i]*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
  {
    Number cached_theta_q_a = trans_rb_sys.eval_theta_q_a(q_a);
    
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      Number cached_theta_q_m = trans_rb_sys.eval_theta_q_m(q_m);
      
      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq += 2.*RB_u_euler_theta[i]*mass_coeffs[j]*
                                         cached_theta_q_a*cached_theta_q_m*
                                         Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  
  if(libmesh_real(residual_norm_sq) < 0)
  {
    libMesh::out << "Warning: Square of residual norm is negative "
                 << "in TransientRBEvaluation::compute_residual_dual_norm()" << std::endl;

    // Sometimes this is negative due to rounding error,
    // but error is on the order of 1.e-10, so shouldn't
    // affect result
//    libmesh_error();
     residual_norm_sq = std::abs(residual_norm_sq);
  }

//   libMesh::out << "slow residual_sq = " << slow_residual_norm_sq
//                << ", fast residual_sq = " << residual_norm_sq << std::endl;

  STOP_LOG("uncached_compute_residual_dual_norm()", "TransientRBEvaluation");

  return libmesh_real(std::sqrt( residual_norm_sq ));
}

void TransientRBEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "TransientRBEvaluation");

  Parent::write_offline_data_to_files(directory_name);
  
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= trans_rb_sys.get_Nmax() );

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Write out the L2 matrix
    std::ofstream RB_L2_matrix_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_L2_matrix.dat";
      RB_L2_matrix_out.open(file_name.str().c_str());
    }
    RB_L2_matrix_out.precision(precision_level);
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        RB_L2_matrix_out << std::scientific << RB_L2_matrix(i,j) << " ";
      }
    }
    RB_L2_matrix_out.close();
          
    // Write out the M_q matrices
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_M_";
      OSSRealzeroright(file_name,3,0,q_m);
      file_name << ".dat";
      std::ofstream RB_M_q_m_out(file_name.str().c_str());

      if ( !RB_M_q_m_out.good() )
      {
        libMesh::err << "Error opening RB_M_" << q_m << ".dat" << std::endl;
        libmesh_error();
      }

      RB_M_q_m_out.precision(precision_level);
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_M_q_m_out << std::scientific << RB_M_q_vector[q_m](i,j) << " ";
        }
      }
      RB_M_q_m_out.close();
    }

    // Write out the initial condition data
    // and the initial L2 error for all N
    std::ofstream initial_conditions_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/initial_conditions.dat";
      initial_conditions_out.open(file_name.str().c_str());
    }
    std::ofstream initial_L2_error_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/initial_L2_error.dat";
      initial_L2_error_out.open(file_name.str().c_str());
    }
    if (!initial_conditions_out.good() || !initial_L2_error_out.good())
    {
      libMesh::err << "Error opening initial conditions output files" << std::endl;
      libmesh_error();
    }

    initial_conditions_out.precision(precision_level);
    initial_L2_error_out.precision(precision_level);

    for(unsigned int i=0; i<n_bfs; i++)
    {
      initial_L2_error_out << initial_L2_error_all_N[i] << " ";
      for(unsigned int j=0; j<=i; j++)
      {
        initial_conditions_out << std::scientific << RB_initial_condition_all_N[i](j) << " ";
      }
    }
    initial_conditions_out.close();
    initial_L2_error_out.close();

    // Next write out the Fq_Mq representor norm data
    std::ofstream RB_Fq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Fq_Mq_norms.dat";
      RB_Fq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Fq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening Fq_Mq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_Mq_norms_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
    {
      for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
      {
        for(unsigned int i=0; i<n_bfs; i++)
        {
          RB_Fq_Mq_norms_out << std::scientific << Fq_Mq_representor_norms[q_f][q_m][i] << " ";
        }
      }
    }
    RB_Fq_Mq_norms_out.close();

    // Next write out the Mq_Mq representor norm data
    std::ofstream RB_Mq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Mq_Mq_norms.dat";
      RB_Mq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Mq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening RB_Mq_Mq_norms_out.dat" << std::endl;
      libmesh_error();
    }
    RB_Mq_Mq_norms_out.precision(precision_level);
    unsigned int Q_m_hat = trans_rb_sys.get_Q_m()*(trans_rb_sys.get_Q_m()+1)/2;
    for(unsigned int q=0; q<Q_m_hat; q++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_Mq_Mq_norms_out << std::scientific << Mq_Mq_representor_norms[q][i][j] << " ";
        }
      }
    }
    RB_Mq_Mq_norms_out.close();

    // Next write out the Aq_Mq representor norm data
    std::ofstream RB_Aq_Mq_norms_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/Aq_Mq_norms.dat";
      RB_Aq_Mq_norms_out.open(file_name.str().c_str());
    }
    if ( !RB_Aq_Mq_norms_out.good() )
    {
      libMesh::err << "Error opening Aq_Mq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Aq_Mq_norms_out.precision(precision_level);
    for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
    {
      for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
      {
        for(unsigned int i=0; i<n_bfs; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
          {
            RB_Aq_Mq_norms_out << std::scientific << Aq_Mq_representor_norms[q_a][q_m][i][j] << " ";
          }
        }
      }
    }
    RB_Aq_Mq_norms_out.close();
  }

  // Write out the residual representors to file if requested
  if (rb_sys.store_representors)
  {
    // Write out the M_q_representors.  These are useful to have when restarting,
    // so you don't have to recompute them all over again.  There should be
    // this->get_n_basis_functions() of these.
    libMesh::out << "Writing out the M_q_representors..." << std::endl;

    std::ostringstream file_name;
    const std::string residual_representor_suffix = (rb_sys.write_binary_residual_representors ? ".xdr" : ".dat");

    // Residual representors written out to their own separate directory
    std::string residual_representors_dir = "residual_representors";

    const unsigned int istop  = this->get_n_basis_functions();
    const unsigned int istart = istop-rb_sys.get_delta_N();

    for (unsigned int q=0; q<M_q_representor.size(); ++q)
      for (unsigned int i=istart; i<istop; ++i)
      {
	libMesh::out << "Writing out M_q_representor[" << q << "][" << i << "]..." << std::endl;
	libmesh_assert(M_q_representor[q][i] != NULL);

	file_name.str(""); // reset filename
	file_name << residual_representors_dir << "/M_q_representor" << i << residual_representor_suffix;

	{
	  // No need to copy!
	  //*solution = *(M_q_representor[q][i]);
	  M_q_representor[q][i]->swap(*rb_sys.solution);

	  Xdr mr_data(file_name.str(),
	              rb_sys.write_binary_residual_representors ? ENCODE : WRITE);

          rb_sys.write_serialized_data(mr_data, false);

	  // Synchronize before moving on
	  Parallel::barrier();

	  // Swap back.
	  M_q_representor[q][i]->swap(*rb_sys.solution);

	  // TODO: bzip the resulting file?  See $LIBMESH_DIR/src/mesh/unstructured_mesh.C
	  // for the system call, be sure to do it only on one processor, etc.
	}
      }
  } // end if store_representors

  STOP_LOG("write_offline_data_to_files()", "TransientRBEvaluation");
}

void TransientRBEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "TransientRBEvaluation");

  Parent::read_offline_data_from_files(directory_name);
  
  TransientRBSystem& trans_rb_sys = libmesh_cast_ref<TransientRBSystem&>(rb_sys);

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();
  
  std::ifstream RB_L2_matrix_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_L2_matrix.dat";
    RB_L2_matrix_in.open(file_name.str().c_str());
  }
  for(unsigned int i=0; i<n_bfs; i++)
  {
    for(unsigned int j=0; j<n_bfs; j++)
    {
      Number value;
      RB_L2_matrix_in >> value;
      RB_L2_matrix(i,j) = value;
    }
  }
  RB_L2_matrix_in.close();

  // Read in the M_q matrices
  for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_M_";
    OSSRealzeroright(file_name,3,0,q_m);
    file_name << ".dat";
    std::ifstream RB_M_q_m_in(file_name.str().c_str());

    if ( !RB_M_q_m_in.good() )
    {
      libMesh::err << "Error opening RB_M_" << q_m << ".dat" << std::endl;
      libmesh_error();
    }

    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        Number value;
        RB_M_q_m_in >> value;
        RB_M_q_vector[q_m](i,j) = value;
      }
    }
    RB_M_q_m_in.close();
  }


  // Read in the initial condition data
  // and the initial L2 error for all N
  std::ifstream initial_conditions_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/initial_conditions.dat";
    initial_conditions_in.open(file_name.str().c_str());
  }
  std::ifstream initial_L2_error_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/initial_L2_error.dat";
    initial_L2_error_in.open(file_name.str().c_str());
  }
  if (!initial_conditions_in.good() || !initial_L2_error_in.good())
  {
    libMesh::err << "Error opening initial conditions output files" << std::endl;
    libmesh_error();
  }

  for(unsigned int i=0; i<n_bfs; i++)
  {
    initial_L2_error_in >> initial_L2_error_all_N[i];
    for(unsigned int j=0; j<=i; j++)
    {
      initial_conditions_in >> RB_initial_condition_all_N[i](j);
    }
  }
  initial_conditions_in.close();
  initial_L2_error_in.close();

  // Next read in the Fq_Mq representor norm data
  std::ifstream RB_Fq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Fq_Mq_norms.dat";
    RB_Fq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Fq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening Fq_Mq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<trans_rb_sys.get_Q_f(); q_f++)
  {
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        RB_Fq_Mq_norms_in >> Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }
  RB_Fq_Mq_norms_in.close();

  // Next read in the Mq_Mq representor norm data
  std::ifstream RB_Mq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Mq_Mq_norms.dat";
    RB_Mq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Mq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening RB_Mq_Mq_norms_in.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_m_hat = trans_rb_sys.get_Q_m()*(trans_rb_sys.get_Q_m()+1)/2;
  for(unsigned int q=0; q<Q_m_hat; q++)
  {
    for(unsigned int i=0; i<n_bfs; i++)
    {
      for(unsigned int j=0; j<n_bfs; j++)
      {
        RB_Mq_Mq_norms_in >> Mq_Mq_representor_norms[q][i][j];
      }
    }
  }
  RB_Mq_Mq_norms_in.close();

  // Next read in the Aq_Mq representor norm data
  std::ifstream RB_Aq_Mq_norms_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/Aq_Mq_norms.dat";
    RB_Aq_Mq_norms_in.open(file_name.str().c_str());
  }
  if ( !RB_Aq_Mq_norms_in.good() )
  {
    libMesh::err << "Error opening Aq_Mq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_a=0; q_a<trans_rb_sys.get_Q_a(); q_a++)
  {
    for(unsigned int q_m=0; q_m<trans_rb_sys.get_Q_m(); q_m++)
    {
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_Aq_Mq_norms_in >> Aq_Mq_representor_norms[q_a][q_m][i][j];
        }
      }
    }
  }
  RB_Aq_Mq_norms_in.close();
  
  // Read in the representors if requested
  if (rb_sys.store_representors)
  {
    const std::string residual_representors_dir = "residual_representors";
    const std::string residual_representor_suffix =
      (rb_sys.read_binary_residual_representors ? ".xdr" : ".dat");

    std::ostringstream file_name;
    struct stat stat_info;
      

    libMesh::out << "Reading in the M_q_representors..." << std::endl;

    // Read in the A_q representors.  The class makes room for [Q_m][Nmax] of these.  We are going to
    // read in [Q_m][this->get_n_basis_functions()].  FIXME:
    // should we be worried about leaks in the locations where we're about to fill entries?
    for (unsigned int i=0; i<M_q_representor.size(); ++i)
      for (unsigned int j=0; j<M_q_representor[i].size(); ++j)
      {
        if (M_q_representor[i][j] != NULL)
        {
          libMesh::out << "Error, must delete existing M_q_representor before reading in from file."
	  	       << std::endl;
          libmesh_error();
        }
      }

    // Now ready to read them in from file!
    for (unsigned int i=0; i<M_q_representor.size(); ++i)
      for (unsigned int j=0; j<this->get_n_basis_functions(); ++j)
      {
        file_name.str(""); // reset filename
        file_name << residual_representors_dir
		  << "/M_q_representor" << i << "_" << j << residual_representor_suffix;

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

	Xdr aqr_data(file_name.str(),
	             rb_sys.read_binary_residual_representors ? DECODE : READ);

	rb_sys.read_serialized_data(aqr_data, false);

	M_q_representor[i][j] = NumericVector<Number>::build().release();
	M_q_representor[i][j]->init (rb_sys.n_dofs(), rb_sys.n_local_dofs(),
	                 	     false, libMeshEnums::PARALLEL);

	// No need to copy, just swap
	//*M_q_representor[i][j] = *solution;
	M_q_representor[i][j]->swap(*rb_sys.solution);
      }
  } // end if (store_representors)

  STOP_LOG("read_offline_data_from_files()", "TransientRBEvaluation");

}

}