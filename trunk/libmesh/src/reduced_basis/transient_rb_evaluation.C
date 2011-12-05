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
#include "transient_rb_theta_expansion.h"
#include "numeric_vector.h"
#include "libmesh_logging.h"
#include "xdr_cxx.h"
#include "parallel.h"

#include "o_string_stream.h"
#include <fstream>
#include <sstream>

namespace libMesh
{

TransientRBEvaluation::TransientRBEvaluation ()
{
  // Indicate that we need to compute the RB
  // inner product matrix in this case
  compute_RB_inner_product = true;
}

TransientRBEvaluation::~TransientRBEvaluation ()
{
  clear();
}

void TransientRBEvaluation::clear()
{
  Parent::clear();

  clear_riesz_representors();
}

void TransientRBEvaluation::clear_riesz_representors()
{
  Parent::clear_riesz_representors();

  // Delete the M_q representors
  for(unsigned int q_m=0; q_m<M_q_representor.size(); q_m++)
  {
    for(unsigned int i=0; i<M_q_representor[q_m].size(); i++)
    {
      if(M_q_representor[q_m][i])
      {
        M_q_representor[q_m][i]->clear();
        delete M_q_representor[q_m][i];
        M_q_representor[q_m][i] = NULL;
      }
    }
  }
}

void TransientRBEvaluation::resize_data_structures(const unsigned int Nmax)
{
  START_LOG("resize_data_structures()", "TransientRBEvaluation");

  Parent::resize_data_structures(Nmax);

  RB_L2_matrix.resize(Nmax,Nmax);

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

  // Allocate dense matrices for RB solves
  RB_M_q_vector.resize(Q_m);
  for(unsigned int q=0; q<Q_m; q++)
  {
    // Initialize the memory for the RB matrices
    RB_M_q_vector[q].resize(Nmax,Nmax);
  }

  // Initialize vectors for the norms of the representors
  Fq_Mq_representor_norms.resize(Q_f);
  for(unsigned int i=0; i<Q_f; i++)
  {
    Fq_Mq_representor_norms[i].resize(Q_m);
    for(unsigned int j=0; j<Q_m; j++)
    {
      Fq_Mq_representor_norms[i][j].resize(Nmax, 0.);
    }
  }

  unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
  Mq_Mq_representor_norms.resize(Q_m_hat);
  for(unsigned int i=0; i<Q_m_hat; i++)
  {
    Mq_Mq_representor_norms[i].resize(Nmax);
    for(unsigned int j=0; j<Nmax; j++)
    {
      Mq_Mq_representor_norms[i][j].resize(Nmax, 0.);
    }
  }

  Aq_Mq_representor_norms.resize(Q_a);
  for(unsigned int i=0; i<Q_a; i++)
  {
    Aq_Mq_representor_norms[i].resize(Q_m);
    for(unsigned int j=0; j<Q_m; j++)
    {
      Aq_Mq_representor_norms[i][j].resize(Nmax);
      for(unsigned int k=0; k<Nmax; k++)
      {
        Aq_Mq_representor_norms[i][j][k].resize(Nmax, 0.);
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

  initial_L2_error_all_N.resize(Nmax, 0.);

  // Resize M_q_representor
  // This is cleared in the call to clear_riesz_representors
  // in Parent::resize_RB_data, so just resize here
  M_q_representor.resize(Q_m);
  for(unsigned int q_m=0; q_m<Q_m; q_m++)
  {
    M_q_representor[q_m].resize(Nmax);
  }

  STOP_LOG("resize_data_structures()", "TransientRBEvaluation");
}

Real TransientRBEvaluation::rb_solve(unsigned int N)
{
  START_LOG("rb_solve()", "TransientRBEvaluation");

  if(N > get_n_basis_functions())
  {
    libMesh::err << "ERROR: N cannot be larger than the number "
                 << "of basis functions in rb_solve" << std::endl;
    libmesh_error();
  }

  const std::vector<Real> mu = get_current_parameters();

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

  const unsigned int n_time_steps = temporal_discretization.get_n_time_steps();
  const Real dt                   = temporal_discretization.get_delta_t();
  const Real euler_theta          = temporal_discretization.get_euler_theta();

  // Resize the RB and error bound vectors
  error_bound_all_k.resize(n_time_steps+1);
  RB_outputs_all_k.resize(trans_theta_expansion.get_n_outputs());
  RB_output_error_bounds_all_k.resize(trans_theta_expansion.get_n_outputs());
  for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
  {
    RB_outputs_all_k[n].resize(n_time_steps+1, 0.);
    RB_output_error_bounds_all_k[n].resize(n_time_steps+1, 0.);
  }

  // First assemble the mass matrix
  DenseMatrix<Number> RB_mass_matrix_N(N,N);
  RB_mass_matrix_N.zero();
  DenseMatrix<Number> RB_M_q_m;
  for(unsigned int q_m=0; q_m<Q_m; q_m++)
  {
    RB_M_q_vector[q_m].get_principal_submatrix(N, RB_M_q_m);
    RB_mass_matrix_N.add(trans_theta_expansion.eval_theta_q_m(q_m, mu), RB_M_q_m);
  }

  DenseMatrix<Number> RB_LHS_matrix(N,N);
  RB_LHS_matrix.zero();

  DenseMatrix<Number> RB_RHS_matrix(N,N);
  RB_RHS_matrix.zero();

  RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);
  RB_RHS_matrix.add(1./dt, RB_mass_matrix_N);

  DenseMatrix<Number> RB_A_q_a;
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    RB_A_q_vector[q_a].get_principal_submatrix(N, RB_A_q_a);

    RB_LHS_matrix.add(       euler_theta*trans_theta_expansion.eval_theta_q_a(q_a,mu), RB_A_q_a);
    RB_RHS_matrix.add( -(1.-euler_theta)*trans_theta_expansion.eval_theta_q_a(q_a,mu), RB_A_q_a);
  }

  // Set system time level to 0
  temporal_discretization.set_time_step(0);

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
  RB_temporal_solution_data.resize(n_time_steps+1);
  for(unsigned int time_level=0; time_level<=n_time_steps; time_level++)
  {
    RB_temporal_solution_data[time_level].resize(N);
  }
  // and load the initial data
  RB_temporal_solution_data[0] = RB_solution;

  // Set outputs at initial time
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
  {
    RB_outputs_all_k[n][0] = 0.;
    for(unsigned int q_l=0; q_l<trans_theta_expansion.get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs_all_k[n][0] += trans_theta_expansion.eval_theta_q_l(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
    }
  }

  // Initialize error bounds, if necessary
  Real error_bound_sum = 0.;
  Real alpha_LB = 0.;
  if(evaluate_RB_error_bound)
  {
    if(N > 0)
    {
      error_bound_sum += pow( initial_L2_error_all_N[N-1], 2.);
    }

    // Set error bound at the initial time
    error_bound_all_k[temporal_discretization.get_time_step()] = std::sqrt(error_bound_sum);

    // Compute the outputs and associated error bounds at the initial time
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
    {
      RB_outputs_all_k[n][0] = 0.;
      for(unsigned int q_l=0; q_l<trans_theta_expansion.get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
        RB_outputs_all_k[n][0] += trans_theta_expansion.eval_theta_q_l(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
      }

      RB_output_error_bounds_all_k[n][0] = error_bound_all_k[0] * eval_output_dual_norm(n,mu);
    }

    alpha_LB = get_stability_lower_bound();

    // Precompute time-invariant parts of the dual norm of the residual.
    cache_online_residual_terms(N);
  }

  for(unsigned int time_level=1; time_level<=n_time_steps; time_level++)
  {
    temporal_discretization.set_time_step(time_level);
    old_RB_solution = RB_solution;

    // Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
    RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

    // Add forcing terms
    DenseVector<Number> RB_F_q_f;
    for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      RB_F_q_vector[q_f].get_principal_subvector(N, RB_F_q_f);
      RB_rhs.add(trans_theta_expansion.eval_theta_q_f(q_f,mu), RB_F_q_f);
    }

    if(N > 0)
    {
      RB_LHS_matrix.lu_solve(RB_rhs, RB_solution);
    }

    // Save RB_solution for current time level
    RB_temporal_solution_data[time_level] = RB_solution;

    // Evaluate outputs
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
    {
      RB_outputs_all_k[n][time_level] = 0.;
      for(unsigned int q_l=0; q_l<trans_theta_expansion.get_Q_l(n); q_l++)
      {
        RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
        RB_outputs_all_k[n][time_level] += trans_theta_expansion.eval_theta_q_l(n,q_l,mu)*
                                           RB_output_vector_N.dot(RB_solution);
      }
    }

    // Calculate RB error bounds
    if(evaluate_RB_error_bound)
    {
      // Evaluate the dual norm of the residual for RB_solution_vector
      // Real epsilon_N = uncached_compute_residual_dual_norm(N);
      Real epsilon_N = compute_residual_dual_norm(N);

      error_bound_sum += residual_scaling_numer(alpha_LB) * pow(epsilon_N, 2.);

      // store error bound at time-level _k
      error_bound_all_k[time_level] = std::sqrt(error_bound_sum/residual_scaling_denom(alpha_LB));

      // Now evaluated output error bounds
      for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
      {
        RB_output_error_bounds_all_k[n][time_level] = error_bound_all_k[time_level] *
                                                      eval_output_dual_norm(n,mu);
      }
    }
  }

  STOP_LOG("rb_solve()", "TransientRBEvaluation");

  if(evaluate_RB_error_bound) // Calculate the error bounds
  {
    return error_bound_all_k[n_time_steps];
  }
  else // Don't calculate the error bounds
  {
    // Just return -1. if we did not compute the error bound
    return -1.;
  }
}

Real TransientRBEvaluation::get_rb_solution_norm()
{
  // Return the L2 norm of RB_solution
  // After an rb_solve, RB_solution will hold the
  // solution vector for the final time level.

  const unsigned int N = RB_solution.size();
  DenseVector<Number> temp(N);
  DenseMatrix<Number> RB_L2_matrix_N;
  RB_L2_matrix.get_principal_submatrix(N, RB_L2_matrix_N);
  RB_L2_matrix_N.vector_mult(temp, RB_solution);

  return libmesh_real(std::sqrt(RB_solution.dot(temp)));
}

Real TransientRBEvaluation::residual_scaling_numer(Real)
{
  return temporal_discretization.get_delta_t();
}

void TransientRBEvaluation::cache_online_residual_terms(const unsigned int N)
{
  START_LOG("cache_online_residual_terms()", "TransientRBEvaluation");

  const std::vector<Real> mu = get_current_parameters();

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

  cached_Fq_term = 0.;
  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<Q_f; q_f1++)
  {
    Number cached_theta_q_f1 = trans_theta_expansion.eval_theta_q_f(q_f1,mu);
    for(unsigned int q_f2=q_f1; q_f2<Q_f; q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      cached_Fq_term += delta*cached_theta_q_f1*trans_theta_expansion.eval_theta_q_f(q_f2,mu) *
                        Fq_representor_norms[q];

      q++;
    }
  }

  cached_Fq_Aq_vector.resize(N);
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
  {
    Number cached_theta_q_f = trans_theta_expansion.eval_theta_q_f(q_f,mu);
    for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      Number cached_theta_q_a = trans_theta_expansion.eval_theta_q_a(q_a,mu);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Aq_vector(i) += 2.*cached_theta_q_f*cached_theta_q_a*
                                  Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  cached_Aq_Aq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_a1=0; q_a1<Q_a; q_a1++)
  {
    Number cached_theta_q_a1 = trans_theta_expansion.eval_theta_q_a(q_a1,mu);
    for(unsigned int q_a2=q_a1; q_a2<Q_a; q_a2++)
    {
      Number cached_theta_q_a2 = trans_theta_expansion.eval_theta_q_a(q_a2,mu);
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
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
  {
    Number cached_theta_q_f = trans_theta_expansion.eval_theta_q_f(q_f,mu);
    for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      Number cached_theta_q_m = trans_theta_expansion.eval_theta_q_m(q_m,mu);
      for(unsigned int i=0; i<N; i++)
      {
        cached_Fq_Mq_vector(i) += 2.*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  cached_Aq_Mq_matrix.resize(N,N);
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    Number cached_theta_q_a = trans_theta_expansion.eval_theta_q_a(q_a,mu);

    for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      Number cached_theta_q_m = trans_theta_expansion.eval_theta_q_m(q_m,mu);

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
  for(unsigned int q_m1=0; q_m1<Q_m; q_m1++)
  {
    Number cached_theta_q_m1 = trans_theta_expansion.eval_theta_q_m(q_m1,mu);
    for(unsigned int q_m2=q_m1; q_m2<Q_m; q_m2++)
    {
      Number cached_theta_q_m2 = trans_theta_expansion.eval_theta_q_m(q_m2,mu);
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
  // and that the rb_solve parameter is constant in time

  const Real dt          = temporal_discretization.get_delta_t();
  const Real euler_theta = temporal_discretization.get_euler_theta();

  DenseVector<Number> RB_u_euler_theta(N);
  DenseVector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta(i)  = euler_theta*RB_solution(i) +
                          (1.-euler_theta)*old_RB_solution(i);
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

  const std::vector<Real> mu = get_current_parameters();

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

  const Real dt          = temporal_discretization.get_delta_t();
  const Real euler_theta = temporal_discretization.get_euler_theta();

  std::vector<Number> RB_u_euler_theta(N);
  std::vector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
  {
    RB_u_euler_theta[i]  = euler_theta*RB_solution(i) +
                          (1.-euler_theta)*old_RB_solution(i);
    mass_coeffs[i] = -(RB_solution(i) - old_RB_solution(i))/dt;
  }

  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<Q_f; q_f1++)
  {
    Number cached_theta_q_f1 = trans_theta_expansion.eval_theta_q_f(q_f1,mu);
    for(unsigned int q_f2=q_f1; q_f2<Q_f; q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      residual_norm_sq += delta*cached_theta_q_f1*trans_theta_expansion.eval_theta_q_f(q_f2,mu) * Fq_representor_norms[q];

      q++;
    }
  }

  for(unsigned int q_f=0; q_f<Q_f; q_f++)
  {
    Number cached_theta_q_f = trans_theta_expansion.eval_theta_q_f(q_f,mu);
    for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      Number cached_theta_q_a = trans_theta_expansion.eval_theta_q_a(q_a,mu);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*RB_u_euler_theta[i]*cached_theta_q_f*cached_theta_q_a*
                               Fq_Aq_representor_norms[q_f][q_a][i];
      }
    }
  }

  q=0;
  for(unsigned int q_a1=0; q_a1<Q_a; q_a1++)
  {
    Number cached_theta_q_a1 = trans_theta_expansion.eval_theta_q_a(q_a1,mu);
    for(unsigned int q_a2=q_a1; q_a2<Q_a; q_a2++)
    {
      Number cached_theta_q_a2 = trans_theta_expansion.eval_theta_q_a(q_a2,mu);
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
  for(unsigned int q_m1=0; q_m1<Q_m; q_m1++)
  {
    Number cached_theta_q_m1 = trans_theta_expansion.eval_theta_q_m(q_m1,mu);
    for(unsigned int q_m2=q_m1; q_m2<Q_m; q_m2++)
    {
      Number cached_theta_q_m2 = trans_theta_expansion.eval_theta_q_m(q_m2,mu);
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

  for(unsigned int q_f=0; q_f<Q_f; q_f++)
  {
    Number cached_theta_q_f = trans_theta_expansion.eval_theta_q_f(q_f,mu);
    for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      Number cached_theta_q_m = trans_theta_expansion.eval_theta_q_m(q_m,mu);
      for(unsigned int i=0; i<N; i++)
      {
        residual_norm_sq += 2.*mass_coeffs[i]*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_norms[q_f][q_m][i];
      }
    }
  }

  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    Number cached_theta_q_a = trans_theta_expansion.eval_theta_q_a(q_a,mu);

    for(unsigned int q_m=0; q_m<Q_m; q_m++)
    {
      Number cached_theta_q_m = trans_theta_expansion.eval_theta_q_m(q_m,mu);

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

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

  const unsigned int n_bfs = get_n_basis_functions();

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
    for(unsigned int q_m=0; q_m<Q_m; q_m++)
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
    for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
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
    unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
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
    for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
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

  STOP_LOG("write_offline_data_to_files()", "TransientRBEvaluation");
}

void TransientRBEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "TransientRBEvaluation");

  Parent::read_offline_data_from_files(directory_name);

  TransientRBThetaExpansion& trans_theta_expansion =
    libmesh_cast_ref<TransientRBThetaExpansion&>(*rb_theta_expansion);
  const unsigned int Q_m = trans_theta_expansion.get_Q_m();
  const unsigned int Q_a = trans_theta_expansion.get_Q_a();
  const unsigned int Q_f = trans_theta_expansion.get_Q_f();

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
  for(unsigned int q_m=0; q_m<Q_m; q_m++)
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
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
  {
    for(unsigned int q_m=0; q_m<Q_m; q_m++)
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
  unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
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
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
  {
    for(unsigned int q_m=0; q_m<Q_m; q_m++)
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

  STOP_LOG("read_offline_data_from_files()", "TransientRBEvaluation");

}

}
