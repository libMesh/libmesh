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
#include "libmesh/transient_rb_evaluation.h"
#include "libmesh/transient_rb_theta_expansion.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/parallel.h"
#include "libmesh/getpot.h"

// C++ includes
#include <fstream>
#include <sstream>
#include <iomanip>

namespace libMesh
{

TransientRBEvaluation::TransientRBEvaluation(const Parallel::Communicator & comm_in) :
  RBEvaluation(comm_in),
  _rb_solve_data_cached(false)
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
  for (std::size_t q_m=0; q_m<M_q_representor.size(); q_m++)
    for (std::size_t i=0; i<M_q_representor[q_m].size(); i++)
      {
        if (M_q_representor[q_m][i])
          {
            M_q_representor[q_m][i]->clear();
            delete M_q_representor[q_m][i];
            M_q_representor[q_m][i] = libmesh_nullptr;
          }
      }
}

void TransientRBEvaluation::resize_data_structures(const unsigned int Nmax,
                                                   bool resize_error_bound_data)
{
  LOG_SCOPE("resize_data_structures()", "TransientRBEvaluation");

  Parent::resize_data_structures(Nmax, resize_error_bound_data);

  RB_L2_matrix.resize(Nmax,Nmax);
  RB_LHS_matrix.resize(Nmax,Nmax);
  RB_RHS_matrix.resize(Nmax,Nmax);
  RB_RHS_save.resize(Nmax);

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  // Allocate dense matrices for RB solves
  RB_M_q_vector.resize(Q_m);
  for(unsigned int q=0; q<Q_m; q++)
    {
      // Initialize the memory for the RB matrices
      RB_M_q_vector[q].resize(Nmax,Nmax);
    }

  // Initialize the initial condition storage
  RB_initial_condition_all_N.resize(Nmax);
  for (std::size_t i=0; i<RB_initial_condition_all_N.size(); i++)
    {
      // The i^th row holds a vector of lenght i+1
      RB_initial_condition_all_N[i].resize(i+1);
    }

  initial_L2_error_all_N.resize(Nmax, 0.);


  if(resize_error_bound_data)
    {
      // Initialize vectors for the norms of the representors
      Fq_Mq_representor_innerprods.resize(Q_f);
      for(unsigned int i=0; i<Q_f; i++)
        {
          Fq_Mq_representor_innerprods[i].resize(Q_m);
          for(unsigned int j=0; j<Q_m; j++)
            {
              Fq_Mq_representor_innerprods[i][j].resize(Nmax, 0.);
            }
        }

      unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
      Mq_Mq_representor_innerprods.resize(Q_m_hat);
      for(unsigned int i=0; i<Q_m_hat; i++)
        {
          Mq_Mq_representor_innerprods[i].resize(Nmax);
          for(unsigned int j=0; j<Nmax; j++)
            {
              Mq_Mq_representor_innerprods[i][j].resize(Nmax, 0.);
            }
        }

      Aq_Mq_representor_innerprods.resize(Q_a);
      for(unsigned int i=0; i<Q_a; i++)
        {
          Aq_Mq_representor_innerprods[i].resize(Q_m);
          for(unsigned int j=0; j<Q_m; j++)
            {
              Aq_Mq_representor_innerprods[i][j].resize(Nmax);
              for(unsigned int k=0; k<Nmax; k++)
                {
                  Aq_Mq_representor_innerprods[i][j][k].resize(Nmax, 0.);
                }
            }
        }

      // Resize M_q_representor
      // This is cleared in the call to clear_riesz_representors
      // in Parent::resize_RB_data, so just resize here
      M_q_representor.resize(Q_m);
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          M_q_representor[q_m].resize(Nmax);
        }
    }
}

Real TransientRBEvaluation::rb_solve(unsigned int N)
{
  LOG_SCOPE("rb_solve()", "TransientRBEvaluation");

  if(N > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  const unsigned int n_time_steps = get_n_time_steps();
  const Real dt                   = get_delta_t();
  const Real euler_theta          = get_euler_theta();

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
      RB_mass_matrix_N.add(trans_theta_expansion.eval_M_theta(q_m, mu), RB_M_q_m);
    }

  RB_LHS_matrix.resize(N,N);
  RB_LHS_matrix.zero();

  RB_RHS_matrix.resize(N,N);
  RB_RHS_matrix.zero();

  RB_LHS_matrix.add(1./dt, RB_mass_matrix_N);
  RB_RHS_matrix.add(1./dt, RB_mass_matrix_N);

  DenseMatrix<Number> RB_Aq_a;
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);

      RB_LHS_matrix.add(       euler_theta*trans_theta_expansion.eval_A_theta(q_a,mu), RB_Aq_a);
      RB_RHS_matrix.add( -(1.-euler_theta)*trans_theta_expansion.eval_A_theta(q_a,mu), RB_Aq_a);
    }

  // Add forcing terms
  DenseVector<Number> RB_Fq_f;
  RB_RHS_save.resize(N);
  RB_RHS_save.zero();
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      RB_Fq_vector[q_f].get_principal_subvector(N, RB_Fq_f);
      RB_RHS_save.add(trans_theta_expansion.eval_F_theta(q_f,mu), RB_Fq_f);
    }

  // Set system time level to 0
  set_time_step(0);

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
  {
    DenseVector<Number> RB_output_vector_N;
    for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
      {
        RB_outputs_all_k[n][0] = 0.;
        for(unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
          {
            RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
            RB_outputs_all_k[n][0] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
          }
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
      error_bound_all_k[get_time_step()] = std::sqrt(error_bound_sum);

      // Compute the outputs and associated error bounds at the initial time
      DenseVector<Number> RB_output_vector_N;
      for(unsigned int n=0; n<trans_theta_expansion.get_n_outputs(); n++)
        {
          RB_outputs_all_k[n][0] = 0.;
          for(unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
            {
              RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
              RB_outputs_all_k[n][0] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
            }

          RB_output_error_bounds_all_k[n][0] = error_bound_all_k[0] * eval_output_dual_norm(n,mu);
        }

      alpha_LB = get_stability_lower_bound();

      // Precompute time-invariant parts of the dual norm of the residual.
      cache_online_residual_terms(N);
    }

  for(unsigned int time_level=1; time_level<=n_time_steps; time_level++)
    {
      set_time_step(time_level);
      old_RB_solution = RB_solution;

      // Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
      RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

      // Add forcing terms
      RB_rhs.add(get_control(time_level), RB_RHS_save);

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
          for(unsigned int q_l=0; q_l<trans_theta_expansion.get_n_output_terms(n); q_l++)
            {
              RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
              RB_outputs_all_k[n][time_level] += trans_theta_expansion.eval_output_theta(n,q_l,mu)*
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

  _rb_solve_data_cached = true ;

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

Real TransientRBEvaluation::rb_solve_again()
{
  libmesh_assert(_rb_solve_data_cached);

  const unsigned int n_time_steps = get_n_time_steps();
  // Set system time level to 0
  set_time_step(0);

  // Resize/clear the solution vector
  const unsigned int N = RB_RHS_save.size();
  RB_solution.resize(N);

  // Load the initial condition into RB_solution
  if(N > 0)
    RB_solution = RB_initial_condition_all_N[N-1];

  // Resize/clear the old solution vector
  old_RB_solution.resize(N);

  // Initialize the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  for (unsigned int time_level=1; time_level<=n_time_steps; time_level++)
    {
      set_time_step(time_level);
      old_RB_solution = RB_solution;

      // Compute RB_rhs, as *RB_lhs_matrix x old_RB_solution
      RB_RHS_matrix.vector_mult(RB_rhs, old_RB_solution);

      // Add forcing terms
      RB_rhs.add(get_control(time_level), RB_RHS_save);

      if (N > 0)
        RB_LHS_matrix.lu_solve(RB_rhs, RB_solution);
    }

  {
    // Just return -1. We did not compute the error bound
    return -1.;
  }
}

Real TransientRBEvaluation::get_error_bound_normalization()
{
  // Just set the normalization factor to 1 in this case.
  // Users can override this method if specific behavior
  // is required.

  return 1.;
}

Real TransientRBEvaluation::residual_scaling_numer(Real)
{
  return get_delta_t();
}

void TransientRBEvaluation::cache_online_residual_terms(const unsigned int N)
{
  LOG_SCOPE("cache_online_residual_terms()", "TransientRBEvaluation");

  const RBParameters mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  cached_Fq_term = 0.;
  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<Q_f; q_f1++)
    {
      Number cached_theta_q_f1 = trans_theta_expansion.eval_F_theta(q_f1,mu);
      for(unsigned int q_f2=q_f1; q_f2<Q_f; q_f2++)
        {
          Real delta = (q_f1==q_f2) ? 1. : 2.;
          cached_Fq_term += delta*cached_theta_q_f1*trans_theta_expansion.eval_F_theta(q_f2,mu) *
            Fq_representor_innerprods[q];

          q++;
        }
    }

  cached_Fq_Aq_vector.resize(N);
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      Number cached_theta_q_f = trans_theta_expansion.eval_F_theta(q_f,mu);
      for(unsigned int q_a=0; q_a<Q_a; q_a++)
        {
          Number cached_theta_q_a = trans_theta_expansion.eval_A_theta(q_a,mu);
          for(unsigned int i=0; i<N; i++)
            {
              cached_Fq_Aq_vector(i) += 2.*cached_theta_q_f*cached_theta_q_a*
                Fq_Aq_representor_innerprods[q_f][q_a][i];
            }
        }
    }

  cached_Aq_Aq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_a1=0; q_a1<Q_a; q_a1++)
    {
      Number cached_theta_q_a1 = trans_theta_expansion.eval_A_theta(q_a1,mu);
      for(unsigned int q_a2=q_a1; q_a2<Q_a; q_a2++)
        {
          Number cached_theta_q_a2 = trans_theta_expansion.eval_A_theta(q_a2,mu);
          Real delta = (q_a1==q_a2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  cached_Aq_Aq_matrix(i,j) += delta*
                    cached_theta_q_a1*cached_theta_q_a2*
                    Aq_Aq_representor_innerprods[q][i][j];
                }
            }
          q++;
        }
    }

  cached_Fq_Mq_vector.resize(N);
  for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      Number cached_theta_q_f = trans_theta_expansion.eval_F_theta(q_f,mu);
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          Number cached_theta_q_m = trans_theta_expansion.eval_M_theta(q_m,mu);
          for(unsigned int i=0; i<N; i++)
            {
              cached_Fq_Mq_vector(i) += 2.*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_innerprods[q_f][q_m][i];
            }
        }
    }

  cached_Aq_Mq_matrix.resize(N,N);
  for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      Number cached_theta_q_a = trans_theta_expansion.eval_A_theta(q_a,mu);

      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          Number cached_theta_q_m = trans_theta_expansion.eval_M_theta(q_m,mu);

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  cached_Aq_Mq_matrix(i,j) += 2.*cached_theta_q_a*cached_theta_q_m*Aq_Mq_representor_innerprods[q_a][q_m][i][j];
                }
            }
        }
    }

  cached_Mq_Mq_matrix.resize(N,N);
  q=0;
  for(unsigned int q_m1=0; q_m1<Q_m; q_m1++)
    {
      Number cached_theta_q_m1 = trans_theta_expansion.eval_M_theta(q_m1,mu);
      for(unsigned int q_m2=q_m1; q_m2<Q_m; q_m2++)
        {
          Number cached_theta_q_m2 = trans_theta_expansion.eval_M_theta(q_m2,mu);
          Real delta = (q_m1==q_m2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  cached_Mq_Mq_matrix(i,j) += delta*
                    cached_theta_q_m1*cached_theta_q_m2*
                    Mq_Mq_representor_innerprods[q][i][j];
                }
            }
          q++;
        }
    }
}

Real TransientRBEvaluation::compute_residual_dual_norm(const unsigned int N)
{
  LOG_SCOPE("compute_residual_dual_norm()", "TransientRBEvaluation");

  // This assembly assumes we have already called cache_online_residual_terms
  // and that the rb_solve parameter is constant in time

  const Real dt          = get_delta_t();
  const Real euler_theta = get_euler_theta();
  const Real current_control = get_control(get_time_step());

  DenseVector<Number> RB_u_euler_theta(N);
  DenseVector<Number> mass_coeffs(N);
  for(unsigned int i=0; i<N; i++)
    {
      RB_u_euler_theta(i)  = euler_theta*RB_solution(i) +
        (1.-euler_theta)*old_RB_solution(i);
      mass_coeffs(i) = -(RB_solution(i) - old_RB_solution(i))/dt;
    }

  Number residual_norm_sq = current_control*current_control*cached_Fq_term;

  residual_norm_sq += current_control*RB_u_euler_theta.dot(cached_Fq_Aq_vector);
  residual_norm_sq += current_control*mass_coeffs.dot(cached_Fq_Mq_vector);

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
      residual_norm_sq = std::abs(residual_norm_sq);
    }

  return libmesh_real(std::sqrt( residual_norm_sq ));
}

Real TransientRBEvaluation::uncached_compute_residual_dual_norm(const unsigned int N)
{
  LOG_SCOPE("uncached_compute_residual_dual_norm()", "TransientRBEvaluation");

  // Use the stored representor inner product values
  // to evaluate the residual norm

  const RBParameters & mu = get_parameters();

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  const Real dt          = get_delta_t();
  const Real euler_theta = get_euler_theta();

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
      Number cached_theta_q_f1 = trans_theta_expansion.eval_F_theta(q_f1,mu);
      for(unsigned int q_f2=q_f1; q_f2<Q_f; q_f2++)
        {
          Real delta = (q_f1==q_f2) ? 1. : 2.;
          residual_norm_sq += delta*cached_theta_q_f1*trans_theta_expansion.eval_F_theta(q_f2,mu) * Fq_representor_innerprods[q];

          q++;
        }
    }

  for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      Number cached_theta_q_f = trans_theta_expansion.eval_F_theta(q_f,mu);
      for(unsigned int q_a=0; q_a<Q_a; q_a++)
        {
          Number cached_theta_q_a = trans_theta_expansion.eval_A_theta(q_a,mu);
          for(unsigned int i=0; i<N; i++)
            {
              residual_norm_sq += 2.*RB_u_euler_theta[i]*cached_theta_q_f*cached_theta_q_a*
                Fq_Aq_representor_innerprods[q_f][q_a][i];
            }
        }
    }

  q=0;
  for(unsigned int q_a1=0; q_a1<Q_a; q_a1++)
    {
      Number cached_theta_q_a1 = trans_theta_expansion.eval_A_theta(q_a1,mu);
      for(unsigned int q_a2=q_a1; q_a2<Q_a; q_a2++)
        {
          Number cached_theta_q_a2 = trans_theta_expansion.eval_A_theta(q_a2,mu);
          Real delta = (q_a1==q_a2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  residual_norm_sq += delta*RB_u_euler_theta[i]*RB_u_euler_theta[j]*
                    cached_theta_q_a1*cached_theta_q_a2*
                    Aq_Aq_representor_innerprods[q][i][j];
                }
            }
          q++;
        }
    }

  // Now add the terms due to the time-derivative
  q=0;
  for(unsigned int q_m1=0; q_m1<Q_m; q_m1++)
    {
      Number cached_theta_q_m1 = trans_theta_expansion.eval_M_theta(q_m1,mu);
      for(unsigned int q_m2=q_m1; q_m2<Q_m; q_m2++)
        {
          Number cached_theta_q_m2 = trans_theta_expansion.eval_M_theta(q_m2,mu);
          Real delta = (q_m1==q_m2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  residual_norm_sq += delta*mass_coeffs[i]*mass_coeffs[j]*
                    cached_theta_q_m1*cached_theta_q_m2*
                    Mq_Mq_representor_innerprods[q][i][j];
                }
            }
          q++;
        }
    }

  for(unsigned int q_f=0; q_f<Q_f; q_f++)
    {
      Number cached_theta_q_f = trans_theta_expansion.eval_F_theta(q_f,mu);
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          Number cached_theta_q_m = trans_theta_expansion.eval_M_theta(q_m,mu);
          for(unsigned int i=0; i<N; i++)
            {
              residual_norm_sq += 2.*mass_coeffs[i]*cached_theta_q_f * cached_theta_q_m * Fq_Mq_representor_innerprods[q_f][q_m][i];
            }
        }
    }

  for(unsigned int q_a=0; q_a<Q_a; q_a++)
    {
      Number cached_theta_q_a = trans_theta_expansion.eval_A_theta(q_a,mu);

      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          Number cached_theta_q_m = trans_theta_expansion.eval_M_theta(q_m,mu);

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  residual_norm_sq += 2.*RB_u_euler_theta[i]*mass_coeffs[j]*
                    cached_theta_q_a*cached_theta_q_m*
                    Aq_Mq_representor_innerprods[q_a][q_m][i][j];
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
      residual_norm_sq = std::abs(residual_norm_sq);
    }

  //   libMesh::out << "slow residual_sq = " << slow_residual_norm_sq
  //                << ", fast residual_sq = " << residual_norm_sq << std::endl;

  return libmesh_real(std::sqrt( residual_norm_sq ));
}

void TransientRBEvaluation::legacy_write_offline_data_to_files(const std::string & directory_name,
                                                               const bool write_binary_data)
{
  LOG_SCOPE("legacy_write_offline_data_to_files()", "TransientRBEvaluation");

  Parent::legacy_write_offline_data_to_files(directory_name);

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  const unsigned int n_bfs = get_n_basis_functions();

  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = write_binary_data ? ENCODE : WRITE;

  // The suffix to use for all the files that are written out
  const std::string suffix = write_binary_data ? ".xdr" : ".dat";

  if(this->processor_id() == 0)
    {
      std::ostringstream file_name;

      // Write out the temporal discretization data
      file_name.str("");
      file_name << directory_name << "/temporal_discretization_data" << suffix;
      Xdr temporal_discretization_data_out(file_name.str(), mode);

      Real real_value; unsigned int int_value;
      real_value = get_delta_t(); temporal_discretization_data_out << real_value;
      real_value = get_euler_theta(); temporal_discretization_data_out << real_value;
      int_value = get_n_time_steps(); temporal_discretization_data_out << int_value;
      int_value = get_time_step(); temporal_discretization_data_out << int_value;
      temporal_discretization_data_out.close();


      // Write out the L2 matrix
      file_name.str("");
      file_name << directory_name << "/RB_L2_matrix" << suffix;
      Xdr RB_L2_matrix_out(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
            {
              RB_L2_matrix_out << RB_L2_matrix(i,j);
            }
        }
      RB_L2_matrix_out.close();

      // Write out the M_q matrices
      for(unsigned int q_m=0; q_m<Q_m; q_m++)
        {
          file_name.str("");
          file_name << directory_name << "/RB_M_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << q_m;
          file_name << suffix;
          Xdr RB_M_q_m_out(file_name.str(), mode);

          for(unsigned int i=0; i<n_bfs; i++)
            {
              for(unsigned int j=0; j<n_bfs; j++)
                {
                  RB_M_q_m_out << RB_M_q_vector[q_m](i,j);
                }
            }
          RB_M_q_m_out.close();
        }

      // Write out the initial condition data
      // and the initial L2 error for all N
      file_name.str("");
      file_name << directory_name << "/initial_conditions" << suffix;
      Xdr initial_conditions_out(file_name.str(), mode);
      file_name.str("");
      file_name << directory_name << "/initial_L2_error" << suffix;
      Xdr initial_L2_error_out(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          initial_L2_error_out << initial_L2_error_all_N[i];
          for(unsigned int j=0; j<=i; j++)
            {
              initial_conditions_out << RB_initial_condition_all_N[i](j);
            }
        }
      initial_conditions_out.close();
      initial_L2_error_out.close();

      // Next write out the Fq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_Mq_terms" << suffix;
      Xdr RB_Fq_Mq_terms_out(file_name.str(), mode);

      for(unsigned int q_f=0; q_f<Q_f; q_f++)
        {
          for(unsigned int q_m=0; q_m<Q_m; q_m++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  RB_Fq_Mq_terms_out << Fq_Mq_representor_innerprods[q_f][q_m][i];
                }
            }
        }
      RB_Fq_Mq_terms_out.close();

      // Next write out the Mq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Mq_Mq_terms" << suffix;
      Xdr RB_Mq_Mq_terms_out(file_name.str(), mode);

      unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
      for(unsigned int q=0; q<Q_m_hat; q++)
        {
          for(unsigned int i=0; i<n_bfs; i++)
            {
              for(unsigned int j=0; j<n_bfs; j++)
                {
                  RB_Mq_Mq_terms_out << Mq_Mq_representor_innerprods[q][i][j];
                }
            }
        }
      RB_Mq_Mq_terms_out.close();

      // Next write out the Aq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Aq_Mq_terms" << suffix;
      Xdr RB_Aq_Mq_terms_out(file_name.str(), mode);

      for(unsigned int q_a=0; q_a<Q_a; q_a++)
        {
          for(unsigned int q_m=0; q_m<Q_m; q_m++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  for(unsigned int j=0; j<n_bfs; j++)
                    {
                      RB_Aq_Mq_terms_out << Aq_Mq_representor_innerprods[q_a][q_m][i][j];
                    }
                }
            }
        }
      RB_Aq_Mq_terms_out.close();
    }
}

void TransientRBEvaluation::legacy_read_offline_data_from_files(const std::string & directory_name,
                                                                bool read_error_bound_data,
                                                                const bool read_binary_data)
{
  LOG_SCOPE("legacy_read_offline_data_from_files()", "TransientRBEvaluation");

  Parent::legacy_read_offline_data_from_files(directory_name);

  TransientRBThetaExpansion & trans_theta_expansion =
    cast_ref<TransientRBThetaExpansion &>(get_rb_theta_expansion());
  const unsigned int Q_m = trans_theta_expansion.get_n_M_terms();
  const unsigned int Q_a = trans_theta_expansion.get_n_A_terms();
  const unsigned int Q_f = trans_theta_expansion.get_n_F_terms();

  // First, find out how many basis functions we had when Greedy terminated
  // This was set in RBSystem::read_offline_data_from_files
  unsigned int n_bfs = this->get_n_basis_functions();

  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  // The string stream we'll use to make the file names
  std::ostringstream file_name;

  // Write out the temporal discretization data
  file_name.str("");
  file_name << directory_name << "/temporal_discretization_data" << suffix;
  assert_file_exists(file_name.str());

  Xdr temporal_discretization_data_in(file_name.str(), mode);

  Real real_value; unsigned int int_value;
  temporal_discretization_data_in >> real_value; set_delta_t(real_value);
  temporal_discretization_data_in >> real_value; set_euler_theta(real_value);
  temporal_discretization_data_in >> int_value; set_n_time_steps(int_value);
  temporal_discretization_data_in >> int_value; set_time_step(int_value);
  temporal_discretization_data_in.close();

  file_name.str("");
  file_name << directory_name << "/RB_L2_matrix" << suffix;
  assert_file_exists(file_name.str());

  Xdr RB_L2_matrix_in(file_name.str(), mode);

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
      file_name.str("");
      file_name << directory_name << "/RB_M_";
      file_name << std::setw(3)
                << std::setprecision(0)
                << std::setfill('0')
                << std::right
                << q_m;

      file_name << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_M_q_m_in(file_name.str(), mode);

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
  file_name.str("");
  file_name << directory_name << "/initial_conditions" << suffix;
  assert_file_exists(file_name.str());

  Xdr initial_conditions_in(file_name.str(), mode);

  file_name.str("");
  file_name << directory_name << "/initial_L2_error" << suffix;
  assert_file_exists(file_name.str());

  Xdr initial_L2_error_in(file_name.str(), mode);

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


  if(read_error_bound_data)
    {
      // Next read in the Fq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_Mq_terms" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Fq_Mq_terms_in(file_name.str(), mode);

      for(unsigned int q_f=0; q_f<Q_f; q_f++)
        {
          for(unsigned int q_m=0; q_m<Q_m; q_m++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  RB_Fq_Mq_terms_in >> Fq_Mq_representor_innerprods[q_f][q_m][i];
                }
            }
        }
      RB_Fq_Mq_terms_in.close();

      // Next read in the Mq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Mq_Mq_terms" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Mq_Mq_terms_in(file_name.str(), mode);

      unsigned int Q_m_hat = Q_m*(Q_m+1)/2;
      for(unsigned int q=0; q<Q_m_hat; q++)
        {
          for(unsigned int i=0; i<n_bfs; i++)
            {
              for(unsigned int j=0; j<n_bfs; j++)
                {
                  RB_Mq_Mq_terms_in >> Mq_Mq_representor_innerprods[q][i][j];
                }
            }
        }
      RB_Mq_Mq_terms_in.close();

      // Next read in the Aq_Mq representor norm data
      file_name.str("");
      file_name << directory_name << "/Aq_Mq_terms" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Aq_Mq_terms_in(file_name.str(), mode);

      for(unsigned int q_a=0; q_a<Q_a; q_a++)
        {
          for(unsigned int q_m=0; q_m<Q_m; q_m++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  for(unsigned int j=0; j<n_bfs; j++)
                    {
                      RB_Aq_Mq_terms_in >> Aq_Mq_representor_innerprods[q_a][q_m][i][j];
                    }
                }
            }
        }
      RB_Aq_Mq_terms_in.close();
    }
}

}
