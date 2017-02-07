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
#include "libmesh/rb_evaluation.h"
#include "libmesh/rb_theta_expansion.h"

// libMesh includes
#include "libmesh/libmesh_version.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/mesh_tools.h"

// C/C++ includes
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace libMesh
{

RBEvaluation::RBEvaluation (const Parallel::Communicator & comm_in)
  :
  ParallelObject(comm_in),
  evaluate_RB_error_bound(true),
  compute_RB_inner_product(false),
  rb_theta_expansion(libmesh_nullptr)
{

}

RBEvaluation::~RBEvaluation()
{
  this->clear();
}

void RBEvaluation::clear()
{
  LOG_SCOPE("clear()", "RBEvaluation");

  // Clear the basis functions
  for (std::size_t i=0; i<basis_functions.size(); i++)
    {
      if (basis_functions[i])
        {
          basis_functions[i]->clear();
          delete basis_functions[i];
          basis_functions[i] = libmesh_nullptr;
        }
    }
  set_n_basis_functions(0);

  clear_riesz_representors();

  // Clear the Greedy param list
  for (std::size_t i=0; i<greedy_param_list.size(); i++)
    greedy_param_list[i].clear();
  greedy_param_list.clear();
}

void RBEvaluation::set_rb_theta_expansion(RBThetaExpansion & rb_theta_expansion_in)
{
  rb_theta_expansion = &rb_theta_expansion_in;
}

RBThetaExpansion & RBEvaluation::get_rb_theta_expansion()
{
  if(!is_rb_theta_expansion_initialized())
    libmesh_error_msg("Error: rb_theta_expansion hasn't been initialized yet");

  return *rb_theta_expansion;
}

bool RBEvaluation::is_rb_theta_expansion_initialized() const
{
  if(rb_theta_expansion)
    {
      return true;
    }
  else
    {
      return false;
    }
}

void RBEvaluation::resize_data_structures(const unsigned int Nmax,
                                          bool resize_error_bound_data)
{
  LOG_SCOPE("resize_data_structures()", "RBEvaluation");

  if(Nmax < this->get_n_basis_functions())
    libmesh_error_msg("Error: Cannot set Nmax to be less than the current number of basis functions.");

  // Resize/clear inner product matrix
  if(compute_RB_inner_product)
    RB_inner_product_matrix.resize(Nmax,Nmax);

  // Allocate dense matrices for RB solves
  RB_Aq_vector.resize(rb_theta_expansion->get_n_A_terms());

  for(unsigned int q=0; q<rb_theta_expansion->get_n_A_terms(); q++)
    {
      // Initialize the memory for the RB matrices
      RB_Aq_vector[q].resize(Nmax,Nmax);
    }

  RB_Fq_vector.resize(rb_theta_expansion->get_n_F_terms());

  for(unsigned int q=0; q<rb_theta_expansion->get_n_F_terms(); q++)
    {
      // Initialize the memory for the RB vectors
      RB_Fq_vector[q].resize(Nmax);
    }


  // Initialize the RB output vectors
  RB_output_vectors.resize(rb_theta_expansion->get_n_outputs());
  for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
    {
      RB_output_vectors[n].resize(rb_theta_expansion->get_n_output_terms(n));
      for(unsigned int q_l=0; q_l<rb_theta_expansion->get_n_output_terms(n); q_l++)
        {
          RB_output_vectors[n][q_l].resize(Nmax);
        }
    }

  // Initialize vectors storing output data
  RB_outputs.resize(rb_theta_expansion->get_n_outputs(), 0.);


  if(resize_error_bound_data)
    {
      // Initialize vectors for the norms of the Fq representors
      unsigned int Q_f_hat = rb_theta_expansion->get_n_F_terms()*(rb_theta_expansion->get_n_F_terms()+1)/2;
      Fq_representor_innerprods.resize(Q_f_hat);

      // Initialize vectors for the norms of the representors
      Fq_Aq_representor_innerprods.resize(rb_theta_expansion->get_n_F_terms());
      for(unsigned int i=0; i<rb_theta_expansion->get_n_F_terms(); i++)
        {
          Fq_Aq_representor_innerprods[i].resize(rb_theta_expansion->get_n_A_terms());
          for(unsigned int j=0; j<rb_theta_expansion->get_n_A_terms(); j++)
            {
              Fq_Aq_representor_innerprods[i][j].resize(Nmax, 0.);
            }
        }

      unsigned int Q_a_hat = rb_theta_expansion->get_n_A_terms()*(rb_theta_expansion->get_n_A_terms()+1)/2;
      Aq_Aq_representor_innerprods.resize(Q_a_hat);
      for(unsigned int i=0; i<Q_a_hat; i++)
        {
          Aq_Aq_representor_innerprods[i].resize(Nmax);
          for(unsigned int j=0; j<Nmax; j++)
            {
              Aq_Aq_representor_innerprods[i][j].resize(Nmax, 0.);
            }
        }

      RB_output_error_bounds.resize(rb_theta_expansion->get_n_outputs(), 0.);

      // Resize the output dual norm vectors
      output_dual_innerprods.resize(rb_theta_expansion->get_n_outputs());
      for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
        {
          unsigned int Q_l_hat = rb_theta_expansion->get_n_output_terms(n)*(rb_theta_expansion->get_n_output_terms(n)+1)/2;
          output_dual_innerprods[n].resize(Q_l_hat);
        }

      // Clear and resize the vector of Aq_representors
      clear_riesz_representors();

      Aq_representor.resize(rb_theta_expansion->get_n_A_terms());
      for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
        {
          Aq_representor[q_a].resize(Nmax);
        }
    }
}

NumericVector<Number> & RBEvaluation::get_basis_function(unsigned int i)
{
  libmesh_assert_less (i, basis_functions.size());

  return *(basis_functions[i]);
}

Real RBEvaluation::rb_solve(unsigned int N)
{
  LOG_SCOPE("rb_solve()", "RBEvaluation");

  if(N > get_n_basis_functions())
    libmesh_error_msg("ERROR: N cannot be larger than the number of basis functions in rb_solve");

  const RBParameters & mu = get_parameters();

  // Resize (and clear) the solution vector
  RB_solution.resize(N);

  // Assemble the RB system
  DenseMatrix<Number> RB_system_matrix(N,N);
  RB_system_matrix.zero();

  DenseMatrix<Number> RB_Aq_a;
  for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
    {
      RB_Aq_vector[q_a].get_principal_submatrix(N, RB_Aq_a);

      RB_system_matrix.add(rb_theta_expansion->eval_A_theta(q_a, mu), RB_Aq_a);
    }

  // Assemble the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  DenseVector<Number> RB_Fq_f;
  for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
    {
      RB_Fq_vector[q_f].get_principal_subvector(N, RB_Fq_f);

      RB_rhs.add(rb_theta_expansion->eval_F_theta(q_f, mu), RB_Fq_f);
    }

  // Solve the linear system
  if(N > 0)
    {
      RB_system_matrix.lu_solve(RB_rhs, RB_solution);
    }

  // Evaluate RB outputs
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
    {
      RB_outputs[n] = 0.;
      for(unsigned int q_l=0; q_l<rb_theta_expansion->get_n_output_terms(n); q_l++)
        {
          RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
          RB_outputs[n] += rb_theta_expansion->eval_output_theta(n,q_l,mu)*RB_output_vector_N.dot(RB_solution);
        }
    }

  if(evaluate_RB_error_bound) // Calculate the error bounds
    {
      // Evaluate the dual norm of the residual for RB_solution_vector
      Real epsilon_N = compute_residual_dual_norm(N);

      // Get lower bound for coercivity constant
      const Real alpha_LB = get_stability_lower_bound();
      // alpha_LB needs to be positive to get a valid error bound
      libmesh_assert_greater ( alpha_LB, 0. );

      // Evaluate the (absolute) error bound
      Real abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

      // Now compute the output error bounds
      for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
        {
          RB_output_error_bounds[n] = abs_error_bound * eval_output_dual_norm(n, mu);
        }

      return abs_error_bound;
    }
  else // Don't calculate the error bounds
    {
      // Just return -1. if we did not compute the error bound
      return -1.;
    }
}

Real RBEvaluation::get_error_bound_normalization()
{
  // Normalize the error based on the error bound in the
  // case of an empty reduced basis. The error bound is based
  // on the residual F - AU, so with an empty basis this gives
  // a value based on the norm of F at the current parameters.

  Real normalization = rb_solve(0);
  return normalization;
}

Real RBEvaluation::compute_residual_dual_norm(const unsigned int N)
{
  LOG_SCOPE("compute_residual_dual_norm()", "RBEvaluation");

  const RBParameters & mu = get_parameters();

  // Use the stored representor inner product values
  // to evaluate the residual norm
  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<rb_theta_expansion->get_n_F_terms(); q_f1++)
    {
      for(unsigned int q_f2=q_f1; q_f2<rb_theta_expansion->get_n_F_terms(); q_f2++)
        {
          Real delta = (q_f1==q_f2) ? 1. : 2.;
          residual_norm_sq += delta * libmesh_real(
                                                   rb_theta_expansion->eval_F_theta(q_f1, mu)
                                                   * libmesh_conj(rb_theta_expansion->eval_F_theta(q_f2, mu)) * Fq_representor_innerprods[q] );

          q++;
        }
    }

  for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
    {
      for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
        {
          for(unsigned int i=0; i<N; i++)
            {
              Real delta = 2.;
              residual_norm_sq +=
                delta * libmesh_real( rb_theta_expansion->eval_F_theta(q_f, mu) *
                                      libmesh_conj(rb_theta_expansion->eval_A_theta(q_a, mu)) *
                                      libmesh_conj(RB_solution(i)) * Fq_Aq_representor_innerprods[q_f][q_a][i] );
            }
        }
    }

  q=0;
  for(unsigned int q_a1=0; q_a1<rb_theta_expansion->get_n_A_terms(); q_a1++)
    {
      for(unsigned int q_a2=q_a1; q_a2<rb_theta_expansion->get_n_A_terms(); q_a2++)
        {
          Real delta = (q_a1==q_a2) ? 1. : 2.;

          for(unsigned int i=0; i<N; i++)
            {
              for(unsigned int j=0; j<N; j++)
                {
                  residual_norm_sq +=
                    delta * libmesh_real( libmesh_conj(rb_theta_expansion->eval_A_theta(q_a1, mu)) *
                                          rb_theta_expansion->eval_A_theta(q_a2, mu) *
                                          libmesh_conj(RB_solution(i)) * RB_solution(j) * Aq_Aq_representor_innerprods[q][i][j] );
                }
            }

          q++;
        }
    }

  if(libmesh_real(residual_norm_sq) < 0.)
    {
      //    libMesh::out << "Warning: Square of residual norm is negative "
      //                 << "in RBSystem::compute_residual_dual_norm()" << std::endl;

      //     Sometimes this is negative due to rounding error,
      //     but when this occurs the error is on the order of 1.e-10,
      //     so shouldn't affect error bound much...
      residual_norm_sq = std::abs(residual_norm_sq);
    }

  return std::sqrt( libmesh_real(residual_norm_sq) );
}

Real RBEvaluation::get_stability_lower_bound()
{
  // Return a default value of 1, this function should
  // be overloaded to specify a problem-dependent stability
  // factor lower bound
  return 1.;
}

Real RBEvaluation::residual_scaling_denom(Real alpha_LB)
{
  // Here we implement the residual scaling for a coercive
  // problem.
  return alpha_LB;
}

Real RBEvaluation::eval_output_dual_norm(unsigned int n, const RBParameters & mu)
{
  Number output_bound_sq = 0.;
  unsigned int q=0;
  for(unsigned int q_l1=0; q_l1<rb_theta_expansion->get_n_output_terms(n); q_l1++)
    {
      for(unsigned int q_l2=q_l1; q_l2<rb_theta_expansion->get_n_output_terms(n); q_l2++)
        {
          Real delta = (q_l1==q_l2) ? 1. : 2.;
          output_bound_sq += delta * libmesh_real(
                                                  libmesh_conj(rb_theta_expansion->eval_output_theta(n,q_l1,mu))*
                                                  rb_theta_expansion->eval_output_theta(n,q_l2,mu) * output_dual_innerprods[n][q] );
          q++;
        }
    }

  return libmesh_real(std::sqrt( output_bound_sq ));
}

void RBEvaluation::clear_riesz_representors()
{

  // Clear the Aq_representors
  for (std::size_t q_a=0; q_a<Aq_representor.size(); q_a++)
    for (std::size_t i=0; i<Aq_representor[q_a].size(); i++)
      {
        delete Aq_representor[q_a][i];
        Aq_representor[q_a][i] = libmesh_nullptr;
      }
}

void RBEvaluation::legacy_write_offline_data_to_files(const std::string & directory_name,
                                                      const bool write_binary_data)
{
  LOG_SCOPE("legacy_write_offline_data_to_files()", "RBEvaluation");

  // Get the number of basis functions
  unsigned int n_bfs = get_n_basis_functions();

  // The writing mode: ENCODE for binary, WRITE for ASCII
  XdrMODE mode = write_binary_data ? ENCODE : WRITE;

  // The suffix to use for all the files that are written out
  const std::string suffix = write_binary_data ? ".xdr" : ".dat";

  if(this->processor_id() == 0)
    {

      // Make a directory to store all the data files
      mkdir(directory_name.c_str(), 0777);
      //    if( mkdir(directory_name.c_str(), 0777) == -1)
      //    {
      //      libMesh::out << "In RBEvaluation::write_offline_data_to_files, directory "
      //                   << directory_name << " already exists, overwriting contents." << std::endl;
      //    }

      // First, write out how many basis functions we have generated
      std::ostringstream file_name;
      {
        file_name << directory_name << "/n_bfs" << suffix;
        Xdr n_bfs_out(file_name.str(), mode);
        n_bfs_out << n_bfs;
        n_bfs_out.close();
      }

      // Write out the parameter ranges
      file_name.str("");
      file_name << directory_name << "/parameter_ranges" << suffix;
      std::string continuous_param_file_name = file_name.str();

      // Write out the discrete parameter values
      file_name.str("");
      file_name << directory_name << "/discrete_parameter_values" << suffix;
      std::string discrete_param_file_name = file_name.str();

      write_parameter_data_to_files(continuous_param_file_name,
                                    discrete_param_file_name,
                                    write_binary_data);

      // Write out Fq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_innerprods" << suffix;
      Xdr RB_Fq_innerprods_out(file_name.str(), mode);
      unsigned int Q_f_hat = rb_theta_expansion->get_n_F_terms()*(rb_theta_expansion->get_n_F_terms()+1)/2;
      for(unsigned int i=0; i<Q_f_hat; i++)
        {
          RB_Fq_innerprods_out << Fq_representor_innerprods[i];
        }
      RB_Fq_innerprods_out.close();

      // Write out output data
      for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
        {
          file_name.str("");
          file_name << directory_name << "/output_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << n;

          file_name << "_dual_innerprods" << suffix;
          Xdr output_dual_innerprods_out(file_name.str(), mode);

          unsigned int Q_l_hat = rb_theta_expansion->get_n_output_terms(n)*(rb_theta_expansion->get_n_output_terms(n)+1)/2;
          for(unsigned int q=0; q<Q_l_hat; q++)
            {
              output_dual_innerprods_out << output_dual_innerprods[n][q];
            }
          output_dual_innerprods_out.close();
        }


      // Write out output data to multiple files
      for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
        {
          for(unsigned int q_l=0; q_l<rb_theta_expansion->get_n_output_terms(n); q_l++)
            {
              file_name.str("");
              file_name << directory_name << "/output_";
              file_name << std::setw(3)
                        << std::setprecision(0)
                        << std::setfill('0')
                        << std::right
                        << n;
              file_name << "_";
              file_name << std::setw(3)
                        << std::setprecision(0)
                        << std::setfill('0')
                        << std::right
                        << q_l;
              file_name << suffix;
              Xdr output_n_out(file_name.str(), mode);

              for(unsigned int j=0; j<n_bfs; j++)
                {
                  output_n_out << RB_output_vectors[n][q_l](j);
                }
              output_n_out.close();
            }
        }

      if(compute_RB_inner_product)
        {
          // Next write out the inner product matrix
          file_name.str("");
          file_name << directory_name << "/RB_inner_product_matrix" << suffix;
          Xdr RB_inner_product_matrix_out(file_name.str(), mode);
          for(unsigned int i=0; i<n_bfs; i++)
            {
              for(unsigned int j=0; j<n_bfs; j++)
                {
                  RB_inner_product_matrix_out << RB_inner_product_matrix(i,j);
                }
            }
          RB_inner_product_matrix_out.close();
        }

      // Next write out the Fq vectors
      for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
        {
          file_name.str("");
          file_name << directory_name << "/RB_F_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << q_f;
          file_name << suffix;
          Xdr RB_Fq_f_out(file_name.str(), mode);

          for(unsigned int i=0; i<n_bfs; i++)
            {
              RB_Fq_f_out << RB_Fq_vector[q_f](i);
            }
          RB_Fq_f_out.close();
        }

      // Next write out the Aq matrices
      for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
        {
          file_name.str("");
          file_name << directory_name << "/RB_A_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << q_a;
          file_name << suffix;
          Xdr RB_Aq_a_out(file_name.str(), mode);

          for(unsigned int i=0; i<n_bfs; i++)
            {
              for(unsigned int j=0; j<n_bfs; j++)
                {
                  RB_Aq_a_out << RB_Aq_vector[q_a](i,j);
                }
            }
          RB_Aq_a_out.close();
        }

      // Next write out Fq_Aq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_Aq_innerprods" << suffix;
      Xdr RB_Fq_Aq_innerprods_out(file_name.str(), mode);

      for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
        {
          for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  RB_Fq_Aq_innerprods_out << Fq_Aq_representor_innerprods[q_f][q_a][i];
                }
            }
        }
      RB_Fq_Aq_innerprods_out.close();

      // Next write out Aq_Aq representor norm data
      file_name.str("");
      file_name << directory_name << "/Aq_Aq_innerprods" << suffix;
      Xdr RB_Aq_Aq_innerprods_out(file_name.str(), mode);

      unsigned int Q_a_hat = rb_theta_expansion->get_n_A_terms()*(rb_theta_expansion->get_n_A_terms()+1)/2;
      for(unsigned int i=0; i<Q_a_hat; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
            {
              for(unsigned int l=0; l<n_bfs; l++)
                {
                  RB_Aq_Aq_innerprods_out << Aq_Aq_representor_innerprods[i][j][l];
                }
            }
        }
      RB_Aq_Aq_innerprods_out.close();

      // Also, write out the greedily selected parameters
      {
        file_name.str("");
        file_name << directory_name << "/greedy_params" << suffix;
        Xdr greedy_params_out(file_name.str(), mode);

        for (std::size_t i=0; i<greedy_param_list.size(); i++)
          {
            RBParameters::const_iterator it     = greedy_param_list[i].begin();
            RBParameters::const_iterator it_end = greedy_param_list[i].end();
            for( ; it != it_end; ++it)
              {
                // Need to make a copy of the value so that it's not const
                // Xdr is not templated on const's
                Real param_value = it->second;
                greedy_params_out << param_value;
              }
          }
        greedy_params_out.close();
      }

    }
}

void RBEvaluation::legacy_read_offline_data_from_files(const std::string & directory_name,
                                                       bool read_error_bound_data,
                                                       const bool read_binary_data)
{
  LOG_SCOPE("legacy_read_offline_data_from_files()", "RBEvaluation");

  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  // The string stream we'll use to make the file names
  std::ostringstream file_name;

  // First, find out how many basis functions we had when Greedy terminated
  unsigned int n_bfs;
  {
    file_name << directory_name << "/n_bfs" << suffix;
    assert_file_exists(file_name.str());

    Xdr n_bfs_in(file_name.str(), mode);
    n_bfs_in >> n_bfs;
    n_bfs_in.close();
  }
  resize_data_structures(n_bfs, read_error_bound_data);

  // Read in the parameter ranges
  file_name.str("");
  file_name << directory_name << "/parameter_ranges" << suffix;
  std::string continuous_param_file_name = file_name.str();

  // Read in the discrete parameter values
  file_name.str("");
  file_name << directory_name << "/discrete_parameter_values" << suffix;
  std::string discrete_param_file_name = file_name.str();
  read_parameter_data_from_files(continuous_param_file_name,
                                 discrete_param_file_name,
                                 read_binary_data);

  // Read in output data in multiple files
  for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
    {
      for(unsigned int q_l=0; q_l<rb_theta_expansion->get_n_output_terms(n); q_l++)
        {
          file_name.str("");
          file_name << directory_name << "/output_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << n;
          file_name << "_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << q_l;
          file_name << suffix;
          assert_file_exists(file_name.str());

          Xdr output_n_in(file_name.str(), mode);

          for(unsigned int j=0; j<n_bfs; j++)
            {
              Number value;
              output_n_in >> value;
              RB_output_vectors[n][q_l](j) = value;
            }
          output_n_in.close();
        }
    }

  if(compute_RB_inner_product)
    {
      // Next read in the inner product matrix
      file_name.str("");
      file_name << directory_name << "/RB_inner_product_matrix" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_inner_product_matrix_in(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
            {
              Number value;
              RB_inner_product_matrix_in >> value;
              RB_inner_product_matrix(i,j) = value;
            }
        }
      RB_inner_product_matrix_in.close();
    }

  // Next read in the Fq vectors
  for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
    {
      file_name.str("");
      file_name << directory_name << "/RB_F_";
      file_name << std::setw(3)
                << std::setprecision(0)
                << std::setfill('0')
                << std::right
                << q_f;
      file_name << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Fq_f_in(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          Number value;
          RB_Fq_f_in >> value;
          RB_Fq_vector[q_f](i) = value;
        }
      RB_Fq_f_in.close();
    }

  // Next read in the Aq matrices
  for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
    {
      file_name.str("");
      file_name << directory_name << "/RB_A_";
      file_name << std::setw(3)
                << std::setprecision(0)
                << std::setfill('0')
                << std::right
                << q_a;
      file_name << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Aq_a_in(file_name.str(), mode);

      for(unsigned int i=0; i<n_bfs; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
            {
              Number  value;
              RB_Aq_a_in >> value;
              RB_Aq_vector[q_a](i,j) = value;
            }
        }
      RB_Aq_a_in.close();
    }


  if(read_error_bound_data)
    {
      // Next read in Fq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_innerprods" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Fq_innerprods_in(file_name.str(), mode);

      unsigned int Q_f_hat = rb_theta_expansion->get_n_F_terms()*(rb_theta_expansion->get_n_F_terms()+1)/2;
      for(unsigned int i=0; i<Q_f_hat; i++)
        {
          RB_Fq_innerprods_in >> Fq_representor_innerprods[i];
        }
      RB_Fq_innerprods_in.close();

      // Read in output data
      for(unsigned int n=0; n<rb_theta_expansion->get_n_outputs(); n++)
        {
          file_name.str("");
          file_name << directory_name << "/output_";
          file_name << std::setw(3)
                    << std::setprecision(0)
                    << std::setfill('0')
                    << std::right
                    << n;
          file_name << "_dual_innerprods" << suffix;
          assert_file_exists(file_name.str());

          Xdr output_dual_innerprods_in(file_name.str(), mode);

          unsigned int Q_l_hat = rb_theta_expansion->get_n_output_terms(n)*(rb_theta_expansion->get_n_output_terms(n)+1)/2;
          for(unsigned int q=0; q<Q_l_hat; q++)
            {
              output_dual_innerprods_in >> output_dual_innerprods[n][q];
            }
          output_dual_innerprods_in.close();
        }


      // Next read in Fq_Aq representor norm data
      file_name.str("");
      file_name << directory_name << "/Fq_Aq_innerprods" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Fq_Aq_innerprods_in(file_name.str(), mode);

      for(unsigned int q_f=0; q_f<rb_theta_expansion->get_n_F_terms(); q_f++)
        {
          for(unsigned int q_a=0; q_a<rb_theta_expansion->get_n_A_terms(); q_a++)
            {
              for(unsigned int i=0; i<n_bfs; i++)
                {
                  RB_Fq_Aq_innerprods_in >> Fq_Aq_representor_innerprods[q_f][q_a][i];
                }
            }
        }
      RB_Fq_Aq_innerprods_in.close();

      // Next read in Aq_Aq representor norm data
      file_name.str("");
      file_name << directory_name << "/Aq_Aq_innerprods" << suffix;
      assert_file_exists(file_name.str());

      Xdr RB_Aq_Aq_innerprods_in(file_name.str(), mode);

      unsigned int Q_a_hat = rb_theta_expansion->get_n_A_terms()*(rb_theta_expansion->get_n_A_terms()+1)/2;
      for(unsigned int i=0; i<Q_a_hat; i++)
        {
          for(unsigned int j=0; j<n_bfs; j++)
            {
              for(unsigned int l=0; l<n_bfs; l++)
                {
                  RB_Aq_Aq_innerprods_in >> Aq_Aq_representor_innerprods[i][j][l];
                }
            }
        }
      RB_Aq_Aq_innerprods_in.close();
    }

  // Resize basis_functions even if we don't read them in so that
  // get_n_bfs() returns the correct value. Initialize the pointers
  // to NULL
  set_n_basis_functions(n_bfs);
  for (std::size_t i=0; i<basis_functions.size(); i++)
    {
      if(basis_functions[i])
        {
          basis_functions[i]->clear();
          delete basis_functions[i];
        }
      basis_functions[i] = libmesh_nullptr;
    }
}

void RBEvaluation::assert_file_exists(const std::string & file_name)
{
  if (!std::ifstream(file_name.c_str()))
    libmesh_error_msg("File missing: " + file_name);
}

void RBEvaluation::write_out_basis_functions(System & sys,
                                             const std::string & directory_name,
                                             const bool write_binary_basis_functions)
{
  LOG_SCOPE("write_out_basis_functions()", "RBEvaluation");

  write_out_vectors(sys,
                    basis_functions,
                    directory_name,
                    "bf",
                    write_binary_basis_functions);
}

void RBEvaluation::write_out_vectors(System & sys,
                                     std::vector<NumericVector<Number> *> & vectors,
                                     const std::string & directory_name,
                                     const std::string & data_name,
                                     const bool write_binary_vectors)
{
  LOG_SCOPE("write_out_vectors()", "RBEvaluation");

  if(this->processor_id() == 0)
    {
      // Make a directory to store all the data files
      mkdir(directory_name.c_str(), 0777);
    }

  // Make sure processors are synced up before we begin
  this->comm().barrier();

  std::ostringstream file_name;
  const std::string basis_function_suffix = (write_binary_vectors ? ".xdr" : ".dat");

  file_name << directory_name << "/" << data_name << "_header" << basis_function_suffix;
  Xdr header_data(file_name.str(),
                  write_binary_vectors ? ENCODE : WRITE);
  sys.write_header(header_data, get_io_version_string(), /*write_additional_data=*/false);

  // Following EquationSystemsIO::write, we use a temporary numbering (node major)
  // before writing out the data
  MeshTools::Private::globally_renumber_nodes_and_elements(sys.get_mesh());

  // // Use System::write_serialized_data to write out the basis functions
  // // by copying them into this->solution one at a time.
  // for (std::size_t i=0; i<vectors.size(); i++)
  // {
  //   // No need to copy, just swap
  //   // *solution = *vectors[i];
  //   vectors[i]->swap(*sys.solution);
  //   file_name.str(""); // reset the string
  //   file_name << directory_name << "/bf" << i << basis_function_suffix;
  //   Xdr bf_data(file_name.str(),
  //               write_binary_vectors ? ENCODE : WRITE);
  //   // set the current version
  //   bf_data.set_version(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION,
  //    LIBMESH_MINOR_VERSION,
  //    LIBMESH_MICRO_VERSION));

  //   sys.write_serialized_data(bf_data, false);

  //   // Synchronize before moving on
  //   this->comm().barrier();
  //   // Swap back
  //   vectors[i]->swap(*sys.solution);
  // }

  file_name.str("");
  file_name << directory_name << "/" << data_name << "_data" << basis_function_suffix;

  Xdr bf_data(file_name.str(),
              write_binary_vectors ? ENCODE : WRITE);

  // Write all vectors at once.
  {
    // Note the API wants pointers to constant vectors, hence this...
    std::vector<const NumericVector<Number> *> bf_out(vectors.begin(),
                                                      vectors.end());
    // for (std::size_t i=0; i<vectors.size(); i++)
    //   bf_out.push_back(vectors[i]);
    sys.write_serialized_vectors (bf_data, bf_out);
  }


  // set the current version
  bf_data.set_version(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION,
                                         LIBMESH_MINOR_VERSION,
                                         LIBMESH_MICRO_VERSION));


  // Undo the temporary renumbering
  sys.get_mesh().fix_broken_node_and_element_numbering();
}

void RBEvaluation::read_in_basis_functions(System & sys,
                                           const std::string & directory_name,
                                           const bool read_binary_basis_functions)
{
  LOG_SCOPE("read_in_basis_functions()", "RBEvaluation");

  read_in_vectors(sys,
                  basis_functions,
                  directory_name,
                  "bf",
                  read_binary_basis_functions);
}

void RBEvaluation::read_in_vectors(System & sys,
                                   std::vector<NumericVector<Number> *> & vectors,
                                   const std::string & directory_name,
                                   const std::string & data_name,
                                   const bool read_binary_vectors)
{
  std::vector<std::vector<NumericVector<Number> *> *> vectors_vec;
  vectors_vec.push_back(&vectors);

  std::vector<std::string> directory_name_vec;
  directory_name_vec.push_back(directory_name);

  std::vector<std::string> data_name_vec;
  data_name_vec.push_back(data_name);

  read_in_vectors_from_multiple_files(sys,
                                      vectors_vec,
                                      directory_name_vec,
                                      data_name_vec,
                                      read_binary_vectors);
}

void RBEvaluation::read_in_vectors_from_multiple_files(System & sys,
                                                       std::vector< std::vector<NumericVector<Number> *> * > multiple_vectors,
                                                       const std::vector<std::string> & multiple_directory_names,
                                                       const std::vector<std::string> & multiple_data_names,
                                                       const bool read_binary_vectors)
{
  LOG_SCOPE("read_in_vectors_from_multiple_files()", "RBEvaluation");

  unsigned int n_files = multiple_vectors.size();
  unsigned int n_directories = multiple_directory_names.size();
  libmesh_assert((n_files == n_directories) && (n_files == multiple_data_names.size()));

  if (n_files == 0)
    return;

  // Make sure processors are synced up before we begin
  this->comm().barrier();

  std::ostringstream file_name;
  const std::string basis_function_suffix = (read_binary_vectors ? ".xdr" : ".dat");
  struct stat stat_info;

  // Assume that all the headers are the same, hence we can just use the first one.
  file_name << multiple_directory_names[0] << "/"
            << multiple_data_names[0] << "_header" << basis_function_suffix;
  assert_file_exists(file_name.str());

  Xdr header_data(file_name.str(),
                  read_binary_vectors ? DECODE : READ);

  // set the version number in header_data from io_version_string
  // (same code as in EquationSystemsIO::_read_impl)
  std::string io_version_string = get_io_version_string();
  std::string::size_type lm_pos = io_version_string.find("libMesh");
  std::istringstream iss(io_version_string.substr(lm_pos + 8));
  int ver_major = 0, ver_minor = 0, ver_patch = 0;
  char dot;
  iss >> ver_major >> dot >> ver_minor >> dot >> ver_patch;
  header_data.set_version(LIBMESH_VERSION_ID(ver_major, ver_minor, ver_patch));

  // We need to call sys.read_header (e.g. to set _written_var_indices properly),
  // but by setting the read_header argument to false, it doesn't reinitialize the system
  sys.read_header(header_data, io_version_string, /*read_header=*/false, /*read_additional_data=*/false);

  // Following EquationSystemsIO::read, we use a temporary numbering (node major)
  // before writing out the data
  MeshTools::Private::globally_renumber_nodes_and_elements(sys.get_mesh());

  for (unsigned int data_index=0; data_index<n_directories; data_index++)
    {
      std::vector<NumericVector<Number> *> & vectors = *multiple_vectors[data_index];

      // Allocate storage for each vector
      for (std::size_t i=0; i<vectors.size(); i++)
        {
          // vectors should all be NULL, otherwise we get a memory leak when
          // we create the new vectors in RBEvaluation::read_in_vectors.
          if (vectors[i])
            libmesh_error_msg("Non-NULL vector passed to read_in_vectors_from_multiple_files");

          vectors[i] = NumericVector<Number>::build(sys.comm()).release();

          vectors[i]->init (sys.n_dofs(),
                            sys.n_local_dofs(),
                            false,
                            PARALLEL);
        }

      file_name.str("");
      file_name << multiple_directory_names[data_index]
                << "/" << multiple_data_names[data_index]
                << "_data" << basis_function_suffix;

      // On processor zero check to be sure the file exists
      if (this->processor_id() == 0)
        {
          int stat_result = stat(file_name.str().c_str(), &stat_info);

          if (stat_result != 0)
            libmesh_error_msg("File does not exist: " << file_name.str());
        }

      assert_file_exists(file_name.str());
      Xdr vector_data(file_name.str(),
                      read_binary_vectors ? DECODE : READ);

      // The vector_data needs to know which version to read.
      vector_data.set_version(LIBMESH_VERSION_ID(ver_major, ver_minor, ver_patch));

      sys.read_serialized_vectors (vector_data, vectors);
    }

  // Undo the temporary renumbering
  sys.get_mesh().fix_broken_node_and_element_numbering();
}

std::string RBEvaluation::get_io_version_string()
{
  std::string retval("libMesh-" + libMesh::get_io_compatibility_version());

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  retval += " with infinite elements";
#endif

  return retval;
}

} // namespace libMesh
