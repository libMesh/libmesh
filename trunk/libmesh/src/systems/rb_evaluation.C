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

#include "rb_evaluation.h"
#include "rb_system.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "libmesh_logging.h"
#include "xdr_cxx.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "o_string_stream.h"
#include <fstream>
#include <sstream>

namespace libMesh
{

RBEvaluation::RBEvaluation (RBSystem& rb_sys_in)
  :
  rb_sys(rb_sys_in)
{
}

RBEvaluation::~RBEvaluation()
{
  this->clear();
}

void RBEvaluation::clear()
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

}

void RBEvaluation::initialize()
{
  const unsigned int Nmax = rb_sys.get_Nmax();

  // Resize/clear inner product matrix
  if(rb_sys.compute_RB_inner_product)
    RB_inner_product_matrix.resize(Nmax,Nmax);

  // Allocate dense matrices for RB solves
  RB_A_q_vector.resize(rb_sys.get_Q_a());

  for(unsigned int q=0; q<rb_sys.get_Q_a(); q++)
  {
    // Initialize the memory for the RB matrices
    RB_A_q_vector[q].resize(Nmax,Nmax);
  }

  RB_F_q_vector.resize(rb_sys.get_Q_f());

  for(unsigned int q=0; q<rb_sys.get_Q_f(); q++)
  {
    // Initialize the memory for the RB vectors
    RB_F_q_vector[q].resize(Nmax);
  }

  // Initialize vectors for the norms of the representors
  unsigned int Q_f_hat = rb_sys.get_Q_f()*(rb_sys.get_Q_f()+1)/2;
  Fq_representor_norms.resize(Q_f_hat);

  Fq_Aq_representor_norms.resize(rb_sys.get_Q_f());
  for(unsigned int i=0; i<rb_sys.get_Q_f(); i++)
  {
    Fq_Aq_representor_norms[i].resize(rb_sys.get_Q_a());
    for(unsigned int j=0; j<rb_sys.get_Q_a(); j++)
    {
      Fq_Aq_representor_norms[i][j].resize(Nmax);
    }
  }

  unsigned int Q_a_hat = rb_sys.get_Q_a()*(rb_sys.get_Q_a()+1)/2;
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
  RB_output_vectors.resize(rb_sys.get_n_outputs());
  for(unsigned int n=0; n<rb_sys.get_n_outputs(); n++)
  {
    RB_output_vectors[n].resize(rb_sys.get_Q_l(n));
    for(unsigned int q_l=0; q_l<rb_sys.get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].resize(Nmax);
    }
  }

  // Initialize vectors storing output data
  RB_outputs.resize(rb_sys.get_n_outputs());
  RB_output_error_bounds.resize(rb_sys.get_n_outputs());
}

Real RBEvaluation::RB_solve(unsigned int N)
{
  START_LOG("RB_solve()", "RBEvaluation");

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

  // Resize (and clear) the solution vector
  RB_solution.resize(N);

  // Assemble the RB system
  DenseMatrix<Number> RB_system_matrix(N,N);
  RB_system_matrix.zero();

  DenseMatrix<Number> RB_A_q_a;
  for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
  {
    RB_A_q_vector[q_a].get_principal_submatrix(N, RB_A_q_a);

    RB_system_matrix.add(rb_sys.eval_theta_q_a(q_a), RB_A_q_a);
  }

  // Assemble the RB rhs
  DenseVector<Number> RB_rhs(N);
  RB_rhs.zero();

  DenseVector<Number> RB_F_q_f;
  for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
  {
    RB_F_q_vector[q_f].get_principal_subvector(N, RB_F_q_f);

    RB_rhs.add(rb_sys.eval_theta_q_f(q_f), RB_F_q_f);
  }
  
  // Solve the linear system
  RB_system_matrix.lu_solve(RB_rhs, RB_solution);

  // Evaluate the dual norm of the residual for RB_solution_vector
  Real epsilon_N = compute_residual_dual_norm(N);

  // Get lower bound for coercivity constant
  const Real alpha_LB = rb_sys.get_SCM_lower_bound();
  // alpha_LB needs to be positive to get a valid error bound
  libmesh_assert( alpha_LB > 0. );

  // Store (absolute) error bound
  Real abs_error_bound = epsilon_N / rb_sys.residual_scaling_denom(alpha_LB);

  // Compute the norm of RB_solution
  Real RB_solution_norm = RB_solution.l2_norm();

  // Now compute the outputs and associated errors
  DenseVector<Number> RB_output_vector_N;
  for(unsigned int n=0; n<rb_sys.get_n_outputs(); n++)
  {
    RB_outputs[n] = 0.;
    for(unsigned int q_l=0; q_l<rb_sys.get_Q_l(n); q_l++)
    {
      RB_output_vectors[n][q_l].get_principal_subvector(N, RB_output_vector_N);
      RB_outputs[n] += libmesh_conj(rb_sys.eval_theta_q_l(n,q_l))*RB_output_vector_N.dot(RB_solution);
    }

    RB_output_error_bounds[n] = abs_error_bound * rb_sys.eval_output_dual_norm(n);
  }

  STOP_LOG("RB_solve()", "RBEvaluation");

  return ( rb_sys.return_rel_error_bound ? abs_error_bound/RB_solution_norm : abs_error_bound );
}

Real RBEvaluation::compute_residual_dual_norm(const unsigned int N)
{
  START_LOG("compute_residual_dual_norm()", "RBEvaluation");

  // Use the stored representor inner product values
  // to evaluate the residual norm
  Number residual_norm_sq = 0.;

  unsigned int q=0;
  for(unsigned int q_f1=0; q_f1<rb_sys.get_Q_f(); q_f1++)
  {
    for(unsigned int q_f2=q_f1; q_f2<rb_sys.get_Q_f(); q_f2++)
    {
      Real delta = (q_f1==q_f2) ? 1. : 2.;
      residual_norm_sq += delta * libmesh_real(
        rb_sys.eval_theta_q_f(q_f1) * libmesh_conj(rb_sys.eval_theta_q_f(q_f2)) * Fq_representor_norms[q] );

      q++;
    }
  }

  for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
  {
    for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
    {
      for(unsigned int i=0; i<N; i++)
      {
        Real delta = 2.;
        residual_norm_sq +=
          delta * libmesh_real( rb_sys.eval_theta_q_f(q_f) * libmesh_conj(rb_sys.eval_theta_q_a(q_a)) *
          libmesh_conj(RB_solution(i)) * Fq_Aq_representor_norms[q_f][q_a][i] );
      }
    }
  }

  q=0;
  for(unsigned int q_a1=0; q_a1<rb_sys.get_Q_a(); q_a1++)
  {
    for(unsigned int q_a2=q_a1; q_a2<rb_sys.get_Q_a(); q_a2++)
    {
      Real delta = (q_a1==q_a2) ? 1. : 2.;

      for(unsigned int i=0; i<N; i++)
      {
        for(unsigned int j=0; j<N; j++)
        {
          residual_norm_sq +=
            delta * libmesh_real( libmesh_conj(rb_sys.eval_theta_q_a(q_a1)) * rb_sys.eval_theta_q_a(q_a2) *
            libmesh_conj(RB_solution(i)) * RB_solution(j) * Aq_Aq_representor_norms[q][i][j] );
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
//     libmesh_error();
    residual_norm_sq = std::abs(residual_norm_sq);
  }

  STOP_LOG("compute_residual_dual_norm()", "RBEvaluation");

  return std::sqrt( libmesh_real(residual_norm_sq) );
}

void RBEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBEvaluation");

  const unsigned int precision_level = 14;

  const unsigned int n_bfs = get_n_basis_functions();
  libmesh_assert( n_bfs <= rb_sys.get_Nmax() );

  if(libMesh::processor_id() == 0)
  {

    // Make a directory to store all the data files
    if( mkdir(directory_name.c_str(), 0777) == -1)
    {
      libMesh::out << "In RBEvaluation::write_offline_data_to_files, directory "
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
        libMesh::err << "Error opening n_bfs.dat" << std::endl;
        libmesh_error();
      }
      n_bfs_out << n_bfs;
      n_bfs_out.close();
    }

    // Write out output data
    for(unsigned int n=0; n<rb_sys.get_n_outputs(); n++)
    {
      for(unsigned int q_l=0; q_l<rb_sys.get_Q_l(n); q_l++)
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
          libMesh::err << "Error opening output file for output " << n << std::endl;
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
    
    if(rb_sys.compute_RB_inner_product)
    {
      // Next write out the inner product matrix
      std::ofstream RB_inner_product_matrix_out;
      {
        OStringStream file_name;
        file_name << directory_name << "/RB_inner_product_matrix.dat";
        RB_inner_product_matrix_out.open(file_name.str().c_str());
      }
      if ( !RB_inner_product_matrix_out.good() )
      {
        libMesh::err << "Error opening RB_inner_product_matrix.dat" << std::endl;
        libmesh_error();
      }
      RB_inner_product_matrix_out.precision(precision_level);
      for(unsigned int i=0; i<n_bfs; i++)
      {
        for(unsigned int j=0; j<n_bfs; j++)
        {
          RB_inner_product_matrix_out << std::scientific
            << RB_inner_product_matrix(i,j) << " ";
        }
      }
      RB_inner_product_matrix_out.close();
    }

    // Next write out the F_q vectors
    for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_F_";
      OSSRealzeroright(file_name,3,0,q_f);
      file_name << ".dat";
      std::ofstream RB_F_q_f_out(file_name.str().c_str());

      if ( !RB_F_q_f_out.good() )
      {
        libMesh::err << "Error opening RB_F_" << q_f << ".dat" << std::endl;
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
    for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_A_";
      OSSRealzeroright(file_name,3,0,q_a);
      file_name << ".dat";
      std::ofstream RB_A_q_a_out(file_name.str().c_str());

      if ( !RB_A_q_a_out.good() )
      {
        libMesh::err << "Error opening RB_A_" << q_a << ".dat" << std::endl;
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
      libMesh::err << "Error opening Fq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_norms_out.precision(precision_level);
    unsigned int Q_f_hat = rb_sys.get_Q_f()*(rb_sys.get_Q_f()+1)/2;
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
      libMesh::err << "Error opening Fq_Aq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Fq_Aq_norms_out.precision(precision_level);
    for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
    {
      for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
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
      libMesh::err << "Error opening Aq_Aq_norms.dat" << std::endl;
      libmesh_error();
    }
    RB_Aq_Aq_norms_out.precision(precision_level);
    unsigned int Q_a_hat = rb_sys.get_Q_a()*(rb_sys.get_Q_a()+1)/2;
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
  write_out_basis_functions(directory_name, precision_level);

  STOP_LOG("write_offline_data_to_files()", "RBEvaluation");
}

void RBEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBEvaluation");

  // First, find out how many basis functions we had when Greedy terminated
  unsigned int n_bfs;
  {
    OStringStream file_name;
    file_name << directory_name << "/n_bfs.dat";
    std::ifstream n_bfs_in(file_name.str().c_str());

    if ( !n_bfs_in.good() )
    {
      libMesh::err << "Error opening n_bfs.dat" << std::endl;
      libmesh_error();
    }

    n_bfs_in >> n_bfs;
    n_bfs_in.close();
  }
  libmesh_assert( n_bfs <= rb_sys.get_Nmax() );

  // Read in output data
  for(unsigned int n=0; n<rb_sys.get_n_outputs(); n++)
  {
    for(unsigned int q_l=0; q_l<rb_sys.get_Q_l(n); q_l++)
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
        libMesh::err << "Error opening input file for output " << n << std::endl;
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
  
  if(rb_sys.compute_RB_inner_product)
  {
    // Next read in the inner product matrix
    std::ifstream RB_inner_product_matrix_in;
    {
      OStringStream file_name;
      file_name << directory_name << "/RB_inner_product_matrix.dat";
      RB_inner_product_matrix_in.open(file_name.str().c_str());
    }
    if ( !RB_inner_product_matrix_in.good() )
    {
      libMesh::err << "Error opening RB_inner_product_matrix.dat" << std::endl;
      libmesh_error();
    }
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

  // Next read in the F_q vectors
  for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_F_";
    OSSRealzeroright(file_name,3,0,q_f);
    file_name << ".dat";
    std::ifstream RB_F_q_f_in(file_name.str().c_str());

    if ( !RB_F_q_f_in.good() )
    {
      libMesh::err << "Error opening RB_F_" << q_f << ".dat" << std::endl;
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
  for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
  {
    OStringStream file_name;
    file_name << directory_name << "/RB_A_";
    OSSRealzeroright(file_name,3,0,q_a);
    file_name << ".dat";
    std::ifstream RB_A_q_a_in(file_name.str().c_str());

    if ( !RB_A_q_a_in.good() )
    {
      libMesh::err << "Error opening RB_A_" << q_a << ".dat" << std::endl;
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
    libMesh::err << "Error opening Fq_norms.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_f_hat = rb_sys.get_Q_f()*(rb_sys.get_Q_f()+1)/2;
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
    libMesh::err << "Error opening Fq_Aq_norms.dat" << std::endl;
    libmesh_error();
  }
  for(unsigned int q_f=0; q_f<rb_sys.get_Q_f(); q_f++)
  {
    for(unsigned int q_a=0; q_a<rb_sys.get_Q_a(); q_a++)
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
    libMesh::err << "Error opening Aq_Aq_norms.dat" << std::endl;
    libmesh_error();
  }
  unsigned int Q_a_hat = rb_sys.get_Q_a()*(rb_sys.get_Q_a()+1)/2;
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
  set_n_basis_functions(n_bfs);
  for(unsigned int i=0; i<basis_functions.size(); i++)
    {
      if(basis_functions[i])
      {
        basis_functions[i]->clear();
        delete basis_functions[i];
      }
      basis_functions[i] = NULL;
    }

  // Read in the basis functions if requested.
  read_in_basis_functions(directory_name);

  STOP_LOG("read_offline_data_from_files()", "RBEvaluation");
}

void RBEvaluation::write_out_basis_functions(const std::string& directory_name,
                                         const unsigned int)
{
  if(rb_sys.store_basis_functions)
  {
    libMesh::out << "Writing out the basis functions..." << std::endl;

    std::ostringstream file_name;
    const std::string basis_function_suffix = (rb_sys.write_binary_basis_functions ? ".xdr" : ".dat");

    // Use System::write_serialized_data to write out the basis functions
    // by copying them into this->solution one at a time.
    for(unsigned int i=0; i<basis_functions.size(); i++)
    {
      // No need to copy, just swap
      // *solution = *basis_functions[i];
      basis_functions[i]->swap(*rb_sys.solution);

      file_name.str(""); // reset the string
      file_name << directory_name << "/bf" << i << basis_function_suffix;

      Xdr bf_data(file_name.str(),
		  rb_sys.write_binary_basis_functions ? ENCODE : WRITE);

      rb_sys.write_serialized_data(bf_data, false);

      // Synchronize before moving on
      Parallel::barrier();

      // Swap back
      basis_functions[i]->swap(*rb_sys.solution);
    }
  }
}

void RBEvaluation::read_in_basis_functions(const std::string& directory_name)
{
  if(rb_sys.store_basis_functions)
  {
    libMesh::out << "Reading in the basis functions..." << std::endl;

    std::ostringstream file_name;
    const std::string basis_function_suffix = (rb_sys.read_binary_basis_functions ? ".xdr" : ".dat");
    struct stat stat_info;

    // Use System::read_serialized_data to read in the basis functions
    // into this->solution and then swap with the appropriate
    // of basis function.
    for(unsigned int i=0; i<basis_functions.size(); i++)
    {
      file_name.str(""); // reset the string
      file_name << directory_name << "/bf" << i << basis_function_suffix;

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

      Xdr bf_data(file_name.str(),
		  rb_sys.read_binary_basis_functions ? DECODE : READ);

      rb_sys.read_serialized_data(bf_data, false);

      basis_functions[i] = NumericVector<Number>::build().release();
      basis_functions[i]->init (rb_sys.n_dofs(), rb_sys.n_local_dofs(), false, libMeshEnums::PARALLEL);

      // No need to copy, just swap
      // *basis_functions[i] = *solution;
      basis_functions[i]->swap(*rb_sys.solution);
    }

    libMesh::out << "Finished reading in the basis functions..." << std::endl;
  }
}

} // namespace libMesh
