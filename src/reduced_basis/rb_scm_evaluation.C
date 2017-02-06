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
#include "libmesh/libmesh_config.h"

// RBSCMEvaluation should only be available
// if SLEPc and GLPK support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

// rbOOmit includes
#include "libmesh/rb_scm_evaluation.h"
#include "libmesh/rb_theta_expansion.h"

// libMesh includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel.h"
#include "libmesh/dof_map.h"
#include "libmesh/xdr_cxx.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// glpk includes
#include <glpk.h>

namespace libMesh
{

RBSCMEvaluation::RBSCMEvaluation (const Parallel::Communicator & comm_in) :
  ParallelObject(comm_in)
{
  // Clear SCM data vectors
  B_min.clear();
  B_max.clear();
  C_J.clear();
  C_J_stability_vector.clear();
  SCM_UB_vectors.clear();
}

RBSCMEvaluation::~RBSCMEvaluation ()
{
}

void RBSCMEvaluation::set_rb_theta_expansion(RBThetaExpansion & rb_theta_expansion_in)
{
  rb_theta_expansion = &rb_theta_expansion_in;
}

RBThetaExpansion & RBSCMEvaluation::get_rb_theta_expansion()
{
  if(!rb_theta_expansion)
    libmesh_error_msg("Error: rb_theta_expansion hasn't been initialized yet");

  return *rb_theta_expansion;
}

void RBSCMEvaluation::set_C_J_stability_constraint(unsigned int j, Real stability_const_in)
{
  if(j >= C_J_stability_vector.size())
    libmesh_error_msg("Error: Input parameter j is too large in set_C_J_stability_constraint.");

  // we assume that C_J_stability_vector is resized elsewhere
  // to be the same size as C_J.
  libmesh_assert_equal_to (C_J_stability_vector.size(), C_J.size());

  C_J_stability_vector[j] = stability_const_in;
}

Real RBSCMEvaluation::get_C_J_stability_constraint(unsigned int j) const
{
  if(j >= C_J_stability_vector.size())
    libmesh_error_msg("Error: Input parameter j is too large in get_C_J_stability_constraint.");

  return C_J_stability_vector[j];
}

void RBSCMEvaluation::set_SCM_UB_vector(unsigned int j, unsigned int q, Real y_q)
{
  // First make sure that j <= J
  if(j >= SCM_UB_vectors.size())
    libmesh_error_msg("Error: We must have j < J in set_SCM_UB_vector.");

  // Next make sure that q <= Q_a or Q_a_hat
  if(q >= SCM_UB_vectors[0].size())
    libmesh_error_msg("Error: q is too large in set_SCM_UB_vector.");

  SCM_UB_vectors[j][q] = y_q;
}

Real RBSCMEvaluation::get_SCM_UB_vector(unsigned int j, unsigned int q)
{
  // First make sure that j <= J
  if(j >= SCM_UB_vectors.size())
    libmesh_error_msg("Error: We must have j < J in get_SCM_UB_vector.");

  if(q >= SCM_UB_vectors[0].size())
    libmesh_error_msg("Error: q is too large in get_SCM_UB_vector.");

  return SCM_UB_vectors[j][q];
}

const RBParameters & RBSCMEvaluation::get_C_J_entry(unsigned int j)
{
  if(j >= C_J.size())
    libmesh_error_msg("Error: Input parameter j is too large in get_C_J.");

  return C_J[j];
}

Real RBSCMEvaluation::get_B_min(unsigned int q) const
{
  if(q >= B_min.size())
    libmesh_error_msg("Error: q is too large in get_B_min.");

  return B_min[q];
}


Real RBSCMEvaluation::get_B_max(unsigned int q) const
{
  if(q >= B_max.size())
    libmesh_error_msg("Error: q is too large in get_B_max.");

  return B_max[q];
}

void RBSCMEvaluation::set_B_min(unsigned int q, Real B_min_val)
{
  if(q >= B_min.size())
    libmesh_error_msg("Error: q is too large in set_B_min.");

  B_min[q] = B_min_val;
}

void RBSCMEvaluation::set_B_max(unsigned int q, Real B_max_val)
{
  if(q >= B_max.size())
    libmesh_error_msg("Error: q is too large in set_B_max.");

  B_max[q] = B_max_val;
}

Real RBSCMEvaluation::get_SCM_LB()
{
  LOG_SCOPE("get_SCM_LB()", "RBSCMEvaluation");

  // Initialize the LP
  glp_prob * lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp,GLP_MIN);

  // Add columns to the LP: corresponds to
  // the variables y_1,...y_Q_a.
  // These are the same for each \mu in the SCM
  // training set, hence can do this up front.
  glp_add_cols(lp,rb_theta_expansion->get_n_A_terms());

  for(unsigned int q=0; q<rb_theta_expansion->get_n_A_terms(); q++)
    {
      if(B_max[q] < B_min[q]) // Invalid bound, set as free variable
        {
          // GLPK indexing is not zero based!
          glp_set_col_bnds(lp, q+1, GLP_FR, 0., 0.);
        }
      else
        {
          // GLPK indexing is not zero based!
          glp_set_col_bnds(lp, q+1, GLP_DB, B_min[q], B_max[q]);
        }

      // If B_max is not defined, just set lower bounds...
      //       glp_set_col_bnds(lp, q+1, GLP_LO, B_min[q], 0.);
    }


  // Add rows to the LP: corresponds to the auxiliary
  // variables that define the constraints at each
  // mu \in C_J_M
  unsigned int n_rows = C_J.size();
  glp_add_rows(lp, n_rows);

  // Now put current_parameters in saved_parameters
  save_current_parameters();

  unsigned int matrix_size = n_rows*rb_theta_expansion->get_n_A_terms();
  std::vector<int> ia(matrix_size+1);
  std::vector<int> ja(matrix_size+1);
  std::vector<double> ar(matrix_size+1);
  unsigned int count=0;
  for(unsigned int m=0; m<n_rows; m++)
    {
      set_current_parameters_from_C_J(m);

      // Set the lower bound on the auxiliary variable
      // due to the stability constant at mu_index
      glp_set_row_bnds(lp, m+1, GLP_LO, C_J_stability_vector[m], 0.);

      // Now define the matrix that relates the y's
      // to the auxiliary variables at the current
      // value of mu.
      for(unsigned int q=0; q<rb_theta_expansion->get_n_A_terms(); q++)
        {
          count++;

          ia[count] = m+1;
          ja[count] = q+1;

          // This can only handle Reals right now
          ar[count] = libmesh_real( rb_theta_expansion->eval_A_theta(q,get_parameters()) );
        }
    }

  // Now load the original parameters back into current_parameters
  // in order to set the coefficients of the objective function
  reload_current_parameters();

  glp_load_matrix(lp, matrix_size, &ia[0], &ja[0], &ar[0]);

  for(unsigned int q=0; q<rb_theta_expansion->get_n_A_terms(); q++)
    {
      glp_set_obj_coef(lp,q+1, libmesh_real( rb_theta_expansion->eval_A_theta(q,get_parameters()) ) );
    }

  // Use this command to initialize the basis for the LP
  // since default behavior is to use the basis from
  // the previous solve, but that might become singular
  // if we switch the order of constraints (as can
  // happen when we generate a new C_J_M)
  //lpx_cpx_basis(lp); //glp_cpx_basis(lp);

  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_ERR;
  parm.meth = GLP_DUAL;


  // use the simplex method and solve the LP
  glp_simplex(lp, &parm);

  Real min_J_obj = glp_get_obj_val(lp);

  //   int simplex_status =  glp_get_status(lp);
  //   if(simplex_status == GLP_UNBND)
  //   {
  //     libMesh::out << "Simplex method gave unbounded solution." << std::endl;
  //     min_J_obj = std::numeric_limits<Real>::quiet_NaN();
  //   }
  //   else
  //   {
  //     min_J_obj = glp_get_obj_val(lp);
  //   }

  // Destroy the LP
  glp_delete_prob(lp);

  return min_J_obj;
}

Real RBSCMEvaluation::get_SCM_UB()
{
  LOG_SCOPE("get_SCM_UB()", "RBSCMEvaluation");

  // Add rows to the LP: corresponds to the auxiliary
  // variables that define the constraints at each
  // mu \in C_J
  unsigned int n_rows = C_J.size();

  // For each mu, we just find the minimum of J_obj over
  // the subset of vectors in SCM_UB_vectors corresponding
  // to C_J_M (SCM_UB_vectors contains vectors for all of
  // C_J).
  Real min_J_obj = 0.;
  for(unsigned int m=0; m<n_rows; m++)
    {
      const std::vector<Real> UB_vector = SCM_UB_vectors[m];

      Real J_obj = 0.;
      for(unsigned int q=0; q<rb_theta_expansion->get_n_A_terms(); q++)
        {
          J_obj += libmesh_real( rb_theta_expansion->eval_A_theta(q,get_parameters()) )*UB_vector[q];
        }

      if( (m==0) || (J_obj < min_J_obj) )
        {
          min_J_obj = J_obj;
        }
    }

  return min_J_obj;
}

void RBSCMEvaluation::set_current_parameters_from_C_J(unsigned int C_J_index)
{
  set_parameters(C_J[C_J_index]);
}

void RBSCMEvaluation::save_current_parameters()
{
  saved_parameters = get_parameters();
}

void RBSCMEvaluation::reload_current_parameters()
{
  set_parameters(saved_parameters);
}

void RBSCMEvaluation::legacy_write_offline_data_to_files(const std::string & directory_name,
                                                         const bool write_binary_data)
{
  LOG_SCOPE("legacy_write_offline_data_to_files()", "RBSCMEvaluation");

  if(this->processor_id() == 0)
    {
      // Make a directory to store all the data files
      if( mkdir(directory_name.c_str(), 0777) == -1)
        {
          libMesh::out << "In RBSCMEvaluation::write_offline_data_to_files, directory "
                       << directory_name << " already exists, overwriting contents." << std::endl;
        }

      // The writing mode: ENCODE for binary, WRITE for ASCII
      XdrMODE mode = write_binary_data ? ENCODE : WRITE;

      // The suffix to use for all the files that are written out
      const std::string suffix = write_binary_data ? ".xdr" : ".dat";

      // Stream for building the file names
      std::ostringstream file_name;

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

      // Write out the bounding box min values
      file_name.str("");
      file_name << directory_name << "/B_min" << suffix;
      Xdr B_min_out(file_name.str(), mode);

      for (std::size_t i=0; i<B_min.size(); i++)
        {
          Real B_min_i = get_B_min(i);
          B_min_out << B_min_i;
        }
      B_min_out.close();


      // Write out the bounding box max values
      file_name.str("");
      file_name << directory_name << "/B_max" << suffix;
      Xdr B_max_out(file_name.str(), mode);

      for (std::size_t i=0; i<B_max.size(); i++)
        {
          Real B_max_i = get_B_max(i);
          B_max_out << B_max_i;
        }
      B_max_out.close();

      // Write out the length of the C_J data
      file_name.str("");
      file_name << directory_name << "/C_J_length" << suffix;
      Xdr C_J_length_out(file_name.str(), mode);

      unsigned int C_J_length = C_J.size();
      C_J_length_out << C_J_length;
      C_J_length_out.close();

      // Write out C_J_stability_vector
      file_name.str("");
      file_name << directory_name << "/C_J_stability_vector" << suffix;
      Xdr C_J_stability_vector_out(file_name.str(), mode);

      for (std::size_t i=0; i<C_J_stability_vector.size(); i++)
        {
          Real C_J_stability_constraint_i = get_C_J_stability_constraint(i);
          C_J_stability_vector_out << C_J_stability_constraint_i;
        }
      C_J_stability_vector_out.close();

      // Write out C_J
      file_name.str("");
      file_name << directory_name << "/C_J" << suffix;
      Xdr C_J_out(file_name.str(), mode);

      for (std::size_t i=0; i<C_J.size(); i++)
        {
          RBParameters::const_iterator it     = C_J[i].begin();
          RBParameters::const_iterator it_end = C_J[i].end();
          for( ; it != it_end; ++it)
            {
              // Need to make a copy of the value so that it's not const
              // Xdr is not templated on const's
              Real param_value = it->second;
              C_J_out << param_value;
            }
        }
      C_J_out.close();

      // Write out SCM_UB_vectors get_SCM_UB_vector
      file_name.str("");
      file_name << directory_name << "/SCM_UB_vectors" << suffix;
      Xdr SCM_UB_vectors_out(file_name.str(), mode);

      for (std::size_t i=0; i<SCM_UB_vectors.size(); i++)
        for (unsigned int j=0; j<rb_theta_expansion->get_n_A_terms(); j++)
          {
            Real SCM_UB_vector_ij = get_SCM_UB_vector(i,j);
            SCM_UB_vectors_out << SCM_UB_vector_ij;
          }
      SCM_UB_vectors_out.close();
    }
}


void RBSCMEvaluation::legacy_read_offline_data_from_files(const std::string & directory_name,
                                                          const bool read_binary_data)
{
  LOG_SCOPE("legacy_read_offline_data_from_files()", "RBSCMEvaluation");

  // The reading mode: DECODE for binary, READ for ASCII
  XdrMODE mode = read_binary_data ? DECODE : READ;

  // The suffix to use for all the files that are written out
  const std::string suffix = read_binary_data ? ".xdr" : ".dat";

  // The string stream we'll use to make the file names
  std::ostringstream file_name;

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

  // Read in the bounding box min values
  // Note that there are Q_a values
  file_name.str("");
  file_name << directory_name << "/B_min" << suffix;
  Xdr B_min_in(file_name.str(), mode);

  B_min.clear();
  for(unsigned int i=0; i<rb_theta_expansion->get_n_A_terms(); i++)
    {
      Real B_min_val;
      B_min_in >> B_min_val;
      B_min.push_back(B_min_val);
    }
  B_min_in.close();


  // Read in the bounding box max values
  // Note that there are Q_a values
  file_name.str("");
  file_name << directory_name << "/B_max" << suffix;
  Xdr B_max_in(file_name.str(), mode);

  B_max.clear();
  for(unsigned int i=0; i<rb_theta_expansion->get_n_A_terms(); i++)
    {
      Real B_max_val;
      B_max_in >> B_max_val;
      B_max.push_back(B_max_val);
    }

  // Read in the length of the C_J data
  file_name.str("");
  file_name << directory_name << "/C_J_length" << suffix;
  Xdr C_J_length_in(file_name.str(), mode);

  unsigned int C_J_length;
  C_J_length_in >> C_J_length;
  C_J_length_in.close();

  // Read in C_J_stability_vector
  file_name.str("");
  file_name << directory_name << "/C_J_stability_vector" << suffix;
  Xdr C_J_stability_vector_in(file_name.str(), mode);

  C_J_stability_vector.clear();
  for(unsigned int i=0; i<C_J_length; i++)
    {
      Real C_J_stability_val;
      C_J_stability_vector_in >> C_J_stability_val;
      C_J_stability_vector.push_back(C_J_stability_val);
    }
  C_J_stability_vector_in.close();

  // Read in C_J
  file_name.str("");
  file_name << directory_name << "/C_J" << suffix;
  Xdr C_J_in(file_name.str(), mode);

  // Resize C_J based on C_J_stability_vector and Q_a
  C_J.resize( C_J_length );
  for (std::size_t i=0; i<C_J.size(); i++)
    {
      RBParameters::const_iterator it     = get_parameters().begin();
      RBParameters::const_iterator it_end = get_parameters().end();
      for( ; it != it_end; ++it)
        {
          std::string param_name = it->first;
          Real param_value;
          C_J_in >> param_value;
          C_J[i].set_value(param_name, param_value);
        }
    }
  C_J_in.close();


  // Read in SCM_UB_vectors get_SCM_UB_vector
  file_name.str("");
  file_name << directory_name << "/SCM_UB_vectors" << suffix;
  Xdr SCM_UB_vectors_in(file_name.str(), mode);

  // Resize SCM_UB_vectors based on C_J_stability_vector and Q_a
  SCM_UB_vectors.resize( C_J_stability_vector.size() );
  for (std::size_t i=0; i<SCM_UB_vectors.size(); i++)
    {
      SCM_UB_vectors[i].resize( rb_theta_expansion->get_n_A_terms() );
      for(unsigned int j=0; j<rb_theta_expansion->get_n_A_terms(); j++)
        {
          SCM_UB_vectors_in >> SCM_UB_vectors[i][j];
        }
    }
  SCM_UB_vectors_in.close();
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
