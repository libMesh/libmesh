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

// Currently, the RBSCMEvaluation should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "rb_scm_evaluation.h"

#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "equation_systems.h"
#include "getpot.h"
#include "parallel.h"
#include "dof_map.h"
// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// glpk includes
#include <glpk.h>

namespace libMesh
{

RBSCMEvaluation::RBSCMEvaluation ()
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

void RBSCMEvaluation::set_C_J_stability_constraint(unsigned int j, Real stability_const_in)
{
  if(j >= C_J_stability_vector.size())
  {
    libMesh::err << "Error: Input parameter j is too large in set_C_J_stability_constraint.";
    libmesh_error();
  }

  // we assume that C_J_stability_vector is resized elsewhere
  // to be the same size as C_J.
  libmesh_assert(C_J_stability_vector.size() == C_J.size());

  C_J_stability_vector[j] = stability_const_in;
}

Real RBSCMEvaluation::get_C_J_stability_constraint(unsigned int j) const
{
  if(j >= C_J_stability_vector.size())
  {
    libMesh::err << "Error: Input parameter j is too large in get_C_J_stability_constraint.";
    libmesh_error();
  }

  return C_J_stability_vector[j];
}

void RBSCMEvaluation::set_SCM_UB_vector(unsigned int j, unsigned int q, Real y_q)
{
  // First make sure that j <= J
  if(j >= SCM_UB_vectors.size())
  {
    libMesh::err << "Error: We must have j < J in set_SCM_UB_vector.";
    libmesh_error();
  }
  // Next make sure that q <= Q_a or Q_a_hat
  if(q >= SCM_UB_vectors[0].size())
  {
    libMesh::err << "Error: q is too large in set_SCM_UB_vector."
                 << std::endl;
    libmesh_error();
  }

  SCM_UB_vectors[j][q] = y_q;
}

Real RBSCMEvaluation::get_SCM_UB_vector(unsigned int j, unsigned int q)
{
  // First make sure that j <= J
  if(j >= SCM_UB_vectors.size())
  {
    libMesh::err << "Error: We must have j < J in get_SCM_UB_vector.";
    libmesh_error();
  }
  if(q >= SCM_UB_vectors[0].size())
  {
    libMesh::err << "Error: q is too large in get_SCM_UB_vector."
                 << std::endl;
    libmesh_error();
  }

  return SCM_UB_vectors[j][q];
}

std::vector<Real> RBSCMEvaluation::get_C_J_entry(unsigned int j)
{
  if(j >= C_J.size())
  {
    libMesh::err << "Error: Input parameter j is too large in get_C_J.";
    libmesh_error();
  }

  return C_J[j];
}

Real RBSCMEvaluation::get_B_min(unsigned int q) const
{
  if(q >= B_min.size())
  {
    libMesh::err << "Error: q is too large in get_B_min."
                 << std::endl;
    libmesh_error();
  }

  return B_min[q];
}


Real RBSCMEvaluation::get_B_max(unsigned int q) const
{
  if(q >= B_max.size())
  {
    libMesh::err << "Error: q is too large in get_B_max."
                 << std::endl;
    libmesh_error();
  }

  return B_max[q];
}

void RBSCMEvaluation::set_B_min(unsigned int q, Real B_min_val)
{
  if(q >= B_min.size())
  {
    libMesh::err << "Error: q is too large in set_B_min."
                 << std::endl;
    libmesh_error();
  }

  B_min[q] = B_min_val;
}

void RBSCMEvaluation::set_B_max(unsigned int q, Real B_max_val)
{
  if(q >= B_max.size())
  {
    libMesh::err << "Error: q is too large in set_B_max."
                 << std::endl;
    libmesh_error();
  }

  B_max[q] = B_max_val;
}

Real RBSCMEvaluation::get_SCM_LB()
{
  START_LOG("get_SCM_LB()", "RBSCMEvaluation");

  // Initialize the LP
  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp,GLP_MIN);

  // Add columns to the LP: corresponds to
  // the variables y_1,...y_Q_a.
  // These are the same for each \mu in the SCM
  // training set, hence can do this up front.
  glp_add_cols(lp,rb_theta_expansion->get_Q_a());

  for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
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

  unsigned int matrix_size = n_rows*rb_theta_expansion->get_Q_a();
  std::vector<int> ia(matrix_size+1);
  std::vector<int> ja(matrix_size+1);
  std::vector<Real> ar(matrix_size+1);
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
    for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
      {
        count++;

        ia[count] = m+1;
        ja[count] = q+1;

        // This can only handle Reals right now
        ar[count] = libmesh_real( rb_theta_expansion->eval_theta_q_a(q,get_current_parameters()) );
      }
  }
  glp_load_matrix(lp, matrix_size, &ia[0], &ja[0], &ar[0]);

  // Now load the original parameters back into current_parameters
  // in order to set the coefficients of the objective function
  reload_current_parameters();
  for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
    {
      glp_set_obj_coef(lp,q+1, libmesh_real( rb_theta_expansion->eval_theta_q_a(q,get_current_parameters()) ) );
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

  STOP_LOG("get_SCM_LB()", "RBSCMEvaluation");

  return min_J_obj;
}

Real RBSCMEvaluation::get_SCM_UB()
{
  START_LOG("get_SCM_UB()", "RBSCMEvaluation");

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
    for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
      {
        J_obj += libmesh_real( rb_theta_expansion->eval_theta_q_a(q,get_current_parameters()) )*UB_vector[q];
      }

    if( (m==0) || (J_obj < min_J_obj) )
    {
      min_J_obj = J_obj;
    }
  }

  STOP_LOG("get_SCM_UB()", "RBSCMEvaluation");

  return min_J_obj;
}

void RBSCMEvaluation::set_current_parameters_from_C_J(unsigned int C_J_index)
{
  set_current_parameters(C_J[C_J_index]);
}

void RBSCMEvaluation::save_current_parameters()
{
  saved_parameters = current_parameters;
}

void RBSCMEvaluation::reload_current_parameters()
{
  current_parameters = saved_parameters;
}

void RBSCMEvaluation::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBSCMEvaluation");

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Make a directory to store all the data files
    if( mkdir(directory_name.c_str(), 0777) == -1)
    {
      libMesh::out << "In RBSCMEvaluation::write_offline_data_to_files, directory "
                   << directory_name << " already exists, overwriting contents." << std::endl;
    }

    // Write out the bounding box min values
    std::ofstream B_min_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/B_min.dat";
      B_min_out.open(file_name.str().c_str());
    }
    if ( !B_min_out.good() )
    {
      libMesh::err << "Error opening B_min.dat" << std::endl;
      libmesh_error();
    }
    B_min_out.precision(precision_level);
    for(unsigned int i=0; i<B_min.size(); i++)
    {
      B_min_out << std::scientific << get_B_min(i) << " ";
    }
    B_min_out << std::endl;
    B_min_out.close();


    // Write out the bounding box max values
    std::ofstream B_max_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/B_max.dat";
      B_max_out.open(file_name.str().c_str());
    }
    if ( !B_max_out.good() )
    {
      libMesh::err << "Error opening B_max.dat" << std::endl;
      libmesh_error();
    }
    B_max_out.precision(precision_level);
    for(unsigned int i=0; i<B_max.size(); i++)
    {
      B_max_out << std::scientific << get_B_max(i) << " ";
    }
    B_max_out << std::endl;
    B_max_out.close();

    // Write out C_J_stability_vector
    std::ofstream C_J_stability_vector_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/C_J_stability_vector.dat";
      C_J_stability_vector_out.open(file_name.str().c_str());
    }
    if ( !C_J_stability_vector_out.good() )
    {
      libMesh::err << "Error opening C_J_stability_vector.dat" << std::endl;
      libmesh_error();
    }
    C_J_stability_vector_out.precision(precision_level);
    for(unsigned int i=0; i<C_J_stability_vector.size(); i++)
    {
      C_J_stability_vector_out << std::scientific << get_C_J_stability_constraint(i) << " ";
    }
    C_J_stability_vector_out << std::endl;
    C_J_stability_vector_out.close();

    // Write out C_J
    std::ofstream C_J_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/C_J.dat";
      C_J_out.open(file_name.str().c_str());
    }
    if ( !C_J_out.good() )
    {
      libMesh::err << "Error opening C_J_out.dat" << std::endl;
      libmesh_error();
    }
    C_J_out.precision(precision_level);
    for(unsigned int i=0; i<C_J.size(); i++)
    {
      for(unsigned int j=0; j<get_n_params(); j++)
      {
        C_J_out << std::scientific << C_J[i][j] << " ";
      }
    }
    C_J_out << std::endl;
    C_J_out.close();

    // Write out SCM_UB_vectors get_SCM_UB_vector
    std::ofstream SCM_UB_vectors_out;
    {
      OStringStream file_name;
      file_name << directory_name << "/SCM_UB_vectors.dat";
      SCM_UB_vectors_out.open(file_name.str().c_str());
    }
    if ( !SCM_UB_vectors_out.good() )
    {
      libMesh::err << "Error opening SCM_UB_vectors.dat" << std::endl;
      libmesh_error();
    }
    SCM_UB_vectors_out.precision(precision_level);
    for(unsigned int i=0; i<SCM_UB_vectors.size(); i++)
    {
      for(unsigned int j=0; j<rb_theta_expansion->get_Q_a(); j++)
      {
        SCM_UB_vectors_out << std::scientific << get_SCM_UB_vector(i,j) << " ";
      }
    }
    SCM_UB_vectors_out << std::endl;
    SCM_UB_vectors_out.close();
  }

  STOP_LOG("write_offline_data_to_files()", "RBSCMEvaluation");
}


void RBSCMEvaluation::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBSCMEvaluation");

  // Read in the bounding box min values
  std::ifstream B_min_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/B_min.dat";
    B_min_in.open(file_name.str().c_str());
  }
  if ( !B_min_in.good() )
  {
    libMesh::err << "Error opening B_min.dat" << std::endl;
    libmesh_error();
  }
  B_min.clear();
  Real B_min_val;
  while( B_min_in >> B_min_val )
  {
    B_min.push_back(B_min_val);
  }
  B_min_in.close();


  // Read in the bounding box max values
  std::ifstream B_max_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/B_max.dat";
    B_max_in.open(file_name.str().c_str());
  }
  if ( !B_max_in.good() )
  {
    libMesh::err << "Error opening B_max.dat" << std::endl;
    libmesh_error();
  }
  B_max.clear();
  Real B_max_val;
  while( B_max_in >> B_max_val )
  {
    B_max.push_back(B_max_val);
  }
  B_max_in.close();


  // Read in C_J_stability_vector
  std::ifstream C_J_stability_vector_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/C_J_stability_vector.dat";
    C_J_stability_vector_in.open(file_name.str().c_str());
  }
  if ( !C_J_stability_vector_in.good() )
  {
    libMesh::err << "Error opening C_J_stability_vector.dat" << std::endl;
    libmesh_error();
  }
  C_J_stability_vector.clear();
  Real C_J_stability_val;
  while( C_J_stability_vector_in >> C_J_stability_val )
  {
    C_J_stability_vector.push_back(C_J_stability_val);
  }
  C_J_stability_vector_in.close();

  // Read in C_J
  std::ifstream C_J_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/C_J.dat";
    C_J_in.open(file_name.str().c_str());
  }
  if ( !C_J_in.good() )
  {
    libMesh::err << "Error opening C_J_in.dat" << std::endl;
    libmesh_error();
  }
  // Resize C_J based on C_J_stability_vector and Q_a
  C_J.resize( C_J_stability_vector.size() );
  for(unsigned int i=0; i<C_J.size(); i++)
  {
    C_J[i].resize(get_n_params());
    for(unsigned int j=0; j<get_n_params(); j++)
    {
      C_J_in >> C_J[i][j];
    }
  }
  C_J_in.close();


  // Read in SCM_UB_vectors get_SCM_UB_vector
  std::ifstream SCM_UB_vectors_in;
  {
    OStringStream file_name;
    file_name << directory_name << "/SCM_UB_vectors.dat";
    SCM_UB_vectors_in.open(file_name.str().c_str());
  }
  if ( !SCM_UB_vectors_in.good() )
  {
    libMesh::err << "Error opening SCM_UB_vectors.dat" << std::endl;
    libmesh_error();
  }
  // Resize SCM_UB_vectors based on C_J_stability_vector and Q_a
  SCM_UB_vectors.resize( C_J_stability_vector.size() );
  for(unsigned int i=0; i<SCM_UB_vectors.size(); i++)
  {
    SCM_UB_vectors[i].resize( rb_theta_expansion->get_Q_a() );
    for(unsigned int j=0; j<rb_theta_expansion->get_Q_a(); j++)
    {
      SCM_UB_vectors_in >> SCM_UB_vectors[i][j];
    }
  }
  SCM_UB_vectors_in.close();

  STOP_LOG("read_offline_data_from_files()", "RBSCMEvaluation");
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
