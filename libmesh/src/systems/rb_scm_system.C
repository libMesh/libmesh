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

// Currently, the RBSCMSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "rb_scm_system.h"

#include "libmesh_logging.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "equation_systems.h"
#include "getpot.h"
#include "rb_system.h"
#include "parallel.h"

// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// glpk includes
#include <glpk.h>

namespace libMesh
{

RBSCMSystem::RBSCMSystem (EquationSystems& es,
                          const std::string& name,
                          const unsigned int number)
  : Parent(es, name, number),
    parameters_filename(""),
    SCM_eps(0.5),
    SCM_M(1),
    RB_system_name("")
{

  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;

  // We symmetrize all operators hence use symmetric solvers
  set_eigenproblem_type(GHEP);

  // Clear SCM data vectors
  B_min.clear();
  B_max.clear();
  C_J.clear();
  C_J_stability_vector.clear();
  SCM_UB_vectors.clear();
}

RBSCMSystem::~RBSCMSystem ()
{
  this->clear();
}


void RBSCMSystem::clear()
{
  Parent::clear();
}

void RBSCMSystem::initialize_SCM_system ()
{
  // First read in data from parameters_filename
  GetPot infile(parameters_filename);

  const unsigned int n_SCM_parameters = infile("n_SCM_parameters",1);
  const unsigned int n_SCM_training_samples = infile("n_SCM_training_samples",1);
  const bool deterministic_training = infile("deterministic_training",false);

  // Set boolean flag to indicate whether or not we initialize
  // mesh dependent matrices and vectors when init_data
  // is called. Default value is true.
  initialize_mesh_dependent_data = infile("initialize_mesh_dependent_data",
                                          initialize_mesh_dependent_data);

  // SCM Greedy termination tolerance
  const Real SCM_eps_in = infile("SCM_eps", SCM_eps);
  set_SCM_eps(SCM_eps_in);

  // The number of SCM stability constraints to impose
  const unsigned int SCM_M_in = infile("SCM_M", SCM_M);
  set_SCM_M(SCM_M_in);

  // Resize the bounding box vectors
  B_min.resize(rb_theta_expansion->get_Q_a());
  B_max.resize(rb_theta_expansion->get_Q_a());

  std::vector<Real> mu_min_vector(n_SCM_parameters);
  std::vector<Real> mu_max_vector(n_SCM_parameters);
  std::vector<bool> log_scaling(n_SCM_parameters);
  for(unsigned int i=0; i<n_SCM_parameters; i++)
  {
    // Read vector-based mu_min values.
    mu_min_vector[i] = infile("SCM_mu_min", mu_min_vector[i], i);

    // Read vector-based mu_max values.
    mu_max_vector[i] = infile("SCM_mu_max", mu_max_vector[i], i);

    // Read vector-based log scaling values.  Note the intermediate conversion to
    // int... this implies log_scaling = '1 1 1...' in the input file.
    log_scaling[i] = static_cast<bool>(infile("SCM_log_scaling", static_cast<int>(log_scaling[i]), i));
  }

  // Make sure this generates training parameters properly!
  initialize_training_parameters(mu_min_vector,
                                 mu_max_vector,
                                 n_SCM_training_samples,
                                 log_scaling,
                                 deterministic_training);

  libMesh::out << std::endl << "RBSCMSystem parameters:" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "Q_a: " << rb_theta_expansion->get_Q_a() << std::endl;
  libMesh::out << "SCM_eps: " << get_SCM_eps() << std::endl;
  libMesh::out << "SCM_M: " << get_SCM_M() << std::endl;
  for(unsigned int i=0; i<n_SCM_parameters; i++)
  {
    libMesh::out <<   "SCM Parameter " << i
                 << ": Min = " << get_parameter_min(i)
                 << ", Max = " << get_parameter_max(i)
                 << ", log scaling = " << log_scaling[i] << std::endl;
  }
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << "using deterministic training samples? " << deterministic_training << std::endl;
  libMesh::out << "initializing mesh-dependent data structures "
               << initialize_mesh_dependent_data << std::endl;
  libMesh::out << std::endl;
}

void RBSCMSystem::resize_SCM_vectors()
{
  // Clear SCM data vectors
  B_min.clear();
  B_max.clear();
  C_J.clear();
  C_J_stability_vector.clear();
  for(unsigned int i=0; i<SCM_UB_vectors.size(); i++)
    SCM_UB_vectors[i].clear();
  SCM_UB_vectors.clear();

  // Resize the bounding box vectors
  B_min.resize(rb_theta_expansion->get_Q_a());
  B_max.resize(rb_theta_expansion->get_Q_a());
}

void RBSCMSystem::add_scaled_symm_Aq(unsigned int q_a, Number scalar)
{
  START_LOG("add_scaled_symm_Aq()", "RBSCMSystem");

  // Load the operators from the RBSystem
  EquationSystems& es = this->get_equation_systems();
  RBSystem& rb_system = es.get_system<RBSystem>(RB_system_name);

  rb_system.add_scaled_Aq(scalar, q_a, matrix_A, true);

  STOP_LOG("add_scaled_symm_Aq()", "RBSCMSystem");
}

void RBSCMSystem::perform_SCM_greedy()
{
  START_LOG("perform_SCM_greedy()", "RBSCMSystem");

  // Copy the local parts of the Dirichlet and
  // non-Dirichlet dofs lists over from
  // the associated RBSystem
  EquationSystems& es = this->get_equation_systems();
  RBSystem& rb_system = es.get_system<RBSystem>(RB_system_name);
  this->initialize_condensed_dofs(rb_system.global_dirichlet_dofs_set);

  load_matrix_B();

  attach_deflation_space();

  compute_SCM_bounding_box();

  // This loads the new parameter into current_parameters
  enrich_C_J(0);

  unsigned int SCM_iter=0;
  while(true)
  {
    evaluate_stability_constant();

    std::pair<unsigned int,Real> SCM_error_pair = compute_SCM_bounds_on_training_set();

    libMesh::out << "SCM iteration " << SCM_iter
                 << ", max_SCM_error = " << SCM_error_pair.second << std::endl;

    if( SCM_error_pair.second < SCM_eps )
    {
      libMesh::out << std::endl << "SCM tolerance of " << SCM_eps << " reached."
                   << std::endl << std::endl;
      break;
    }

    // If we need another SCM iteration, then enrich C_J
    enrich_C_J(SCM_error_pair.first);

    libMesh::out << std::endl << "-----------------------------------" << std::endl << std::endl;

    SCM_iter++;
  }

  STOP_LOG("perform_SCM_greedy()", "RBSCMSystem");
}

void RBSCMSystem::compute_SCM_bounding_box()
{
  START_LOG("compute_SCM_bounding_box()", "RBSCMSystem");

  for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
  {
    matrix_A->zero();
    add_scaled_symm_Aq(q, 1.);

    // Compute B_min(q)
    eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
    set_eigensolver_properties(q);

    solve();
    unsigned int nconv = get_n_converged();
    if (nconv != 0)
    {
      std::pair<Real, Real> eval = get_eigenpair(0);

      // ensure that the eigenvalue is real
      libmesh_assert(eval.second < TOLERANCE);

      set_B_min(q, eval.first);
      libMesh::out << std::endl << "B_min("<<q<<") = " << get_B_min(q) << std::endl;
    }
    else
    {
      libMesh::err << "Eigen solver for computing B_min did not converge" << std::endl;
      libmesh_error();
    }

    // Compute B_max(q)
    eigen_solver->set_position_of_spectrum(LARGEST_REAL);
    set_eigensolver_properties(q);

    solve();
    nconv = get_n_converged();
    if (nconv != 0)
    {
      std::pair<Real, Real> eval = get_eigenpair(0);

      // ensure that the eigenvalue is real
      libmesh_assert(eval.second < TOLERANCE);

      set_B_max(q,eval.first);
      libMesh::out << "B_max("<<q<<") = " << get_B_max(q) << std::endl;
    }
    else
    {
      libMesh::err << "Eigen solver for computing B_max did not converge" << std::endl;
      libmesh_error();
    }
  }

  STOP_LOG("compute_SCM_bounding_box()", "RBSCMSystem");
}

void RBSCMSystem::evaluate_stability_constant()
{
  START_LOG("evaluate_stability_constant()", "RBSCMSystem");

  // Get current index of C_J
  const unsigned int j = C_J.size()-1;

  eigen_solver->set_position_of_spectrum(SMALLEST_REAL);

  // We assume B is set in system assembly
  // For coercive problems, B is set to the inner product matrix
  // For non-coercive time-dependent problems, B is set to the mass matrix

  // Set matrix A corresponding to mu_star
  matrix_A->zero();
  for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
  {
    add_scaled_symm_Aq(q, rb_theta_expansion->eval_theta_q_a(q));
  }

  set_eigensolver_properties(-1);
  solve();
  unsigned int nconv = get_n_converged();
  if (nconv != 0)
  {
    std::pair<Real, Real> eval = get_eigenpair(0);

    // ensure that the eigenvalue is real
    libmesh_assert(eval.second < TOLERANCE);

    // Store the coercivity constant corresponding to mu_star
    set_C_J_stability_constraint(j,eval.first);
    libMesh::out << std::endl << "Stability constant for C_J("<<j<<") = "
                 << get_C_J_stability_constraint(j) << std::endl << std::endl;

    // Compute and store the vector y = (y_1, \ldots, y_Q) for the
    // eigenvector currently stored in eigen_system.solution.
    // We use this later to compute the SCM upper bounds.
    Real norm_B2 = libmesh_real( B_inner_product(*solution, *solution) );

    for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
    {
      Real norm_Aq2 = libmesh_real( Aq_inner_product(q, *solution, *solution) );

      set_SCM_UB_vector(j,q,norm_Aq2/norm_B2);
    }
  }
  else
  {
    libMesh::err << "Error: Eigensolver did not converge in evaluate_stability_constant"
                 << std::endl;
    libmesh_error();
  }

  STOP_LOG("evaluate_stability_constant()", "RBSCMSystem");
}

Number RBSCMSystem::B_inner_product(const NumericVector<Number>& v, const NumericVector<Number>& w) const
{
  matrix_B->vector_mult(*inner_product_storage_vector, w);

  return v.dot(*inner_product_storage_vector);
}

Number RBSCMSystem::Aq_inner_product(unsigned int q,
                                    const NumericVector<Number>& v,
                                    const NumericVector<Number>& w)
{
  if(q >= rb_theta_expansion->get_Q_a())
  {
    libMesh::err << "Error: We must have q < Q_a in Aq_inner_product."
                 << std::endl;
    libmesh_error();
  }

  matrix_A->zero();
  add_scaled_symm_Aq(q, 1.);
  matrix_A->vector_mult(*inner_product_storage_vector, w);

  return v.dot(*inner_product_storage_vector);
}


void RBSCMSystem::set_C_J_stability_constraint(unsigned int j, Real stability_const_in)
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

Real RBSCMSystem::get_C_J_stability_constraint(unsigned int j) const
{
  if(j >= C_J_stability_vector.size())
  {
    libMesh::err << "Error: Input parameter j is too large in get_C_J_stability_constraint.";
    libmesh_error();
  }

  return C_J_stability_vector[j];
}

void RBSCMSystem::set_SCM_UB_vector(unsigned int j, unsigned int q, Real y_q)
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

Real RBSCMSystem::get_SCM_UB_vector(unsigned int j, unsigned int q)
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

void RBSCMSystem::enrich_C_J(unsigned int new_C_J_index)
{
  START_LOG("enrich_C_J()", "RBSCMSystem");

  load_training_parameter_globally(new_C_J_index);

  C_J.push_back(get_current_parameters());

  libMesh::out << std::endl << "SCM: Added mu = (";
  for(unsigned int i=0; i<get_n_params(); i++)
  {
    libMesh::out << C_J[C_J.size()-1][i];
    if(i < (get_n_params()-1)) libMesh::out << ",";
  }
  libMesh::out << ")" << std::endl;

  // Finally, resize C_J_stability_vector and SCM_UB_vectors
  C_J_stability_vector.push_back(0.);

  std::vector<Real> zero_vector(rb_theta_expansion->get_Q_a());
  SCM_UB_vectors.push_back(zero_vector);

  STOP_LOG("enrich_C_J()", "RBSCMSystem");
}

std::vector<Real> RBSCMSystem::get_C_J_entry(unsigned int j)
{
  if(j >= C_J.size())
  {
    libMesh::err << "Error: Input parameter j is too large in get_C_J.";
    libmesh_error();
  }

  return C_J[j];
}

Real RBSCMSystem::get_B_min(unsigned int q) const
{
  if(q >= B_min.size())
  {
    libMesh::err << "Error: q is too large in get_B_min."
                 << std::endl;
    libmesh_error();
  }

  return B_min[q];
}


Real RBSCMSystem::get_B_max(unsigned int q) const
{
  if(q >= B_max.size())
  {
    libMesh::err << "Error: q is too large in get_B_max."
                 << std::endl;
    libmesh_error();
  }

  return B_max[q];
}

void RBSCMSystem::set_B_min(unsigned int q, Real B_min_val)
{
  if(q >= B_min.size())
  {
    libMesh::err << "Error: q is too large in set_B_min."
                 << std::endl;
    libmesh_error();
  }

  B_min[q] = B_min_val;
}

void RBSCMSystem::set_B_max(unsigned int q, Real B_max_val)
{
  if(q >= B_max.size())
  {
    libMesh::err << "Error: q is too large in set_B_max."
                 << std::endl;
    libmesh_error();
  }

  B_max[q] = B_max_val;
}

std::pair<unsigned int,Real> RBSCMSystem::compute_SCM_bounds_on_training_set()
{
  START_LOG("compute_SCM_bounds_on_training_set()", "RBSCMSystem");

  // Now compute the maximum bound error over training_parameters
  unsigned int new_C_J_index = 0;
  Real max_SCM_error = 0.;

  unsigned int first_index = get_first_local_training_index();
  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
  {
    load_training_parameter_locally(first_index+i);

    Real LB = get_SCM_LB();
    Real UB = get_SCM_UB();

    Real error_i = SCM_greedy_error_indicator(LB, UB);

    if( error_i > max_SCM_error )
    {
      max_SCM_error = error_i;
      new_C_J_index = i;
    }
  }

  unsigned int global_index = first_index + new_C_J_index;
  std::pair<unsigned int,Real> error_pair(global_index, max_SCM_error);
  get_global_max_error_pair(error_pair);

  STOP_LOG("compute_SCM_bounds_on_training_set()", "RBSCMSystem");

  return error_pair;
}

Real RBSCMSystem::get_SCM_LB()
{
  START_LOG("get_SCM_LB()", "RBSCMSystem");

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
  unsigned int n_rows = ( SCM_M < C_J.size()) ? SCM_M : C_J.size();
  glp_add_rows(lp, n_rows);

  // Find the distance between mu and each element of C_J
  const unsigned int J = C_J.size();
  std::vector< std::pair<Real,unsigned int> > dist_from_mu(J);
  for(unsigned int j=0; j<J; j++)
  {
    dist_from_mu[j].first  = param_dist(get_current_parameters(), C_J[j]);
    dist_from_mu[j].second = j;
  }

  // Now put current_parameters in saved_parameters
  save_current_parameters();

  // Now sort dist_from_mu based on the first element of each pair,
  // the first M elements of this sorted list give us C_J_M.
  // We use the default "less than" operator for std::pair, which
  // compares based on the first element (which is the dist for us)
  // and if they're the same, then compares the second elements.
  std::sort( dist_from_mu.begin(), dist_from_mu.end() );


  unsigned int matrix_size = n_rows*rb_theta_expansion->get_Q_a();
  std::vector<int> ia(matrix_size+1);
  std::vector<int> ja(matrix_size+1);
  std::vector<Real> ar(matrix_size+1);
  unsigned int count=0;
  for(unsigned int m=0; m<n_rows; m++)
  {
    unsigned int mu_index = dist_from_mu[m].second;
    set_current_parameters_from_C_J(mu_index);

    // Set the lower bound on the auxiliary variable
    // due to the stability constant at mu_index
    glp_set_row_bnds(lp, m+1, GLP_LO, C_J_stability_vector[mu_index], 0.);

    // Now define the matrix that relates the y's
    // to the auxiliary variables at the current
    // value of mu.
    for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
      {
        count++;

        ia[count] = m+1;
        ja[count] = q+1;
        ar[count] = libmesh_real( rb_theta_expansion->eval_theta_q_a(q) ); // This can only handle Reals right now
      }
  }
  glp_load_matrix(lp, matrix_size, &ia[0], &ja[0], &ar[0]);

  // Now load the original parameters back into current_parameters
  // in order to set the coefficients of the objective function
  reload_current_parameters();
  for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
    {
      glp_set_obj_coef(lp,q+1, libmesh_real( rb_theta_expansion->eval_theta_q_a(q) ) );
    }

  // Use this command to initialize the basis for the LP
  // since default behavior is to use the basis from
  // the previous solve, but that might become singular
  // if we switch the order of constraints (as can
  // happen when we generate a new C_J_M)
  lpx_cpx_basis(lp); //glp_cpx_basis(lp);

  glp_simplex(lp, NULL);

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

  STOP_LOG("get_SCM_LB()", "RBSCMSystem");

  return min_J_obj;
}

Real RBSCMSystem::get_SCM_UB()
{
  START_LOG("get_SCM_UB()", "RBSCMSystem");

  // Add rows to the LP: corresponds to the auxiliary
  // variables that define the constraints at each
  // mu \in C_J_M
  unsigned int n_rows = ( SCM_M < C_J.size()) ? SCM_M : C_J.size();

  // Find the distance between mu and each element of C_J
  const unsigned int J = C_J.size();
  std::vector< std::pair<Real,unsigned int> > dist_from_mu(J);
  for(unsigned int j=0; j<J; j++)
  {
    dist_from_mu[j].first  = param_dist(get_current_parameters(), C_J[j]);
    dist_from_mu[j].second = j;
  }

  // Now sort dist_from_mu based on the first element of each pair,
  // the first M elements of this sorted list give us C_J_M.
  // We use the default "less than" operator for std::pair, which
  // compares based on the first element (which is the dist for us)
  // and if they're the same, then compares the second elements.
  std::sort( dist_from_mu.begin(), dist_from_mu.end() );

  // For each mu, we just find the minimum of J_obj over
  // the subset of vectors in SCM_UB_vectors corresponding
  // to C_J_M (SCM_UB_vectors contains vectors for all of
  // C_J).
  Real min_J_obj = 0.;
  for(unsigned int m=0; m<n_rows; m++)
  {
    // This is the index from C_J of the \mu in C_J_M
    unsigned int mu_index = dist_from_mu[m].second;

    const std::vector<Real> UB_vector = SCM_UB_vectors[mu_index];

    Real J_obj = 0.;
    for(unsigned int q=0; q<rb_theta_expansion->get_Q_a(); q++)
      {
        J_obj += libmesh_real( rb_theta_expansion->eval_theta_q_a(q) )*UB_vector[q];
      }

    if( (m==0) || (J_obj < min_J_obj) )
    {
      min_J_obj = J_obj;
    }
  }

  STOP_LOG("get_SCM_UB()", "RBSCMSystem");

  return min_J_obj;
}

void RBSCMSystem::set_current_parameters_from_C_J(unsigned int C_J_index)
{
  set_current_parameters(C_J[C_J_index]);
}

Real RBSCMSystem::param_dist(const std::vector<Real>& mu_1, const std::vector<Real>& mu_2)
{
  // Default distance is Euclidean norm
  Real sum = 0.;

  for(unsigned int i=0; i<get_n_params(); i++)
  {
    sum += pow(mu_1[i] - mu_2[i],2.);
  }

  return sqrt(sum);
}

void RBSCMSystem::save_current_parameters()
{
  saved_parameters = current_parameters;
}

void RBSCMSystem::reload_current_parameters()
{
  current_parameters = saved_parameters;
}

void RBSCMSystem::write_offline_data_to_files(const std::string& directory_name)
{
  START_LOG("write_offline_data_to_files()", "RBSCMSystem");

  const unsigned int precision_level = 14;

  if(libMesh::processor_id() == 0)
  {
    // Make a directory to store all the data files
    if( mkdir(directory_name.c_str(), 0777) == -1)
    {
      libMesh::out << "In RBSCMSystem::write_offline_data_to_files, directory "
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

  STOP_LOG("write_offline_data_to_files()", "RBSCMSystem");
}

void RBSCMSystem::read_offline_data_from_files(const std::string& directory_name)
{
  START_LOG("read_offline_data_from_files()", "RBSCMSystem");

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
  for(unsigned int i=0; i<B_min.size(); i++)
  {
    B_min_in >> B_min[i];
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
  for(unsigned int i=0; i<B_max.size(); i++)
  {
    B_max_in >> B_max[i];
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

  STOP_LOG("read_offline_data_from_files()", "RBSCMSystem");
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
