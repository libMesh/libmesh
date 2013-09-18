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

// Currently, the RBSCMConstruction should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC) && (LIBMESH_HAVE_GLPK)

#include "libmesh/rb_scm_construction.h"
#include "libmesh/rb_construction.h"

#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"
#include "libmesh/parallel.h"
#include "libmesh/dof_map.h"
// For creating a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

namespace libMesh
{

RBSCMConstruction::RBSCMConstruction (EquationSystems& es,
                          const std::string& name,
                          const unsigned int number)
  : Parent(es, name, number),
    SCM_training_tolerance(0.5),
    RB_system_name(""),
    rb_scm_eval(NULL)
{

  // set assemble_before_solve flag to false
  // so that we control matrix assembly.
  assemble_before_solve = false;

  // We symmetrize all operators hence use symmetric solvers
  set_eigenproblem_type(GHEP);

}

RBSCMConstruction::~RBSCMConstruction ()
{
  this->clear();
}


void RBSCMConstruction::clear()
{
  Parent::clear();
}

void RBSCMConstruction::set_rb_scm_evaluation(RBSCMEvaluation& rb_scm_eval_in)
{
  rb_scm_eval = &rb_scm_eval_in;
}

RBSCMEvaluation& RBSCMConstruction::get_rb_scm_evaluation()
{
  if(!rb_scm_eval)
  {
    libMesh::out << "Error: RBSCMEvaluation object hasn't been initialized yet" << std::endl;
    libmesh_error();
  }

  return *rb_scm_eval;
}

RBThetaExpansion& RBSCMConstruction::get_rb_theta_expansion()
{
  return get_rb_scm_evaluation().get_rb_theta_expansion();
}

void RBSCMConstruction::process_parameters_file(const std::string& parameters_filename)
{
  // First read in data from parameters_filename
  GetPot infile(parameters_filename);
  const unsigned int n_training_samples = infile("n_training_samples",1);
  const bool deterministic_training     = infile("deterministic_training",false);

  // Read in training_parameters_random_seed value.  This is used to
  // seed the RNG when picking the training parameters.  By default the
  // value is -1, which means use std::time to seed the RNG.
  unsigned int training_parameters_random_seed_in = static_cast<int>(-1);
  training_parameters_random_seed_in = infile("training_parameters_random_seed",
					   training_parameters_random_seed_in);
  set_training_random_seed(training_parameters_random_seed_in);

  // SCM Greedy termination tolerance
  const Real SCM_training_tolerance_in = infile("SCM_training_tolerance", SCM_training_tolerance);
  set_SCM_training_tolerance(SCM_training_tolerance_in);

  // Initialize the parameter ranges and the parameters themselves
  const unsigned int n_parameters = infile("n_parameters",1);
  RBParameters mu_min_in;
  RBParameters mu_max_in;
  RBParameters initial_mu_in;
  for(unsigned int i=0; i<n_parameters; i++)
  {
    // Read in the parameter names
    std::string param_name = infile("parameter_names", "NONE", i);

    for(unsigned int j=0; j<3; j++)
    {
      if(j==0)
      {
        Real min_val = infile(param_name, 0., j);
        mu_min_in.set_value(param_name, min_val);
      }
      else if(j==1)
      {
        Real max_val = infile(param_name, 0., j);
        mu_max_in.set_value(param_name, max_val);
      }
      else
      {
        Real init_val = infile(param_name, 0., j);
        initial_mu_in.set_value(param_name, init_val);
      }
    }
  }
  initialize_parameters(mu_min_in, mu_max_in, initial_mu_in);

  std::map<std::string,bool> log_scaling;
  const RBParameters& mu = get_parameters();
  RBParameters::const_iterator it     = mu.begin();
  RBParameters::const_iterator it_end = mu.end();
  unsigned int i=0;
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    log_scaling[param_name] = static_cast<bool>(infile("log_scaling", 0, i));
    i++;
  }

  initialize_training_parameters(this->get_parameters_min(),
                                 this->get_parameters_max(),
                                 n_training_samples,
                                 log_scaling,
                                 deterministic_training);   // use deterministic parameters
}

void RBSCMConstruction::print_info()
{
  // Print out info that describes the current setup
  libMesh::out << std::endl << "RBSCMConstruction parameters:" << std::endl;
  libMesh::out << "system name: " << this->name() << std::endl;
  libMesh::out << "SCM Greedy tolerance: " << get_SCM_training_tolerance() << std::endl;
  if(rb_scm_eval)
  {
    libMesh::out << "A_q operators attached: " << get_rb_theta_expansion().get_n_A_terms() << std::endl;
    libMesh::out << "Number of parameters: "   << get_n_params() << std::endl;
  }
  else
  {
    libMesh::out << "RBThetaExpansion member is not set yet" << std::endl;
  }
  RBParameters::const_iterator it     = get_parameters().begin();
  RBParameters::const_iterator it_end = get_parameters().end();
  for( ; it != it_end; ++it)
  {
    std::string param_name = it->first;
    libMesh::out <<   "Parameter " << param_name
                 << ": Min = " << get_parameter_min(param_name)
                 << ", Max = " << get_parameter_max(param_name)
                 << ", value = " << get_parameters().get_value(param_name) << std::endl;
  }
  libMesh::out << "n_training_samples: " << get_n_training_samples() << std::endl;
  libMesh::out << std::endl;
}

void RBSCMConstruction::resize_SCM_vectors()
{
  // Clear SCM data vectors
  rb_scm_eval->B_min.clear();
  rb_scm_eval->B_max.clear();
  rb_scm_eval->C_J.clear();
  rb_scm_eval->C_J_stability_vector.clear();
  for(unsigned int i=0; i<rb_scm_eval->SCM_UB_vectors.size(); i++)
    rb_scm_eval->SCM_UB_vectors[i].clear();
  rb_scm_eval->SCM_UB_vectors.clear();

  // Resize the bounding box vectors
  rb_scm_eval->B_min.resize(get_rb_theta_expansion().get_n_A_terms());
  rb_scm_eval->B_max.resize(get_rb_theta_expansion().get_n_A_terms());
}

void RBSCMConstruction::add_scaled_symm_Aq(unsigned int q_a, Number scalar)
{
  START_LOG("add_scaled_symm_Aq()", "RBSCMConstruction");
  // Load the operators from the RBConstruction
  EquationSystems& es = this->get_equation_systems();
  RBConstruction& rb_system = es.get_system<RBConstruction>(RB_system_name);
  rb_system.add_scaled_Aq(scalar, q_a, matrix_A, true);
  STOP_LOG("add_scaled_symm_Aq()", "RBSCMConstruction");
}

void RBSCMConstruction::load_matrix_B()
{
  // Load the operators from the RBConstruction
  EquationSystems& es = this->get_equation_systems();
  RBConstruction& rb_system = es.get_system<RBConstruction>(RB_system_name);

  matrix_B->zero();
  matrix_B->close();
  matrix_B->add(1.,*rb_system.get_inner_product_matrix());
}

void RBSCMConstruction::perform_SCM_greedy()
{
  START_LOG("perform_SCM_greedy()", "RBSCMConstruction");

  // initialize rb_scm_eval's parameters
  rb_scm_eval->initialize_parameters(*this);

  // Get a list of constrained dofs from rb_system
  std::set<unsigned int> constrained_dofs_set;
  EquationSystems& es = this->get_equation_systems();
  RBConstruction& rb_system = es.get_system<RBConstruction>(RB_system_name);

  for(unsigned int i=0; i<rb_system.n_dofs(); i++)
  {
    if( rb_system.get_dof_map().is_constrained_dof(i) )
    {
      constrained_dofs_set.insert(i);
    }
  }

  // Use these constrained dofs to identify which dofs we want to "get rid of"
  // (i.e. condense) in our eigenproblems.
  this->initialize_condensed_dofs(constrained_dofs_set);

  // Copy the inner product matrix over from rb_system to be used as matrix_B
  load_matrix_B();

  attach_deflation_space();

  compute_SCM_bounding_box();
  // This loads the new parameter into current_parameters
  enrich_C_J(0);

  unsigned int SCM_iter=0;
  while(true)
  {
    // matrix_A is reinitialized for the current parameters
    // on each call to evaluate_stability_constant
    evaluate_stability_constant();

    std::pair<unsigned int,Real> SCM_error_pair = compute_SCM_bounds_on_training_set();

    libMesh::out << "SCM iteration " << SCM_iter
                 << ", max_SCM_error = " << SCM_error_pair.second << std::endl;

    if( SCM_error_pair.second < SCM_training_tolerance )
    {
      libMesh::out << std::endl << "SCM tolerance of " << SCM_training_tolerance << " reached."
                   << std::endl << std::endl;
      break;
    }

    // If we need another SCM iteration, then enrich C_J
    enrich_C_J(SCM_error_pair.first);

    libMesh::out << std::endl << "-----------------------------------" << std::endl << std::endl;

    SCM_iter++;
  }

  STOP_LOG("perform_SCM_greedy()", "RBSCMConstruction");
}

void RBSCMConstruction::compute_SCM_bounding_box()
{
  START_LOG("compute_SCM_bounding_box()", "RBSCMConstruction");

  // Resize the bounding box vectors
  rb_scm_eval->B_min.resize(get_rb_theta_expansion().get_n_A_terms());
  rb_scm_eval->B_max.resize(get_rb_theta_expansion().get_n_A_terms());

  for(unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
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
      libmesh_assert_less (eval.second, TOLERANCE);

      rb_scm_eval->set_B_min(q, eval.first);
      libMesh::out << std::endl << "B_min("<<q<<") = " << rb_scm_eval->get_B_min(q) << std::endl;
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
      libmesh_assert_less (eval.second, TOLERANCE);

      rb_scm_eval->set_B_max(q,eval.first);
      libMesh::out << "B_max("<<q<<") = " << rb_scm_eval->get_B_max(q) << std::endl;
    }
    else
    {
      libMesh::err << "Eigen solver for computing B_max did not converge" << std::endl;
      libmesh_error();
    }
  }

  STOP_LOG("compute_SCM_bounding_box()", "RBSCMConstruction");
}

void RBSCMConstruction::evaluate_stability_constant()
{
  START_LOG("evaluate_stability_constant()", "RBSCMConstruction");

  // Get current index of C_J
  const unsigned int j = rb_scm_eval->C_J.size()-1;

  eigen_solver->set_position_of_spectrum(SMALLEST_REAL);

  // We assume B is set in system assembly
  // For coercive problems, B is set to the inner product matrix
  // For non-coercive time-dependent problems, B is set to the mass matrix

  // Set matrix A corresponding to mu_star
  matrix_A->zero();
  for(unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
  {
    add_scaled_symm_Aq(q, get_rb_theta_expansion().eval_A_theta(q,get_parameters()));
  }

  set_eigensolver_properties(-1);
  solve();
  unsigned int nconv = get_n_converged();
  if (nconv != 0)
  {
    std::pair<Real, Real> eval = get_eigenpair(0);

    // ensure that the eigenvalue is real
    libmesh_assert_less (eval.second, TOLERANCE);

    // Store the coercivity constant corresponding to mu_star
    rb_scm_eval->set_C_J_stability_constraint(j,eval.first);
    libMesh::out << std::endl << "Stability constant for C_J("<<j<<") = "
                 << rb_scm_eval->get_C_J_stability_constraint(j) << std::endl << std::endl;

    // Compute and store the vector y = (y_1, \ldots, y_Q) for the
    // eigenvector currently stored in eigen_system.solution.
    // We use this later to compute the SCM upper bounds.
    Real norm_B2 = libmesh_real( B_inner_product(*solution, *solution) );

    for(unsigned int q=0; q<get_rb_theta_expansion().get_n_A_terms(); q++)
    {
      Real norm_Aq2 = libmesh_real( Aq_inner_product(q, *solution, *solution) );

      rb_scm_eval->set_SCM_UB_vector(j,q,norm_Aq2/norm_B2);
    }
  }
  else
  {
    libMesh::err << "Error: Eigensolver did not converge in evaluate_stability_constant"
                 << std::endl;
    libmesh_error();
  }

  STOP_LOG("evaluate_stability_constant()", "RBSCMConstruction");
}

Number RBSCMConstruction::B_inner_product(const NumericVector<Number>& v, const NumericVector<Number>& w) const
{
  matrix_B->vector_mult(*inner_product_storage_vector, w);

  return v.dot(*inner_product_storage_vector);
}

Number RBSCMConstruction::Aq_inner_product(unsigned int q,
                                    const NumericVector<Number>& v,
                                    const NumericVector<Number>& w)
{
  if(q >= get_rb_theta_expansion().get_n_A_terms())
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

std::pair<unsigned int,Real> RBSCMConstruction::compute_SCM_bounds_on_training_set()
{
  START_LOG("compute_SCM_bounds_on_training_set()", "RBSCMConstruction");

  // Now compute the maximum bound error over training_parameters
  unsigned int new_C_J_index = 0;
  Real max_SCM_error = 0.;

  unsigned int first_index = get_first_local_training_index();
  for(unsigned int i=0; i<get_local_n_training_samples(); i++)
  {
    set_params_from_training_set(first_index+i);
    rb_scm_eval->set_parameters( get_parameters() );
    Real LB = rb_scm_eval->get_SCM_LB();
    Real UB = rb_scm_eval->get_SCM_UB();

    Real error_i = SCM_greedy_error_indicator(LB, UB);

    if( error_i > max_SCM_error )
    {
      max_SCM_error = error_i;
      new_C_J_index = i;
    }
  }

  unsigned int global_index = first_index + new_C_J_index;
  std::pair<unsigned int,Real> error_pair(global_index, max_SCM_error);
  get_global_max_error_pair(this->comm(),error_pair);

  STOP_LOG("compute_SCM_bounds_on_training_set()", "RBSCMConstruction");

  return error_pair;
}

void RBSCMConstruction::enrich_C_J(unsigned int new_C_J_index)
{
  START_LOG("enrich_C_J()", "RBSCMConstruction");

  set_params_from_training_set_and_broadcast(new_C_J_index);

  rb_scm_eval->C_J.push_back(get_parameters());

  libMesh::out << std::endl << "SCM: Added mu = (";

  RBParameters::const_iterator it     = get_parameters().begin();
  RBParameters::const_iterator it_end = get_parameters().end();
  for( ; it != it_end; ++it)
  {
    if(it != get_parameters().begin()) libMesh::out << ",";
    std::string param_name = it->first;
    RBParameters C_J_params = rb_scm_eval->C_J[rb_scm_eval->C_J.size()-1];
    libMesh::out << C_J_params.get_value(param_name);
  }
  libMesh::out << ")" << std::endl;

  // Finally, resize C_J_stability_vector and SCM_UB_vectors
  rb_scm_eval->C_J_stability_vector.push_back(0.);

  std::vector<Real> zero_vector(get_rb_theta_expansion().get_n_A_terms());
  rb_scm_eval->SCM_UB_vectors.push_back(zero_vector);

  STOP_LOG("enrich_C_J()", "RBSCMConstruction");
}


} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC && LIBMESH_HAVE_GLPK
