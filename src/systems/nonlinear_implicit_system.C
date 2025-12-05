// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/diff_solver.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/static_condensation.h"
#include "libmesh/static_condensation_preconditioner.h"

namespace libMesh
{

NonlinearImplicitSystem::NonlinearImplicitSystem (EquationSystems & es,
                                                  const std::string & name_in,
                                                  const unsigned int number_in) :

  Parent                    (es, name_in, number_in),
  nonlinear_solver          (NonlinearSolver<Number>::build(*this)),
  diff_solver               (),
  _n_nonlinear_iterations   (0),
  _final_nonlinear_residual (1.e20)
{
  // Set default parameters
  // These were chosen to match the Petsc defaults
  es.parameters.set<Real>        ("linear solver tolerance") = 1e-5;
  es.parameters.set<Real>        ("linear solver minimum tolerance") = 1e-5;
  es.parameters.set<unsigned int>("linear solver maximum iterations") = 10000;

  es.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 50;
  es.parameters.set<unsigned int>("nonlinear solver maximum function evaluations") = 10000;

  es.parameters.set<Real>("nonlinear solver absolute residual tolerance") = 1e-35;
  es.parameters.set<Real>("nonlinear solver relative residual tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver divergence tolerance") = 1e+4;
  es.parameters.set<Real>("nonlinear solver absolute step tolerance") = 1e-8;
  es.parameters.set<Real>("nonlinear solver relative step tolerance") = 1e-8;

  es.parameters.set<bool>("reuse preconditioner") = false;
  es.parameters.set<unsigned int>("reuse preconditioner maximum linear iterations") = 1;

  if (this->has_static_condensation())
    this->setup_static_condensation_preconditioner(*nonlinear_solver);
}



NonlinearImplicitSystem::~NonlinearImplicitSystem () = default;



void NonlinearImplicitSystem::create_static_condensation()
{
  Parent::create_static_condensation();
  this->setup_static_condensation_preconditioner(*nonlinear_solver);
}



void NonlinearImplicitSystem::clear ()
{
  // clear the nonlinear solver
  nonlinear_solver->clear();

  // FIXME - this is necessary for petsc_auto_fieldsplit
  // nonlinear_solver->init_names(*this);

  // clear the parent data
  Parent::clear();
}



void NonlinearImplicitSystem::reinit ()
{
  // re-initialize the nonlinear solver interface
  nonlinear_solver->clear();

  // force the solver to get a new preconditioner, in
  // case reuse was set
  nonlinear_solver->force_new_preconditioner();

  // FIXME - this is necessary for petsc_auto_fieldsplit
  // nonlinear_solver->init_names(*this);

  if (diff_solver.get())
    diff_solver->reinit();

  // initialize parent data
  Parent::reinit();
}



void NonlinearImplicitSystem::set_solver_parameters ()
{
  // Get a reference to the EquationSystems
  const EquationSystems & es =
    this->get_equation_systems();

  // Get the user-specified nonlinear solver tolerances
  const unsigned int maxits = parameters.have_parameter<unsigned int>("nonlinear solver maximum iterations") ?
    parameters.get<unsigned int>("nonlinear solver maximum iterations") :
    es.parameters.get<unsigned int>("nonlinear solver maximum iterations");

  const unsigned int maxfuncs = parameters.have_parameter<unsigned int>("nonlinear solver maximum function evaluations") ?
    parameters.get<unsigned int>("nonlinear solver maximum function evaluations") :
    es.parameters.get<unsigned int>("nonlinear solver maximum function evaluations");

  const double abs_resid_tol = parameters.have_parameter<Real>("nonlinear solver absolute residual tolerance") ?
    double(parameters.get<Real>("nonlinear solver absolute residual tolerance")) :
    double(es.parameters.get<Real>("nonlinear solver absolute residual tolerance"));

  const double rel_resid_tol = parameters.have_parameter<Real>("nonlinear solver relative residual tolerance") ?
    double(parameters.get<Real>("nonlinear solver relative residual tolerance")) :
    double(es.parameters.get<Real>("nonlinear solver relative residual tolerance"));

  const double div_tol = parameters.have_parameter<Real>("nonlinear solver divergence tolerance") ?
    double(parameters.get<Real>("nonlinear solver divergence tolerance")) :
    double(es.parameters.get<Real>("nonlinear solver divergence tolerance"));

  const double abs_step_tol = parameters.have_parameter<Real>("nonlinear solver absolute step tolerance") ?
    double(parameters.get<Real>("nonlinear solver absolute step tolerance")) :
    double(es.parameters.get<Real>("nonlinear solver absolute step tolerance"));

  const double rel_step_tol = parameters.have_parameter<Real>("nonlinear solver relative step tolerance")?
    double(parameters.get<Real>("nonlinear solver relative step tolerance")) :
    double(es.parameters.get<Real>("nonlinear solver relative step tolerance"));

  // Get the user-specified linear solver tolerances
  const auto [maxlinearits, linear_tol] = this->Parent::get_linear_solve_parameters();

  const double linear_min_tol = parameters.have_parameter<Real>("linear solver minimum tolerance") ?
    double(parameters.get<Real>("linear solver minimum tolerance")) :
    double(es.parameters.get<Real>("linear solver minimum tolerance"));

  const bool reuse_preconditioner = parameters.have_parameter<unsigned int>("reuse preconditioner") ?
    parameters.get<unsigned int>("reuse preconditioner") :
      es.parameters.get<bool>("reuse preconditioner");
  const unsigned int reuse_preconditioner_max_linear_its =
    parameters.have_parameter<unsigned int>("reuse preconditioner maximum linear iterations") ?
    parameters.get<unsigned int>("reuse preconditioner maximum linear iterations") :
      es.parameters.get<unsigned int>("reuse preconditioner maximum linear iterations");

  // Set all the parameters on the NonlinearSolver
  nonlinear_solver->max_nonlinear_iterations = maxits;
  nonlinear_solver->max_function_evaluations = maxfuncs;
  nonlinear_solver->absolute_residual_tolerance = abs_resid_tol;
  nonlinear_solver->relative_residual_tolerance = rel_resid_tol;
  nonlinear_solver->divergence_tolerance = div_tol;
  nonlinear_solver->absolute_step_tolerance = abs_step_tol;
  nonlinear_solver->relative_step_tolerance = rel_step_tol;
  nonlinear_solver->max_linear_iterations = maxlinearits;
  nonlinear_solver->initial_linear_tolerance = linear_tol;
  nonlinear_solver->minimum_linear_tolerance = linear_min_tol;
  nonlinear_solver->set_reuse_preconditioner(reuse_preconditioner);
  nonlinear_solver->set_reuse_preconditioner_max_linear_its(reuse_preconditioner_max_linear_its);

  if (diff_solver.get())
    {
      diff_solver->max_nonlinear_iterations = maxits;
      diff_solver->absolute_residual_tolerance = abs_resid_tol;
      diff_solver->relative_residual_tolerance = rel_resid_tol;
      diff_solver->absolute_step_tolerance = abs_step_tol;
      diff_solver->relative_step_tolerance = rel_step_tol;
      diff_solver->max_linear_iterations = maxlinearits;
      diff_solver->initial_linear_tolerance = linear_tol;
      diff_solver->minimum_linear_tolerance = linear_min_tol;
    }
}



void NonlinearImplicitSystem::solve ()
{
  // Log how long the nonlinear solve takes.
  LOG_SCOPE("solve()", "System");

  this->set_solver_parameters();

  if (diff_solver.get())
    {
      diff_solver->solve();

      // Store the number of nonlinear iterations required to
      // solve and the final residual.
      _n_nonlinear_iterations   = diff_solver->total_outer_iterations();
      _final_nonlinear_residual = 0.; // FIXME - support this!
    }
  else
    {
      if (this->prefix_with_name())
        nonlinear_solver->init(this->prefix().c_str());
      else
        nonlinear_solver->init();

      // FIXME - this is necessary for petsc_auto_fieldsplit
      // nonlinear_solver->init_names(*this);

      // Solve the nonlinear system.
      // Store the number of nonlinear iterations required to
      // solve and the final residual.
      std::tie(_n_nonlinear_iterations, _final_nonlinear_residual) =
        nonlinear_solver->solve (*matrix, *solution, *rhs,
                                 nonlinear_solver->relative_residual_tolerance,
                                 nonlinear_solver->max_linear_iterations);
    }

  // Update the system after the solve
  this->update();
}



std::pair<unsigned int, Real> NonlinearImplicitSystem::get_linear_solve_parameters() const
{
  if (diff_solver.get())
    return std::make_pair(this->diff_solver->max_linear_iterations,
                          this->diff_solver->relative_residual_tolerance);
  return std::make_pair(this->nonlinear_solver->max_linear_iterations,
                        this->nonlinear_solver->relative_residual_tolerance);
}



void NonlinearImplicitSystem::assembly(bool get_residual,
                                       bool get_jacobian,
                                       bool /*apply_heterogeneous_constraints*/,
                                       bool /*apply_no_constraints*/)
{
  libmesh_assert(this->get_mesh().is_prepared());
#ifdef DEBUG
  libmesh_assert(MeshTools::valid_is_prepared(this->get_mesh()));
#endif

  // Get current_local_solution in sync
  this->update();

  //-----------------------------------------------------------------------------
  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  libmesh_error_msg_if(nonlinear_solver->jacobian && nonlinear_solver->jacobian_object,
                       "ERROR: cannot specify both a function and object to compute the Jacobian!");

  libmesh_error_msg_if(nonlinear_solver->residual && nonlinear_solver->residual_object,
                       "ERROR: cannot specify both a function and object to compute the Residual!");

  libmesh_error_msg_if(nonlinear_solver->matvec && nonlinear_solver->residual_and_jacobian_object,
                       "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!");


  if (get_jacobian)
    {
      if (nonlinear_solver->jacobian != nullptr)
        nonlinear_solver->jacobian (*current_local_solution.get(), *matrix, *this);

      else if (nonlinear_solver->jacobian_object != nullptr)
        nonlinear_solver->jacobian_object->jacobian (*current_local_solution.get(), *matrix, *this);

      else if (nonlinear_solver->matvec != nullptr)
        nonlinear_solver->matvec (*current_local_solution.get(), get_residual ? rhs : nullptr, matrix, *this);

      else if (nonlinear_solver->residual_and_jacobian_object != nullptr)
        nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian (*current_local_solution.get(), get_residual ? rhs : nullptr, matrix, *this);

      else
        libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");
    }

  if (get_residual)
    {
      if (nonlinear_solver->residual != nullptr)
        nonlinear_solver->residual (*current_local_solution.get(), *rhs, *this);

      else if (nonlinear_solver->residual_object != nullptr)
        nonlinear_solver->residual_object->residual (*current_local_solution.get(), *rhs, *this);

      else if (nonlinear_solver->matvec != nullptr)
        {
          // we might have already grabbed the residual and jacobian together
          if (!get_jacobian)
            nonlinear_solver->matvec (*current_local_solution.get(), rhs, nullptr, *this);
        }

      else if (nonlinear_solver->residual_and_jacobian_object != nullptr)
        {
          // we might have already grabbed the residual and jacobian together
          if (!get_jacobian)
            nonlinear_solver->residual_and_jacobian_object->residual_and_jacobian (*current_local_solution.get(), rhs, nullptr, *this);
        }

      else
        libmesh_error_msg("Error! Unable to compute residual and/or Jacobian!");
    }
  else
    libmesh_assert(get_jacobian);  // I can't believe you really wanted to assemble *nothing*
}




unsigned NonlinearImplicitSystem::get_current_nonlinear_iteration_number() const
{
  return nonlinear_solver->get_current_nonlinear_iteration_number();
}



} // namespace libMesh
