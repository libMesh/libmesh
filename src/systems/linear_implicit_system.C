// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/linear_solver.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h" // for parameter sensitivity calcs
//#include "libmesh/parameter_vector.h"
#include "libmesh/sparse_matrix.h" // for get_transpose
#include "libmesh/system_subset.h"

namespace libMesh
{


// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (EquationSystems & es,
                                            const std::string & name_in,
                                            const unsigned int number_in) :

  Parent                 (es, name_in, number_in),
  linear_solver          (LinearSolver<Number>::build(es.comm())),
  _n_linear_iterations   (0),
  _final_linear_residual (1.e20),
  _shell_matrix(libmesh_nullptr),
  _subset(libmesh_nullptr),
  _subset_solve_mode(SUBSET_ZERO)
{
}



LinearImplicitSystem::~LinearImplicitSystem ()
{
  // Clear data
  this->clear();
}



void LinearImplicitSystem::clear ()
{
  // clear the linear solver
  linear_solver->clear();

  this->restrict_solve_to(libmesh_nullptr);

  // clear the parent data
  Parent::clear();
}



void LinearImplicitSystem::init_data ()
{
  // initialize parent data
  Parent::init_data();

  // re-initialize the linear solver interface
  linear_solver->clear();
}



void LinearImplicitSystem::reinit ()
{
  // re-initialize the linear solver interface
  linear_solver->clear();

  // initialize parent data
  Parent::reinit();
}



void LinearImplicitSystem::restrict_solve_to (const SystemSubset * subset,
                                              const SubsetSolveMode subset_solve_mode)
{
  _subset = subset;
  _subset_solve_mode = subset_solve_mode;

  if (subset != libmesh_nullptr)
    libmesh_assert_equal_to (&subset->get_system(), this);
}



void LinearImplicitSystem::solve ()
{
  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assemble ();

  // Get a reference to the EquationSystems
  const EquationSystems & es =
    this->get_equation_systems();

  // If the linear solver hasn't been initialized, we do so here.
  if (libMesh::on_command_line("--solver_system_names"))
    linear_solver->init((this->name()+"_").c_str());
  else
    linear_solver->init();

  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  if (_subset != libmesh_nullptr)
    linear_solver->restrict_solve_to(&_subset->dof_ids(),_subset_solve_mode);

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
  if(_shell_matrix)
    // 1.) Shell matrix with or without user-supplied preconditioner.
    rval = linear_solver->solve(*_shell_matrix, this->request_matrix("Preconditioner"), *solution, *rhs, tol, maxits);
  else
    // 2.) No shell matrix, with or without user-supplied preconditioner
    rval = linear_solver->solve (*matrix, this->request_matrix("Preconditioner"), *solution, *rhs, tol, maxits);

  if (_subset != libmesh_nullptr)
    linear_solver->restrict_solve_to(libmesh_nullptr);

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;

  // Update the system after the solve
  this->update();
}



void LinearImplicitSystem::attach_shell_matrix (ShellMatrix<Number> * shell_matrix)
{
  _shell_matrix = shell_matrix;
}


/*
  void LinearImplicitSystem::sensitivity_solve (const ParameterVector & parameters)
  {
  if (this->assemble_before_solve)
  {
  // Assemble the linear system
  this->assemble ();

  // But now assemble right hand sides with the residual's
  // parameter derivatives
  this->assemble_residual_derivatives(parameters);
  }

  // Get a reference to the EquationSystems
  const EquationSystems & es =
  this->get_equation_systems();

  // Get the user-specifiied linear solver tolerance
  const Real tol            =
  es.parameters.get<Real>("sensitivity solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
  es.parameters.get<unsigned int>("sensitivity solver maximum iterations");

  // Our iteration counts and residuals will be sums of the individual
  // results
  _n_linear_iterations = 0;
  _final_linear_residual = 0.0;
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);

  // Solve the linear system.
  SparseMatrix<Number> * pc = this->request_matrix("Preconditioner");
  for (std::size_t p=0; p != parameters.size(); ++p)
  {
  rval = linear_solver->solve (*matrix, pc,
  this->get_sensitivity_solution(p),
  this->get_sensitivity_rhs(p), tol, maxits);
  _n_linear_iterations   += rval.first;
  _final_linear_residual += rval.second;
  }

  // Our matrix is the *negative* of the Jacobian for b-A*u, so our
  // solutions are all inverted
  for (std::size_t p=0; p != parameters.size(); ++p)
  {
  this->get_sensitivity_solution(p) *= -1.0;
  }
  }



  void LinearImplicitSystem::adjoint_solve (const QoISet & qoi_indices)
  {
  const unsigned int Nq = this->qoi.size();

  // We currently don't support adjoint solves of shell matrices
  // FIXME - we should let shell matrices support
  // vector_transpose_mult so that we can use them here.
  if(_shell_matrix!=libmesh_nullptr)
  libmesh_not_implemented();

  if (this->assemble_before_solve)
  {
  // Assemble the linear system
  this->assemble ();

  // And take the adjoint
  matrix->get_transpose(*matrix);

  // Including of any separate preconditioner
  SparseMatrix<Number> * pc = this->request_matrix("Preconditioner");
  if(pc)
  pc->get_transpose(*pc);

  // But now replace the right hand sides with the quantity of
  // interest functional's derivative
  this->assemble_qoi_derivative(qoi_indices);
  }

  // Get a reference to the EquationSystems
  const EquationSystems & es =
  this->get_equation_systems();

  // Get the user-specifiied linear solver tolerance
  const Real tol            =
  es.parameters.get<Real>("adjoint solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
  es.parameters.get<unsigned int>("adjoint solver maximum iterations");

  // Our iteration counts and residuals will be sums of the individual
  // results
  _n_linear_iterations = 0;
  _final_linear_residual = 0.0;
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);

  // Solve the linear system.
  SparseMatrix<Number> * pc = this->request_matrix("Preconditioner");
  for (unsigned int i=0; i != Nq; ++i)
  if (qoi_indices.has_index(i))
  {
  rval = linear_solver->solve (*matrix, pc,
  this->add_adjoint_solution(i),
  this->get_adjoint_rhs(i), tol, maxits);
  _n_linear_iterations   += rval.first;
  _final_linear_residual += rval.second;
  }
  }



  void LinearImplicitSystem::forward_qoi_parameter_sensitivity
  (const QoISet &          qoi_indices,
  const ParameterVector & parameters,
  SensitivityData &       sensitivities)
  {
  const unsigned int Np = parameters.size();
  const unsigned int Nq = this->qoi.size();

  // An introduction to the problem:
  //
  // A(p)*u(p) = b(p), where x is determined implicitly.
  // Residual R(u(p),p) := b(p) - A(p)*u(p)
  // partial R / partial u = -A
  //
  // This implies that:
  // d/dp(R) = 0
  // (partial b / partial p) -
  // (partial A / partial p) * u -
  // A * (partial u / partial p) = 0
  // A * (partial u / partial p) = (partial R / partial p)
  //   = (partial b / partial p) - (partial A / partial p) * u

  // We first solve for (partial u / partial p) for each parameter:
  // -A * (partial u / partial p) = - (partial R / partial p)

  this->sensitivity_solve(parameters);

  // Get ready to fill in senstivities:
  sensitivities.allocate_data(qoi_indices, *this, parameters);

  // We use the identity:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)

  // We get (partial q / partial u) from the user
  this->assemble_qoi_derivative(qoi_indices);

  for (unsigned int j=0; j != Np; ++j)
  {
  // We currently get partial derivatives via central differencing
  Number delta_p = 1e-6;

  // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)

  Number old_parameter = *parameters[j];

  *parameters[j] = old_parameter - delta_p;
  this->assemble_qoi(qoi_indices);
  std::vector<Number> qoi_minus = this->qoi;

  *parameters[j] = old_parameter + delta_p;
  this->assemble_qoi(qoi_indices);
  std::vector<Number> & qoi_plus = this->qoi;
  std::vector<Number> partialq_partialp(Nq, 0);
  for (unsigned int i=0; i != Nq; ++i)
  if (qoi_indices.has_index(i))
  partialq_partialp[i] = (qoi_plus[i] - qoi_minus[i]) / (2.*delta_p);

  for (unsigned int i=0; i != Nq; ++i)
  if (qoi_indices.has_index(i))
  sensitivities[i][j] = partialq_partialp[i] +
  this->get_adjoint_rhs(i).dot(this->get_sensitivity_solution(i));
  }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assemble();
  this->rhs->close();
  this->matrix->close();
  this->assemble_qoi(qoi_indices);
  }
*/



LinearSolver<Number> * LinearImplicitSystem::get_linear_solver() const
{
  return linear_solver.get();
}



void LinearImplicitSystem::release_linear_solver(LinearSolver<Number> *) const
{
}



void LinearImplicitSystem::assembly(bool,
                                    bool,
                                    bool)
{
  // Residual R(u(p),p) := A(p)*u(p) - b(p)
  // partial R / partial u = A

  this->assemble();
  this->rhs->close();
  this->matrix->close();

  *(this->rhs) *= -1.0;
  this->rhs->add_vector(*this->solution, *this->matrix);
}

} // namespace libMesh
