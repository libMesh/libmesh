// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "linear_implicit_system.h"
#include "linear_solver.h"
#include "equation_systems.h"
#include "libmesh_logging.h"
#include "numeric_vector.h" // for parameter sensitivity calcs
#include "sparse_matrix.h" // for get_transpose


// ------------------------------------------------------------
// LinearImplicitSystem implementation
LinearImplicitSystem::LinearImplicitSystem (EquationSystems& es,
					    const std::string& name,
					    const unsigned int number) :
  
  Parent                 (es, name, number),
  linear_solver          (LinearSolver<Number>::build()),
  _n_linear_iterations   (0),
  _final_linear_residual (1.e20),
  _shell_matrix(NULL)
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
  
  // clear the parent data
  Parent::clear();
}



void LinearImplicitSystem::reinit ()
{
  // re-initialize the linear solver interface
  linear_solver->clear();
  
  // initialize parent data
  Parent::reinit();  
}



void LinearImplicitSystem::solve ()
{
  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assemble (); 

  // Log how long the linear solve takes.
  // This gets done by the LinearSolver classes now [RHS]
  // START_LOG("solve()", "System");

  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
  if(_shell_matrix)
    // 1.) Shell matrix with or without user-supplied preconditioner.
    rval = linear_solver->solve(*_shell_matrix, this->request_matrix("Preconditioner"), *solution, *rhs, tol, maxits);
  else
    // 2.) No shell matrix, with or without user-supplied preconditioner
    rval = linear_solver->solve (*matrix, this->request_matrix("Preconditioner"), *solution, *rhs, tol, maxits);

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  // This gets done by the LinearSolver classes now [RHS]
  // STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



void LinearImplicitSystem::adjoint_solve ()
{
  // The adjoint_solve API (and all APIs using it) are about to see a
  // series of non-backwards-compatible changes, primarily to add
  // multiple-QoI support
  libmesh_experimental();

  // We currently don't support adjoint solves of shell matrices
  // FIXME - we should let shell matrices support
  // vector_transpose_mult so that we can use them here.
  if(_shell_matrix!=NULL)
    libmesh_not_implemented();

  // Adding an adjoint_solution vector, allocate an adjoint_solution if it doesn't already exist 

  NumericVector<Number> & adjoint_solution = System::add_vector("adjoint_solution");     

  if (this->assemble_before_solve)
    {
      // Assemble the linear system
      this->assemble (); 

      // And take the adjoint
      matrix->get_transpose(*matrix);

      // Including of any separate preconditioner
      if(this->have_matrix("Preconditioner"))
        {
          SparseMatrix<Number> &pre = this->get_matrix("Preconditioner");
          pre.get_transpose(pre);
        }

      // But now replace the right hand side with the quantity of
      // interest functional's derivative
      this->assemble_qoi_derivative();
    }

  // Get a reference to the EquationSystems
  const EquationSystems& es =
    this->get_equation_systems();
  
  // Get the user-specifiied linear solver tolerance
  const Real tol            =
    es.parameters.get<Real>("adjoint solver tolerance");

  // Get the user-specified maximum # of linear solver iterations
  const unsigned int maxits =
    es.parameters.get<unsigned int>("adjoint solver maximum iterations");

  // Solve the linear system.  Several cases:
  std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
  if(this->have_matrix("Preconditioner"))
    {
      // 3.) No shell matrix, but with user-supplied preconditioner
      rval = linear_solver->solve (*matrix, this->get_matrix("Preconditioner"), adjoint_solution, *rhs, tol, maxits);
    }
  else
    {
      // 4.) No shell matrix, and use system matrix for the preconditioner
      rval = linear_solver->solve (*matrix, adjoint_solution, *rhs, tol, maxits);
    }

  // Store the number of linear iterations required to
  // solve and the final residual.
  _n_linear_iterations   = rval.first;
  _final_linear_residual = rval.second;
    
  // Stop logging the linear solve
  // This gets done by the LinearSolver classes now [RHS]
  // STOP_LOG("solve()", "System");

  // Update the system after the solve
  this->update();  
}



void LinearImplicitSystem::qoi_parameter_sensitivity
  (std::vector<Number *>& parameters,
   std::vector<Number>&   sensitivities)
{
  // Get ready to fill in senstivities:
  sensitivities.clear();
  sensitivities.resize(parameters.size(), 0);

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

  // We first do an adjoint solve:
  // A^T * z = (partial q / partial u)

  this->adjoint_solve();

  // We use the identities:
  // dq/dp = (partial q / partial p) + (partial q / partial u) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + (A^T * z) *
  //         (partial u / partial p)
  // dq/dp = (partial q / partial p) + z * A *
  //         (partial u / partial p)
 
  // Leading to our final formula:
  // dq/dp = (partial q / partial p) - z * (partial R / partial p)

  for (unsigned int i=0; i != parameters.size(); ++i)
    {
      // We currently get partial derivatives via central differencing
      Number delta_p = 1e-6;

      // (partial q / partial p) ~= (q(p+dp)-q(p-dp))/(2*dp)

      Number old_parameter = *parameters[i];
      // Number old_qoi = this->qoi;

      *parameters[i] = old_parameter - delta_p;
      this->assemble_qoi();
      Number qoi_minus = this->qoi;

      *parameters[i] = old_parameter + delta_p;
      this->assemble_qoi();
      Number qoi_plus = this->qoi;
      Number partialq_partialp = (qoi_plus - qoi_minus) / (2.*delta_p);

      // (partial R / partial p) = (partial b / partial p) - 
      //                           (partial A / partial p) * u
      //   ~= (b(p+dp)-A(p+dp)*u-b(p-dp)+A(p-dp)*u)/(2*dp)

      // We're still at p+delta_p, so start with b(p+dp)
      this->assemble();
      this->rhs->close();
      this->matrix->close();
      AutoPtr<NumericVector<Number> > partialR_partialp = this->rhs->clone();

      // PETSc doesn't implement SGEMX, so neither does NumericVector,
      // so here's a hackish workaround for v-=A*u:
      *partialR_partialp *= -1;
      partialR_partialp->add_vector(*this->solution, *this->matrix);
      *partialR_partialp *= -1;
      
      *parameters[i] = old_parameter - delta_p;
      this->assemble();
      this->rhs->close();
      this->matrix->close();
      *partialR_partialp -= *this->rhs;
      partialR_partialp->add_vector(*this->solution, *this->matrix);

      *partialR_partialp /= (2.*delta_p);

      // Don't leave the parameter changed
      *parameters[i] = old_parameter;

      sensitivities[i] = partialq_partialp -
			 partialR_partialp->dot(this->get_adjoint_solution());
    }

  // All parameters have been reset.
  // Don't leave the qoi or system changed - principle of least
  // surprise.
  this->assemble();
  this->assemble_qoi();
  this->rhs->close();
  this->matrix->close();
}



void LinearImplicitSystem::attach_shell_matrix (ShellMatrix<Number>* shell_matrix)
{
  _shell_matrix = shell_matrix;
}



