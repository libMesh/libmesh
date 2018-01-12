// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

#include "libmesh/diff_system.h"
#include "libmesh/eigen_time_solver.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

// Constructor
EigenTimeSolver::EigenTimeSolver (sys_type & s)
  : Parent(s),
    eigen_solver     (EigenSolver<Number>::build(s.comm())),
    tol(pow(TOLERANCE, 5./3.)),
    maxits(1000),
    n_eigenpairs_to_compute(5),
    n_basis_vectors_to_use(3*n_eigenpairs_to_compute),
    n_converged_eigenpairs(0),
    n_iterations_reqd(0)
{
  libmesh_experimental();
  eigen_solver->set_eigenproblem_type(GHEP);//or GNHEP
  eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
}

EigenTimeSolver::~EigenTimeSolver ()
{
}

void EigenTimeSolver::reinit ()
{
  // empty...
}

void EigenTimeSolver::init ()
{
  // Add matrix "B" to _system if not already there.
  // The user may have already added a matrix "B" before
  // calling the System initialization.  This would be
  // necessary if e.g. the System originally started life
  // with a different type of TimeSolver and only later
  // had its TimeSolver changed to an EigenTimeSolver.
  if (!_system.have_matrix("B"))
    _system.add_matrix("B");
}

void EigenTimeSolver::solve ()
{
  // The standard implementation is basically to call:
  // _diff_solver->solve();
  // which internally assembles (when necessary) matrices and vectors
  // and calls linear solver software while also doing Newton steps (see newton_solver.C)
  //
  // The element_residual and side_residual functions below control
  // what happens in the interior of the element assembly loops.
  // We have a system reference, so it's possible to call _system.assembly()
  // ourselves if we want to...
  //
  // Interestingly, for the EigenSolver we don't need residuals...just Jacobians.
  // The Jacobian should therefore always be requested, and always return
  // jacobian_computed as being true.

  // The basic plan of attack is:
  // .) Construct the Jacobian using _system.assembly(true,true) as we
  //    would for a steady system.  Use a flag in this class to
  //    control behavior in element_residual and side_residual
  // .) Swap _system.matrix to matrix "B" (be sure to add this extra matrix during init)
  // .) Call _system.assembly(true,true) again, using the flag in element_residual
  //    and side_residual to only get the mass matrix terms.
  // .) Send A and B to Steffen's EigenSolver interface.

  // Assemble the spatial part (matrix A) of the operator
  if (!this->quiet)
    libMesh::out << "Assembling matrix A." << std::endl;
  _system.matrix =   &( _system.get_matrix ("System Matrix") );
  this->now_assembling = Matrix_A;
  _system.assembly(true, true);
  //_system.matrix->print_matlab("matrix_A.m");

  // Point the system's matrix at B, call assembly again.
  if (!this->quiet)
    libMesh::out << "Assembling matrix B." << std::endl;
  _system.matrix =   &( _system.get_matrix ("B") );
  this->now_assembling = Matrix_B;
  _system.assembly(true, true);
  //_system.matrix->print_matlab("matrix_B.m");

  // Send matrices A, B to Steffen's SlepcEigenSolver interface
  //libmesh_here();
  if (!this->quiet)
    libMesh::out << "Calling the EigenSolver." << std::endl;
  std::pair<unsigned int, unsigned int> solve_data =
    eigen_solver->solve_generalized (_system.get_matrix ("System Matrix"),
                                     _system.get_matrix ("B"),
                                     n_eigenpairs_to_compute,
                                     n_basis_vectors_to_use,
                                     tol,
                                     maxits);

  this->n_converged_eigenpairs = solve_data.first;
  this->n_iterations_reqd      = solve_data.second;
}



bool EigenTimeSolver::element_residual(bool request_jacobian,
                                       DiffContext & context)
{
  // The EigenTimeSolver always computes jacobians!
  libmesh_assert (request_jacobian);

  context.elem_solution_rate_derivative = 1;
  context.elem_solution_derivative = 1;

  // Assemble the operator for the spatial part.
  if (now_assembling == Matrix_A)
    {
      bool jacobian_computed =
        _system.get_physics()->element_time_derivative(request_jacobian, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert(request_jacobian || !jacobian_computed);

      bool jacobian_computed2 =
        _system.get_physics()->element_constraint(jacobian_computed, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert (jacobian_computed || !jacobian_computed2);

      return jacobian_computed && jacobian_computed2;

    }

  // Assemble the mass matrix operator
  else if (now_assembling == Matrix_B)
    {
      bool mass_jacobian_computed =
        _system.get_physics()->mass_residual(request_jacobian, context);

      // Scale Jacobian by -1 to get positive matrix from new negative
      // mass_residual convention
      context.get_elem_jacobian() *= -1.0;

      return mass_jacobian_computed;
    }

  else
    libmesh_error_msg("Unrecognized value now_assembling = " << now_assembling);
}



bool EigenTimeSolver::side_residual(bool request_jacobian,
                                    DiffContext & context)
{
  // The EigenTimeSolver always requests jacobians?
  //libmesh_assert (request_jacobian);

  context.elem_solution_rate_derivative = 1;
  context.elem_solution_derivative = 1;

  // Assemble the operator for the spatial part.
  if (now_assembling == Matrix_A)
    {
      bool jacobian_computed =
        _system.get_physics()->side_time_derivative(request_jacobian, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert (request_jacobian || !jacobian_computed);

      bool jacobian_computed2 =
        _system.get_physics()->side_constraint(jacobian_computed, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert (jacobian_computed || !jacobian_computed2);

      return jacobian_computed && jacobian_computed2;

    }

  // There is now a "side" equivalent for the mass matrix
  else if (now_assembling == Matrix_B)
    {
      bool mass_jacobian_computed =
        _system.get_physics()->side_mass_residual(request_jacobian, context);

      // Scale Jacobian by -1 to get positive matrix from new negative
      // mass_residual convention
      context.get_elem_jacobian() *= -1.0;

      return mass_jacobian_computed;
    }

  else
    libmesh_error_msg("Unrecognized value now_assembling = " << now_assembling);
}



bool EigenTimeSolver::nonlocal_residual(bool request_jacobian,
                                        DiffContext & context)
{
  // The EigenTimeSolver always requests jacobians?
  //libmesh_assert (request_jacobian);

  // Assemble the operator for the spatial part.
  if (now_assembling == Matrix_A)
    {
      bool jacobian_computed =
        _system.get_physics()->nonlocal_time_derivative(request_jacobian, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert (request_jacobian || !jacobian_computed);

      bool jacobian_computed2 =
        _system.get_physics()->nonlocal_constraint(jacobian_computed, context);

      // The user shouldn't compute a jacobian unless requested
      libmesh_assert (jacobian_computed || !jacobian_computed2);

      return jacobian_computed && jacobian_computed2;

    }

  // There is now a "side" equivalent for the mass matrix
  else if (now_assembling == Matrix_B)
    {
      bool mass_jacobian_computed =
        _system.get_physics()->nonlocal_mass_residual(request_jacobian, context);

      // Scale Jacobian by -1 to get positive matrix from new negative
      // mass_residual convention
      context.get_elem_jacobian() *= -1.0;

      return mass_jacobian_computed;
    }

  else
    libmesh_error_msg("Unrecognized value now_assembling = " << now_assembling);
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC
