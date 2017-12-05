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


#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// C++ includes

// Local includes
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_base.h"

namespace libMesh
{


// ------------------------------------------------------------
// EigenSystem implementation
EigenSystem::EigenSystem (EquationSystems & es,
                          const std::string & name_in,
                          const unsigned int number_in
                          ) :
  Parent           (es, name_in, number_in),
  eigen_solver     (EigenSolver<Number>::build(es.comm())),
  _n_converged_eigenpairs (0),
  _n_iterations           (0),
  _is_generalized_eigenproblem (false),
  _eigen_problem_type (NHEP)
{
}



EigenSystem::~EigenSystem ()
{
  // clear data
  this->clear();
}



void EigenSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // Clean up the matrices
  matrix_A.reset();
  matrix_B.reset();

  // clear the solver
  eigen_solver->clear();
}


void EigenSystem::set_eigenproblem_type (EigenProblemType ept)
{
  _eigen_problem_type = ept;

  eigen_solver->set_eigenproblem_type(ept);

  // libMesh::out << "The Problem type is set to be: " << std::endl;

  switch (_eigen_problem_type)
    {
    case HEP:
      // libMesh::out << "Hermitian" << std::endl;
      break;

    case NHEP:
      // libMesh::out << "Non-Hermitian" << std::endl;
      break;

    case GHEP:
      // libMesh::out << "Generalized Hermitian" << std::endl;
      break;

    case GNHEP:
      // libMesh::out << "Generalized Non-Hermitian" << std::endl;
      break;

    case GHIEP:
      // libMesh::out << "Generalized indefinite Hermitian" << std::endl;
      break;

    default:
      // libMesh::out << "not properly specified" << std::endl;
      libmesh_error_msg("Unrecognized _eigen_problem_type = " << _eigen_problem_type);
    }
}




void EigenSystem::init_data ()
{
  // initialize parent data
  Parent::init_data();

  // define the type of eigenproblem
  if (_eigen_problem_type == GNHEP ||
      _eigen_problem_type == GHEP  ||
      _eigen_problem_type == GHIEP)
    _is_generalized_eigenproblem = true;

  // build the system matrix
  matrix_A = SparseMatrix<Number>::build(this->comm());

  this->init_matrices();
}



void EigenSystem::init_matrices ()
{
  DofMap & dof_map = this->get_dof_map();

  dof_map.attach_matrix(*matrix_A);

  // build matrix_B only in case of a
  // generalized problem
  if (_is_generalized_eigenproblem)
    {
      matrix_B = SparseMatrix<Number>::build(this->comm());
      dof_map.attach_matrix(*matrix_B);
    }

  dof_map.compute_sparsity(this->get_mesh());

  // initialize and zero system matrix
  matrix_A->init();
  matrix_A->zero();

  // eventually initialize and zero system matrix_B
  if (_is_generalized_eigenproblem)
    {
      matrix_B->init();
      matrix_B->zero();
    }
}



void EigenSystem::reinit ()
{
  // initialize parent data
  Parent::reinit();

  // Clear the matrices
  matrix_A->clear();

  if (_is_generalized_eigenproblem)
    matrix_B->clear();

  DofMap & dof_map = this->get_dof_map();

  // Clear the sparsity pattern
  dof_map.clear_sparsity();

  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // both matrices, \p DofMap now knows them
  dof_map.compute_sparsity(this->get_mesh());

  matrix_A->init();
  matrix_A->zero();

  if (_is_generalized_eigenproblem)
    {
      matrix_B->init();
      matrix_B->zero();
    }
}



void EigenSystem::solve ()
{

  // A reference to the EquationSystems
  EquationSystems & es = this->get_equation_systems();

  // check that necessary parameters have been set
  libmesh_assert (es.parameters.have_parameter<unsigned int>("eigenpairs"));
  libmesh_assert (es.parameters.have_parameter<unsigned int>("basis vectors"));

  if (this->assemble_before_solve)
    // Assemble the linear system
    this->assemble ();

  // Get the tolerance for the solver and the maximum
  // number of iterations. Here, we simply adopt the linear solver
  // specific parameters.
  const Real tol            =
    es.parameters.get<Real>("linear solver tolerance");

  const unsigned int maxits =
    es.parameters.get<unsigned int>("linear solver maximum iterations");

  const unsigned int nev    =
    es.parameters.get<unsigned int>("eigenpairs");

  const unsigned int ncv    =
    es.parameters.get<unsigned int>("basis vectors");

  std::pair<unsigned int, unsigned int> solve_data;

  // call the solver depending on the type of eigenproblem

  // Generalized eigenproblem
  if (_is_generalized_eigenproblem)
    solve_data = eigen_solver->solve_generalized (*matrix_A, *matrix_B, nev, ncv, tol, maxits);

  // Standard eigenproblem
  else
    {
      libmesh_assert (!matrix_B);
      solve_data = eigen_solver->solve_standard (*matrix_A, nev, ncv, tol, maxits);
    }

  this->_n_converged_eigenpairs = solve_data.first;
  this->_n_iterations           = solve_data.second;
}


void EigenSystem::assemble ()
{

  // Assemble the linear system
  Parent::assemble ();

}


std::pair<Real, Real> EigenSystem::get_eigenpair (dof_id_type i)
{
  // call the eigen_solver get_eigenpair method
  return eigen_solver->get_eigenpair (i, *solution);
}

} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC
