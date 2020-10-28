// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_base.h"
#include "libmesh/enum_eigen_solver_type.h"

namespace libMesh
{


EigenSystem::EigenSystem (EquationSystems & es,
                          const std::string & name_in,
                          const unsigned int number_in
                          ) :
  Parent           (es, name_in, number_in),
  matrix_A         (nullptr),
  matrix_B         (nullptr),
  precond_matrix   (nullptr),
  eigen_solver     (EigenSolver<Number>::build(es.comm())),
  _n_converged_eigenpairs (0),
  _n_iterations           (0),
  _eigen_problem_type (NHEP),
  _use_shell_matrices (false),
  _use_shell_precond_matrix(false)
{
}



EigenSystem::~EigenSystem () = default;



void EigenSystem::clear ()
{
  // Clear the parent data
  Parent::clear();

  // The SparseMatrices are contained in _matrices
  // and are cleared with the above call
  matrix_A = nullptr;
  matrix_B = nullptr;
  precond_matrix = nullptr;

  // Operators
  shell_matrix_A.reset();
  shell_matrix_B.reset();

  // Preconditioning matrices
  shell_precond_matrix.reset();

  // clear the solver
  eigen_solver->clear();
}


void EigenSystem::set_eigenproblem_type (EigenProblemType ept)
{
  if (!can_add_matrices())
    libmesh_error_msg("ERROR: Cannot change eigen problem type after system initialization");

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



void EigenSystem::add_matrices ()
{
  Parent::add_matrices();

  if (_use_shell_matrices)
  {
    if (!shell_matrix_A)
      shell_matrix_A = ShellMatrix<Number>::build(this->comm());

    if (generalized() && !shell_matrix_B)
      shell_matrix_B = ShellMatrix<Number>::build(this->comm());

    if (_use_shell_precond_matrix)
      {
        if (!shell_precond_matrix)
          shell_precond_matrix = ShellMatrix<Number>::build(this->comm());
      }
    else if (!precond_matrix)
      precond_matrix = &(this->add_matrix("Eigen Preconditioner"));
  }
  else
  {
    if (!matrix_A)
      matrix_A = &(this->add_matrix("Eigen Matrix A"));

    if (generalized() && !matrix_B)
      matrix_B = &(this->add_matrix("Eigen Matrix B"));
  }
}



void EigenSystem::init_matrices ()
{
  Parent::init_matrices();

  if (shell_matrix_A)
    {
      shell_matrix_A->attach_dof_map(this->get_dof_map());
      shell_matrix_A->init();
    }

  if (shell_matrix_B)
    {
      shell_matrix_B->attach_dof_map(this->get_dof_map());
      shell_matrix_B->init();
    }

  if (shell_precond_matrix)
    {
      shell_precond_matrix->attach_dof_map(this->get_dof_map());
      shell_precond_matrix->init();
    }
}


void EigenSystem::reinit ()
{
  // initialize parent data
  // this calls reinit on matrix_A, matrix_B, and precond_matrix (if any)
  Parent::reinit();

  if (shell_matrix_A)
    {
      shell_matrix_A->clear();
      shell_matrix_A->init();
    }

  if (shell_matrix_B)
    {
      shell_matrix_B->clear();
      shell_matrix_B->init();
    }

  if (shell_precond_matrix)
    {
      shell_precond_matrix->clear();
      shell_precond_matrix->init();
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

  if (_use_shell_matrices)
  {
    // Generalized eigenproblem
    if (generalized())
      // Shell preconditioning matrix
      if (_use_shell_precond_matrix)
        solve_data = eigen_solver->solve_generalized (*shell_matrix_A, *shell_matrix_B,*shell_precond_matrix, nev, ncv, tol, maxits);
      else
        solve_data = eigen_solver->solve_generalized (*shell_matrix_A, *shell_matrix_B,*precond_matrix, nev, ncv, tol, maxits);

    // Standard eigenproblem
    else
      {
        libmesh_assert (!shell_matrix_B);
        // Shell preconditioning matrix
        if (_use_shell_precond_matrix)
          solve_data = eigen_solver->solve_standard (*shell_matrix_A,*shell_precond_matrix, nev, ncv, tol, maxits);
        else
          solve_data = eigen_solver->solve_standard (*shell_matrix_A,*precond_matrix, nev, ncv, tol, maxits);
      }
  }
  else
  {
    // Generalized eigenproblem
    if (generalized())
      solve_data = eigen_solver->solve_generalized (*matrix_A, *matrix_B, nev, ncv, tol, maxits);

    // Standard eigenproblem
    else
      {
        libmesh_assert (!matrix_B);
        solve_data = eigen_solver->solve_standard (*matrix_A, nev, ncv, tol, maxits);
      }
  }

  this->_n_converged_eigenpairs = solve_data.first;
  this->_n_iterations           = solve_data.second;
}

std::pair<Real, Real> EigenSystem::get_eigenpair (dof_id_type i)
{
  // call the eigen_solver get_eigenpair method
  return eigen_solver->get_eigenpair (i, *solution);
}

std::pair<Real, Real> EigenSystem::get_eigenvalue (dof_id_type i)
{
  return eigen_solver->get_eigenvalue (i);
}

void EigenSystem::set_initial_space (NumericVector<Number> & initial_space_in)
{
  eigen_solver->set_initial_space (initial_space_in);
}


bool EigenSystem::generalized () const
{
  return _eigen_problem_type == GNHEP ||
         _eigen_problem_type == GHEP ||
         _eigen_problem_type == GHIEP;
}


const SparseMatrix<Number> & EigenSystem::get_matrix_A() const
{
  libmesh_assert(matrix_A);
  libmesh_assert(!_use_shell_matrices);
  libmesh_assert(has_matrix_A());
  return *matrix_A;
}


SparseMatrix<Number> & EigenSystem::get_matrix_A()
{
  libmesh_assert(matrix_A);
  libmesh_assert(!_use_shell_matrices);
  libmesh_assert(has_matrix_A());
  return *matrix_A;
}


const SparseMatrix<Number> & EigenSystem::get_matrix_B() const
{
  libmesh_assert(matrix_B);
  libmesh_assert(!_use_shell_matrices);
  libmesh_assert(generalized());
  libmesh_assert(has_matrix_B());
  return *matrix_B;
}


SparseMatrix<Number> & EigenSystem::get_matrix_B()
{
  libmesh_assert(matrix_B);
  libmesh_assert(!_use_shell_matrices);
  libmesh_assert(generalized());
  libmesh_assert(has_matrix_B());
  return *matrix_B;
}


const SparseMatrix<Number> & EigenSystem::get_precond_matrix() const
{
  libmesh_assert(precond_matrix);
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(!_use_shell_precond_matrix);
  libmesh_assert(has_precond_matrix());
  return *precond_matrix;
}


SparseMatrix<Number> & EigenSystem::get_precond_matrix()
{
  libmesh_assert(precond_matrix);
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(!_use_shell_precond_matrix);
  libmesh_assert(has_precond_matrix());
  return *precond_matrix;
}


const ShellMatrix<Number> & EigenSystem::get_shell_matrix_A() const
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(has_shell_matrix_A());
  return *shell_matrix_A;
}


ShellMatrix<Number> & EigenSystem::get_shell_matrix_A()
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(has_shell_matrix_A());
  return *shell_matrix_A;
}


const ShellMatrix<Number> & EigenSystem::get_shell_matrix_B() const
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(generalized());
  libmesh_assert(has_shell_matrix_B());
  return *shell_matrix_B;
}


ShellMatrix<Number> & EigenSystem::get_shell_matrix_B()
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(generalized());
  libmesh_assert(has_shell_matrix_B());
  return *shell_matrix_B;
}


const ShellMatrix<Number> & EigenSystem::get_shell_precond_matrix() const
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(_use_shell_precond_matrix);
  libmesh_assert(has_shell_precond_matrix());
  return *shell_precond_matrix;
}


ShellMatrix<Number> & EigenSystem::get_shell_precond_matrix()
{
  libmesh_assert(_use_shell_matrices);
  libmesh_assert(_use_shell_precond_matrix);
  libmesh_assert(has_shell_precond_matrix());
  return *shell_precond_matrix;
}


bool EigenSystem::has_matrix_A() const
{
  libmesh_assert_equal_to(request_matrix("Eigen Matrix A"), matrix_A);
  return matrix_A;
}


bool EigenSystem::has_matrix_B() const
{
  libmesh_assert_equal_to(request_matrix("Eigen Matrix B"), matrix_B);
  return matrix_B;
}


bool EigenSystem::has_precond_matrix() const
{
  libmesh_assert_equal_to(request_matrix("Eigen Preconditioner"), precond_matrix);
  return precond_matrix;
}


bool EigenSystem::has_shell_matrix_A() const
{
  return shell_matrix_A.get();
}


bool EigenSystem::has_shell_matrix_B() const
{
  return shell_matrix_B.get();
}


bool EigenSystem::has_shell_precond_matrix() const
{
  return shell_precond_matrix.get();
}


const EigenSolver<Number> & EigenSystem::get_eigen_solver() const
{
  libmesh_assert(eigen_solver);
  return *eigen_solver;
}


EigenSolver<Number> & EigenSystem::get_eigen_solver()
{
  libmesh_assert(eigen_solver);
  return *eigen_solver;
}


} // namespace libMesh

#endif // LIBMESH_HAVE_SLEPC
