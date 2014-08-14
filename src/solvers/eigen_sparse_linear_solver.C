// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN


// C++ includes

// Local Includes
#include "libmesh/eigen_sparse_linear_solver.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

/*----------------------- functions ----------------------------------*/
template <typename T>
void EigenSparseLinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      this->_solver_type         = BICGSTAB;
      this->_preconditioner_type = ILU_PRECOND;
    }
}



template <typename T>
void EigenSparseLinearSolver<T>::init (const char* /*name*/)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;
    }
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (SparseMatrix<T> &matrix_in,
                                   NumericVector<T> &solution_in,
                                   NumericVector<T> &rhs_in,
                                   const double tol,
                                   const unsigned int m_its)
{
  START_LOG("solve()", "EigenSparseLinearSolver");
  this->init ();

  // Make sure the data passed in are really Eigen types
  EigenSparseMatrix<T>& matrix   = cast_ref<EigenSparseMatrix<T>&>(matrix_in);
  EigenSparseVector<T>& solution = cast_ref<EigenSparseVector<T>&>(solution_in);
  EigenSparseVector<T>& rhs      = cast_ref<EigenSparseVector<T>&>(rhs_in);

  // Close the matrix and vectors in case this wasn't already done.
  matrix.close();
  solution.close();
  rhs.close();

  std::pair<unsigned int, Real> retval(0,0.);

  // Solve the linear system
  switch (this->_solver_type)
    {
      // Conjugate-Gradient
    case CG:
      {
        Eigen::ConjugateGradient<EigenSM> solver (matrix._mat);
        solver.setMaxIterations(m_its);
        solver.setTolerance(tol);
        solution._vec = solver.solveWithGuess(rhs._vec,solution._vec);
        libMesh::out << "#iterations: " << solver.iterations() << std::endl;
        libMesh::out << "estimated error: " << solver.error() << std::endl;
        retval = std::make_pair(solver.iterations(), solver.error());
        break;
      }

      // Bi-Conjugate Gradient Stabilized
    case BICGSTAB:
      {
        Eigen::BiCGSTAB<EigenSM> solver (matrix._mat);
        solver.setMaxIterations(m_its);
        solver.setTolerance(tol);
        solution._vec = solver.solveWithGuess(rhs._vec,solution._vec);
        libMesh::out << "#iterations: " << solver.iterations() << std::endl;
        libMesh::out << "estimated error: " << solver.error() << std::endl;
        retval = std::make_pair(solver.iterations(), solver.error());
        break;
      }

      //   // Generalized Minimum Residual
      // case GMRES:
      //   {
      // libmesh_not_implemented();
      // break;
      //   }

      // Unknown solver, use BICGSTAB
    default:
      {
        libMesh::err << "ERROR:  Unsupported Eigen Solver: "
                     << Utility::enum_to_string(this->_solver_type) << std::endl
                     << "Continuing with BICGSTAB" << std::endl;

        this->_solver_type = BICGSTAB;

        STOP_LOG("solve()", "EigenSparseLinearSolver");

        return this->solve (matrix,
                            solution,
                            rhs,
                            tol,
                            m_its);
      }
    }

  STOP_LOG("solve()", "EigenSparseLinearSolver");
  return retval;
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::adjoint_solve (SparseMatrix<T> &matrix_in,
                                           NumericVector<T> &solution_in,
                                           NumericVector<T> &rhs_in,
                                           const double tol,
                                           const unsigned int m_its)
{

  START_LOG("adjoint_solve()", "EigenSparseLinearSolver");

  libmesh_experimental();
  EigenSparseMatrix<T> mat_trans(this->comm());
  matrix_in.get_transpose(mat_trans);

  std::pair<unsigned int, Real> retval = this->solve (mat_trans,
                                                      solution_in,
                                                      rhs_in,
                                                      tol,
                                                      m_its);

  STOP_LOG("adjoint_solve()", "EigenSparseLinearSolver");

  return retval;
}




template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (const ShellMatrix<T>& /*shell_matrix*/,
                                   NumericVector<T>& /*solution_in*/,
                                   NumericVector<T>& /*rhs_in*/,
                                   const double /*tol*/,
                                   const unsigned int /*m_its*/)
{
  libmesh_not_implemented();
  return std::make_pair(0,0.0);
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (const ShellMatrix<T>& /*shell_matrix*/,
                                   const SparseMatrix<T>& /*precond_matrix*/,
                                   NumericVector<T>& /*solution_in*/,
                                   NumericVector<T>& /*rhs_in*/,
                                   const double /*tol*/,
                                   const unsigned int /*m_its*/)
{
  libmesh_not_implemented();
  return std::make_pair(0,0.0);
}



template <typename T>
void EigenSparseLinearSolver<T>::set_eigen_preconditioner_type ()
{
  libmesh_not_implemented();

  // switch (this->_preconditioner_type)
  //   {
  //   case IDENTITY_PRECOND:
  //     _precond_type = NULL; return;

  //   case ILU_PRECOND:
  //     _precond_type = ILUPrecond; return;

  //   case JACOBI_PRECOND:
  //     _precond_type = JacobiPrecond; return;

  //   case SSOR_PRECOND:
  //     _precond_type = SSORPrecond; return;


  //   default:
  //     libMesh::err << "ERROR:  Unsupported LASPACK Preconditioner: "
  //     << this->_preconditioner_type << std::endl
  //     << "Continuing with ILU"      << std::endl;
  //     this->_preconditioner_type = ILU_PRECOND;
  //     this->set_laspack_preconditioner_type();
  //   }
}



template <typename T>
void EigenSparseLinearSolver<T>::print_converged_reason() const
{
  libMesh::out << "print_converged_reason() is currently only supported"
               << "with Petsc 2.3.1 and later." << std::endl;
}



template <typename T>
LinearConvergenceReason EigenSparseLinearSolver<T>::get_converged_reason() const
{
  libmesh_not_implemented();
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseLinearSolver<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
