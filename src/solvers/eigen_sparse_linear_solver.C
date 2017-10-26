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



#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_EIGEN


// C++ includes

// Local Includes
#include "libmesh/eigen_sparse_linear_solver.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/solver_configuration.h"

// GMRES is an "unsupported" iterative solver in Eigen.
#include <unsupported/Eigen/IterativeSolvers>

namespace libMesh
{

template <typename T>
EigenSparseLinearSolver<T>::
EigenSparseLinearSolver(const Parallel::Communicator & comm_in) :
  LinearSolver<T>(comm_in),
  _comp_info(Eigen::Success)
{
  // The GMRES _solver_type can be used in EigenSparseLinearSolver,
  // however, the GMRES iterative solver is currently in the Eigen
  // "unsupported" directory, so we use BICGSTAB as our default.
  this->_solver_type = BICGSTAB;
}



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
void EigenSparseLinearSolver<T>::init (const char * /*name*/)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;
    }
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (SparseMatrix<T> & matrix_in,
                                   NumericVector<T> & solution_in,
                                   NumericVector<T> & rhs_in,
                                   const double tol,
                                   const unsigned int m_its)
{
  LOG_SCOPE("solve()", "EigenSparseLinearSolver");
  this->init ();

  // Make sure the data passed in are really Eigen types
  EigenSparseMatrix<T> & matrix   = cast_ref<EigenSparseMatrix<T> &>(matrix_in);
  EigenSparseVector<T> & solution = cast_ref<EigenSparseVector<T> &>(solution_in);
  EigenSparseVector<T> & rhs      = cast_ref<EigenSparseVector<T> &>(rhs_in);

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
        _comp_info = solver.info();
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
        _comp_info = solver.info();
        break;
      }

      // Generalized Minimum Residual
    case GMRES:
      {
        Eigen::GMRES<EigenSM> solver (matrix._mat);
        solver.setMaxIterations(m_its);
        solver.setTolerance(tol);

        // If there is an int parameter called "gmres_restart" in the
        // SolverConfiguration object, pass it to the Eigen GMRES
        // solver.
        if (this->_solver_configuration)
          {
            std::map<std::string, int>::iterator it =
              this->_solver_configuration->int_valued_data.find("gmres_restart");

            if (it != this->_solver_configuration->int_valued_data.end())
              solver.set_restart(it->second);
          }

        libMesh::out << "Eigen GMRES solver, restart = " << solver.get_restart() << std::endl;
        solution._vec = solver.solveWithGuess(rhs._vec, solution._vec);
        libMesh::out << "#iterations: " << solver.iterations() << std::endl;
        libMesh::out << "estimated error: " << solver.error() << std::endl;
        retval = std::make_pair(solver.iterations(), solver.error());
        _comp_info = solver.info();
        break;
      }

    case SPARSELU:
      {
        // SparseLU solver code adapted from:
        // http://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseLU.html
        //
        // From Eigen docs:
        // The input matrix A should be in a compressed and
        // column-major form. Otherwise an expensive copy will be
        // made. You can call the inexpensive makeCompressed() to get
        // a compressed matrix.
        //
        // Note: we don't have a column-major storage format here, so
        // I think a copy must be made in order to use SparseLU.  It
        // appears that we also have to call makeCompressed(),
        // otherwise you get a segfault.
        matrix._mat.makeCompressed();

        // Build the SparseLU solver object.  Note, there is one other
        // sparse direct solver available in Eigen:
        //
        // Eigen::SparseQR<EigenSM, Eigen::AMDOrdering<int>> solver;
        //
        // I've tested it, and it works, but it is much slower than
        // SparseLU.  The main benefit of SparseQR is that it can
        // handle non-square matrices, but we don't allow non-square
        // sparse matrices to be built in libmesh...
        Eigen::SparseLU<EigenSM> solver;

        // Compute the ordering permutation vector from the structural pattern of the matrix.
        solver.analyzePattern(matrix._mat);

        // Compute the numerical factorization
        solver.factorize(matrix._mat);

        // Use the factors to solve the linear system
        solution._vec = solver.solve(rhs._vec);

        // Set up the return value.  The SparseLU solver doesn't
        // support asking for the number of iterations or the final
        // error, so we'll just report back 1 and 0, respectively.
        retval = std::make_pair(/*n. iterations=*/1, /*error=*/0);

        // Store the success/failure reason and break out.
        _comp_info = solver.info();
        break;
      }

      // Unknown solver, use BICGSTAB
    default:
      {
        libMesh::err << "ERROR:  Unsupported Eigen Solver: "
                     << Utility::enum_to_string(this->_solver_type) << std::endl
                     << "Continuing with BICGSTAB" << std::endl;

        this->_solver_type = BICGSTAB;

        return this->solve (matrix,
                            solution,
                            rhs,
                            tol,
                            m_its);
      }
    }

  return retval;
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::adjoint_solve (SparseMatrix<T> & matrix_in,
                                           NumericVector<T> & solution_in,
                                           NumericVector<T> & rhs_in,
                                           const double tol,
                                           const unsigned int m_its)
{
  LOG_SCOPE("adjoint_solve()", "EigenSparseLinearSolver");

  libmesh_experimental();
  EigenSparseMatrix<T> mat_trans(this->comm());
  matrix_in.get_transpose(mat_trans);

  std::pair<unsigned int, Real> retval = this->solve (mat_trans,
                                                      solution_in,
                                                      rhs_in,
                                                      tol,
                                                      m_its);

  return retval;
}




template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (const ShellMatrix<T> & /*shell_matrix*/,
                                   NumericVector<T> & /*solution_in*/,
                                   NumericVector<T> & /*rhs_in*/,
                                   const double /*tol*/,
                                   const unsigned int /*m_its*/)
{
  libmesh_not_implemented();
  return std::make_pair(0,0.0);
}



template <typename T>
std::pair<unsigned int, Real>
EigenSparseLinearSolver<T>::solve (const ShellMatrix<T> & /*shell_matrix*/,
                                   const SparseMatrix<T> & /*precond_matrix*/,
                                   NumericVector<T> & /*solution_in*/,
                                   NumericVector<T> & /*rhs_in*/,
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
  //     _precond_type = libmesh_nullptr; return;

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
LinearConvergenceReason EigenSparseLinearSolver<T>::get_converged_reason() const
{
  std::map<Eigen::ComputationInfo, LinearConvergenceReason>::iterator it =
    _convergence_reasons.find(_comp_info);

  // If later versions of Eigen start returning new enumerations,
  // we'll need to add them to the map...
  if (it == _convergence_reasons.end())
    {
      libmesh_warning("Warning: unknown Eigen::ComputationInfo: " \
                      << _comp_info \
                      << " returning CONVERGED_ITS." \
                      << std::endl);
      return CONVERGED_ITS;
    }
  else
    return it->second;
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseLinearSolver<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
