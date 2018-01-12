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



#include "libmesh/libmesh_common.h"

#if defined(LIBMESH_HAVE_LASPACK)


// C++ includes

// Local Includes
#include "libmesh/laspack_linear_solver.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/string_to_enum.h"

namespace libMesh
{

// #ifndef LIBMESH_USE_COMPLEX_NUMBERS
// extern "C"
// {
// #endif

//   void print_iter_accuracy(int Iter,
//    _LPReal rNorm,
//    _LPReal bNorm,
//    IterIdType IterId)
//     /* put out accuracy reached after each solver iteration */
//   {

//     //FILE * out = fopen("residual.hist", "a");
//     static int icall=0;

//     if (!icall)
//       {
// printf("Iter   ||r||/||f||\n");
// printf("------------------\n");
// icall=1;
//       }

//     if ( Iter%1==0 && (IterId == CGIterId ||
//        IterId == CGNIterId ||
//        IterId == GMRESIterId ||
//        IterId == BiCGIterId ||
//        IterId == QMRIterId ||
//        IterId == CGSIterId ||
//        IterId == BiCGSTABIterId)  )
//       {
// if (!_LPIsZeroReal(bNorm))
//   printf("%d    \t %g\n", Iter, rNorm/bNorm);
// else
//   printf("%d     (fnorm == 0)\n", Iter);
//       }

//     //fclose(out);
//   }

// #ifndef LIBMESH_USE_COMPLEX_NUMBERS
// }
// #endif

/*----------------------- functions ----------------------------------*/
template <typename T>
void LaspackLinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      this->_solver_type         = GMRES;
      this->_preconditioner_type = ILU_PRECOND;
    }
}



template <typename T>
void LaspackLinearSolver<T>::init (const char * /* name */)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;
    }

  // SetRTCAuxProc (print_iter_accuracy);
}



template <typename T>
std::pair<unsigned int, Real>
LaspackLinearSolver<T>::solve (SparseMatrix<T> & matrix_in,
                               NumericVector<T> & solution_in,
                               NumericVector<T> & rhs_in,
                               const double tol,
                               const unsigned int m_its)
{
  LOG_SCOPE("solve()", "LaspackLinearSolver");
  this->init ();

  // Make sure the data passed in are really in Laspack types
  LaspackMatrix<T> * matrix   = cast_ptr<LaspackMatrix<T> *>(&matrix_in);
  LaspackVector<T> * solution = cast_ptr<LaspackVector<T> *>(&solution_in);
  LaspackVector<T> * rhs      = cast_ptr<LaspackVector<T> *>(&rhs_in);

  // Zero-out the solution to prevent the solver from exiting in 0
  // iterations (?)
  //TODO:[BSK] Why does Laspack do this?  Comment out this and try ex13...
  solution->zero();

  // Close the matrix and vectors in case this wasn't already done.
  matrix->close ();
  solution->close ();
  rhs->close ();

  // Set the preconditioner type
  this->set_laspack_preconditioner_type ();

  // Set the solver tolerance
  SetRTCAccuracy (tol);

  // Solve the linear system
  switch (this->_solver_type)
    {
      // Conjugate-Gradient
    case CG:
      {
        CGIter (&matrix->_QMat,
                &solution->_vec,
                &rhs->_vec,
                m_its,
                _precond_type,
                1.);
        break;
      }

      // Conjugate-Gradient Normalized
    case CGN:
      {
        CGNIter (&matrix->_QMat,
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Conjugate-Gradient Squared
    case CGS:
      {
        CGSIter (&matrix->_QMat,
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Bi-Conjugate Gradient
    case BICG:
      {
        BiCGIter (&matrix->_QMat,
                  &solution->_vec,
                  &rhs->_vec,
                  m_its,
                  _precond_type,
                  1.);
        break;
      }

      // Bi-Conjugate Gradient Stabilized
    case BICGSTAB:
      {
        BiCGSTABIter (&matrix->_QMat,
                      &solution->_vec,
                      &rhs->_vec,
                      m_its,
                      _precond_type,
                      1.);
        break;
      }

      // Quasi-Minimum Residual
    case QMR:
      {
        QMRIter (&matrix->_QMat,
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Symmetric over-relaxation
    case SSOR:
      {
        SSORIter (&matrix->_QMat,
                  &solution->_vec,
                  &rhs->_vec,
                  m_its,
                  _precond_type,
                  1.);
        break;
      }

      // Jacobi Relaxation
    case JACOBI:
      {
        JacobiIter (&matrix->_QMat,
                    &solution->_vec,
                    &rhs->_vec,
                    m_its,
                    _precond_type,
                    1.);
        break;
      }

      // Generalized Minimum Residual
    case GMRES:
      {
        SetGMRESRestart (30);
        GMRESIter (&matrix->_QMat,
                   &solution->_vec,
                   &rhs->_vec,
                   m_its,
                   _precond_type,
                   1.);
        break;
      }

      // Unknown solver, use GMRES
    default:
      {
        libMesh::err << "ERROR:  Unsupported LASPACK Solver: "
                     << Utility::enum_to_string(this->_solver_type) << std::endl
                     << "Continuing with GMRES" << std::endl;

        this->_solver_type = GMRES;

        return this->solve (*matrix,
                            *solution,
                            *rhs,
                            tol,
                            m_its);
      }
    }

  // Check for an error
  if (LASResult() != LASOK)
    {
      WriteLASErrDescr(stdout);
      libmesh_error_msg("Exiting after LASPACK Error!");
    }

  // Get the convergence step # and residual
  return std::make_pair(GetLastNoIter(), GetLastAccuracy());
}



template <typename T>
std::pair<unsigned int, Real>
LaspackLinearSolver<T>::adjoint_solve (SparseMatrix<T> & matrix_in,
                                       NumericVector<T> & solution_in,
                                       NumericVector<T> & rhs_in,
                                       const double tol,
                                       const unsigned int m_its)
{
  LOG_SCOPE("adjoint_solve()", "LaspackLinearSolver");
  this->init ();

  // Make sure the data passed in are really in Laspack types
  LaspackMatrix<T> * matrix   = cast_ptr<LaspackMatrix<T> *>(&matrix_in);
  LaspackVector<T> * solution = cast_ptr<LaspackVector<T> *>(&solution_in);
  LaspackVector<T> * rhs      = cast_ptr<LaspackVector<T> *>(&rhs_in);

  // Zero-out the solution to prevent the solver from exiting in 0
  // iterations (?)
  //TODO:[BSK] Why does Laspack do this?  Comment out this and try ex13...
  solution->zero();

  // Close the matrix and vectors in case this wasn't already done.
  matrix->close ();
  solution->close ();
  rhs->close ();

  // Set the preconditioner type
  this->set_laspack_preconditioner_type ();

  // Set the solver tolerance
  SetRTCAccuracy (tol);

  // Solve the linear system
  switch (this->_solver_type)
    {
      // Conjugate-Gradient
    case CG:
      {
        CGIter (Transp_Q(&matrix->_QMat),
                &solution->_vec,
                &rhs->_vec,
                m_its,
                _precond_type,
                1.);
        break;
      }

      // Conjugate-Gradient Normalized
    case CGN:
      {
        CGNIter (Transp_Q(&matrix->_QMat),
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Conjugate-Gradient Squared
    case CGS:
      {
        CGSIter (Transp_Q(&matrix->_QMat),
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Bi-Conjugate Gradient
    case BICG:
      {
        BiCGIter (Transp_Q(&matrix->_QMat),
                  &solution->_vec,
                  &rhs->_vec,
                  m_its,
                  _precond_type,
                  1.);
        break;
      }

      // Bi-Conjugate Gradient Stabilized
    case BICGSTAB:
      {
        BiCGSTABIter (Transp_Q(&matrix->_QMat),
                      &solution->_vec,
                      &rhs->_vec,
                      m_its,
                      _precond_type,
                      1.);
        break;
      }

      // Quasi-Minimum Residual
    case QMR:
      {
        QMRIter (Transp_Q(&matrix->_QMat),
                 &solution->_vec,
                 &rhs->_vec,
                 m_its,
                 _precond_type,
                 1.);
        break;
      }

      // Symmetric over-relaxation
    case SSOR:
      {
        SSORIter (Transp_Q(&matrix->_QMat),
                  &solution->_vec,
                  &rhs->_vec,
                  m_its,
                  _precond_type,
                  1.);
        break;
      }

      // Jacobi Relaxation
    case JACOBI:
      {
        JacobiIter (Transp_Q(&matrix->_QMat),
                    &solution->_vec,
                    &rhs->_vec,
                    m_its,
                    _precond_type,
                    1.);
        break;
      }

      // Generalized Minimum Residual
    case GMRES:
      {
        SetGMRESRestart (30);
        GMRESIter (Transp_Q(&matrix->_QMat),
                   &solution->_vec,
                   &rhs->_vec,
                   m_its,
                   _precond_type,
                   1.);
        break;
      }

      // Unknown solver, use GMRES
    default:
      {
        libMesh::err << "ERROR:  Unsupported LASPACK Solver: "
                     << Utility::enum_to_string(this->_solver_type) << std::endl
                     << "Continuing with GMRES" << std::endl;

        this->_solver_type = GMRES;

        return this->solve (*matrix,
                            *solution,
                            *rhs,
                            tol,
                            m_its);
      }
    }

  // Check for an error
  if (LASResult() != LASOK)
    {
      WriteLASErrDescr(stdout);
      libmesh_error_msg("Exiting after LASPACK Error!");
    }

  // Get the convergence step # and residual
  return std::make_pair(GetLastNoIter(), GetLastAccuracy());
}




template <typename T>
std::pair<unsigned int, Real>
LaspackLinearSolver<T>::solve (const ShellMatrix<T> & /*shell_matrix*/,
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
LaspackLinearSolver<T>::solve (const ShellMatrix<T> & /*shell_matrix*/,
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
void LaspackLinearSolver<T>::set_laspack_preconditioner_type ()
{
  switch (this->_preconditioner_type)
    {
    case IDENTITY_PRECOND:
      _precond_type = libmesh_nullptr; return;

    case ILU_PRECOND:
      _precond_type = ILUPrecond; return;

    case JACOBI_PRECOND:
      _precond_type = JacobiPrecond; return;

    case SSOR_PRECOND:
      _precond_type = SSORPrecond; return;


    default:
      libMesh::err << "ERROR:  Unsupported LASPACK Preconditioner: "
                   << this->_preconditioner_type << std::endl
                   << "Continuing with ILU"      << std::endl;
      this->_preconditioner_type = ILU_PRECOND;
      this->set_laspack_preconditioner_type();
    }
}



template <typename T>
void LaspackLinearSolver<T>::print_converged_reason() const
{
  switch (LASResult())
    {
    case LASOK :
      libMesh::out << "Laspack converged.\n";
      break;
    default    :
      libMesh::out << "Laspack diverged.\n";
    }
  libMesh::out << "Detailed reporting is currently only supported"
               << "with Petsc 2.3.1 and later." << std::endl;
}



template <typename T>
LinearConvergenceReason LaspackLinearSolver<T>::get_converged_reason() const
{
  switch (LASResult())
    {
    case LASOK : return CONVERGED_RTOL_NORMAL;
    default    : return DIVERGED_NULL;
    }
}



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackLinearSolver<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_LASPACK
