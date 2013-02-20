// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

namespace libMesh
{

/*----------------------- functions ----------------------------------*/
template <typename T>
void EigenSparseLinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      this->_solver_type         = GMRES;
      this->_preconditioner_type = ILU_PRECOND;
    }
}



template <typename T>
void EigenSparseLinearSolver<T>::init ()
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

  // Make sure the data passed in are really in Laspack types
  EigenSparseMatrix<T>* matrix   = libmesh_cast_ptr<EigenSparseMatrix<T>*>(&matrix_in);
  EigenSparseVector<T>* solution = libmesh_cast_ptr<EigenSparseVector<T>*>(&solution_in);
  EigenSparseVector<T>* rhs      = libmesh_cast_ptr<EigenSparseVector<T>*>(&rhs_in);

  libmesh_not_implemented();

  return std::make_pair(0,0.);
  
  // // Zero-out the solution to prevent the solver from exiting in 0
  // // iterations (?)
  // //TODO:[BSK] Why does Laspack do this?  Comment out this and try ex13...
  // solution->zero();

  // // Close the matrix and vectors in case this wasn't already done.
  // matrix->close ();
  // solution->close ();
  // rhs->close ();

  // // Set the preconditioner type
  // this->set_laspack_preconditioner_type ();

  // // Set the solver tolerance
  // SetRTCAccuracy (tol);

  // // Solve the linear system
  // switch (this->_solver_type)
  //   {
  //     // Conjugate-Gradient
  //   case CG:
  //     {
  // 	CGIter (&matrix->_QMat,
  // 		&solution->_vec,
  // 		&rhs->_vec,
  // 		m_its,
  // 		_precond_type,
  // 		1.);
  // 	break;
  //     }

  //     // Conjugate-Gradient Normalized
  //   case CGN:
  //     {
  // 	CGNIter (&matrix->_QMat,
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Conjugate-Gradient Squared
  //   case CGS:
  //     {
  // 	CGSIter (&matrix->_QMat,
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Bi-Conjugate Gradient
  //   case BICG:
  //     {
  // 	BiCGIter (&matrix->_QMat,
  // 		  &solution->_vec,
  // 		  &rhs->_vec,
  // 		  m_its,
  // 		  _precond_type,
  // 		  1.);
  // 	break;
  //     }

  //     // Bi-Conjugate Gradient Stabilized
  //   case BICGSTAB:
  //     {
  // 	BiCGSTABIter (&matrix->_QMat,
  // 		      &solution->_vec,
  // 		      &rhs->_vec,
  // 		      m_its,
  // 		      _precond_type,
  // 		      1.);
  // 	break;
  //     }

  //     // Quasi-Minimum Residual
  //   case QMR:
  //     {
  // 	QMRIter (&matrix->_QMat,
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Symmetric over-relaxation
  //   case SSOR:
  //     {
  // 	SSORIter (&matrix->_QMat,
  // 		  &solution->_vec,
  // 		  &rhs->_vec,
  // 		  m_its,
  // 		  _precond_type,
  // 		  1.);
  // 	break;
  //     }

  //     // Jacobi Relaxation
  //   case JACOBI:
  //     {
  // 	JacobiIter (&matrix->_QMat,
  // 		    &solution->_vec,
  // 		    &rhs->_vec,
  // 		    m_its,
  // 		    _precond_type,
  // 		    1.);
  // 	break;
  //     }

  //     // Generalized Minimum Residual
  //   case GMRES:
  //     {
  // 	SetGMRESRestart (30);
  // 	GMRESIter (&matrix->_QMat,
  // 		   &solution->_vec,
  // 		   &rhs->_vec,
  // 		   m_its,
  // 		   _precond_type,
  // 		   1.);
  // 	break;
  //     }

  //     // Unknown solver, use GMRES
  //   default:
  //     {
  // 	libMesh::err << "ERROR:  Unsupported LASPACK Solver: "
  // 		      << this->_solver_type      << std::endl
  // 		      << "Continuing with GMRES" << std::endl;

  // 	this->_solver_type = GMRES;

  // 	return this->solve (*matrix,
  // 			    *solution,
  // 			    *rhs,
  // 			    tol,
  // 			    m_its);
  //     }
  //   }

  // // Check for an error
  // if (LASResult() != LASOK)
  //   {
  //     libMesh::err << "ERROR:  LASPACK Error: " << std::endl;
  //     WriteLASErrDescr(stdout);
  //     libmesh_error();
  //   }

  // STOP_LOG("solve()", "EigenSparseLinearSolver");
  // // Get the convergence step # and residual
  // return std::make_pair(GetLastNoIter(), GetLastAccuracy());
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
  this->init ();

  // Make sure the data passed in are really in Laspack types
  EigenSparseMatrix<T>* matrix   = libmesh_cast_ptr<EigenSparseMatrix<T>*>(&matrix_in);
  EigenSparseVector<T>* solution = libmesh_cast_ptr<EigenSparseVector<T>*>(&solution_in);
  EigenSparseVector<T>* rhs      = libmesh_cast_ptr<EigenSparseVector<T>*>(&rhs_in);

  libmesh_not_implemented();
  return std::make_pair(0,0.);
  
  // // Zero-out the solution to prevent the solver from exiting in 0
  // // iterations (?)
  // //TODO:[BSK] Why does Laspack do this?  Comment out this and try ex13...
  // solution->zero();

  // // Close the matrix and vectors in case this wasn't already done.
  // matrix->close ();
  // solution->close ();
  // rhs->close ();

  // // Set the preconditioner type
  // this->set_laspack_preconditioner_type ();

  // // Set the solver tolerance
  // SetRTCAccuracy (tol);

  // // Solve the linear system
  // switch (this->_solver_type)
  //   {
  //     // Conjugate-Gradient
  //   case CG:
  //     {
  // 	CGIter (Transp_Q(&matrix->_QMat),
  // 		&solution->_vec,
  // 		&rhs->_vec,
  // 		m_its,
  // 		_precond_type,
  // 		1.);
  // 	break;
  //     }

  //     // Conjugate-Gradient Normalized
  //   case CGN:
  //     {
  // 	CGNIter (Transp_Q(&matrix->_QMat),
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Conjugate-Gradient Squared
  //   case CGS:
  //     {
  // 	CGSIter (Transp_Q(&matrix->_QMat),
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Bi-Conjugate Gradient
  //   case BICG:
  //     {
  // 	BiCGIter (Transp_Q(&matrix->_QMat),
  // 		  &solution->_vec,
  // 		  &rhs->_vec,
  // 		  m_its,
  // 		  _precond_type,
  // 		  1.);
  // 	break;
  //     }

  //     // Bi-Conjugate Gradient Stabilized
  //   case BICGSTAB:
  //     {
  // 	BiCGSTABIter (Transp_Q(&matrix->_QMat),
  // 		      &solution->_vec,
  // 		      &rhs->_vec,
  // 		      m_its,
  // 		      _precond_type,
  // 		      1.);
  // 	break;
  //     }

  //     // Quasi-Minimum Residual
  //   case QMR:
  //     {
  // 	QMRIter (Transp_Q(&matrix->_QMat),
  // 		 &solution->_vec,
  // 		 &rhs->_vec,
  // 		 m_its,
  // 		 _precond_type,
  // 		 1.);
  // 	break;
  //     }

  //     // Symmetric over-relaxation
  //   case SSOR:
  //     {
  // 	SSORIter (Transp_Q(&matrix->_QMat),
  // 		  &solution->_vec,
  // 		  &rhs->_vec,
  // 		  m_its,
  // 		  _precond_type,
  // 		  1.);
  // 	break;
  //     }

  //     // Jacobi Relaxation
  //   case JACOBI:
  //     {
  // 	JacobiIter (Transp_Q(&matrix->_QMat),
  // 		    &solution->_vec,
  // 		    &rhs->_vec,
  // 		    m_its,
  // 		    _precond_type,
  // 		    1.);
  // 	break;
  //     }

  //     // Generalized Minimum Residual
  //   case GMRES:
  //     {
  // 	SetGMRESRestart (30);
  // 	GMRESIter (Transp_Q(&matrix->_QMat),
  // 		   &solution->_vec,
  // 		   &rhs->_vec,
  // 		   m_its,
  // 		   _precond_type,
  // 		   1.);
  // 	break;
  //     }

  //     // Unknown solver, use GMRES
  //   default:
  //     {
  // 	libMesh::err << "ERROR:  Unsupported LASPACK Solver: "
  // 		      << this->_solver_type      << std::endl
  // 		      << "Continuing with GMRES" << std::endl;

  // 	this->_solver_type = GMRES;

  // 	return this->solve (*matrix,
  // 			    *solution,
  // 			    *rhs,
  // 			    tol,
  // 			    m_its);
  //     }
  //   }

  // // Check for an error
  // if (LASResult() != LASOK)
  //   {
  //     libMesh::err << "ERROR:  LASPACK Error: " << std::endl;
  //     WriteLASErrDescr(stdout);
  //     libmesh_error();
  //   }

  // STOP_LOG("adjoint_solve()", "EigenSparseLinearSolver");
  // // Get the convergence step # and residual
  // return std::make_pair(GetLastNoIter(), GetLastAccuracy());
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
  // 		    << this->_preconditioner_type << std::endl
  // 		    << "Continuing with ILU"      << std::endl;
  //     this->_preconditioner_type = ILU_PRECOND;
  //     this->set_laspack_preconditioner_type();
  //   }
}


  
template <typename T>
void EigenSparseLinearSolver<T>::print_converged_reason()
{
  libMesh::out << "print_converged_reason() is currently only supported"
	       << "with Petsc 2.3.1 and later." << std::endl;
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSparseLinearSolver<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_EIGEN
