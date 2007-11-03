// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#include "libmesh_common.h"

#if defined(HAVE_LASPACK)


// C++ includes

// Local Includes
#include "laspack_linear_solver.h"
#include "libmesh_logging.h"

// #ifndef USE_COMPLEX_NUMBERS
// extern "C"
// {
// #endif

//   void print_iter_accuracy(int Iter,
// 			   _LPReal rNorm,
// 			   _LPReal bNorm,
// 			   IterIdType IterId)
//     /* put out accuracy reached after each solver iteration */
//   {
    
//     //FILE* out = fopen("residual.hist", "a");
//     static int icall=0;
    
//     if (!icall)
//       {
// 	printf("Iter   ||r||/||f||\n");
// 	printf("------------------\n");
// 	icall=1;
//       }
    
//     if ( Iter%1==0 && (IterId == CGIterId ||
// 		       IterId == CGNIterId ||
// 		       IterId == GMRESIterId ||
// 		       IterId == BiCGIterId ||
// 		       IterId == QMRIterId ||
// 		       IterId == CGSIterId ||
// 		       IterId == BiCGSTABIterId)  )
//       {
// 	if (!_LPIsZeroReal(bNorm))
// 	  printf("%d    \t %g\n", Iter, rNorm/bNorm);
// 	else
// 	  printf("%d     (fnorm == 0)\n", Iter);
//       }
    
//     //fclose(out);
//   }

// #ifndef USE_COMPLEX_NUMBERS
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
void LaspackLinearSolver<T>::init ()
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
LaspackLinearSolver<T>::solve (SparseMatrix<T> &matrix_in,
			       NumericVector<T> &solution_in,
			       NumericVector<T> &rhs_in,
			       const double tol,
			       const unsigned int m_its)
{
  START_LOG("solve()", "LaspackLinearSolver");
  this->init ();

  LaspackMatrix<T>* matrix   = dynamic_cast<LaspackMatrix<T>*>(&matrix_in);
  LaspackVector<T>* solution = dynamic_cast<LaspackVector<T>*>(&solution_in);
  LaspackVector<T>* rhs      = dynamic_cast<LaspackVector<T>*>(&rhs_in);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(matrix   != NULL);
  assert(solution != NULL);
  assert(rhs      != NULL);
  
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
	std::cerr << "ERROR:  Unsupported LASPACK Solver: "
		  << this->_solver_type      << std::endl
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
      std::cerr << "ERROR:  LASPACK Error: " << std::endl;
      WriteLASErrDescr(stdout);
      error();
    }

  STOP_LOG("solve()", "LaspackLinearSolver");
  // Get the convergence step # and residual 
  return std::make_pair(GetLastNoIter(), GetLastAccuracy());
}



template <typename T>
void LaspackLinearSolver<T>::set_laspack_preconditioner_type ()
{
  switch (this->_preconditioner_type)
    {
    case IDENTITY_PRECOND:
      _precond_type = NULL; return;

    case ILU_PRECOND:
      _precond_type = ILUPrecond; return;

    case JACOBI_PRECOND:
      _precond_type = JacobiPrecond; return;

    case SSOR_PRECOND:
      _precond_type = SSORPrecond; return;


    default:
      std::cerr << "ERROR:  Unsupported LASPACK Preconditioner: "
		<< this->_preconditioner_type << std::endl
		<< "Continuing with ILU"      << std::endl;
      this->_preconditioner_type = ILU_PRECOND;
      this->set_laspack_preconditioner_type();      
    }
}



template <typename T>
void LaspackLinearSolver<T>::print_converged_reason()
{
  std::cout << "print_converged_reason() is currently only supported"
            << "with Petsc 2.3.1 and later." << std::endl;
}



//------------------------------------------------------------------
// Explicit instantiations
template class LaspackLinearSolver<Number>;
 

#endif // #ifdef HAVE_LASPACK
