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



#include "libmesh_common.h"

#ifdef LIBMESH_HAVE_TRILINOS


// C++ includes

// Local Includes
#include "libmesh_logging.h"
#include "trilinos_aztec_linear_solver.h"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_vector.h"

namespace libMesh
{


/*----------------------- functions ----------------------------------*/
template <typename T>
void AztecLinearSolver<T>::clear ()
{
  if (this->initialized())
  {
    this->_is_initialized = false;
	     
    // Mimic PETSc default solver and preconditioner
    this->_solver_type           = GMRES;

    if (libMesh::n_processors() == 1)
      this->_preconditioner_type = ILU_PRECOND;
    else
      this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
  }
}



template <typename T>
void AztecLinearSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
  {
    this->_is_initialized = true;
 
    _linear_solver = new AztecOO();

    switch(this->_preconditioner_type)
    {
    case ILU_PRECOND:
      _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
      break;
    case BLOCK_JACOBI_PRECOND:
      _linear_solver->SetAztecOption(AZ_precond,AZ_Jacobi);
      break;
    default:
      _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
    }
  }
}




template <typename T>
std::pair<unsigned int, Real> 
AztecLinearSolver<T>::solve (SparseMatrix<T>& matrix_in,
			     SparseMatrix<T>& precond_in,
			     NumericVector<T>& solution_in,
			     NumericVector<T>& rhs_in,
			     const double tol,
			     const unsigned int m_its)
{
  START_LOG("solve()", "AztecLinearSolver");  

  // Make sure the data passed in are really of Epetra types
  EpetraMatrix<T>* matrix   = libmesh_cast_ptr<EpetraMatrix<T>*>(&matrix_in);
  EpetraMatrix<T>* precond  = libmesh_cast_ptr<EpetraMatrix<T>*>(&precond_in);
  EpetraVector<T>* solution = libmesh_cast_ptr<EpetraVector<T>*>(&solution_in);
  EpetraVector<T>* rhs      = libmesh_cast_ptr<EpetraVector<T>*>(&rhs_in);

  this->init();

  // Close the matrices and vectors in case this wasn't already done.
  matrix->close ();
  precond->close ();
  solution->close ();
  rhs->close ();
  
  _linear_solver->SetAztecOption(AZ_max_iter,m_its);
  _linear_solver->SetAztecParam(AZ_tol,tol);

  Epetra_FECrsMatrix * emat = matrix->mat();
  Epetra_Vector * esol = solution->vec();
  Epetra_Vector * erhs = rhs->vec();

  _linear_solver->Iterate(emat, esol, erhs, m_its, tol);

  STOP_LOG("solve()", "AztecLinearSolver");
  
  // return the # of its. and the final residual norm.
  return std::make_pair(_linear_solver->NumIters(), _linear_solver->TrueResidual());
}



template <typename T>
std::pair<unsigned int, Real> 
AztecLinearSolver<T>::solve (const ShellMatrix<T>&,
			     NumericVector<T>&,
			     NumericVector<T>&,
			     const double,
			     const unsigned int)
//AztecLinearSolver<T>::solve (const ShellMatrix<T>& shell_matrix,
//			     NumericVector<T>& solution_in,
//			     NumericVector<T>& rhs_in,
//			     const double tol,
//			     const unsigned int m_its)
{
  libmesh_not_implemented();
}



template <typename T>
std::pair<unsigned int, Real> 
AztecLinearSolver<T>::solve (const ShellMatrix<T>&,
			     const SparseMatrix<T>&,
			     NumericVector<T> &,
			     NumericVector<T> &,
			     const double,
			     const unsigned int)
//AztecLinearSolver<T>::solve (const ShellMatrix<T>& shell_matrix,
//			     const SparseMatrix<T>& precond_matrix,
//			     NumericVector<T> &solution_in,
//			     NumericVector<T> &rhs_in,
//			     const double tol,
//			     const unsigned int m_its)
{
  libmesh_not_implemented();
}



template <typename T>
void AztecLinearSolver<T>::get_residual_history(std::vector<double>& /* hist */)
{
  libmesh_not_implemented();

//   int ierr = 0;
//   int its  = 0;

//   // Fill the residual history vector with the residual norms
//   // Note that GetResidualHistory() does not copy any values, it
//   // simply sets the pointer p.  Note that for some Krylov subspace
//   // methods, the number of residuals returned in the history
//   // vector may be different from what you are expecting.  For
//   // example, TFQMR returns two residual values per iteration step.
//   PetscReal* p;
//   ierr = KSPGetResidualHistory(_ksp, &p, &its);
//   CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   // Check for early return
//   if (its == 0) return;
  
//   // Create space to store the result
//   hist.resize(its);

//   // Copy history into the vector provided by the user.
//   for (int i=0; i<its; ++i)
//     {
//       hist[i] = *p;
//       p++;
//     }
}




template <typename T>
Real AztecLinearSolver<T>::get_initial_residual()
{
  return _linear_solver->TrueResidual();
}



template <typename T>
void AztecLinearSolver<T>::print_converged_reason()
{
  libmesh_not_implemented();

// #if PETSC_VERSION_LESS_THAN(2,3,1)
//   libMesh::out << "This method is currently not supported "
// 	    << "(but may work!) for Petsc 2.3.0 and earlier." << std::endl;
// #else
//   KSPConvergedReason reason;
//   KSPGetConvergedReason(_ksp, &reason);
  
//   //  KSP_CONVERGED_RTOL (residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side)
//   //  KSP_CONVERGED_ATOL (residual 2-norm less than abstol)
//   //  KSP_CONVERGED_ITS (used by the preonly preconditioner that always uses ONE iteration) 
//   //  KSP_CONVERGED_STEP_LENGTH
//   //  KSP_DIVERGED_ITS  (required more than its to reach convergence)
//   //  KSP_DIVERGED_DTOL (residual norm increased by a factor of divtol)
//   //  KSP_DIVERGED_NAN (residual norm became Not-a-number likely do to 0/0)
//   //  KSP_DIVERGED_BREAKDOWN (generic breakdown in method)

//   switch (reason)
//     {
//     case KSP_CONVERGED_RTOL:
//        {
// 	libMesh::out << "Linear solver converged, relative tolerance reached." << std::endl;
// 	break;
//        }
//     case KSP_CONVERGED_ATOL:
//        {
// 	 libMesh::out << "Linear solver converged, absolute tolerance reached." << std::endl;
// 	 break;
//        }

//       // Divergence
//     case KSP_DIVERGED_ITS:
//        {
// 	 libMesh::out << "Linear solver diverged, max no. of iterations reached." << std::endl;
// 	 break;
//        }
//     case KSP_DIVERGED_DTOL:
//        {
// 	 libMesh::out << "Linear solver diverged, residual norm increase by dtol (default 1.e5)." << std::endl;
// 	 break;
//        }
//     case KSP_DIVERGED_NAN:
//        {
// 	 libMesh::out << "Linear solver diverged, residual norm is NaN." << std::endl;
// 	 break;
//        }
//     case KSP_DIVERGED_BREAKDOWN:
//        {
// 	 libMesh::out << "Linear solver diverged, generic breakdown in the method." << std::endl;
// 	 break;
//        }
//     default:
//       {
// 	libMesh::out << "Unknown/unsupported con(di)vergence reason: " << reason << std::endl;
//       }
//     }
// #endif
}


//------------------------------------------------------------------
// Explicit instantiations
template class AztecLinearSolver<Number>;

} // namespace libMesh
 


#endif // #ifdef LIBMESH_HAVE_TRILINOS
