// $Id: laspack_interface.C,v 1.1 2003-02-10 03:55:51 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#include "mesh_common.h"

#if defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)


// C++ includes
#include <math.h>

// Local Includes
#include "laspack_interface.h"

extern "C"
{
  void print_iter_accuracy(int Iter,
			   double rNorm,
			   double bNorm,
			   Laspack::IterIdType IterId)
    /* put out accuracy reached after each solver iteration */
  {
    
    //FILE* out = fopen("residual.hist", "a");
    static int icall=0;
    
    if (!icall)
      {
	printf("Iter   ||r||/||f||\n");
	printf("------------------\n");
	icall=1;
      }
    
    if ( Iter%1==0 && (IterId == Laspack::CGIterId ||
		       IterId == Laspack::CGNIterId ||
		       IterId == Laspack::GMRESIterId ||
		       IterId == Laspack::BiCGIterId ||
		       IterId == Laspack::QMRIterId ||
		       IterId == Laspack::CGSIterId ||
		       IterId == Laspack::BiCGSTABIterId)  )
      {
	if (!IsZero(bNorm))
	  printf("%d    \t %g\n", Iter, rNorm/bNorm);
	else
	  printf("%d     (fnorm == 0)\n", Iter);
      }
    
    //fclose(out);
  }
}
  

/*----------------------- functions ----------------------------------*/
void LaspackInterface::clear ()
{
  if (initialized())
    {
      _is_initialized = false;
      
      _solver_type = GMRES;
      _preconditioner_type = ILU_PRECOND;
    };
};



void LaspackInterface::init ()
{
  // Initialize the data structures if not done so already.
  if (!initialized())
    {
      _is_initialized = true;
    };

  Laspack::SetRTCAuxProc (print_iter_accuracy);
};



std::pair<unsigned int, Real> 
LaspackInterface::solve (SparseMatrix &matrix_in,
			 NumericVector &solution_in,
			 NumericVector &rhs_in,
			 const double tol,
			 const unsigned int m_its)
{
  init ();

  LaspackMatrix& matrix   = reinterpret_cast<LaspackMatrix&>(matrix_in);
  LaspackVector& solution = reinterpret_cast<LaspackVector&>(solution_in);
  LaspackVector& rhs      = reinterpret_cast<LaspackVector&>(rhs_in);
  
  // Close the matrix and vectors in case this wasn't already done.
  matrix.close ();
  solution.close ();
  rhs.close ();

  // Set the preconditioner type
  set_preconditioner_type(JACOBI_PRECOND);
  set_solver_type(CG);
  set_laspack_preconditioner_type ();


  // Solve the linear system
  switch (_solver_type)
    {
    case CG:
      Laspack::CGIter (matrix._QMat, solution._vec, rhs._vec,
		       m_its, _precond_type, 1.); break;

    case CGN:
      Laspack::CGNIter (matrix._QMat, solution._vec, rhs._vec,
			m_its, _precond_type, 1.); break;
      
    case CGS:
      Laspack::CGSIter (matrix._QMat, solution._vec, rhs._vec,
			m_its, _precond_type, 1.); break;

    case BICG:
      Laspack::BiCGIter (matrix._QMat, solution._vec, rhs._vec,
			 m_its, _precond_type, 1.); break;

    case BICGSTAB:
      Laspack::BiCGSTABIter (matrix._QMat, solution._vec, rhs._vec,
			     m_its, _precond_type, 1.); break;

    case QMR:
      Laspack::QMRIter (matrix._QMat, solution._vec, rhs._vec,
			m_its, _precond_type, 1.); break;

    case SSOR:
      Laspack::SSORIter (matrix._QMat, solution._vec, rhs._vec,
			 m_its, _precond_type, 1.); break;

    case JACOBI:
      Laspack::JacobiIter (matrix._QMat, solution._vec, rhs._vec,
			   m_its, _precond_type, 1.); break;

     case GMRES:
       Laspack::SetGMRESRestart (30);
       Laspack::GMRESIter (matrix._QMat, solution._vec, rhs._vec,
			   m_its, _precond_type, 1.); break;
      
    default:
      std::cerr << "ERROR:  Unsupported LASPACK Solver: "
		<< _solver_type            << std::endl
		<< "Continuing with GMRES" << std::endl;
      _solver_type = GMRES;
      return solve (matrix, solution, rhs, tol, m_its);
    };

  here();

  if (Laspack::LASResult() != Laspack::LASOK)
    {
      std::cerr << "ERROR:  LASPACK Error: " << std::endl;
      Laspack::WriteLASErrDescr(stdout);
      error();
    }
        
  std::pair<unsigned int, Real> p (Laspack::GetLastNoIter(),
				   Laspack::GetLastAccuracy());
  here();
  return p;
};



void LaspackInterface::set_laspack_preconditioner_type ()
{
  switch (_preconditioner_type)
    {
    case IDENTITY_PRECOND:
      _precond_type = NULL; return;

    case ILU_PRECOND:
      _precond_type = Laspack::ILUPrecond; return;

    case JACOBI_PRECOND:
      _precond_type = Laspack::JacobiPrecond; return;

    case SSOR_PRECOND:
      _precond_type = Laspack::SSORPrecond; return;


    default:
      std::cerr << "ERROR:  Unsupported LASPACK Preconditioner: "
		<< _preconditioner_type  << std::endl
		<< "Continuing with ILU" << std::endl;
      _preconditioner_type = ILU_PRECOND;
      set_laspack_preconditioner_type();      
    };
};



#endif // #ifdef HAVE_LASPACK
