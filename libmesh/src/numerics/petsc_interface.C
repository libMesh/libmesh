//    $Id: petsc_interface.C,v 1.3 2003-01-20 17:06:46 jwpeterson Exp $

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

#ifdef HAVE_PETSC


// C++ includes



// Local Includes
#include "petsc_interface.h"


/*----------------------- functions ----------------------------------*/
void PetscInterface::clear ()
{
  if (initialized)
    {
      initialized = false;

      int ierr=0;

      ierr = SLESDestroy(sles); CHKERRQ(ierr);
    };
};



void PetscInterface::init ()
{
  int ierr=0;
  
  // Initialize the data structures if not done so already.
  if (!initialized)
    {
      initialized = true;

      // Create the linear solver context
      ierr = SLESCreate (PETSC_COMM_WORLD, &sles); CHKERRQ(ierr);
      
      // Create the Krylov subspace & preconditioner contexts
      ierr = SLESGetKSP       (sles, &ksp); CHKERRQ(ierr);
      ierr = SLESGetPC        (sles, &pc);  CHKERRQ(ierr);
      
      // Have the Krylov subspace method use our good initial guess rather than 0
      // 2.1.0 - style
      //ierr = KSPSetInitialGuessNonzero (ksp); CHKERRQ(ierr);
      // 2.1.1 - style
      ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE); CHKERRQ(ierr);
      
      // Set the options from user-input
      // Set runtime options, e.g.,
      //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      //  These options will override those specified above as long as
      //  SLESSetFromOptions() is called _after_ any other customization
      //  routines.
      
      ierr = SLESSetFromOptions (sles);                   CHKERRQ(ierr);
    };
};



std::pair<unsigned int, real> 
PetscInterface::solve (PetscMatrix &matrix,
		       PetscVector &solution,
		       PetscVector &rhs,
		       const double tol,
		       const unsigned int m_its)
{
  init ();
  
  int ierr=0;
  int its=0, max_its = static_cast<int>(m_its);
  PetscReal final_resid=0.;

  // Close the matrix and vectors in case this wasn't already done.
  matrix.close ();
  solution.close ();
  rhs.close ();

  
  // Set operators. The input matrix works as the preconditioning matrix
  ierr = SLESSetOperators(sles, matrix.mat, matrix.mat,
			  SAME_NONZERO_PATTERN);             CHKERRQ(ierr);

  
  ierr = KSPSetTolerances (ksp, tol, tol,
			   PETSC_DEFAULT,
			   max_its);                         CHKERRQ(ierr);

  
  // Solve the linear system
  ierr = SLESSolve (sles, rhs.vec, solution.vec, &its);      CHKERRQ(ierr);
  
  
  ierr = KSPGetResidualNorm (ksp, &final_resid);             CHKERRQ(ierr);
  
  //ierr = KSPTrueMonitor (ksp, its, final_resid, PETSC_NULL); CHKERRQ(ierr);
  
  std::pair<unsigned int, real> p;

  p.first  = static_cast<unsigned int>(its);
  p.second = static_cast<real>(final_resid);

  return p;
};



#endif
