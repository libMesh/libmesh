// $Id: petsc_nonlinear_solver.C,v 1.1 2005-01-03 20:14:41 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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

#ifdef HAVE_PETSC


// C++ includes

// Local Includes
#include "petsc_nonlinear_solver.h"



//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods,
extern "C"
{
  //-------------------------------------------------------------------
  // this function is called by PETSc at the end of each nonlinear step  
  PetscErrorCode
  __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *ctx)
  {
    int ierr=0;
    
    assert (ctx != NULL);
    
    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);
    
    std::cout << "  NL conv: step " << its
	      << std::scientific
	      << ", |u|_oo = "      << 0
	      << ", |resid|_2 = "   << fnorm
              << std::endl;

    return ierr;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __libmesh_petsc_snes_function (SNES,Vec,Vec,void*)
  {
    int ierr=0;

    return ierr;
  }


  
  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __libmesh_petsc_snes_jacobian (SNES, Vec, Mat*, Mat*, MatStructure*, void*)
  {
    int ierr=0;

    return ierr;
  }
    
} // end extern "C"


/*----------------------- functions ----------------------------------*/
template <typename T>
void PetscNonlinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

      ierr = SNESDestroy(_snes);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
}



template <typename T>
void PetscNonlinearSolver<T>::init ()
{  
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;
      
      int ierr=0;
      
      ierr = SNESCreate(PETSC_COMM_WORLD,&_snes);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);

      ierr = SNESSetMonitor (_snes, __libmesh_petsc_snes_monitor,
			     this, PETSC_NULL);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
	     
      ierr = SNESSetFromOptions(_snes);
             CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
}








template <typename T>
std::pair<unsigned int, Real> 
PetscNonlinearSolver<T>::solve () 
{
  this->init ();
  

  // return the # of its. and the final residual norm.
  return std::make_pair(0, 0.);
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscNonlinearSolver<Number>;
 


#endif // #ifdef HAVE_PETSC
