// $Id: petsc_nonlinear_solver.C,v 1.17 2007-07-16 15:12:36 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "petsc_vector.h"
#include "petsc_matrix.h"



//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
// 
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  // Older versions of PETSc do not have the different int typedefs.
  // On 64-bit machines, PetscInt may actually be a long long int.
  // This change occurred in Petsc-2.2.1.
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif
  
  //-------------------------------------------------------------------
  // this function is called by PETSc at the end of each nonlinear step  
  PetscErrorCode
  __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *)
  {
    //int ierr=0;

    //if (its > 0)
      std::cout << "  NL step " << its
		<< std::scientific
		<< ", |residual|_2 = " << fnorm
		<< std::endl;

    //return ierr;
    return 0;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __libmesh_petsc_snes_residual (SNES, Vec x, Vec r, void *ctx)
  {
    int ierr=0;

    assert (x   != NULL);
    assert (r   != NULL);
    assert (ctx != NULL);
    
    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);
    
    PetscVector<Number> X_global(x), R(r);
    PetscVector<Number> X_local(X_global.size());

    X_global.localize (X_local);
  
    if (solver->residual != NULL) solver->residual (X_local, R);
    if (solver->matvec   != NULL) solver->matvec   (X_local, &R, NULL);

    R.close();
    
    return ierr;
  }


  
  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __libmesh_petsc_snes_jacobian (SNES, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    int ierr=0;
    
    assert (ctx != NULL);
    
    PetscNonlinearSolver<Number>* solver =
      static_cast<PetscNonlinearSolver<Number>*> (ctx);
    
    PetscMatrix<Number> PC(*pc);
    PetscMatrix<Number> Jac(*jac);
    PetscVector<Number> X_global(x);
    PetscVector<Number> X_local (X_global.size());

    X_global.localize (X_local);

    if (solver->jacobian != NULL) solver->jacobian (X_local, PC);
    if (solver->matvec   != NULL) solver->matvec   (X_local, NULL, &PC);
    
    PC.close();
    Jac.close();
    
    *msflag = SAME_NONZERO_PATTERN;
    
    return ierr;
  }
    
} // end extern "C"
//---------------------------------------------------------------------



//---------------------------------------------------------------------
// PetscNonlinearSolver<> methods
template <typename T>
void PetscNonlinearSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      int ierr=0;

      ierr = SNESDestroy(_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
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

#if PETSC_VERSION_LESS_THAN(2,1,2)
      // At least until Petsc 2.1.1, the SNESCreate had a different calling syntax.
      // The second argument was of type SNESProblemType, and could have a value of
      // either SNES_NONLINEAR_EQUATIONS or SNES_UNCONSTRAINED_MINIMIZATION.
      ierr = SNESCreate(libMesh::COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#else

      ierr = SNESCreate(libMesh::COMM_WORLD,&_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);

#endif	     

#if PETSC_VERSION_LESS_THAN(2,3,3)
      ierr = SNESSetMonitor (_snes, __libmesh_petsc_snes_monitor,
			     this, PETSC_NULL);
#else
      // API name change in PETSc 2.3.3
      ierr = SNESMonitorSet (_snes, __libmesh_petsc_snes_monitor,
			     this, PETSC_NULL);
#endif
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
	     
      ierr = SNESSetFromOptions(_snes);
             CHKERRABORT(libMesh::COMM_WORLD,ierr);
    }
}



template <typename T>
std::pair<unsigned int, Real> 
PetscNonlinearSolver<T>::solve (SparseMatrix<T>&  jac_in,  // System Jacobian Matrix
				NumericVector<T>& x_in,    // Solution vector
				NumericVector<T>& r_in,    // Residual vector
				const double,              // Stopping tolerance
				const unsigned int) 
{
  this->init ();
  
  PetscMatrix<T>* jac = dynamic_cast<PetscMatrix<T>*>(&jac_in);
  PetscVector<T>* x   = dynamic_cast<PetscVector<T>*>(&x_in);
  PetscVector<T>* r   = dynamic_cast<PetscVector<T>*>(&r_in);

  // We cast to pointers so we can be sure that they succeeded
  // by comparing the result against NULL.
  assert(jac != NULL); assert(jac->mat() != NULL);
  assert(x   != NULL); assert(x->vec()   != NULL);
  assert(r   != NULL); assert(r->vec()   != NULL);
  
  int ierr=0;
  int n_iterations =0;

  ierr = SNESSetFunction (_snes, r->vec(), __libmesh_petsc_snes_residual, this);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);

  ierr = SNESSetJacobian (_snes, jac->mat(), jac->mat(), __libmesh_petsc_snes_jacobian, this);
         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	 
// Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
// the last one being a pointer to an int to hold the number of iterations required.
# if PETSC_VERSION_LESS_THAN(2,2,0)

 ierr = SNESSolve (_snes, x->vec(), &n_iterations);
        CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.2.x style	
#elif PETSC_VERSION_LESS_THAN(2,3,0)
	
 ierr = SNESSolve (_snes, x->vec());
        CHKERRABORT(libMesh::COMM_WORLD,ierr);

// 2.3.x & newer style	
#else
	
 ierr = SNESSolve (_snes, PETSC_NULL, x->vec());
        CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  	
#endif

  this->clear();
		 
  // return the # of its. and the final residual norm.  Note that
  // n_iterations may be zero for PETSc versions 2.2.x and greater.
  return std::make_pair(n_iterations, 0.);
}




//------------------------------------------------------------------
// Explicit instantiations
template class PetscNonlinearSolver<Number>;
 


#endif // #ifdef HAVE_PETSC
