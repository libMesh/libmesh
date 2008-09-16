// $Id: petsc_nonlinear_solver.C 2891 2008-06-25 15:25:47Z friedmud $

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

#ifdef HAVE_NOX


// C++ includes

// Local Includes
#include "nonlinear_implicit_system.h"
#include "trilinos_nox_nonlinear_solver.h"
#include "system.h"
#include "trilinos_epetra_vector.h"

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"	// class definition

// ---------- Forward Declarations ----------

class Problem_Interface : public NOX::Epetra::Interface::Required,
			  public NOX::Epetra::Interface::Jacobian,
			  public NOX::Epetra::Interface::Preconditioner
{
public:
  Problem_Interface(NoxNonlinearSolver<Number> * solver);
  ~Problem_Interface();
  
  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
		NOX::Epetra::Interface::Required::FillType fillType);
	
  //! Compute an explicit Jacobian
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);
  
  //! Compute the Epetra_RowMatrix M, that will be used by the Aztec preconditioner instead of the Jacobian.  This is used when there is no explicit Jacobian present (i.e. Matrix-Free Newton-Krylov).  This MUST BE an Epetra_RowMatrix since the Aztec preconditioners need to know the sparsity pattern of the matrix.  Returns true if computation was successful.
  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M);
  
  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& Prec,
			     Teuchos::ParameterList* p);
  
  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */ 
  NoxNonlinearSolver<Number> * _solver;
};

// ----------   Includes   ----------
//#include <iostream>
//#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
//#include "Brusselator.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(NoxNonlinearSolver<Number> * solver) :
  _solver(solver)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& r,
                       NOX::Epetra::Interface::Required::FillType fillType)
{
  libmesh_assert(xfe != NULL);
  libmesh_assert(rfe != NULL);
  
  NonlinearImplicitSystem &sys = _solver->system();

  EpetraVector<Number> X_global(x), R(r);
  EpetraVector<Number>& X_sys = *dynamic_cast<EpetraVector<Number>*>(sys.solution.get());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);
  sys.update();
  X_global.swap(X_sys);
  
  R.zero();
  
  if( _solver->residual == NULL )
    return false;
  _solver->residual (*sys.current_local_solution.get(), R);
  return true;
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                                        Epetra_Operator& Jac)
{
  throw 1;
  
  // return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
//   cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
   throw 1;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                                              Epetra_Operator& Prec,
                                              Teuchos::ParameterList* p)
{
//   cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
   throw 1;
}

// //--------------------------------------------------------------------
// // Functions with C linkage to pass to PETSc.  PETSc will call these
// // methods as needed.
// // 
// // Since they must have C linkage they have no knowledge of a namespace.
// // Give them an obscure name to avoid namespace pollution.
// extern "C"
// {
//   // Older versions of PETSc do not have the different int typedefs.
//   // On 64-bit machines, PetscInt may actually be a long long int.
//   // This change occurred in Petsc-2.2.1.
// #if PETSC_VERSION_LESS_THAN(2,2,1)
//   typedef int PetscErrorCode;
//   typedef int PetscInt;
// #endif
  
//   //-------------------------------------------------------------------
//   // this function is called by PETSc at the end of each nonlinear step  
//   PetscErrorCode
//   __libmesh_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *)
//   {
//     //int ierr=0;

//     //if (its > 0)
//       std::cout << "  NL step " << its
// 		<< std::scientific
// 		<< ", |residual|_2 = " << fnorm
// 		<< std::endl;

//     //return ierr;
//     return 0;
//   }



//   //---------------------------------------------------------------
//   // this function is called by PETSc to evaluate the residual at X
//   PetscErrorCode
//   __libmesh_petsc_snes_residual (SNES, Vec x, Vec r, void *ctx)
//   {
//     int ierr=0;

//     libmesh_assert (x   != NULL);
//     libmesh_assert (r   != NULL);
//     libmesh_assert (ctx != NULL);
    
//     PetscNonlinearSolver<Number>* solver =
//       static_cast<PetscNonlinearSolver<Number>*> (ctx);
    
//     NonlinearImplicitSystem &sys = solver->system();

//     PetscVector<Number> X_global(x), R(r);
//     PetscVector<Number>& X_sys = *dynamic_cast<PetscVector<Number>*>(sys.solution.get());

//     // Use the systems update() to get a good local version of the parallel solution
//     X_global.swap(X_sys);
//     sys.update();
//     X_global.swap(X_sys);
  
//     R.zero();

//     if      (solver->residual != NULL) solver->residual (*sys.current_local_solution.get(), R);
//     else if (solver->matvec   != NULL) solver->matvec   (*sys.current_local_solution.get(), &R, NULL);


//     R.close();

    
//     return ierr;
//   }


  
//   //---------------------------------------------------------------
//   // this function is called by PETSc to evaluate the Jacobian at X
//   PetscErrorCode
//   __libmesh_petsc_snes_jacobian (SNES, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
//   {
//     int ierr=0;
    
//     libmesh_assert (ctx != NULL);
    
//     PetscNonlinearSolver<Number>* solver =
//       static_cast<PetscNonlinearSolver<Number>*> (ctx);

//     NonlinearImplicitSystem &sys = solver->system();
    
//     PetscMatrix<Number> PC(*pc);
//     PetscMatrix<Number> Jac(*jac);
//     PetscVector<Number> X_global(x);
//     PetscVector<Number>& X_sys = *dynamic_cast<PetscVector<Number>*>(sys.solution.get());

//     // Use the systems update() to get a good local version of the parallel solution
//     X_global.swap(X_sys);
//     sys.update();
//     X_global.swap(X_sys);

//     PC.zero();

//     if      (solver->jacobian != NULL) solver->jacobian (*sys.current_local_solution.get(), PC);
//     else if (solver->matvec   != NULL) solver->matvec   (*sys.current_local_solution.get(), NULL, &PC);
    
//     PC.close();
//     Jac.close();
    
//     *msflag = SAME_NONZERO_PATTERN;
    
//     return ierr;
//   }
    
// } // end extern "C"
// //---------------------------------------------------------------------



//---------------------------------------------------------------------
// NoxNonlinearSolver<> methods
template <typename T>
void NoxNonlinearSolver<T>::clear ()
{
 //  if (this->initialized())
//     {
//       this->_is_initialized = false;

//       int ierr=0;

//       ierr = SNESDestroy(_snes);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }
}



template <typename T>
void NoxNonlinearSolver<T>::init ()
{
  if (!this->initialized())
    _interface = new Problem_Interface(this);
}



template <typename T>
std::pair<unsigned int, Real> 
NoxNonlinearSolver<T>::solve (SparseMatrix<T>&  jac_in,  // System Jacobian Matrix
				NumericVector<T>& x_in,    // Solution vector
				NumericVector<T>& r_in,    // Residual vector
				const double,              // Stopping tolerance
				const unsigned int) 
{
  this->init ();

  EpetraVector<T> * x_epetra = dynamic_cast<EpetraVector<T>*>(&x_in);
  EpetraVector<T> * r_epetra = dynamic_cast<EpetraVector<T>*>(&r_in);

  libmesh_assert(x_epetra != NULL);
  libmesh_assert(r_epetra != NULL);

  Teuchos::RCP<Epetra_Vector> x_t(x_epetra->vec());
  Teuchos::RCP<Epetra_Vector> r_t(r_epetra->vec());

  NOX::Epetra::Vector x(x_t, NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector r(r_t, NOX::Epetra::Vector::CreateView);
  
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());
  nlParams.set("Nonlinear Solver", "Line Search Based");

  //print params	
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  
  //create linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq(_interface);	
  Teuchos::RCP<NOX::Epetra::MatrixFree> MF = 
    Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams,
                                             iReq,
                                             x));
  	
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");
  
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES"); 
  lsParams.set("Max Iterations", 800); 
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 100);	 
//  lsParams.set("Preconditioner", "AztecOO");
  
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = MF;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, MF,
                                                      x));
  //create group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
					iReq,
					x,
					linSys)); 
  NOX::Epetra::Group& grp = *(grpPtr.get());
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;
  status = solver->solve(); 






//   this->init ();
  
//   NoxMatrix<T>* jac = dynamic_cast<NoxMatrix<T>*>(&jac_in);
//   NoxVector<T>* x   = dynamic_cast<NoxVector<T>*>(&x_in);
//   NoxVector<T>* r   = dynamic_cast<NoxVector<T>*>(&r_in);

//   // We cast to pointers so we can be sure that they succeeded
//   // by comparing the result against NULL.
//   libmesh_assert(jac != NULL); libmesh_assert(jac->mat() != NULL);
//   libmesh_assert(x   != NULL); libmesh_assert(x->vec()   != NULL);
//   libmesh_assert(r   != NULL); libmesh_assert(r->vec()   != NULL);
  
//   int ierr=0;
//   int n_iterations =0;

//   ierr = SNESSetFunction (_snes, r->vec(), __libmesh_nox_snes_residual, this);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   ierr = SNESSetJacobian (_snes, jac->mat(), jac->mat(), __libmesh_nox_snes_jacobian, this);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

//    // Have the Krylov subspace method use our good initial guess rather than 0
//    KSP ksp;	 
//    ierr = SNESGetKSP (_snes, &ksp);
//           CHKERRABORT(libMesh::COMM_WORLD,ierr);
		      
// //    ierr = KSPSetInitialGuessNonzero (ksp, NOX_TRUE);
// //           CHKERRABORT(libMesh::COMM_WORLD,ierr);
      	 
// // Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
// // the last one being a pointer to an int to hold the number of iterations required.
// # if NOX_VERSION_LESS_THAN(2,2,0)

//  ierr = SNESSolve (_snes, x->vec(), &n_iterations);
//         CHKERRABORT(libMesh::COMM_WORLD,ierr);

// // 2.2.x style	
// #elif NOX_VERSION_LESS_THAN(2,3,0)
	
//  ierr = SNESSolve (_snes, x->vec());
//         CHKERRABORT(libMesh::COMM_WORLD,ierr);

// // 2.3.x & newer style	
// #else
	
//  ierr = SNESSolve (_snes, NOX_NULL, x->vec());
//         CHKERRABORT(libMesh::COMM_WORLD,ierr);
	  	
// #endif

//   this->clear();
		 
//   // return the # of its. and the final residual norm.  Note that
//   // n_iterations may be zero for Nox versions 2.2.x and greater.

  return std::make_pair(1, 0.);
}




//------------------------------------------------------------------
// Explicit instantiations
template class NoxNonlinearSolver<Number>;
 


#endif // #ifdef HAVE_NOX
