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

#ifdef LIBMESH_HAVE_NOX


// C++ includes

// Local Includes
#include "dof_map.h"
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
#include "NOX_Epetra_Vector.H"

namespace libMesh
{

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
  
  NoxNonlinearSolver<Number> * _solver;
};

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(NoxNonlinearSolver<Number> * solver) :
  _solver(solver)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& r,
                       NOX::Epetra::Interface::Required::FillType)
//bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& r,
//                       NOX::Epetra::Interface::Required::FillType fillType)
{
  NonlinearImplicitSystem &sys = _solver->system();

  EpetraVector<Number> X_global(*const_cast<Epetra_Vector *>(&x)), R(r);
  EpetraVector<Number>& X_sys = *libmesh_cast_ptr<EpetraVector<Number>*>(sys.solution.get());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);

  sys.get_dof_map().enforce_constraints_exactly(sys);
  sys.update();

  X_global.swap(X_sys);
  
  R.zero();
  
  //-----------------------------------------------------------------------------
   // if the user has provided both function pointers and objects only the pointer
   // will be used, so catch that as an error

    if (_solver->residual && _solver->residual_object)
      {
	libMesh::err << "ERROR: cannot specifiy both a function and object to compute the Residual!" << std::endl;
	libmesh_error();
      }
    
    if (_solver->matvec && _solver->residual_and_jacobian_object)
      {
	libMesh::err << "ERROR: cannot specifiy both a function and object to compute the combined Residual & Jacobian!" << std::endl;
	libmesh_error();
      }
    //-----------------------------------------------------------------------------

    if      (_solver->residual != NULL)                     _solver->residual                                            (*sys.current_local_solution.get(), R, sys);
    else if (_solver->residual_object != NULL)              _solver->residual_object->residual                           (*sys.current_local_solution.get(), R, sys);
    else if (_solver->matvec   != NULL)                     _solver->matvec                                              (*sys.current_local_solution.get(), &R, NULL, sys);
    else if (_solver->residual_and_jacobian_object != NULL) _solver->residual_and_jacobian_object->residual_and_jacobian (*sys.current_local_solution.get(), &R, NULL, sys);
    else return false;

  R.close();
  X_global.close();
  return true;
}

bool Problem_Interface::computeJacobian(const Epetra_Vector&,
                                        Epetra_Operator&)
//bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
//                                        Epetra_Operator& Jac)
{
  throw 1;
  
  // return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0);
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector&, Epetra_RowMatrix&)
//bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
//   cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
   throw 1;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector&,
                                              Epetra_Operator&,
                                              Teuchos::ParameterList*)
//bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
//                                              Epetra_Operator& Prec,
//                                              Teuchos::ParameterList* p)
{
//   cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
   throw 1;
}


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

  EpetraVector<T> * x_epetra = libmesh_cast_ptr<EpetraVector<T>*>(&x_in);
  EpetraVector<T> * r_epetra = libmesh_cast_ptr<EpetraVector<T>*>(&r_in);

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
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
                  NOX::Utils::OuterIteration +
                  NOX::Utils::OuterIterationStatusTest +
                  NOX::Utils::InnerIteration +
                  NOX::Utils::LinearSolverDetails +
                  NOX::Utils::Parameters +
                  NOX::Utils::Details +
                  NOX::Utils::Warning);
  
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
  lsParams.set("Max Iterations", static_cast<int>(this->max_linear_iterations));
  lsParams.set("Tolerance", this->initial_linear_tolerance);  lsParams.set("Output Frequency", 1);	 
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

  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(this->absolute_residual_tolerance, NOX::StatusTest::NormF::Unscaled));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, this->relative_residual_tolerance));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(this->max_nonlinear_iterations));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteval =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue());
  Teuchos::RCP<NOX::StatusTest::NormUpdate> normupdate =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(this->absolute_step_tolerance));
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(absresid);
  combo->addStatusTest(relresid);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteval);
  combo->addStatusTest(normupdate);
  
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;
  status = solver->solve();
  
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const NOX::Abstract::Vector& finalSolution = finalGroup.getX();
//  const Epetra_Vector& finalSolution =
//    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  x = finalSolution;

  x_in.close();

  return std::make_pair(1, 0.);
}




//------------------------------------------------------------------
// Explicit instantiations
template class NoxNonlinearSolver<Number>;

} // namespace libMesh
 


#endif // #ifdef LIBMESH_HAVE_NOX
