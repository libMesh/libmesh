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

#ifdef LIBMESH_HAVE_NOX


// C++ includes

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/trilinos_nox_nonlinear_solver.h"
#include "libmesh/system.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/trilinos_preconditioner.h"

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
  explicit
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
                                 NOX::Epetra::Interface::Required::FillType /*fillType*/)
{
  START_LOG("residual()", "TrilinosNoxNonlinearSolver");

  NonlinearImplicitSystem &sys = _solver->system();

  EpetraVector<Number> X_global(*const_cast<Epetra_Vector *>(&x), sys.communicator()), R(r, sys.communicator());
  EpetraVector<Number>& X_sys = *libmesh_cast_ptr<EpetraVector<Number>*>(sys.solution.get());
  EpetraVector<Number>& R_sys = *libmesh_cast_ptr<EpetraVector<Number>*>(sys.rhs);

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);
  R.swap(R_sys);

  sys.get_dof_map().enforce_constraints_exactly(sys);
  sys.update();

  // Swap back
  X_global.swap(X_sys);
  R.swap(R_sys);

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

  STOP_LOG("residual()", "TrilinosNoxNonlinearSolver");

  return true;
}

bool Problem_Interface::computeJacobian(const Epetra_Vector & x,
                                        Epetra_Operator & jac)
{
  START_LOG("jacobian()", "TrilinosNoxNonlinearSolver");

  NonlinearImplicitSystem &sys = _solver->system();

  EpetraMatrix<Number> Jac(&dynamic_cast<Epetra_FECrsMatrix &>(jac));
  EpetraVector<Number>& X_sys = *libmesh_cast_ptr<EpetraVector<Number>*>(sys.solution.get());
  EpetraVector<Number> X_global(*const_cast<Epetra_Vector *>(&x), sys.communicator());

  // Set the dof maps
  Jac.attach_dof_map(sys.get_dof_map());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);

  sys.get_dof_map().enforce_constraints_exactly(sys);
  sys.update();

  X_global.swap(X_sys);

  //-----------------------------------------------------------------------------
  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  if (_solver->jacobian && _solver->jacobian_object)
  {
    libMesh::err << "ERROR: cannot specify both a function and object to compute the Jacobian!" << std::endl;
    libmesh_error();
  }

  if (_solver->matvec && _solver->residual_and_jacobian_object)
  {
    libMesh::err << "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!" << std::endl;
    libmesh_error();
  }
  //-----------------------------------------------------------------------------

  if      (_solver->jacobian != NULL)                     _solver->jacobian                                            (*sys.current_local_solution.get(), Jac, sys);
  else if (_solver->jacobian_object != NULL)              _solver->jacobian_object->jacobian                           (*sys.current_local_solution.get(), Jac, sys);
  else if (_solver->matvec   != NULL)                     _solver->matvec                                              (*sys.current_local_solution.get(), NULL, &Jac, sys);
  else if (_solver->residual_and_jacobian_object != NULL) _solver->residual_and_jacobian_object->residual_and_jacobian (*sys.current_local_solution.get(), NULL, &Jac, sys);
  else libmesh_error();

  Jac.close();
  X_global.close();

  STOP_LOG("jacobian()", "TrilinosNoxNonlinearSolver");

  return true;
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector & /*x*/, Epetra_RowMatrix & /*M*/)
{
//   cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
   throw 1;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector & x,
                                              Epetra_Operator & prec,
                                              Teuchos::ParameterList * /*p*/)
{
  START_LOG("preconditioner()", "TrilinosNoxNonlinearSolver");

  NonlinearImplicitSystem &sys = _solver->system();
  TrilinosPreconditioner<Number> & tpc = dynamic_cast<TrilinosPreconditioner<Number> &>(prec);

  EpetraMatrix<Number> Jac(dynamic_cast<Epetra_FECrsMatrix *>(tpc.mat()));
  EpetraVector<Number>& X_sys = *libmesh_cast_ptr<EpetraVector<Number>*>(sys.solution.get());
  EpetraVector<Number> X_global(*const_cast<Epetra_Vector *>(&x), sys.communicator());

  // Set the dof maps
  Jac.attach_dof_map(sys.get_dof_map());

  // Use the systems update() to get a good local version of the parallel solution
  X_global.swap(X_sys);

  sys.get_dof_map().enforce_constraints_exactly(sys);
  sys.update();

  X_global.swap(X_sys);

  //-----------------------------------------------------------------------------
  // if the user has provided both function pointers and objects only the pointer
  // will be used, so catch that as an error
  if (_solver->jacobian && _solver->jacobian_object)
  {
    libMesh::err << "ERROR: cannot specify both a function and object to compute the Jacobian!" << std::endl;
    libmesh_error();
  }

  if (_solver->matvec && _solver->residual_and_jacobian_object)
  {
    libMesh::err << "ERROR: cannot specify both a function and object to compute the combined Residual & Jacobian!" << std::endl;
    libmesh_error();
  }
  //-----------------------------------------------------------------------------

  if      (_solver->jacobian != NULL)                     _solver->jacobian                                            (*sys.current_local_solution.get(), Jac, sys);
  else if (_solver->jacobian_object != NULL)              _solver->jacobian_object->jacobian                           (*sys.current_local_solution.get(), Jac, sys);
  else if (_solver->matvec   != NULL)                     _solver->matvec                                              (*sys.current_local_solution.get(), NULL, &Jac, sys);
  else if (_solver->residual_and_jacobian_object != NULL) _solver->residual_and_jacobian_object->residual_and_jacobian (*sys.current_local_solution.get(), NULL, &Jac, sys);
  else libmesh_error();

  Jac.close();
  X_global.close();

  tpc.compute();

  STOP_LOG("preconditioner()", "TrilinosNoxNonlinearSolver");

  return true;
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

//       ierr = LibMeshSNESDestroy(_snes);
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
NoxNonlinearSolver<T>::solve (SparseMatrix<T>&  /* jac_in */,  // System Jacobian Matrix
                              NumericVector<T>& x_in,          // Solution vector
                              NumericVector<T>& /* r_in */,    // Residual vector
                              const double,                    // Stopping tolerance
                              const unsigned int)
{
  this->init ();

  if (this->user_presolve)
    this->user_presolve(this->system());

  EpetraVector<T> * x_epetra = libmesh_cast_ptr<EpetraVector<T>*>(&x_in);
  // Creating a Teuchos::RCP as they do in NOX examples does not work here - we get some invalid memory references
  // thus we make a local copy
  NOX::Epetra::Vector x(*x_epetra->vec());

  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
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

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", static_cast<int>(this->max_linear_iterations));
  lsParams.set("Tolerance", this->initial_linear_tolerance);
  lsParams.set("Output Frequency", 1);
  lsParams.set("Size of Krylov Subspace", 1000);

  //create linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq(_interface);
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys;
  Teuchos::RCP<Epetra_Operator> pc;

  if (this->jacobian || this->jacobian_object || this->residual_and_jacobian_object)
  {
    if(this->_preconditioner)
    {
      // PJNFK
      lsParams.set("Preconditioner", "User Defined");

      TrilinosPreconditioner<Number> * trilinos_pc = libmesh_cast_ptr<TrilinosPreconditioner<Number> *>(this->_preconditioner);
      pc = Teuchos::rcp(trilinos_pc);

      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec(_interface);
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iPrec, pc, x));
    }
    else
    {
      lsParams.set("Preconditioner", "None");
//      lsParams.set("Preconditioner", "Ifpack");
//      lsParams.set("Preconditioner", "AztecOO");

      // full jacobian
      NonlinearImplicitSystem & sys = _interface->_solver->system();
      EpetraMatrix<Number> & jacSys = *libmesh_cast_ptr<EpetraMatrix<Number>*>(sys.matrix);
      Teuchos::RCP<Epetra_RowMatrix> jacMat = Teuchos::rcp(jacSys.mat());

      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac(_interface);
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iJac, jacMat, x));
    }
  }
  else
  {
    // matrix free
    Teuchos::RCP<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, iReq, x));

    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac(MF);
    linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iJac, MF, x));
  }

  //create group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, x, linSys));
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

  Teuchos::RCP<Teuchos::ParameterList> finalPars = nlParamsPtr;

  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = solver->solve();
  this->converged = (status == NOX::StatusTest::Converged);

  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const NOX::Epetra::Vector& noxFinalSln = dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX());

  *x_epetra->vec() = noxFinalSln.getEpetraVector();
  x_in.close();

  Real residual_norm = finalGroup.getNormF();
  unsigned int total_iters = solver->getNumIterations();
  _n_linear_iterations = finalPars->sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations", -1);

  // do not let Trilinos to deallocate what we allocated
  pc.release();
  iReq.release();

  return std::make_pair(total_iters, residual_norm);
}

template <typename T>
int
NoxNonlinearSolver<T>::get_total_linear_iterations()
{
  return _n_linear_iterations;
}



//------------------------------------------------------------------
// Explicit instantiations
template class NoxNonlinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_HAVE_NOX
