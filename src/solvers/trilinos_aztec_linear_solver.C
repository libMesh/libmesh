// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifdef LIBMESH_TRILINOS_HAVE_AZTECOO


// C++ includes

// Local Includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/trilinos_aztec_linear_solver.h"
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/trilinos_epetra_vector.h"

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

      if (this->n_processors() == 1)
        this->_preconditioner_type = ILU_PRECOND;
      else
        this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
}



template <typename T>
void AztecLinearSolver<T>::init (const char * /*name*/)
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      _linear_solver = new AztecOO();

      set_solver_type();

      switch(this->_preconditioner_type)
        {
        case ILU_PRECOND:
          _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
          _linear_solver->SetAztecOption(AZ_subdomain_solve,AZ_ilu);
          break;

        case BLOCK_JACOBI_PRECOND:
          _linear_solver->SetAztecOption(AZ_precond,AZ_Jacobi);
          break;

        case ICC_PRECOND:
          _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
          _linear_solver->SetAztecOption(AZ_subdomain_solve,AZ_icc);
          break;

        case LU_PRECOND:
          _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
          _linear_solver->SetAztecOption(AZ_subdomain_solve,AZ_lu);
          break;

        default:
          _linear_solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
          _linear_solver->SetAztecOption(AZ_subdomain_solve,AZ_ilu);
        }
    }
}




template <typename T>
std::pair<unsigned int, Real>
AztecLinearSolver<T>::solve (SparseMatrix<T> & matrix_in,
                             SparseMatrix<T> & precond_in,
                             NumericVector<T> & solution_in,
                             NumericVector<T> & rhs_in,
                             const double tol,
                             const unsigned int m_its)
{
  LOG_SCOPE("solve()", "AztecLinearSolver");

  // Make sure the data passed in are really of Epetra types
  EpetraMatrix<T> * matrix   = cast_ptr<EpetraMatrix<T> *>(&matrix_in);
  EpetraMatrix<T> * precond  = cast_ptr<EpetraMatrix<T> *>(&precond_in);
  EpetraVector<T> * solution = cast_ptr<EpetraVector<T> *>(&solution_in);
  EpetraVector<T> * rhs      = cast_ptr<EpetraVector<T> *>(&rhs_in);

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

  // return the # of its. and the final residual norm.
  return std::make_pair(_linear_solver->NumIters(), _linear_solver->TrueResidual());
}



template <typename T>
std::pair<unsigned int, Real>
AztecLinearSolver<T>::solve (const ShellMatrix<T> &,
                             NumericVector<T> &,
                             NumericVector<T> &,
                             const double,
                             const unsigned int)
{
  libmesh_not_implemented();
}



template <typename T>
std::pair<unsigned int, Real>
AztecLinearSolver<T>::solve (const ShellMatrix<T> &,
                             const SparseMatrix<T> &,
                             NumericVector<T> &,
                             NumericVector<T> &,
                             const double,
                             const unsigned int)
{
  libmesh_not_implemented();
}



template <typename T>
void AztecLinearSolver<T>::get_residual_history(std::vector<double> & /* hist */)
{
  libmesh_not_implemented();
}




template <typename T>
Real AztecLinearSolver<T>::get_initial_residual()
{
  return _linear_solver->TrueResidual();
}



template <typename T>
void AztecLinearSolver<T>::print_converged_reason() const
{
  const double *status = _linear_solver->GetAztecStatus();

  switch (static_cast<int>(status[AZ_why]))
    {
    case AZ_normal :
      libMesh::out << "AztecOO converged.\n";
      break;
    case AZ_maxits :
      libMesh::out << "AztecOO failed to converge within maximum iterations.\n";
      break;
    case AZ_param :
      libMesh::out << "AztecOO failed to support a user-requested parameter.\n";
      break;
    case AZ_breakdown :
      libMesh::out << "AztecOO encountered numerical breakdown.\n";
      break;
    case AZ_loss :
      libMesh::out << "AztecOO encountered numerical loss of precision.\n";
      break;
    case AZ_ill_cond :
      libMesh::out << "AztecOO encountered an ill-conditioned GMRES Hessian.\n";
      break;
    default:
      libMesh::out << "AztecOO reported an unrecognized condition.\n";
      break;
    }
}



template <typename T>
LinearConvergenceReason AztecLinearSolver<T>::get_converged_reason() const
{
  const double *status = _linear_solver->GetAztecStatus();

  switch (static_cast<int>(status[AZ_why]))
    {
    case AZ_normal :
      return CONVERGED_RTOL_NORMAL;
    case AZ_maxits :
      return DIVERGED_ITS;
    }
  return DIVERGED_NULL;
}



template <typename T>
void AztecLinearSolver<T>::set_solver_type()
{
  switch (this->_solver_type)
    {
    case CG:
      _linear_solver->SetAztecOption(AZ_solver, AZ_cg); return;

    case CGS:
      _linear_solver->SetAztecOption(AZ_solver, AZ_cgs); return;

    case TFQMR:
      _linear_solver->SetAztecOption(AZ_solver, AZ_tfqmr); return;

    case BICGSTAB:
      _linear_solver->SetAztecOption(AZ_solver, AZ_bicgstab); return;

    case GMRES:
      _linear_solver->SetAztecOption(AZ_solver, AZ_gmres); return;

    default:
      libMesh::err << "ERROR:  Unsupported AztecOO Solver: "
                   << Utility::enum_to_string(this->_solver_type) << std::endl
                   << "Continuing with AztecOO defaults" << std::endl;
    }
}

//------------------------------------------------------------------
// Explicit instantiations
template class AztecLinearSolver<Number>;

} // namespace libMesh



#endif // #ifdef LIBMESH_TRILINOS_HAVE_AZTECOO
