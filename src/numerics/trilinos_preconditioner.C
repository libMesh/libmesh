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

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA

// Local Includes
#include "libmesh/trilinos_preconditioner.h"
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_preconditioner_type.h"

#include "libmesh/ignore_warnings.h"
#ifdef LIBMESH_TRILINOS_HAVE_IFPACK
#include "Ifpack.h"
#include "Ifpack_DiagPreconditioner.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#endif

#ifdef LIBMESH_TRILINOS_HAVE_ML
#include "ml_MultiLevelPreconditioner.h"
#endif
#include "libmesh/restore_warnings.h"

namespace libMesh
{

template <typename T>
void TrilinosPreconditioner<T>::apply(const NumericVector<T> & /* x */,
                                      NumericVector<T> & /* y */ )
{
}




template <typename T>
void TrilinosPreconditioner<T>::init ()
{
  if (!this->_matrix)
    libmesh_error_msg("ERROR: No matrix set for PetscPreconditioner, but init() called");

  // Clear the preconditioner in case it has been created in the past
  if (!this->_is_initialized)
    {
      EpetraMatrix<T> * matrix = cast_ptr<EpetraMatrix<T> *, SparseMatrix<T>>(this->_matrix);
      _mat = matrix->mat();
    }

  set_preconditioner_type(this->_preconditioner_type);

  this->_is_initialized = true;
}


template <typename T>
void
TrilinosPreconditioner<T>::set_params(Teuchos::ParameterList & list)
{
  _param_list = list;
}


template <typename T>
void
TrilinosPreconditioner<T>::compute()
{
#ifdef LIBMESH_TRILINOS_HAVE_IFPACK
  Ifpack_Preconditioner * ifpack = libmesh_nullptr;
#endif

#ifdef LIBMESH_TRILINOS_HAVE_ML
  ML_Epetra::MultiLevelPreconditioner * ml = libmesh_nullptr;
#endif

  switch (this->_preconditioner_type)
    {
#ifdef LIBMESH_TRILINOS_HAVE_IFPACK
      // IFPACK preconditioners
    case ILU_PRECOND:
    case SOR_PRECOND:
      ifpack = dynamic_cast<Ifpack_Preconditioner *>(_prec);
      ifpack->Compute();
      break;
#endif

#ifdef LIBMESH_TRILINOS_HAVE_ML
      // ML preconditioners
    case AMG_PRECOND:
      ml = dynamic_cast<ML_Epetra::MultiLevelPreconditioner *>(_prec);
      ml->ComputePreconditioner();
      break;
#endif

    default:
      // If we made it here, there were no TrilinosPreconditioners
      // active, so that's probably an error.
      libmesh_error_msg("ERROR: No valid TrilinosPreconditioners available!");
      break;
    }
}


template <typename T>
void
TrilinosPreconditioner<T>::set_preconditioner_type (const PreconditionerType & preconditioner_type)
{
#ifdef LIBMESH_TRILINOS_HAVE_IFPACK
  Ifpack_Preconditioner * pc = libmesh_nullptr;
#endif

#ifdef LIBMESH_TRILINOS_HAVE_ML
  ML_Epetra::MultiLevelPreconditioner * ml = libmesh_nullptr;
#endif

  switch (preconditioner_type)
    {
    case IDENTITY_PRECOND:
      //    pc = new Ifpack_DiagPreconditioner();
      break;

    case CHOLESKY_PRECOND:
      break;

    case ICC_PRECOND:
      break;

#ifdef LIBMESH_TRILINOS_HAVE_IFPACK
    case ILU_PRECOND:
      pc = new Ifpack_ILU(_mat);
      pc->SetParameters(_param_list);
      pc->Initialize();
      _prec = pc;
      break;
#endif

    case LU_PRECOND:
      break;

    case ASM_PRECOND:
      break;

    case JACOBI_PRECOND:
      break;

    case BLOCK_JACOBI_PRECOND:
      break;

    case SOR_PRECOND:
      break;

    case EISENSTAT_PRECOND:
      break;

#ifdef LIBMESH_TRILINOS_HAVE_ML
    case AMG_PRECOND:
      ml = new ML_Epetra::MultiLevelPreconditioner(*_mat, _param_list, false);;
      _prec = ml;
      break;
#endif

    default:
      libmesh_error_msg("ERROR:  Unsupported Trilinos Preconditioner: " << preconditioner_type << "\nContinuing with Trilinos defaults");
    }

}


template <typename T>
int
TrilinosPreconditioner<T>::SetUseTranspose(bool UseTranspose)
{
  return _prec->SetUseTranspose(UseTranspose);
}

template <typename T>
int
TrilinosPreconditioner<T>::Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
  return _prec->Apply(X, Y);
}

template <typename T>
int
TrilinosPreconditioner<T>::ApplyInverse(const Epetra_MultiVector & r, Epetra_MultiVector & z) const
{
  return _prec->ApplyInverse(r, z);
}

template <typename T>
double
TrilinosPreconditioner<T>::NormInf() const
{
  return _prec->NormInf();
}

template <typename T>
const char *
TrilinosPreconditioner<T>::Label() const
{
  return _prec->Label();
}

template <typename T>
bool
TrilinosPreconditioner<T>::UseTranspose() const
{
  return _prec->UseTranspose();
}

template <typename T>
bool
TrilinosPreconditioner<T>::HasNormInf() const
{
  return _prec->HasNormInf();
}

template <typename T>
const Epetra_Comm &
TrilinosPreconditioner<T>::Comm() const
{
  return _prec->Comm();
}

template <typename T>
const Epetra_Map &
TrilinosPreconditioner<T>::OperatorDomainMap() const
{
  return _prec->OperatorDomainMap();
}

template <typename T>
const Epetra_Map &
TrilinosPreconditioner<T>::OperatorRangeMap() const
{
  return _prec->OperatorRangeMap();
}

//------------------------------------------------------------------
// Explicit instantiations
template class TrilinosPreconditioner<Number>;

} // namespace libMesh

#endif // LIBMESH_TRILINOS_HAVE_EPETRA
