// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifdef LIBMESH_HAVE_PETSC

// Local Includes
#include "libmesh/petsc_preconditioner.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_preconditioner_type.h"

namespace libMesh
{

template <typename T>
PetscPreconditioner<T>::PetscPreconditioner (const libMesh::Parallel::Communicator & comm_in) :
  Preconditioner<T>(comm_in)
{}



template <typename T>
void PetscPreconditioner<T>::apply(const NumericVector<T> & x, NumericVector<T> & y)
{
  PetscVector<T> & x_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(x));
  PetscVector<T> & y_pvec = cast_ref<PetscVector<T> &>(const_cast<NumericVector<T> &>(y));

  Vec x_vec = x_pvec.vec();
  Vec y_vec = y_pvec.vec();

  LibmeshPetscCall(PCApply(_pc, x_vec, y_vec));
}




template <typename T>
void PetscPreconditioner<T>::init ()
{
  libmesh_error_msg_if(!this->_matrix, "ERROR: No matrix set for PetscPreconditioner, but init() called");

  // Clear the preconditioner in case it has been created in the past
  if (!this->_is_initialized)
    {
      // Should probably use PCReset(), but it's not working at the moment so we'll destroy instead
      if (_pc)
        _pc.destroy();

      LibmeshPetscCall(PCCreate(this->comm().get(), _pc.get()));

      auto pmatrix = cast_ptr<PetscMatrixBase<T> *>(this->_matrix);
      _mat = pmatrix->mat();
    }

  LibmeshPetscCall(PCSetOperators(_pc, _mat, _mat));

  // Set the PCType.  Note: this used to be done *before* the call to
  // PCSetOperators(), and only when !_is_initialized, but
  // 1.) Some preconditioners (those employing sub-preconditioners,
  // for example) have to call PCSetUp(), and can only do this after
  // the operators have been set.
  // 2.) It should be safe to call set_petsc_preconditioner_type()
  // multiple times.
  set_petsc_preconditioner_type(this->_preconditioner_type, *_pc);

  this->_is_initialized = true;
}



template <typename T>
void PetscPreconditioner<T>::clear()
{
  // Calls custom deleter
  _pc.destroy();
}



template <typename T>
PC PetscPreconditioner<T>::pc()
{
  return _pc;
}



template <typename T>
void PetscPreconditioner<T>::set_petsc_preconditioner_type (const PreconditionerType & preconditioner_type, PC & pc)
{
  // get the communicator from the PETSc object
  Parallel::communicator comm;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)pc, & comm);
  if (ierr != LIBMESH_PETSC_SUCCESS)
    libmesh_error_msg("Error retrieving communicator");
  Parallel::Communicator communicator(comm);

  switch (preconditioner_type)
    {
    case IDENTITY_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCNONE)));
      break;

    case CHOLESKY_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCCHOLESKY)));
      break;

    case ICC_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCICC)));
      break;

    case ILU_PRECOND:
      {
        // In serial, just set the ILU preconditioner type
        if (communicator.size())
          LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCILU)));
        else
          {
            // But PETSc has no truly parallel ILU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign ILU sub-preconditioners.
            LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCBJACOBI)));

            // Set ILU as the sub preconditioner type
            set_petsc_subpreconditioner_type(PCILU, pc);
          }
        break;
      }

    case LU_PRECOND:
      {
        // In serial, just set the LU preconditioner type
        if (communicator.size())
          LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCLU)));
        else
          {
            // But PETSc has no truly parallel LU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign LU sub-preconditioners.
            LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCBJACOBI)));

            // Set ILU as the sub preconditioner type
            set_petsc_subpreconditioner_type(PCLU, pc);
          }
        break;
      }

    case ASM_PRECOND:
      {
        // In parallel, I think ASM uses ILU by default as the sub-preconditioner...
        // I tried setting a different sub-preconditioner here, but apparently the matrix
        // is not in the correct state (at this point) to call PCSetUp().
        LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCASM)));
        break;
      }

    case JACOBI_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCJACOBI)));
      break;

    case BLOCK_JACOBI_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCBJACOBI)));
      break;

    case SOR_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCSOR)));
      break;

    case EISENSTAT_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCEISENSTAT)));
      break;

    case AMG_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCHYPRE)));
      break;

    case SVD_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCSVD)));
      break;

    case USER_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCMAT)));
      break;

    case SHELL_PRECOND:
      LibmeshPetscCall2(communicator, PCSetType (pc, const_cast<KSPType>(PCSHELL)));
      break;

    default:
      libMesh::err << "ERROR:  Unsupported PETSC Preconditioner: "
                   << preconditioner_type       << std::endl
                   << "Continuing with PETSC defaults" << std::endl;
    }

  // Set additional options if we are doing AMG and
  // HYPRE is available
#ifdef LIBMESH_HAVE_PETSC_HYPRE
  if (preconditioner_type == AMG_PRECOND)
    LibmeshPetscCall2(communicator, PCHYPRESetType(pc, "boomeramg"));
#endif

  // Let the commandline override stuff
  LibmeshPetscCall2(communicator, PCSetFromOptions(pc));
}


template <typename T>
void PetscPreconditioner<T>::set_petsc_subpreconditioner_type(const PCType type, PC & pc)
{
  // get the communicator from the PETSc object
  Parallel::communicator comm;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)pc, & comm);
  if (ierr != LIBMESH_PETSC_SUCCESS)
    libmesh_error_msg("Error retrieving communicator");
  Parallel::Communicator communicator(comm);

  // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
  // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
  //
  // "Object is in wrong state!"
  // "Matrix must be set first."
  //
  // error messages...
  LibmeshPetscCall2(communicator, PCSetUp(pc));

  // To store array of local KSP contexts on this processor
  KSP * subksps;

  // the number of blocks on this processor
  PetscInt n_local;

  // The global number of the first block on this processor.
  // This is not used, so we just pass null instead.
  // int first_local;

  // Fill array of local KSP contexts
  LibmeshPetscCall2(communicator, PCBJacobiGetSubKSP(pc, &n_local, LIBMESH_PETSC_NULLPTR,
                                                     &subksps));

  // Loop over sub-ksp objects, set ILU preconditioner
  for (PetscInt i=0; i<n_local; ++i)
    {
      // Get pointer to sub KSP object's PC
      PC subpc;
      LibmeshPetscCall2(communicator, KSPGetPC(subksps[i], &subpc));

      // Set requested type on the sub PC
      LibmeshPetscCall2(communicator, PCSetType(subpc, type));
    }
}




//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT PetscPreconditioner<Number>;

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC
