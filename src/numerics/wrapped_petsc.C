// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// libMesh includes
#include "libmesh/wrapped_petsc.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_solver_exception.h"

// PETSc includes
#ifdef I
#define LIBMESH_SAW_I
#endif

#include <petscksp.h>
#include <petscvec.h>
#include <petscis.h>
#include <petscmat.h>
#include <petscviewer.h>
#include <petscpc.h>
#include <petscsnes.h>
#include <petscdm.h>
#include <petscsf.h>
// PetscSection got its own header file in 3.12, prior to that it was
// defined in petscis.h
#if !PETSC_VERSION_LESS_THAN(3, 12, 0)
#include <petscsection.h>
#endif

#ifndef LIBMESH_SAW_I
#undef I // Avoid complex.h contamination
#endif

namespace libMesh
{
// Specializations of the destroy() method for different PETSc classes
// This could also be macro-ized, but maybe it's not necessary?
template <>
void
WrappedPetsc<Vec>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, VecDestroy(&obj));
}

template <>
void
WrappedPetsc<KSP>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, KSPDestroy(&obj));
}

template <>
void
WrappedPetsc<IS>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, ISDestroy(&obj));
}

template <>
void
WrappedPetsc<Mat>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, MatDestroy(&obj));
}

template <>
void
WrappedPetsc<VecScatter>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, VecScatterDestroy(&obj));
}

template <>
void
WrappedPetsc<PetscViewer>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, PetscViewerDestroy(&obj));
}

template <>
void
WrappedPetsc<MatNullSpace>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, MatNullSpaceDestroy(&obj));
}

template <>
void
WrappedPetsc<DM>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, DMDestroy(&obj));
}

template <>
void
WrappedPetsc<MatPartitioning>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, MatPartitioningDestroy(&obj));
}

template <>
void
WrappedPetsc<SNES>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, SNESDestroy(&obj));
}

template <>
void
WrappedPetsc<PC>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, PCDestroy(&obj));
}

template <>
void
WrappedPetsc<PetscSection>::destroy()
{
  LibmeshPetscCallA(PETSC_COMM_SELF, PetscSectionDestroy(&obj));
}
}

#endif // LIBMESH_HAVE_PETSC
