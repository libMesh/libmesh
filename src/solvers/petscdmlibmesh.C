// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/petsc_macro.h"
// This only works with petsc-3.3 and above.
#if !PETSC_VERSION_LESS_THAN(3,3,0)

#if !PETSC_RELEASE_LESS_THAN(3,6,0)
# include <petsc/private/petscimpl.h>
#else
# include <petsc-private/petscimpl.h>
#endif

#include "libmesh/petscdmlibmesh.h"

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshSetSystem"
PetscErrorCode DMlibMeshSetSystem(DM dm, libMesh::NonlinearImplicitSystem & sys)
{
  PetscErrorCode (*f)(DM,libMesh::NonlinearImplicitSystem &) = nullptr;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
#if PETSC_RELEASE_LESS_THAN(3,4,0)
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshSetSystem_C",(PetscVoidFunction*)&f);CHKERRQ(ierr);
#else
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshSetSystem_C",&f);CHKERRQ(ierr);
#endif
  if (!f) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, "DM has no implementation for DMlibMeshSetSystem");
  ierr = (*f)(dm,sys);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetSystem"
PetscErrorCode DMlibMeshGetSystem(DM dm, libMesh::NonlinearImplicitSystem *& sys)
{
  PetscErrorCode (*f)(DM,libMesh::NonlinearImplicitSystem *&) = nullptr;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
#if PETSC_RELEASE_LESS_THAN(3,4,0)
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshGetSystem_C",(PetscVoidFunction*)&f);CHKERRQ(ierr);
#else
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshGetSystem_C",&f);CHKERRQ(ierr);
#endif
  if (!f) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, "DM has no implementation for DMlibMeshGetSystem");
  ierr = (*f)(dm,sys);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
