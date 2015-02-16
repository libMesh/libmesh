#include "libmesh/petsc_macro.h"
// This only works with petsc-3.3 and above.
#if !PETSC_VERSION_LESS_THAN(3,3,0)

#include <petsc-private/petscimpl.h>
#include "libmesh/petscdmlibmesh.h"

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshSetSystem"
PetscErrorCode DMlibMeshSetSystem(DM dm, libMesh::NonlinearImplicitSystem& sys)
{
  PetscErrorCode (*f)(DM,libMesh::NonlinearImplicitSystem&) = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
#if PETSC_RELEASE_LESS_THAN(3,4,0)
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshSetSystem_C",(PetscVoidFunction*)&f);CHKERRQ(ierr);
#else
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshSetSystem_C",&f);CHKERRQ(ierr);
#endif
  if(!f) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, "DM has no implementation for DMlibMeshSetSystem");
  ierr = (*f)(dm,sys);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "DMlibMeshGetSystem"
PetscErrorCode DMlibMeshGetSystem(DM dm, libMesh::NonlinearImplicitSystem*& sys)
{
  PetscErrorCode (*f)(DM,libMesh::NonlinearImplicitSystem*&) = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
#if PETSC_RELEASE_LESS_THAN(3,4,0)
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshGetSystem_C",(PetscVoidFunction*)&f);CHKERRQ(ierr);
#else
  ierr = PetscObjectQueryFunction((PetscObject)dm,"DMlibMeshGetSystem_C",&f);CHKERRQ(ierr);
#endif
  if(!f) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, "DM has no implementation for DMlibMeshGetSystem");
  ierr = (*f)(dm,sys);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
