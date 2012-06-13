#ifndef _petscdmlibmesh
#define _petscdmlibmesh

#include "petsc_macro.h"
// This only works with petsc-3.3 and above.

#if !PETSC_VERSION_LESS_THAN(3,3,0)
// libMesh include
#include <nonlinear_implicit_system.h>

#include <petscdm.h>
#define DMLIBMESH "libmesh"

extern PetscErrorCode DMCreateLibMesh(MPI_Comm,NonlinearImplicitSystem&,DM*);
extern PetscErrorCode DMLibMeshSetSystem(DM,NonlinearImplicitSystem&);
extern PetscErrorCode DMLibMeshGetSystem(DM,NonlinearImplicitSystem*&);
extern PetscErrorCode DMLibMeshGetBlocks(DM,PetscInt*,char***);
extern PetscErrorCode DMLibMeshGetVariables(DM,PetscInt*,char***);
extern PetscErrorCode DMLibMeshCreateFieldDecompositionDM(DM,PetscInt,PetscInt*,const char***,DM*);
extern PetscErrorCode DMLibMeshCreateDomainDecompositionDM(DM,PetscInt,PetscInt*,const char***,DM*);

#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
#endif // #ifdef _petscdmlibmesh
