#ifndef _petscdmlibmesh
#define _petscdmlibmesh

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


#endif // #ifdef _petscdmlibmesh
