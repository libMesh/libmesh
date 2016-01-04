#ifndef PETSCDMLIBMESH_H
#define PETSCDMLIBMESH_H
// This is intended as a kind of public header for DMlibMesh
// (does it get installed?  should it? is this public, but
// developer-facing only?)

#include "libmesh/petsc_macro.h"
// This only works with petsc-3.3 and above.
#if !PETSC_VERSION_LESS_THAN(3,3,0)

#include <petscdm.h>
#include "libmesh/nonlinear_implicit_system.h"

/**
 * Any functional implementation of the DMlibMesh API must compose the following functions with the DM object.
 * (See PETSc documentation on PetscObjectComposeFunction(), a polymorphism mechanism.)
 * The following functions are called in PetscNonlinear Solver (others can be called by users):
 * DMlibMeshSetSystem(), DMlibMeshGetSystem()
 *
 * Any implementation needs to register its creation routine, DMCreate_libMesh, with PETSc using DMRegister().
 */
PETSC_EXTERN PetscErrorCode DMlibMeshSetSystem(DM,libMesh::NonlinearImplicitSystem &);
PETSC_EXTERN PetscErrorCode DMlibMeshGetSystem(DM,libMesh::NonlinearImplicitSystem *&);


#define DMLIBMESH "libmesh"
EXTERN_C_BEGIN
PETSC_EXTERN PetscErrorCode DMCreate_libMesh(DM);
EXTERN_C_END

#endif // #if !PETSC_VERSION_LESS_THAN(3,3,0)
#endif // #ifdef PETSCDMLIBEMSH_H
