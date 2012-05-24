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



#ifndef __petsc_dm_nonlinear_solver_h__
#define __petsc_dm_nonlinear_solver_h__

#include "petsc_macro.h"

// This only works with a recent petsc-dev (post petsc-3.2).
// Replace with a test for petsc-3.3 after it's released.
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,2,0) && !PETSC_VERSION_RELEASE

// Petsc include files.
EXTERN_C_BEGIN
#include <petscsnes.h>
#include <petsc-private/dmimpl.h>

#define DMLIBMESH "libmesh"
//typedef struct
//{
//  NonlinearImplicitSystem* sys;
//} DM_libMesh;

struct DM_libMesh
{
  NonlinearImplicitSystem* sys;
};

EXTERN_C_END



// Local includes
#include "petsc_nonlinear_solver.h"

namespace libMesh
{
  // Register the DM_libMesh constructor with PETSc so that we can do DMSetType(dm, DMLIBMESH);
  void PetscDMRegister();

  // Set a libMesh NonlinearImplicitSystem on a DM_libMesh
  void PetscDMSetSystem(DM, NonlinearImplicitSystem&);

  // Get a libMesh NonlinearImplicitSystem from a DM_libMesh, or throw error if none has been set
  void PetscDMGetSystem(DM, NonlinearImplicitSystem*&);





/**
 * This class provides an interface to PETSc DM-based
 * iterative solvers that is compatible with the \p libMesh
 * \p NonlinearSolver<>
 *
 * @author Dmitry Karpeev, 2012
 */
template <typename T>
class PetscDMNonlinearSolver : public PetscNonlinearSolver<T>
{
public:
  /**
   * The type of system
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   *  Constructor. Initializes PETSc data structures
   */
  explicit
  PetscDMNonlinearSolver (sys_type& system);

  /**
   * Destructor.
   */
  ~PetscDMNonlinearSolver ();


  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();

  /**
   * Call the Petsc solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real> solve (SparseMatrix<T>&,    // System Jacobian Matrix
					       NumericVector<T>&,   // Solution vector
					       NumericVector<T>&,   // Residual vector
					       const double,        // Stopping tolerance
					       const unsigned int); // N. Iterations

};




} // namespace libMesh


#endif // #if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,2,0) && !PETSC_VERSION_RELEASE
#endif // #ifdef __petsc_dm_nonlinear_solver_h__
