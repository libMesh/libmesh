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



#ifndef LIBMESH_PETSC_DIFF_SOLVER_H
#define LIBMESH_PETSC_DIFF_SOLVER_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/diff_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_dm_wrapper.h"

// PETSc includes
# include <petscsnes.h>

namespace libMesh
{

/**
 * This class defines a solver which uses a PETSc SNES
 * context to handle a DifferentiableSystem
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2008
 */
class PetscDiffSolver : public DiffSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  PetscDiffSolver (sys_type & system);

  /**
   * Destructor.
   */
  virtual ~PetscDiffSolver ();

  typedef DiffSolver Parent;

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh.
   */
  virtual void reinit () override;

  /**
   * The initialization function.
   */
  virtual void init () override;

  /**
   * The clear function. Called if we reinit or when we destroy this object.
   */
  void clear ();

  /**
   * This method performs a solve.  What occurs in
   * this method will depend on the PETSc SNES settings.
   * See the PETSc documentation for more details.
   */
  virtual unsigned int solve () override;

protected:

  /**
   * Nonlinear solver context
   */
  SNES _snes;

  /**
   * Wrapper object for interacting with PetscDM
   */
  PetscDMWrapper _dm_wrapper;

private:

  /**
   * Common helper function to setup PETSc data structures
   */
  void setup_petsc_data();
};

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_PETSC

#endif // LIBMESH_PETSC_DIFF_SOLVER_H
