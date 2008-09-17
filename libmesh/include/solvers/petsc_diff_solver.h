// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __petsc_diff_solver_h__
#define __petsc_diff_solver_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "diff_solver.h"

#ifdef HAVE_PETSC

#include "petsc_macro.h"

EXTERN_C_FOR_PETSC_BEGIN
# include <petscsnes.h>
EXTERN_C_FOR_PETSC_END

/**
 * This class defines a solver which uses a PETSc SNES 
 * context to handle a DifferentiableSystem
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2008
 */

// ------------------------------------------------------------
// Solver class definition
class PetscDiffSolver : public DiffSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  PetscDiffSolver (sys_type& system);
  
  /**
   * Destructor.
   */
  virtual ~PetscDiffSolver ();

  typedef DiffSolver Parent;

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh.
   */
  virtual void reinit ();

  /**
   * This method performs a solve.  What occurs in
   * this method will depend on the type of solver.  See
   * the PETSc documentation for more details.
   */
  virtual unsigned int solve ();

protected:

  /**
   * Nonlinear solver context
   */
  SNES _snes;
};

#endif // #ifdef HAVE_PETSC

#endif // #define __petsc_diff_solver_h__
