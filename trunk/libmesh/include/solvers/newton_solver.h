// $Id: newton_solver.h,v 1.1 2006-06-05 00:32:23 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __newton_solver_h__
#define __newton_solver_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "linear_solver.h"
#include "reference_counted_object.h"
#include "diff_solver.h"

/**
 * This class defines a solver which uses the default
 * libMesh linear solver in a quasiNewtonto handle a 
 * DifferentiableSystem
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// Solver class definition
class NewtonSolver : public DiffSolver
{
public:
  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  NewtonSolver (sys_type& system);
  
  /**
   * Destructor.
   */
  virtual ~NewtonSolver ();

  typedef DiffSolver Parent;

  /**
   * This method performs a solve.  What occurs in
   * this method will depend on the type of solver.  See
   * the subclasses for more details.
   */
  virtual void solve ();

protected:

  /**
   * The \p LinearSolver defines the interface used to
   * solve the linear_implicit system.  This class handles all the
   * details of interfacing with various linear algebra packages
   * like PETSc or LASPACK.
   */
  AutoPtr<LinearSolver<Number> > linear_solver;

};



#endif // #define __newton_solver_h__
