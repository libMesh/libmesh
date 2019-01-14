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



#ifndef LIBMESH_TRILINOS_NOX_NONLINEAR_SOLVER_H
#define LIBMESH_TRILINOS_NOX_NONLINEAR_SOLVER_H

// It is possible to build Trilinos with NOX, but without Epetra.  The
// way that libmesh uses NOX requires Epetra, however, so we have to
// check for both.
#if defined(LIBMESH_TRILINOS_HAVE_NOX) && defined(LIBMESH_TRILINOS_HAVE_EPETRA)

// Local includes
#include "libmesh/nonlinear_solver.h"

// Trilinos includes
#include "libmesh/ignore_warnings.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#include "NOX.H"
#include "libmesh/restore_warnings.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
class Problem_Interface;


/**
 * This class provides an interface to nox
 * iterative solvers that is compatible with the \p libMesh
 * \p NonlinearSolver<>
 *
 * \author Chris Newman
 * \date 2008
 */
template <typename T>
class NoxNonlinearSolver : public NonlinearSolver<T>
{
public:
  /**
   * The type of system
   */
  typedef NonlinearImplicitSystem sys_type;

  /**
   *  Constructor. Initializes Nox data structures
   */
  explicit
  NoxNonlinearSolver (sys_type & system);

  /**
   * Destructor.
   */
  virtual ~NoxNonlinearSolver ();

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () override;

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init (const char * name = nullptr) override;

  /**
   * Call the Nox solver.  It calls the method below, using the
   * same matrix for the system and preconditioner matrices.
   */
  virtual std::pair<unsigned int, Real>
  solve (SparseMatrix<T> &,                     // System Jacobian Matrix
         NumericVector<T> &,                    // Solution vector
         NumericVector<T> &,                    // Residual vector
         const double,                          // Stopping tolerance
         const unsigned int) override;  // N. Iterations
  /**
   * Get the total number of linear iterations done in the last solve
   */
  virtual int get_total_linear_iterations() override;

  /**
   * \returns The current nonlinear iteration number if called
   * *during* the solve(), for example by the user-specified residual
   * or Jacobian function.
   *
   * \note Not currently implemented for the NoxNonlinearSolver.
   */
  virtual unsigned get_current_nonlinear_iteration_number() const override
  { libmesh_not_implemented(); return 0; }

private:

  /**
   * Nonlinear solver context
   */
  NOX::Solver::Generic * _solver;

  /**
   * Solver interface
   */
  Problem_Interface * _interface;

  /**
   * Stores the total number of linear iterations from the last solve.
   */
  int _n_linear_iterations;
};


/*----------------------- functions ----------------------------------*/
template <typename T>
inline
NoxNonlinearSolver<T>::NoxNonlinearSolver (sys_type & system) :
  NonlinearSolver<T>(system),
  _solver(nullptr),
  _interface(nullptr),
  _n_linear_iterations(0)
{
}



template <typename T>
inline
NoxNonlinearSolver<T>::~NoxNonlinearSolver ()
{
  this->clear ();
}


} // namespace libMesh


#endif // LIBMESH_TRILINOS_HAVE_NOX && LIBMESH_TRILINOS_HAVE_EPETRA
#endif // LIBMESH_TRILINOS_NOX_NONLINEAR_SOLVER_H
