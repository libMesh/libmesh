// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_LINEAR_IMPLICIT_SYSTEM_H
#define LIBMESH_LINEAR_IMPLICIT_SYSTEM_H

// Local Includes
#include "libmesh/implicit_system.h"

// C++ includes
#include <cstddef>

namespace libMesh
{


// Forward Declarations
template <typename T> class LinearSolver;
template <typename T> class ShellMatrix;


/**
 * \brief Manages consistently variables, degrees of freedom, coefficient
 * vectors, matrices and linear solvers for implicit systems.
 *
 * An implicit system is a system that requires the solution of a
 * system of equations. This class has the ability to create and use a
 * linear solver to solve the system.
 *
 * The matrix LinearImplicitSystem::matrix and the vector
 * LinearImplicitSystem::rhs should be filled during assembly.
 *
 * \note Additional vectors/matrices can be added via parent class
 * interfaces.
 *
 * \author Benjamin Kirk
 * \date 2005
 */
class LinearImplicitSystem : public ImplicitSystem
{
public:

  /**
   * Constructor.
   */
  LinearImplicitSystem (EquationSystems & es,
                        const std::string & name,
                        const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - The destructor is defaulted out-of-line.
   */
  LinearImplicitSystem (const LinearImplicitSystem &) = delete;
  LinearImplicitSystem & operator= (const LinearImplicitSystem &) = delete;
  LinearImplicitSystem (LinearImplicitSystem &&) = default;
  LinearImplicitSystem & operator= (LinearImplicitSystem &&) = delete;
  virtual ~LinearImplicitSystem ();

  /**
   * The type of system.
   */
  typedef LinearImplicitSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef ImplicitSystem Parent;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () override;

  /**
   * Initializes new data members of the system
   */
  virtual void init_data () override;

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () override;

  /**
   * Prepares \p matrix and \p _dof_map for matrix assembly.
   * Does not actually assemble anything.  For matrix assembly,
   * use the \p assemble() in derived classes.
   * Should be overridden in derived classes.
   */
  virtual void assemble () override { ImplicitSystem::assemble(); }

  /**
   * After calling this method, any solve will be limited to the given
   * subset.  To disable this mode, call this method with \p subset
   * being a \p nullptr.
   */
  virtual void restrict_solve_to (const SystemSubset * subset,
                                  const SubsetSolveMode subset_solve_mode=SUBSET_ZERO) override;

  /**
   * Assembles & solves the linear system A*x=b.
   */
  virtual void solve () override;

  /**
   * \returns A pointer to a linear solver appropriate for use in
   * adjoint and/or sensitivity solves
   */
  virtual LinearSolver<Number> * get_linear_solver() const override;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   */
  virtual void assembly(bool get_residual,
                        bool get_jacobian,
                        bool apply_heterogeneous_constraints = false,
                        bool apply_no_constraints = false) override;

  /**
   * \returns \p "LinearImplicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const override { return "LinearImplicit"; }

  /**
   * \returns The number of iterations
   * taken for the most recent linear solve.
   */
  unsigned int n_linear_iterations() const { return _n_linear_iterations; }

  /**
   * \returns The final residual for the linear system solve.
   */
  Real final_linear_residual() const { return _final_linear_residual; }

  /**
   * This function enables the user to provide a shell matrix, i.e. a
   * matrix that is not stored element-wise, but as a function.  When
   * you register your shell matrix using this function, calling \p
   * solve() will no longer use the \p matrix member but the
   * registered shell matrix instead.  You can reset this behaviour to
   * its original state by supplying a \p nullptr to this
   * function.
   */
  void attach_shell_matrix (ShellMatrix<Number> * shell_matrix);

  /**
   * Detaches a shell matrix.  Same as \p attach_shell_matrix(nullptr).
   */
  void detach_shell_matrix () { attach_shell_matrix(nullptr); }

  /**
   * \returns A pointer to the currently attached shell matrix, if any,
   * otherwise \p nullptr.
   */
  ShellMatrix<Number> * get_shell_matrix() { return _shell_matrix; }

protected:

  /**
   * The number of linear iterations required to solve the linear
   * system Ax=b.
   */
  unsigned int _n_linear_iterations;

  /**
   * The final residual for the linear system Ax=b.
   */
  Real _final_linear_residual;

  /**
   * User supplies shell matrix or \p nullptr if no shell matrix is used.
   */
  ShellMatrix<Number> * _shell_matrix;

  /**
   * The current subset on which to solve (or \p nullptr if none).
   */
  const SystemSubset * _subset;

  /**
   * If restrict-solve-to-subset mode is active, this member decides
   * what happens with the dofs outside the subset.
   */
  SubsetSolveMode _subset_solve_mode;
};

} // namespace libMesh

#endif // LIBMESH_LINEAR_IMPLICIT_SYSTEM_H
