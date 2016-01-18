// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_EIGEN_TIME_SOLVER_H
#define LIBMESH_EIGEN_TIME_SOLVER_H

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// Local includes
#include "libmesh/time_solver.h"

// C++ includes

namespace libMesh
{

// Forward declarations
template <typename T> class EigenSolver;


/**
 * The name of this class is confusing...it's meant to refer to the
 * base class (TimeSolver) while still telling one that it's for solving
 * (generalized) EigenValue problems that arise from finite element
 * discretizations.  For a time-dependent problem du/dt=F(u), with a
 * steady solution 0=F(u_0), we look at the time evolution of a small
 * perturbation, p=u-u_0, for which the (linearized) governing equation is
 *
 * dp/dt = F'(u_0)p
 *
 * where F'(u_0) is the Jacobian.  The generalized eigenvalue problem arises
 * by considering perturbations of the general form p = exp(lambda*t)x, which
 * leads to
 *
 *  Ax = lambda*Bx
 *
 * where A is the (discretized by FEM) Jacobian matrix and B is the
 * (discretized by FEM) mass matrix.
 *
 * The EigenSystem class (by Steffen Petersen) is related but does not
 * fall under the FEMSystem paradigm invented by Roy Stogner.  The EigenSolver
 * class (also by Steffen) is meant to provide a generic "linear solver"
 * interface for EigenValue software.  The only current concrete implementation
 * is a SLEPc-based eigensolver class, which we make use of here as well.
 *
 * \author John W. Peterson
 * \date 2007
 */
class EigenTimeSolver : public TimeSolver
{
public:
  /**
   * The type of system
   */
  typedef DifferentiableSystem sys_type;

  /**
   * The parent class
   */
  typedef TimeSolver Parent;

  /**
   * Constructor. Requires a reference to the system
   * to be solved.
   */
  explicit
  EigenTimeSolver (sys_type & s);

  /**
   * Destructor.
   */
  virtual ~EigenTimeSolver ();

  /**
   * The initialization function.  This method is used to
   * initialize internal data structures before a simulation begins.
   */
  virtual void init () libmesh_override;

  /**
   * The reinitialization function.  This method is used after
   * changes in the mesh
   */
  virtual void reinit () libmesh_override;

  /**
   * Implements the assembly of both matrices A and B, and calls
   * the EigenSolver to compute the eigenvalues.
   */
  virtual void solve () libmesh_override;

  /**
   * It doesn't make sense to advance the timestep, so we shouldn't call this.
   */
  virtual void advance_timestep () libmesh_override {}

  /**
   * error convergence order against deltat is
   * not applicable to an eigenvalue problem.
   */
  Real error_order() const { return 0.; }

  /**
   * Forms either the spatial (Jacobian) or mass matrix part of the
   * operator, depending on which is requested.
   */
  virtual bool element_residual (bool get_jacobian,
                                 DiffContext &) libmesh_override;

  /**
   * Forms the jacobian of the boundary terms.
   */
  virtual bool side_residual (bool get_jacobian,
                              DiffContext &) libmesh_override;

  /**
   * Forms the jacobian of the nonlocal terms.
   */
  virtual bool nonlocal_residual (bool get_jacobian,
                                  DiffContext &) libmesh_override;

  /**
   * Nominally computes the size of the difference between
   * successive solution iterates ||u^{n+1} - u^{n}|| in some norm,
   * but for this class just returns 0.
   */
  virtual Real du (const SystemNorm &) const libmesh_override { return 0.; }

  /**
   * This is effectively a steady-state solver.
   */
  virtual bool is_steady() const libmesh_override { return true; }

  /**
   * The EigenSolver object.  This is what actually
   * makes the calls to SLEPc.
   */
  UniquePtr<EigenSolver<Number> > eigen_solver;

  /**
   * The linear solver tolerance to be used when solving the
   * eigenvalue problem. FIXME: need more info...
   */
  Real tol;

  /**
   * The maximum number of iterations allowed to solve the problem.
   */
  unsigned int maxits;

  /**
   * The number of eigenvectors/values to be computed.
   */
  unsigned int n_eigenpairs_to_compute;

  /**
   * The number of basis vectors to use in the computation.  According
   * to ex16, the number of basis vectors must be >= the number of
   * eigenpairs requested, and ncv >= 2*nev is recommended.
   * Increasing this number, even by a little bit, can _greatly_
   * reduce the number of (EigenSolver) iterations required to compute
   * the desired number of eigenpairs, but the _cost per iteration_
   * goes up drastically as well.
   */
  unsigned int n_basis_vectors_to_use;

  /**
   * After a solve, holds the number of eigenpairs successfully
   * converged.
   */
  unsigned int n_converged_eigenpairs;

  /**
   * After a solve, holds the number of iterations required to converge
   * the requested number of eigenpairs.
   */
  unsigned int n_iterations_reqd;

private:

  enum NowAssembling {
    /**
     * The matrix associated with the spatial part of the operator.
     */
    Matrix_A,

    /**
     * The matrix associated with the time derivative (mass matrix).
     */
    Matrix_B,

    /**
     * The enum is in an invalid state.
     */
    Invalid_Matrix
  };

  /**
   * Flag which controls the internals of element_residual() and side_residual().
   */
  NowAssembling now_assembling;
};

} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC
#endif // LIBMESH_EIGEN_TIME_SOLVER_H
