// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIFF_SYSTEM_H
#define LIBMESH_DIFF_SYSTEM_H

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/diff_context.h"
#include "libmesh/diff_physics.h"
#include "libmesh/diff_qoi.h"
#include "libmesh/implicit_system.h"
#include "libmesh/time_solver.h"

// C++ includes

namespace libMesh
{

// Forward Declarations
class TimeSolver;

template <typename T> class NumericVector;

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// DifferentiableSystem class definition

class DifferentiableSystem : public ImplicitSystem,
                             public virtual DifferentiablePhysics,
                             public virtual DifferentiableQoI
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DifferentiableSystem (EquationSystems& es,
			const std::string& name,
			const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~DifferentiableSystem ();

  /**
   * The type of system.
   */
  typedef DifferentiableSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef ImplicitSystem Parent;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();

  /**
   * Prepares \p matrix and \p rhs for matrix assembly.
   * Users should not reimplement this
   */
  virtual void assemble ();

  /**
   * Returns a pointer to a linear solver appropriate for use in
   * adjoint and/or sensitivity solves
   */
  virtual LinearSolver<Number> *get_linear_solver() const;

  /**
   * Returns an integer corresponding to the upper iteration count
   * limit and a Real corresponding to the convergence tolerance to
   * be used in linear adjoint and/or sensitivity solves
   */
  virtual std::pair<unsigned int, Real>
    get_linear_solve_parameters() const;

  /**
   * Releases a pointer to a linear solver acquired by
   * \p this->get_linear_solver()
   */
  virtual void release_linear_solver(LinearSolver<Number> *) const;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   */
  virtual void assembly (bool get_residual, bool get_jacobian) = 0;

  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   */
  virtual void solve ();

  /**
   * This function sets the _is_adjoint boolean member of TimeSolver to
   * true and then calls the adjoint_solve in implicit system
   */
  virtual std::pair<unsigned int, Real> adjoint_solve (const QoISet& qoi_indices = QoISet());

  /**
   * We don't allow systems to be attached to each other
   */
  virtual AutoPtr<DifferentiablePhysics> clone_physics()
  { libmesh_error();
    // dummy
    return AutoPtr<DifferentiablePhysics>(this); }

  /**
   * We don't allow systems to be attached to each other
   */
  virtual AutoPtr<DifferentiableQoI> clone()
  { libmesh_error();
    // dummy
    return AutoPtr<DifferentiableQoI>(this); }

  /**
   * Returns const reference to DifferentiablePhysics object. Note
   * that if no external Physics object is attached, the default is
   * this.
   */
  const DifferentiablePhysics* get_physics() const
  { return this->_diff_physics; }

  /**
   * Returns reference to DifferentiablePhysics object. Note that if
   * no external Physics object is attached, the default is this.
   */
  DifferentiablePhysics* get_physics()
  { return this->_diff_physics; }

  /**
   * Attach external Physics object.
   */
  void attach_physics( DifferentiablePhysics* physics_in )
  { this->_diff_physics = (physics_in->clone_physics()).release();
    this->_diff_physics->init_physics(*this);}

  /**
   * Returns const reference to DifferentiableQoI object. Note that if no external
   * QoI object is attached, the default is this.
   */
  const DifferentiableQoI* get_qoi() const
  { return this->diff_qoi; }

  /**
   * Returns reference to DifferentiableQoI object. Note that if no external
   * QoI object is attached, the default is this.
   */
  DifferentiableQoI* get_qoi()
  { return this->diff_qoi; }

  /**
   * Attach external QoI object.
   */
  void attach_qoi( DifferentiableQoI* qoi_in )
  { this->diff_qoi = (qoi_in->clone()).release();
    // User needs to resize qoi system qoi accordingly
    this->diff_qoi->init_qoi( this->qoi );}

  /**
   * A pointer to the solver object we're going to use.
   * This must be instantiated by the user before solving!
   */
  AutoPtr<TimeSolver> time_solver;

  /**
   * Sets the time_solver
   */
  void set_time_solver(AutoPtr<TimeSolver> _time_solver)
  {
    time_solver = _time_solver;
  }

  /**
   * Returns a pointer to the time solver attached to the calling system
   */
  TimeSolver& get_time_solver();

  /**
   * Non-const version of the above
   */
  const TimeSolver& get_time_solver() const;

  /**
   * For time-dependent problems, this is the amount delta t to advance the
   * solution in time.
   */
  Real deltat;

  /**
   * Builds a DiffContext object with enough information to do
   * evaluations on each element.
   *
   * For most problems, the default "Let FEMSystem build an * FEMContext"
   * reimplementation is correct; users who subclass FEMContext will need to
   * also reimplement this method to build it.
   */
  virtual AutoPtr<DiffContext> build_context();

  /**
   * Executes a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess (){}

  /**
   * Does any work that needs to be done on \p elem in a postprocessing loop.
   */
  virtual void element_postprocess (DiffContext &) {}

  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * postprocessing loop.
   */
  virtual void side_postprocess (DiffContext &) {}

  /**
   * If \p postprocess_sides is true (it is false by default), the
   * postprocessing loop will loop over all sides as well as all elements.
   */
  bool postprocess_sides;

  /**
   * Set print_residual_norms to true to print |U| whenever it is
   * used in an assembly() call
   */
  bool print_solution_norms;

  /**
   * Set print_solutions to true to print U whenever it is used in an
   * assembly() call
   */
  bool print_solutions;

  /**
   * Set print_residual_norms to true to print |F| whenever it is assembled.
   */
  bool print_residual_norms;

  /**
   * Set print_residuals to true to print F whenever it is assembled.
   */
  bool print_residuals;

  /**
   * Set print_jacobian_norms to true to print |J| whenever it is assembled.
   */
  bool print_jacobian_norms;

  /**
   * Set print_jacobians to true to print J whenever it is assembled.
   */
  bool print_jacobians;

  /**
   * Set print_element_jacobians to true to print each J_elem contribution.
   */
  bool print_element_jacobians;

protected:

  /**
   * Pointer to object to use for physics assembly evaluations.
   * Defaults to \p this for backwards compatibility; in the future
   * users should create separate physics objects.
   */
  DifferentiablePhysics *_diff_physics;

  /**
   * Pointer to object to use for quantity of interest assembly
   * evaluations.  Defaults to \p this for backwards compatibility; in
   * the future users should create separate physics objects.
   */
  DifferentiableQoI* diff_qoi;

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();
};

// --------------------------------------------------------------
// DifferentiableSystem inline methods
inline
  TimeSolver& DifferentiableSystem::get_time_solver()
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return *time_solver;
}

inline
 const TimeSolver& DifferentiableSystem::get_time_solver() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return *time_solver;
}

} // namespace libMesh


#endif // LIBMESH_DIFF_SYSTEM_H
