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



#ifndef LIBMESH_DIFF_SYSTEM_H
#define LIBMESH_DIFF_SYSTEM_H

// Local Includes
#include "libmesh/diff_context.h"
#include "libmesh/diff_physics.h"
#include "libmesh/diff_qoi.h"
#include "libmesh/implicit_system.h"
#include "libmesh/time_solver.h"

// C++ includes
#include <memory>

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
 * \author Roy H. Stogner
 * \date 2006
 */
class DifferentiableSystem : public ImplicitSystem,
                             public virtual DifferentiablePhysics,
                             public virtual DifferentiableQoI
{
public:

  /**
   * Constructor.
   */
  DifferentiableSystem (EquationSystems & es,
                        const std::string & name,
                        const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions as its base class.
   * - This class also can't be default move constructed because it
   *   manually manages the memory of the _diff_physics and diff_qoi
   *   objects.
   * - The destructor potentially cleans up the _diff_physics and
   *   diff_qoi objects.
   */
  DifferentiableSystem (const DifferentiableSystem &) = delete;
  DifferentiableSystem & operator= (const DifferentiableSystem &) = delete;
  DifferentiableSystem (DifferentiableSystem &&) = delete;
  DifferentiableSystem & operator= (DifferentiableSystem &&) = delete;
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
  virtual void clear () override;

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () override;

  /**
   * Prepares \p matrix and \p rhs for matrix assembly.
   * Users should not reimplement this.
   * Note that in some cases only
   * \link current_local_solution \endlink is used during assembly,
   * and, therefore, if \link solution \endlink has been altered
   * without \link update() \endlink being called, then the
   * user must call \link update() \endlink before calling
   * this function.
   */
  virtual void assemble () override;

  /**
   * \returns A pointer to a linear solver appropriate for use in
   * adjoint and/or sensitivity solves
   */
  virtual LinearSolver<Number> * get_linear_solver() const override;

  /**
   * \returns An integer corresponding to the upper iteration count
   * limit and a Real corresponding to the convergence tolerance to
   * be used in linear adjoint and/or sensitivity solves
   */
  virtual std::pair<unsigned int, Real>
  get_linear_solve_parameters() const override;

  /**
   * Assembles a residual in \p rhs and/or a jacobian in \p matrix,
   * as requested.
   * Note that in some cases only
   * \link current_local_solution \endlink is used during assembly,
   * and, therefore, if \link solution \endlink has been altered
   * without \link update() \endlink being called, then the
   * user must call \link update() \endlink before calling
   * this function.
   */
  virtual void assembly (bool get_residual,
                         bool get_jacobian,
                         bool apply_heterogeneous_constraints = false,
                         bool apply_no_constraints = false) override = 0;

  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   */
  virtual void solve () override;

  /**
   * This function sets the _is_adjoint boolean member of TimeSolver to
   * true and then calls the adjoint_solve in implicit system
   */
  virtual std::pair<unsigned int, Real>
  adjoint_solve (const QoISet & qoi_indices = QoISet()) override;

  /**
   * We don't allow systems to be attached to each other
   */
  virtual std::unique_ptr<DifferentiablePhysics> clone_physics() override
  {
    libmesh_not_implemented();
    // dummy to avoid compiler warnings, not a real implementation
    return std::unique_ptr<DifferentiablePhysics>(nullptr);
  }

  /**
   * We don't allow systems to be attached to each other
   */
  virtual std::unique_ptr<DifferentiableQoI> clone() override
  {
    libmesh_not_implemented();
    // dummy to avoid compiler warnings, not a real implementation
    return std::unique_ptr<DifferentiableQoI>(nullptr);
  }

  /**
   * \returns A const reference to a DifferentiablePhysics object.
   *
   * \note If no external Physics object is attached, the default is
   * \p this.
   */
  const DifferentiablePhysics * get_physics() const
  { return this->_diff_physics; }

  /**
   * \returns A reference to a DifferentiablePhysics object.
   *
   * \note If no external Physics object is attached, the default is
   * \p this.
   */
  DifferentiablePhysics * get_physics()
  { return this->_diff_physics; }

  /**
   * Attach external Physics object.
   */
  void attach_physics( DifferentiablePhysics * physics_in )
  { this->_diff_physics = (physics_in->clone_physics()).release();
    this->_diff_physics->init_physics(*this);}

  /**
   * Swap current physics object with external object
   */
  void swap_physics ( DifferentiablePhysics * & swap_physics );

  /**
   * \returns A const reference to a DifferentiableQoI object.
   *
   * \note If no external QoI object is attached, the default is \p this.
   */
  const DifferentiableQoI * get_qoi() const
  { return this->diff_qoi; }

  /**
   * \returns A reference to a DifferentiableQoI object.
   *
   * \note If no external QoI object is attached, the default is this.
   */
  DifferentiableQoI * get_qoi()
  { return this->diff_qoi; }

  /**
   * Attach external QoI object.
   */
  void attach_qoi( DifferentiableQoI * qoi_in );

  /**
   * A pointer to the solver object we're going to use.
   * This must be instantiated by the user before solving!
   */
  std::unique_ptr<TimeSolver> time_solver;

  /**
   * Sets the time_solver
   * FIXME: This code is a little dangerous as it transfers ownership
   * from the TimeSolver creator to this class.  The user must no longer
   * access his original TimeSolver object after calling this function.
   */
  void set_time_solver(std::unique_ptr<TimeSolver> _time_solver)
  {
    time_solver.reset(_time_solver.release());
  }

  /**
   * \returns A pointer to the time solver attached to the calling system
   */
  TimeSolver & get_time_solver();

  /**
   * Non-const version of the above
   */
  const TimeSolver & get_time_solver() const;

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
  virtual std::unique_ptr<DiffContext> build_context();

  /**
   * Executes a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess () {}

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
   * For a given second order (in time) variable var, this method will return
   * the index to the corresponding "dot" variable. For FirstOrderUnsteadySolver
   * classes, the "dot" variable would automatically be added and the returned
   * index will correspond to that variable. For SecondOrderUnsteadySolver classes,
   * this method will return var as there this is no "dot" variable per se, but
   * having this function allows one to use the interface to treat both
   * FirstOrderUnsteadySolver and SecondOrderUnsteadySolver simultaneously.
   */
  unsigned int get_second_order_dot_var( unsigned int var ) const;

  /**
   * Check for any first order vars that are also belong to FEFamily::SCALAR
   */
  bool have_first_order_scalar_vars() const;

  /**
   * Check for any second order vars that are also belong to FEFamily::SCALAR
   */
  bool have_second_order_scalar_vars() const;

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
   * Set print_element_solutions to true to print each U_elem input.
   */
  bool print_element_solutions;

  /**
   * Set print_element_residuals to true to print each R_elem contribution.
   */
  bool print_element_residuals;

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
  DifferentiablePhysics * _diff_physics;

  /**
   * Pointer to object to use for quantity of interest assembly
   * evaluations.  Defaults to \p this for backwards compatibility; in
   * the future users should create separate physics objects.
   */
  DifferentiableQoI * diff_qoi;

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   *
   * If the TimeSolver is a FirstOrderUnsteadySolver,
   * we check for second order in time variables and then add them
   * to the System as dot_<varname>. Then, during assembly, the
   * TimeSolver will populate the elem_accel vectors with the
   * dot_<varname> values so the user's element assembly function
   * can still treat the variable as a second order in time variable.
   */
  virtual void init_data () override;

  /**
   * Helper function to add "velocity" variables that are cousins to
   * second order-in-time variables in the DifferentiableSystem. This
   * function is only called if the TimeSolver is a FirstOrderUnsteadySolver.
   */
  void add_second_order_dot_vars();

#ifdef LIBMESH_ENABLE_DIRICHLET
  /**
   * Helper function to and Dirichlet boundary conditions to "dot" variable
   * cousins of second order variables in the system. The function takes the
   * second order variable index, it's corresponding "dot" variable index and
   * then searches for DirichletBoundary objects for var_idx and then adds a
   * DirichletBoundary object for dot_var_idx using the same boundary ids and
   * functors for the var_idx DirichletBoundary.
   */
  void add_dot_var_dirichlet_bcs( unsigned int var_idx, unsigned int dot_var_idx);
#endif

};

// --------------------------------------------------------------
// DifferentiableSystem inline methods
inline
TimeSolver & DifferentiableSystem::get_time_solver()
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return *time_solver;
}

inline
const TimeSolver & DifferentiableSystem::get_time_solver() const
{
  libmesh_assert(time_solver.get());
  libmesh_assert_equal_to (&(time_solver->system()), this);
  return *time_solver;
}

} // namespace libMesh


#endif // LIBMESH_DIFF_SYSTEM_H
