
// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIFF_PHYSICS_H
#define LIBMESH_DIFF_PHYSICS_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h"

// C++ includes
#include <vector>

namespace libMesh
{

// Forward Declarations
class System;
class DiffContext;

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 * For first order (in time) systems, the (nonlinear)
 * residual computed at each time step is
 * \f[ F(u) + G(u) - M(u,\dot{u})\dot{u} = 0 \f]
 * for unsteady TimeSolver and \f$ F(u) = 0\f$ for steady
 * TimeSolver. \f$F(u)\f$ is computed by element/side_time_derivative,
 * \f$ G(u) \f$ is computed using element/side_constraint, and
 * \f$ M(u,\dot{u})\dot{u} \f$ is computed using the mass_residual
 * methods.
 *
 * For second order (in time) systems, the (nonlinear)
 * residual computed at each time step is
 * \f[ M(u,\ddot{u})\ddot{u} + C(u,\dot{u})\dot{u} + F(u) + G(u) = 0 \f]
 * for unsteady TimeSolver and \f$ F(u) = 0\f$ for steady
 * TimeSolver. \f$F(u)\f$ is computed by element/side_time_derivative,
 * \f$G(u)\f$ is computed using element/side_constraint,
 * \f$C(u,\dot{u})\dot{u}\f$ is computed using the dampling_residual methods
 * and \f$ -M(u,\ddot{u})\ddot{u}\f$ is computed using the mass_residual
 * methods. This is the sign convention used by the default implementation;
 * if the method is overridden, the user can choose any self-consistent sign
 * convention they wish.
 *
 * FEMContext provides methods for querying values of the solution \f${u}\f$,
 * its "rate" \f$\dot{u}\f$ and its "acceleration" \f$\ddot{u}\f$. Furthermore,
 * derivatives of each of these w.r.t the nonlinear iteration unknown (e.g. in
 * EulerSolver, the solution at the next time step \f$ u_{n+1} \f$) are provided
 * through DiffContext::get_elem_solution_derivative(),
 * DiffContext::get_elem_solution_rate_derivative(), and
 * DiffContext::get_elem_solution_accel_derivative(). The should be incorporated
 * into the Jacobian evaluations, if the Jacobian is being provided.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class DifferentiablePhysics
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DifferentiablePhysics () :
    compute_internal_sides (false),
    _mesh_sys              (libmesh_nullptr),
    _mesh_x_var            (libMesh::invalid_uint),
    _mesh_y_var            (libMesh::invalid_uint),
    _mesh_z_var            (libMesh::invalid_uint)
  {}

  /**
   * Destructor.
   */
  virtual ~DifferentiablePhysics ();

  /**
   * Copy of this object. User should override to copy any needed state.
   */
  virtual UniquePtr<DifferentiablePhysics> clone_physics() = 0;

  /**
   * Clear any data structures associated with the physics.
   */
  virtual void clear_physics ();

  /**
   * Initialize any data structures associated with the physics.
   */
  virtual void init_physics (const System & sys);

  /**
   * Adds the time derivative contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Users need to reimplement this for their particular PDE.
   *
   * To implement the physics model du/dt = F(u), the user should
   * examine u = elem_solution and add (F(u), phi_i) to elem_residual
   * in elem_time_derivative().
   */
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext &) {
    return request_jacobian;
  }

  /**
   * Adds the constraint contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Users may need to reimplement this for their particular PDE.
   *
   * To implement the constraint 0 = G(u), the user should
   * examine u = elem_solution and add (G(u), phi_i) to elem_residual
   * in elem_constraint().
   */
  virtual bool element_constraint (bool request_jacobian,
                                   DiffContext &) {
    return request_jacobian;
  }

  /**
   * \p compute_internal_sides is false by default, indicating that
   * side_* computations will only be done on boundary sides.  If
   * compute_internal_sides is true, computations will be done
   * on sides between elements as well.
   */
  bool compute_internal_sides;

  /**
   * Adds the time derivative contribution on \p side of \p elem to
   * elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Users may need to reimplement this for their particular PDE
   * depending on the boundary conditions.
   *
   * To implement a weak form of the source term du/dt = F(u) on
   * sides, such as might arise in a flux boundary condition, the user
   * should examine u = elem_solution and add (F(u), phi_i) boundary
   * integral contributions to elem_residual in side_constraint().
   */
  virtual bool side_time_derivative (bool request_jacobian,
                                     DiffContext &) {
    return request_jacobian;
  }

  /**
   * Adds the constraint contribution on \p side of \p elem to
   * elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Users may need to reimplement this for their particular PDE
   * depending on the boundary conditions.
   *
   * To implement a weak form of the constraint 0 = G(u), the user
   * should examine u = elem_solution and add (G(u), phi_i) boundary
   * integral contributions to elem_residual in side_constraint().
   */
  virtual bool side_constraint (bool request_jacobian,
                                DiffContext &) {
    return request_jacobian;
  }

  /**
   * Adds any nonlocal time derivative contributions (e.g. some
   * components of time derivatives in scalar variable equations) to
   * elem_residual
   *
   * If this method receives request_jacobian = true, then it
   * should also modify elem_jacobian and return true if possible.  If
   * the Jacobian changes have not been computed then the method
   * should return false.
   *
   * Users may need to reimplement this for PDEs on systems to which
   * SCALAR variables have been added.
   */
  virtual bool nonlocal_time_derivative (bool request_jacobian,
                                         DiffContext &) {
    return request_jacobian;
  }

  /**
   * Adds any nonlocal constraint contributions (e.g. some
   * components of constraints in scalar variable equations) to
   * elem_residual
   *
   * If this method receives request_jacobian = true, then it
   * should also modify elem_jacobian and return true if possible.  If
   * the Jacobian changes have not been computed then the method
   * should return false.
   *
   * Users may need to reimplement this for PDEs on systems to which
   * SCALAR variables with non-tranient equations have been added.
   */
  virtual bool nonlocal_constraint (bool request_jacobian,
                                    DiffContext &) {
    return request_jacobian;
  }


  /**
   * Tells the DiffSystem that variable var is evolving with
   * respect to time.  In general, the user's init() function
   * should call time_evolving() for any variables which
   * behave like du/dt = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   *
   * Most derived systems will not have to reimplment this function; however
   * any system which reimplements mass_residual() may have to reimplement
   * time_evolving() to prepare data structures.
   *
   * This method is deprecated. Instead, use the time_evolving override
   * and specify the order-in-time of the variable, either 1 or 2. This method
   * assumes the variable is first order for backward compatibility.
   */
  virtual void time_evolving (unsigned int var)
  {
    libmesh_deprecated();
    this->time_evolving(var,1);
  }

  /**
   * Tells the DiffSystem that variable var is evolving with
   * respect to time.  In general, the user's init() function
   * should call time_evolving() with order 1 for any variables which
   * behave like du/dt = F(u), with order 2 for any variables that
   * behave like d^2u/dt^2 = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   *
   * Most derived systems will not have to reimplment this function; however
   * any system which reimplements mass_residual() may have to reimplement
   * time_evolving() to prepare data structures.
   */
  virtual void time_evolving (unsigned int var, unsigned int order);

  /**
   * Returns true iff variable \p var is evolving with
   * respect to time.  In general, the user's init() function
   * should have set time_evolving() for any variables which
   * behave like du/dt = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   */
  bool is_time_evolving (unsigned int var) const
  {
    libmesh_assert_less(var,_time_evolving.size());
    libmesh_assert( _time_evolving[var] == 0 ||
                    _time_evolving[var] == 1 ||
                    _time_evolving[var] == 2 );
    return _time_evolving[var];
  }

  /**
   * Adds a pseudo-convection contribution on \p elem to
   * elem_residual, if the nodes of \p elem are being translated by a
   * moving mesh.
   *
   * The library provides a basic implementation in
   * FEMPhysics::eulerian_residual()
   */
  virtual bool eulerian_residual (bool request_jacobian,
                                  DiffContext &) {
    return request_jacobian;
  }

  /**
   * Subtracts a mass vector contribution on \p elem from
   * elem_residual. For first-order-in-time problems, this
   * is the \f$ M(u,\dot{u})\dot{u} \f$ term. For
   * second-order-in-time problems, this is the
   * \f$ M(u,\ddot{u})\ddot{u} \f$ term. This method is only
   * called for UnsteadySolver-based TimeSolvers.
   *
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Many first-order-in-time problems can use the reimplementation in
   * FEMPhysics::mass_residual which subtracts (du/dt,v) for each
   * transient variable u; users with more complicated transient
   * problems or second-order-in-time problems will need to reimplement
   * this themselves.
   */
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext &) {
    return request_jacobian;
  }

  /**
   * Subtracts a mass vector contribution on \p side of \p elem from
   * elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * For most problems, the default implementation of "do nothing"
   * is correct; users with boundary conditions including time
   * derivatives may need to reimplement this themselves.
   */
  virtual bool side_mass_residual (bool request_jacobian,
                                   DiffContext &) {
    return request_jacobian;
  }

  /**
   * Subtracts any nonlocal mass vector contributions (e.g. any time
   * derivative coefficients in scalar variable equations) from
   * elem_residual
   *
   * If this method receives request_jacobian = true, then it
   * should also modify elem_jacobian and return true if possible.  If
   * the Jacobian changes have not been computed then the method
   * should return false.
   *
   * Many problems can use the reimplementation in
   * FEMPhysics::mass_residual which subtracts (du/dt,v) for each
   * transient scalar variable u; users with more complicated
   * transient scalar variable equations will need to reimplement this
   * themselves.
   */
  virtual bool nonlocal_mass_residual (bool request_jacobian,
                                       DiffContext & c);

  /**
   * Subtracts a damping vector contribution on \p elem from
   * elem_residual. This method is not used in first-order-in-time
   * problems. For second-order-in-time problems, this is the
   * \f$ C(u,\ddot{u})\ddot{u} \f$ term. This method is only
   * called for UnsteadySolver-based TimeSolvers.
   *
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * If the problem has no damping, the default "do-nothing" is correct.
   * Otherwise, this must be reimplemented.
   */
  virtual bool damping_residual (bool request_jacobian,
                                 DiffContext &) {
    return request_jacobian;
  }

  /**
   * Subtracts a damping vector contribution on \p side of \p elem from
   * elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * For most problems, the default implementation of "do nothing"
   * is correct; users with boundary conditions including first time
   * derivatives may need to reimplement this themselves.
   */
  virtual bool side_damping_residual (bool request_jacobian,
                                      DiffContext &) {
    return request_jacobian;
  }

  /**
   * Subtracts any nonlocal damping vector contributions (e.g. any
   * first time derivative coefficients in scalar variable equations) from
   * elem_residual
   *
   * If this method receives request_jacobian = true, then it
   * should also modify elem_jacobian and return true if possible.  If
   * the Jacobian changes have not been computed then the method
   * should return false.
   */
  virtual bool nonlocal_damping_residual (bool request_jacobian,
                                          DiffContext &) {
    return request_jacobian;
  }

  /*
   * Prepares the result of a build_context() call for use.
   *
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular physics requires.
   */
  virtual void init_context(DiffContext &) {}

  /**
   * Tells the DifferentiablePhysics that system \p sys contains the
   * isoparametric Lagrangian variables which correspond to the
   * coordinates of mesh nodes, in problems where the mesh itself is
   * expected to move in time.
   *
   * The system with mesh coordinate data (which may be \p this system
   * itself, for fully coupled moving mesh problems) is currently
   * assumed to have new (end of time step) mesh coordinates stored in
   * solution, old (beginning of time step) mesh coordinates stored in
   * _old_nonlinear_solution, and constant velocity motion during each
   * time step.
   *
   * Activating this function ensures that local (but not neighbor!) element
   * geometry is correctly repositioned when evaluating element residuals.
   *
   * Currently \p sys must be \p *this for a tightly coupled moving
   * mesh problem or NULL to stop mesh movement; loosely coupled
   * moving mesh problems are not implemented.
   *
   * This code is experimental.  "Trust but verify, and not in that
   * order"
   */
  virtual void set_mesh_system(System * sys);

  /**
   * Returns a const reference to the system with variables corresponding to
   * mesh nodal coordinates, or NULL if the mesh is fixed.
   * Useful for ALE calculations.
   */
  const System * get_mesh_system() const;

  /**
   * Returns a reference to the system with variables corresponding to
   * mesh nodal coordinates, or NULL if the mesh is fixed.
   */
  System * get_mesh_system();

  /**
   * Tells the DifferentiablePhysics that variable \p var from the mesh system
   * should be used to update the x coordinate of mesh nodes, in problems where
   * the mesh itself is expected to move in time.
   *
   * The system with mesh coordinate data (which may be this system itself, for
   * fully coupled moving mesh problems) is currently assumed to have new (end
   * of time step) mesh coordinates stored in solution, old (beginning of time
   * step) mesh coordinates stored in _old_nonlinear_solution, and constant
   * velocity motion during each time step.
   *
   * Activating this function ensures that local (but not neighbor!) element
   * geometry is correctly repositioned when evaluating element residuals.
   */
  virtual void set_mesh_x_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh x coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_x_var() const;

  /**
   * Tells the DifferentiablePhysics that variable \p var from the mesh system
   * should be used to update the y coordinate of mesh nodes.
   */
  virtual void set_mesh_y_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh y coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_y_var() const;

  /**
   * Tells the DifferentiablePhysics that variable \p var from the mesh system
   * should be used to update the z coordinate of mesh nodes.
   */
  virtual void set_mesh_z_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh z coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_z_var() const;

  /**
   * This method simply combines element_time_derivative() and
   * eulerian_residual(), which makes its address useful as a
   * pointer-to-member-function when refactoring.
   */
  bool _eulerian_time_deriv (bool request_jacobian,
                             DiffContext &);

  bool have_first_order_vars() const
  { return !_first_order_vars.empty(); }

  /**
   * Returns the set of first order in time variable indices. May be empty.
   */
  const std::set<unsigned int> & get_first_order_vars() const
  { return _first_order_vars; }

  bool is_first_order_var( unsigned int var ) const
  { return _first_order_vars.find(var) != _first_order_vars.end(); }


  bool have_second_order_vars() const
  { return !_second_order_vars.empty(); }

  /**
   * Returns the set of second order in time variable indices. May be empty.
   */
  const std::set<unsigned int> & get_second_order_vars() const
  { return _second_order_vars; }

  bool is_second_order_var( unsigned int var ) const
  { return _second_order_vars.find(var) != _second_order_vars.end(); }


protected:

  /**
   * System from which to acquire moving mesh information
   */
  System * _mesh_sys;

  /**
   * Variables from which to acquire moving mesh information
   */
  unsigned int _mesh_x_var, _mesh_y_var, _mesh_z_var;

  /**
   * Stores unsigned int to tell us which variables are evolving
   * as first order in time (1), second order in time (2), or are
   * not time evolving (0).
   */
  std::vector<unsigned int> _time_evolving;

  /**
   * Variable indices for those variables that are first order in time.
   */
  std::set<unsigned int> _first_order_vars;

  /**
   * Variable indices for those variables that are second order in time.
   */
  std::set<unsigned int> _second_order_vars;

  /**
   * If the user adds any second order variables, then we need to also
   * cache the map to their corresponding dot variable that will
   * be added by this TimeSolver class.
   */
  std::map<unsigned int,unsigned int> _second_order_dot_vars;

};

// ------------------------------------------------------------
// DifferentiablePhysics inline methods


inline
void DifferentiablePhysics::set_mesh_system(System * sys)
{
  // For now we assume that we're doing fully coupled mesh motion
  //  if (sys && sys != this)
  //    libmesh_not_implemented();

  // For the foreseeable future we'll assume that we keep these
  // Systems in the same EquationSystems
  // libmesh_assert_equal_to (&this->get_equation_systems(),
  //                          &sys->get_equation_systems());

  // And for the immediate future this code may not even work
  libmesh_experimental();

  _mesh_sys = sys;
}



inline
void DifferentiablePhysics::set_mesh_x_var (unsigned int var)
{
  _mesh_x_var = var;
}



inline
void DifferentiablePhysics::set_mesh_y_var (unsigned int var)
{
  _mesh_y_var = var;
}



inline
void DifferentiablePhysics::set_mesh_z_var (unsigned int var)
{
  _mesh_z_var = var;
}



inline
const System * DifferentiablePhysics::get_mesh_system() const
{
  return _mesh_sys;
}

inline
System * DifferentiablePhysics::get_mesh_system()
{
  return _mesh_sys;
}

inline
unsigned int DifferentiablePhysics::get_mesh_x_var() const
{
  return _mesh_x_var;
}

inline
unsigned int DifferentiablePhysics::get_mesh_y_var() const
{
  return _mesh_y_var;
}

inline
unsigned int DifferentiablePhysics::get_mesh_z_var() const
{
  return _mesh_z_var;
}



} // namespace libMesh


#endif // LIBMESH_DIFF_PHYSICS_H
