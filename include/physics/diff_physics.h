
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



#ifndef LIBMESH_DIFF_PHYSICS_H
#define LIBMESH_DIFF_PHYSICS_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/diff_context.h"

// C++ includes
#include <vector>

namespace libMesh
{

// Forward Declarations
class System;
class TimeSolver;

template <typename T> class NumericVector;

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// DifferentiablePhysics class definition

class DifferentiablePhysics
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DifferentiablePhysics () :
    compute_internal_sides (false),
    _mesh_sys              (NULL),
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
  virtual AutoPtr<DifferentiablePhysics> clone_physics() = 0;

  /**
   * Clear any data structures associated with the physics.
   */
  virtual void clear_physics ();

  /**
   * Initialize any data structures associated with the physics.
   */
  virtual void init_physics (const System& sys);

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
   */
  virtual void time_evolving (unsigned int var) {
    if (_time_evolving.size() <= var)
      _time_evolving.resize(var+1, false);
    _time_evolving[var] = true;
  }

  /**
   * Returns true iff variable \p var is evolving with
   * respect to time.  In general, the user's init() function
   * should have set time_evolving() for any variables which
   * behave like du/dt = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   */
  bool is_time_evolving (unsigned int var) const {
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
   * Adds a mass vector contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Most problems can use the reimplementation in
   * FEMPhysics::mass_residual; few users will need to reimplement
   * this themselves.
   */
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext &) {
    return request_jacobian;
  }

  /**
   * Adds a mass vector contribution on \p side of \p elem to
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
   * Adds any nonlocal mass vector contributions (e.g. any time
   * derivative coefficients in scalar variable equations) to
   * elem_residual
   *
   * If this method receives request_jacobian = true, then it
   * should also modify elem_jacobian and return true if possible.  If
   * the Jacobian changes have not been computed then the method
   * should return false.
   *
   * Many problems can use the implementation in
   * DifferentiablePhysics::nonlocal_mass_residual, but users trying
   * to solve SCALAR variable equations which have time derivative
   * terms with mass coefficients != 1.0 will need to reimplement this
   * themselves.
   */
  virtual bool nonlocal_mass_residual (bool request_jacobian,
                                       DiffContext &c);

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
  virtual void set_mesh_system(System* sys);

  /**
   * Returns a const reference to the system with variables corresponding to
   * mesh nodal coordinates, or NULL if the mesh is fixed.
   * Useful for ALE calculations.
   */
  const System* get_mesh_system() const;

  /**
   * Returns a reference to the system with variables corresponding to
   * mesh nodal coordinates, or NULL if the mesh is fixed.
   */
  System* get_mesh_system();

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
                             DiffContext&);


protected:

  /**
   * System from which to acquire moving mesh information
   */
  System *_mesh_sys;

  /**
   * Variables from which to acquire moving mesh information
   */
  unsigned int _mesh_x_var, _mesh_y_var, _mesh_z_var;

  /**
   * Stores bools to tell us which variables are evolving
   * in time and which are just constraints
   */
  std::vector<bool> _time_evolving;
};



// ------------------------------------------------------------
// DifferentiablePhysics inline methods


inline
void DifferentiablePhysics::set_mesh_system(System* sys)
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
const System* DifferentiablePhysics::get_mesh_system() const
{
  return _mesh_sys;
}

inline
System* DifferentiablePhysics::get_mesh_system()
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
