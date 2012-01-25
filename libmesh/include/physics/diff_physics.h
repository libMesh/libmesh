
// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __diff_physics_h__
#define __diff_physics_h__

// C++ includes
#include <vector>

// Local Includes
#include "auto_ptr.h"
#include "diff_context.h"

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
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// DifferentiableSystem class definition

class DifferentiablePhysics
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DifferentiablePhysics () : compute_internal_sides(false) {}

  /**
   * Destructor.
   */
  virtual ~DifferentiablePhysics ();

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
   * To implement the physics model du/dt = F(u), the user should
   * examine u = elem_solution and add (F(u), phi_i) to 
   * elem_residual.
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
   * To implement the constraint 0 = G(u), the user should
   * examine u = elem_solution and add (G(u), phi_i) to 
   * elem_residual.
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
   */
  virtual bool side_time_derivative (bool request_jacobian,
                                     DiffContext &) {
    return request_jacobian;
  }

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
   */
  virtual bool side_constraint (bool request_jacobian,
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
    libmesh_assert(_time_evolving.size() > var);
    _time_evolving[var] = true;
  }

  /**
   * Adds a pseudo-convection contribution on \p elem to
   * elem_residual, if the nodes of \p elem are being translated by a
   * moving mesh.
   *
   * The library provides a basic implementation in
   * FEMSystem::eulerian_residual()
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
   * FEMSystem::mass_residual; few users will need to reimplement
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

  /*
   * Prepares the result of a build_context() call for use.
   * 
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular physics requires.
   */
  virtual void init_context(DiffContext &) {}
 
protected:

  /**
   * Stores bools to tell us which variables are evolving
   * in time and which are just constraints
   */
  std::vector<bool> _time_evolving;
};

} // namespace libMesh


#endif
