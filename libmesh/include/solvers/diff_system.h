
// $Id: diff_system.h,v 1.1 2006-06-05 00:32:23 roystgnr Exp $

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



#ifndef __diff_system_h__
#define __diff_system_h__

// C++ includes

// Local Includes
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "transient_system.h"

// Forward Declarations
class FEBase;
class QBase;
class TimeSolver;

template <typename T> class NumericVector;

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 */

// ------------------------------------------------------------
// DifferentiableSystem class definition

class DifferentiableSystem : public ImplicitSystem
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
   * Prepares \p matrix and \p rhs for matrix assembly.
   * Users should not reimplement this
   */
  virtual void assemble ();

  /**
   * Evaluates \p matrix and/or \p rhs at the current iterate
   */
  virtual void assembly (bool get_residual, bool get_jacobian) = 0;

  /**
   * Adds the time derivative contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Users need to reimplement this for their particular PDE.
   * Given the physics model du/dt = F(u), the user should
   * examine u = _nonlinear_solution and add (F(u), phi_i) to 
   * elem_residual.
   */
  virtual bool element_time_derivative (bool request_jacobian) {
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
   * Given the constrain 0 = G(u), the user should
   * examine u = _nonlinear_solution and add (G(u), phi_i) to 
   * elem_residual.
   */
  virtual bool element_constraint (bool request_jacobian) {
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
  virtual bool side_time_derivative (bool request_jacobian) {
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
  virtual bool side_constraint (bool request_jacobian) {
    return request_jacobian;
  }
 
  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   */
  virtual void solve ();
 
  /**
   * @returns \p "Differentiable".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Differentiable"; }

  /**
   * A pointer to the solver object we're going to use.
   * This must be instantiated by the user before solving!
   */
  AutoPtr<TimeSolver> time_solver;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();
};


#endif
