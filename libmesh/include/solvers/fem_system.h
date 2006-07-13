
// $Id: fem_system.h,v 1.7 2006-07-13 05:18:30 roystgnr Exp $

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



#ifndef __fem_system_h__
#define __fem_system_h__

// C++ includes

// Local Includes
#include "diff_system.h"
#include "elem.h"

#ifdef ENABLE_SECOND_DERIVATIVES
#include "tensor_value.h"
#endif

// Forward Declarations

template <typename T> class NumericVector;

/**
 * This class provides a specific system class.  It aims
 * at nonlinear implicit systems, requiring only a
 * cell residual calculation from the user.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the class \p ExplicitSystem.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// FEMSystem class definition

class FEMSystem : public DifferentiableSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FEMSystem (EquationSystems& es,
	         const std::string& name,
	         const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~FEMSystem ();

  /**
   * The type of system.
   */
  typedef FEMSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef DifferentiableSystem Parent;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Prepares \p matrix or \p rhs for matrix assembly.
   * Users should not need to reimplement this
   */
  void assembly (bool get_residual, bool get_jacobian);

  /**
   * Tells the FEMSystem that variable var is evolving with
   * respect to time.  In general, the user's init() function
   * should call time_evolving() for any variables which
   * behave like du/dt = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   *
   * Most derived systems will not have to reimplment this function; however
   * any system which reimplements mass_residual() may have to reimplement
   * time_evolving() to prepare data structures.
   */
  virtual void time_evolving (unsigned int var);

  /**
   * Adds a mass vector contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Most problems can use the FEMSystem::mass_residual implementation,
   * which calculates the residual (u, phi_i) and jacobian (phi_i, phi_j);
   * few users will need to reimplement this themselves.  Using a custom
   * mass matrix (e.g. for divergence-free elements or mass lumping)
   * requires reimplementing mass_residual().
   */
  virtual bool mass_residual (bool request_jacobian);

  /**
   * Returns the value of the solution variable \p at the quadrature
   * point \qp on the current element interior
   */
  Number interior_value(unsigned int var, unsigned int qp);

  /**
   * Returns the value of the solution variable \p at the quadrature
   * point \qp on the current element side
   */
  Number side_value(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p at the quadrature
   * point \qp on the current element interior
   */
  Gradient interior_gradient(unsigned int var, unsigned int qp);

  /**
   * Returns the gradient of the solution variable \p at the quadrature
   * point \qp on the current element side
   */
  Gradient side_gradient(unsigned int var, unsigned int qp);

#ifdef ENABLE_SECOND_DERIVATIVES
  /**
   * Returns the hessian of the solution variable \p at the quadrature
   * point \qp on the current element interior
   */
  Tensor interior_hessian(unsigned int var, unsigned int qp);

  /**
   * Returns the hessian of the solution variable \p at the quadrature
   * point \qp on the current element side
   */
  Tensor side_hessian(unsigned int var, unsigned int qp);
#endif // ENABLE_SECOND_DERIVATIVES

  /**
   * @returns \p "General".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "PDE"; }

  /**
   * If calculating numeric jacobians is required, the FEMSystem
   * will perturb each solution vector entry by numerical_jacobian_h
   * when calculating finite differences.
   */ 
  Real numerical_jacobian_h;

  /**
   * If verify_analytic_jacobian is equal to zero (as it is by
   * default), no numeric jacobians will be calculated unless
   * an overloaded element_time_derivative(), element_constraint(),
   * side_time_derivative(), or side_constraint() function cannot
   * provide an analytic jacobian upon request.
   * 
   * If verify_analytic_jacobian is equal to the positive value tol,
   * then any time a full analytic element jacobian can be calculated
   * it will be tested against a numerical jacobian on the same element,
   * and the program will abort if the relative error (in matrix l1 norms)
   * exceeds tol.
   */
  Real verify_analytic_jacobians;

protected:

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Uses the results of multiple element_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element.
   */
  void numerical_elem_jacobian ();

  /**
   * Uses the results of multiple side_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element's side.
   */
  void numerical_side_jacobian ();

  /**
   * Finite element objects for each variable's interior and sides.
   */
  std::map<FEType, FEBase *> element_fe;
  std::map<FEType, FEBase *> side_fe;

  /**
   * Quadrature rules for element interior and sides.
   * The FEM system will try to find a quadrature rule that
   * correctly integrates all variables
   */
  QBase *element_qrule;
  QBase *side_qrule;

  /**
   * Current element for element_* to examine
   */
  Elem *elem;

  /**
   * Current element for side_* to examine
   */
  unsigned int side;
};


#endif
