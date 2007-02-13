
// $Id: diff_system.h,v 1.16 2007-02-13 21:29:02 roystgnr Exp $

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
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
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
   * To implement the physics model du/dt = F(u), the user should
   * examine u = elem_solution and add (F(u), phi_i) to 
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
   * To implement the constraint 0 = G(u), the user should
   * examine u = elem_solution and add (G(u), phi_i) to 
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
   * Executes a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess () = 0;

  /**
   * If \p postprocess_sides is true (it is false by default), the postprocessing
   * loop will loop over all sides as well as all elements.
   */
  bool postprocess_sides;

  /**
   * Does any work that needs to be done on \p elem in a postprocessing loop.
   */
  virtual void element_postprocess () {}
 
  /**
   * Does any work that needs to be done on \p side of \p elem in a
   * postprocessing loop.
   */
  virtual void side_postprocess () {}
 
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
    assert(_time_evolving.size() > var);
    _time_evolving[var] = true;
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
  virtual bool mass_residual (bool request_jacobian) {
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
//  virtual std::string system_type () const { return "Differentiable"; }

  /**
   * A pointer to the solver object we're going to use.
   * This must be instantiated by the user before solving!
   */
  AutoPtr<TimeSolver> time_solver;

  /**
   * For time-dependent problems, this is the time t for which the current
   * nonlinear_solution is defined.
   */
  Real time;

  /**
   * For time-dependent problems, this is the amount delta t to advance the
   * solution in time.
   */
  Real deltat;

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

  /**
   * Element by element components of nonlinear_solution
   */
  DenseVector<Number> elem_solution;
  std::vector<DenseSubVector<Number> *> elem_subsolutions;

  /**
   * Element residual vector and Jacobian matrix
   */
  DenseVector<Number> elem_residual;
  DenseMatrix<Number> elem_jacobian;

  /**
   * Element residual subvectors and Jacobian submatrices
   */
  std::vector<DenseSubVector<Number> *> elem_subresiduals;
  std::vector<std::vector<DenseSubMatrix<Number> *> > elem_subjacobians;

  /** 
   * Global Degree of freedom index lists
   */
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<unsigned int> > dof_indices_var;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Stores bools to tell us which variables are evolving
   * in time and which are just constraints
   */
  std::vector<bool> _time_evolving;
};


#endif
