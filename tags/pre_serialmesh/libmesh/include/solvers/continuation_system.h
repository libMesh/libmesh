// $Id: continuation_system.h,v 1.1 2007-08-02 17:27:18 jwpeterson Exp $

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



#ifndef __continuation_system_h__
#define __continuation_system_h__

// C++ includes

// Local Includes
#include "fem_system.h"

// Forward Declarations
template <typename T> class LinearSolver;
class NewtonSolver;





/**
 * This class inherits from the FEMSystem.  It can be
 * used to do arclength continuation.  Most of the ideas and
 * the notation here come from HB Keller's 1977 paper:
 *
 * @InProceedings{Kell-1977,
 *  author = {H.~B.~Keller},
 *  title = {{Numerical solution of bifurcation and nonlinear eigenvalue problems}},
 *  booktitle = {Applications of Bifurcation Theory, P.~H.~Rabinowitz (ed.)},
 *  year = 1977,
 *  publisher = {Academic Press},
 *  pages = {359--389},
 *  notes = {QA 3 U45 No.\ 38 (PMA)}
 * }
 *
 * @author John W. Peterson 2007
 */
class ContinuationSystem : public FEMSystem
{
public:
  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ContinuationSystem (EquationSystems& es,
		      const std::string& name,
		      const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~ContinuationSystem ();

  /**
   * The type of system.
   */
  typedef ContinuationSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef FEMSystem Parent;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Perform a standard "solve" of the system, without doing continuation.
   */
  virtual void solve();
  
  /**
   * Perform a continuation solve of the system.  In general, you can only
   * begin the continuation solves after either reading in or solving for two
   * previous values of the control parameter.  The prior two solutions are
   * required for starting up the continuation method.
   */
  void continuation_solve();
  
  /**
   * The continuation parameter must be a member variable of the
   * derived class, and the "continuation_parameter" pointer defined
   * here must be a pointer to that variable.  This is how the
   * continuation system updates the derived class's continuation
   * parameter.
   *
   * Also sometimes referred to as "lambda" in the code comments.
   */
  Real* continuation_parameter;

  /**
   * If quiet==false, the System prints extra information about what
   * it is doing.
   */
  bool quiet;

  /**
   * The arclength step size.
   */
  Real ds;

  /**
   * How tightly should the Newton iterations attempt to converge delta_lambda.
   * Defaults to 1.e-6.
   */
  Real continuation_parameter_tolerance;

  /**
   * How tightly should the Newton iterations attempt to converge ||delta_u||
   * Defaults to 1.e-6.
   */
  Real solution_tolerance;

  /**
   * Stores the current solution and continuation parameter
   * (as "previous_u" and "old_continuation_paramter") for later referral.
   * You may need to call this e.g. after the first regular solve, in order
   * to store the first solution, before computing a second solution and
   * beginning arclength continuation.
   */
  void save_current_solution();
  
protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * There are (so far) two different vectors which may be assembled using the assembly routine:
   * 1.) The Residual = the normal PDE weighted residual
   * 2.) G_Lambda     = the derivative wrt the parameter lambda of the PDE weighted residual
   *
   * It is up to the derived class to handle writing separate assembly code for the different cases.
   * Usually something like:
   * switch (rhs_mode)
   * {
   *   case Residual:
   *   {
   *     Fu(i) +=  ... // normal PDE residual
   *     break;
   *   }
   * 	      
   * case G_Lambda:
   *   {
   *     Fu(i) += ... // derivative wrt control parameter
   *     break;
   *   }
   */
  enum RHS_Mode {Residual,
		 G_Lambda};
  
  RHS_Mode rhs_mode;



private:
  /**
   * Before starting arclength continuation, we need at least 2 prior solutions
   * (both solution and u_previous should be filled with meaningful values)
   * And we need to initialize the tangent vector.  This only needs to be
   * called once.
   */
  void initialize_tangent();

  /**
   * Special solve algorithm for solving the tangent system.
   */
  void solve_tangent();

  /**
   * This function (which assumes the most recent tangent vector has been
   * computed) updates the solution and the control parameter with the initial
   * guess for the next point on the continuation path.
   */
  void update_solution();

  /**
   * Extra work vectors used by the continuation algorithm.
   * These are added to the system by the init_data() routine.
   *
   *
   * The "solution" tangent vector du/ds.
   */
  NumericVector<Number>* du_ds;

  /**
   * The solution at the previous value of the continuation variable.
   */ 
  NumericVector<Number>* previous_u;

  /**
   * Temporary vector "y" ... the solution of Ay=G_{\lambda}.
   */
  NumericVector<Number>* y;

  /**
   * Temporary vector "z" ... the solution of Az = -G
   */
  NumericVector<Number>* z;

  /**
   * Temporary vector "delta u" ... the Newton step update in our custom
   * augmented PDE solve.
   */
  NumericVector<Number>* delta_u;

  /**
   * We maintain our own linear solver interface, for solving 
   * custom systems of equations and/or things which do not require
   * a full-blown NewtonSolver.
   */ 
  AutoPtr<LinearSolver<Number> > linear_solver;
  
  /**
   * False until initialize_tangent() is called
   */
  bool tangent_initialized;

  /**
   * A pointer to the underlying Newton solver used by the DiffSystem.
   * From this pointer, we can get access to all the parameters and options
   * which are available to the "normal" Newton solver.
   */
  NewtonSolver* newton_solver;
  
  /**
   * The system also keeps track of the old value of the continuation parameter.
   */
  Real old_continuation_parameter;

  /**
   * The most recent value of the derivative of the continuation parameter
   * with respect to s.  We use "lambda" here for shortness of notation, lambda
   * always refers to the continuation parameter.
   */ 
  Real dlambda_ds;
};

#endif
