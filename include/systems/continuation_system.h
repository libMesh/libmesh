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



#ifndef LIBMESH_CONTINUATION_SYSTEM_H
#define LIBMESH_CONTINUATION_SYSTEM_H

// Local Includes
#include "libmesh/fem_system.h"

// C++ includes

namespace libMesh
{

// Forward Declarations
template <typename T> class LinearSolver;
class NewtonSolver;

/**
 * This class inherits from the FEMSystem.  It can be
 * used to do arclength continuation.  Most of the ideas and
 * the notation here come from HB Keller's 1977 paper:
 *
 * \verbatim
 * @InProceedings{Kell-1977,
 *   author    = {H.~B.~Keller},
 *   title     = {{Numerical solution of bifurcation and nonlinear eigenvalue problems}},
 *   booktitle = {Applications of Bifurcation Theory, P.~H.~Rabinowitz (ed.)},
 *   year      = 1977,
 *   publisher = {Academic Press},
 *   pages     = {359--389},
 *   notes     = {QA 3 U45 No.\ 38 (PMA)}
 * }
 * \endverbatim
 *
 * \author John W. Peterson
 * \date 2007
 */
class ContinuationSystem : public FEMSystem
{
public:
  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ContinuationSystem (EquationSystems & es,
                      const std::string & name,
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
  virtual void clear () libmesh_override;

  /**
   * Perform a standard "solve" of the system, without doing continuation.
   */
  virtual void solve () libmesh_override;

  /**
   * Perform a continuation solve of the system.  In general, you can only
   * begin the continuation solves after either reading in or solving for two
   * previous values of the control parameter.  The prior two solutions are
   * required for starting up the continuation method.
   */
  void continuation_solve();

  /**
   * Call this function after a continuation solve to compute the tangent and
   * get the next initial guess.
   */
  void advance_arcstep();

  /**
   * The continuation parameter must be a member variable of the
   * derived class, and the "continuation_parameter" pointer defined
   * here must be a pointer to that variable.  This is how the
   * continuation system updates the derived class's continuation
   * parameter.
   *
   * Also sometimes referred to as "lambda" in the code comments.
   */
  Real * continuation_parameter;

  /**
   * If quiet==false, the System prints extra information about what
   * it is doing.
   */
  bool quiet;

  /**
   * Sets (initializes) the max-allowable ds value and the current ds value.
   * Call this before beginning arclength continuation.  The default max stepsize
   * is 0.1
   */
  void set_max_arclength_stepsize(Real maxds) { ds=maxds; ds_current=maxds; }

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
   * How much to try to reduce the residual by at the first (inexact) Newton step.
   * This is frequently something large like 1/2 in an inexact Newton method, to
   * prevent oversolving.
   */
  Real initial_newton_tolerance;

  /**
   * Stores the current solution and continuation parameter
   * (as "previous_u" and "old_continuation_paramter") for later referral.
   * You may need to call this e.g. after the first regular solve, in order
   * to store the first solution, before computing a second solution and
   * beginning arclength continuation.
   */
  void save_current_solution();

  /**
   * The system also keeps track of the old value of the continuation parameter.
   */
  Real old_continuation_parameter;

  /**
   * The minimum-allowable value of the continuation parameter.  The Newton iterations
   * will quit if the continuation parameter falls below this value.
   */
  Real min_continuation_parameter;

  /**
   * The maximum-allowable value of the continuation parameter.  The Newton iterations
   * will quit if the continuation parameter goes above this value.  If this value is zero,
   * there is no maximum value for the continuation parameter.
   */
  Real max_continuation_parameter;

  /**
   * Arclength normalization parameter.  Defaults to 1.0 (no normalization).  Used to
   * ensure that one term in the arclength contstraint equation does not wash out all
   * the others.
   */
  Real Theta;

  /**
   * Another normalization parameter, which is described in the LOCA manual.
   * This one is designed to maintain a relatively "fixed" value of d(lambda)/ds.
   * It is initially 1.0 and is updated after each solve.
   */
  Real Theta_LOCA;

  /**
   * Another scaling parameter suggested by the LOCA people.  This one attempts
   * to shrink the stepsize ds whenever the angle between the previous two
   * tangent vectors gets large.
   */
  //Real tau;

  /**
   * Number of (Newton) backtracking steps to attempt if a Newton step does not
   * reduce the residual.  This is backtracking within a *single* Newton step,
   * if you want to try a smaller arcstep, set n_arclength_reductions > 0.
   */
  unsigned int n_backtrack_steps;

  /**
   * Number of arclength reductions to try when Newton fails to reduce
   * the residual.  For each arclength reduction, the arcstep size is
   * cut in half.
   */
  unsigned int n_arclength_reductions;

  /**
   * The minimum-allowed steplength, defaults to 1.e-8.
   */
  Real ds_min;

  /**
   * The code provides the ability to select from different predictor
   * schemes for getting the initial guess for the solution at the next
   * point on the solution arc.
   */
  enum Predictor {
    /**
     * First-order Euler predictor
     */
    Euler,

    /**
     * Second-order explicit Adams-Bashforth predictor
     */
    AB2,

    /**
     * Invalid predictor
     */
    Invalid_Predictor
  };

  Predictor predictor;

  /**
   * A measure of how rapidly one should attempt to grow the arclength
   * stepsize based on the number of Newton iterations required to solve
   * the problem. Default value is 1.0, if set to zero, will not try to
   * grow or shrink the arclength stepsize based on the number of Newton
   * iterations required.
   */
  Real newton_stepgrowth_aggressiveness;

  /**
   * True by default, the Newton progress check allows the Newton loop
   * to exit if half the allowed iterations have elapsed without a reduction
   * in the *initial* residual.  In our experience this usually means the
   * Newton steps are going to fail eventually and we could save some time
   * by quitting early.
   */
  bool newton_progress_check;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data () libmesh_override;

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
   * A centralized function for setting the normalization parameter theta
   */
  void set_Theta();

  /**
   * A centralized function for setting the other normalization parameter, i.e.
   * the one suggested by the LOCA developers.
   */
  void set_Theta_LOCA();

  /**
   * Applies the predictor to the current solution to get a guess for the
   * next solution.
   */
  void apply_predictor();

  /**
   * Extra work vectors used by the continuation algorithm.
   * These are added to the system by the init_data() routine.
   *
   *
   * The "solution" tangent vector du/ds.
   */
  NumericVector<Number> * du_ds;

  /**
   * The value of du_ds from the previous solution
   */
  NumericVector<Number> * previous_du_ds;

  /**
   * The solution at the previous value of the continuation variable.
   */
  NumericVector<Number> * previous_u;

  /**
   * Temporary vector "y" ... stores -du/dlambda, the solution of \f$ Ay=G_{\lambda} \f$.
   */
  NumericVector<Number> * y;

  /**
   * Temporary vector "y_old" ... stores the previous value of -du/dlambda,
   * which is the solution of \f$ Ay=G_{\lambda} \f$.
   */
  NumericVector<Number> * y_old;

  /**
   * Temporary vector "z" ... the solution of \f$ Az = -G \f$
   */
  NumericVector<Number> * z;

  /**
   * Temporary vector "delta u" ... the Newton step update in our custom
   * augmented PDE solve.
   */
  NumericVector<Number> * delta_u;

  /**
   * We maintain our own linear solver interface, for solving
   * custom systems of equations and/or things which do not require
   * a full-blown NewtonSolver.
   */
  UniquePtr<LinearSolver<Number> > linear_solver;

  /**
   * False until initialize_tangent() is called
   */
  bool tangent_initialized;

  /**
   * A pointer to the underlying Newton solver used by the \p DiffSystem.
   * From this pointer, we can get access to all the parameters and options
   * which are available to the "normal" Newton solver.
   */
  NewtonSolver * newton_solver;

  /**
   * The most recent value of the derivative of the continuation parameter
   * with respect to s.  We use "lambda" here for shortness of notation, lambda
   * always refers to the continuation parameter.
   */
  Real dlambda_ds;

  /**
   * The initial arclength stepsize, selected by the user.  This is
   * the max-allowable arclength stepsize, but the algorithm may frequently
   * reduce ds near turning points.
   */
  Real ds;

  /**
   * Value of stepsize currently in use.  Will not exceed user-provided maximum
   * arclength stepize ds.
   */
  Real ds_current;

  /**
   * The old parameter tangent value.
   */
  Real previous_dlambda_ds;

  /**
   * The previous arcstep length used.
   */
  Real previous_ds;

  /**
   * Loop counter for nonlinear (Newton) iteration loop.
   */
  unsigned int newton_step;
};

} // namespace libMesh

#endif // LIBMESH_CONTINUATION_SYSTEM_H
