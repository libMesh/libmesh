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



#include "libmesh/libmesh_common.h"

#if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)


// C++ includes

// Local Includes
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/tao_optimization_solver.h"
#include "libmesh/equation_systems.h"

namespace libMesh
{

//--------------------------------------------------------------------
// Functions with C linkage to pass to Tao. Tao will call these
// methods as needed.
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{

  //---------------------------------------------------------------
  // This function is called by Tao to evaluate objective function at x
  PetscErrorCode
  __libmesh_tao_objective (Tao /*tao*/, Vec x, PetscReal * objective, void * ctx)
  {
    LOG_SCOPE("objective()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(objective);
    libmesh_assert(ctx);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    if (solver->objective_object != libmesh_nullptr)
      (*objective) = solver->objective_object->objective(*(sys.current_local_solution), sys);
    else
      libmesh_error_msg("Objective function not defined in __libmesh_tao_objective");

    return ierr;
  }



  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the gradient at x
  PetscErrorCode
  __libmesh_tao_gradient(Tao /*tao*/, Vec x, Vec g, void * ctx)
  {
    LOG_SCOPE("gradient()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(g);
    libmesh_assert(ctx);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // We'll also pass the gradient in to the assembly routine
    // so let's make a PETSc vector for that too.
    PetscVector<Number> gradient(g, sys.comm());

    // Clear the gradient prior to assembly
    gradient.zero();

    if (solver->gradient_object != libmesh_nullptr)
      solver->gradient_object->gradient(*(sys.current_local_solution), gradient, sys);
    else
      libmesh_error_msg("Gradient function not defined in __libmesh_tao_gradient");

    gradient.close();

    return ierr;
  }

  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the Hessian at x
  PetscErrorCode
  __libmesh_tao_hessian(Tao /*tao*/, Vec x, Mat h, Mat pc, void * ctx)
  {
    LOG_SCOPE("hessian()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(h);
    libmesh_assert(pc);
    libmesh_assert(ctx);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // Let's also wrap pc and h in PetscMatrix objects for convenience
    PetscMatrix<Number> PC(pc, sys.comm());
    PetscMatrix<Number> hessian(h, sys.comm());
    PC.attach_dof_map(sys.get_dof_map());
    hessian.attach_dof_map(sys.get_dof_map());

    if (solver->hessian_object != libmesh_nullptr)
      {
        // Following PetscNonlinearSolver by passing in PC. It's not clear
        // why we pass in PC and not hessian though?
        solver->hessian_object->hessian(*(sys.current_local_solution), PC, sys);
      }
    else
      libmesh_error_msg("Hessian function not defined in __libmesh_tao_hessian");

    PC.close();
    hessian.close();

    return ierr;
  }


  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the equality constraints at x
  PetscErrorCode
  __libmesh_tao_equality_constraints(Tao /*tao*/, Vec x, Vec ce, void * ctx)
  {
    LOG_SCOPE("equality_constraints()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(ce);
    libmesh_assert(ctx);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // We'll also pass the constraints vector ce into the assembly routine
    // so let's make a PETSc vector for that too.
    PetscVector<Number> eq_constraints(ce, sys.comm());

    // Clear the gradient prior to assembly
    eq_constraints.zero();

    if (solver->equality_constraints_object != libmesh_nullptr)
      solver->equality_constraints_object->equality_constraints(*(sys.current_local_solution), eq_constraints, sys);
    else
      libmesh_error_msg("Constraints function not defined in __libmesh_tao_equality_constraints");

    eq_constraints.close();

    return ierr;
  }

  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the Jacobian of the
  // equality constraints at x
  PetscErrorCode
  __libmesh_tao_equality_constraints_jacobian(Tao /*tao*/, Vec x, Mat J, Mat Jpre, void * ctx)
  {
    LOG_SCOPE("equality_constraints_jacobian()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(J);
    libmesh_assert(Jpre);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // Let's also wrap J and Jpre in PetscMatrix objects for convenience
    PetscMatrix<Number> J_petsc(J, sys.comm());
    PetscMatrix<Number> Jpre_petsc(Jpre, sys.comm());

    if (solver->equality_constraints_jacobian_object != libmesh_nullptr)
      solver->equality_constraints_jacobian_object->equality_constraints_jacobian(*(sys.current_local_solution), J_petsc, sys);
    else
      libmesh_error_msg("Constraints function not defined in __libmesh_tao_equality_constraints_jacobian");

    J_petsc.close();
    Jpre_petsc.close();

    return ierr;
  }

  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the inequality constraints at x
  PetscErrorCode
  __libmesh_tao_inequality_constraints(Tao /*tao*/, Vec x, Vec cineq, void * ctx)
  {
    LOG_SCOPE("inequality_constraints()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(cineq);
    libmesh_assert(ctx);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // We'll also pass the constraints vector ce into the assembly routine
    // so let's make a PETSc vector for that too.
    PetscVector<Number> ineq_constraints(cineq, sys.comm());

    // Clear the gradient prior to assembly
    ineq_constraints.zero();

    if (solver->inequality_constraints_object != libmesh_nullptr)
      solver->inequality_constraints_object->inequality_constraints(*(sys.current_local_solution), ineq_constraints, sys);
    else
      libmesh_error_msg("Constraints function not defined in __libmesh_tao_inequality_constraints");

    ineq_constraints.close();

    return ierr;
  }

  //---------------------------------------------------------------
  // This function is called by Tao to evaluate the Jacobian of the
  // equality constraints at x
  PetscErrorCode
  __libmesh_tao_inequality_constraints_jacobian(Tao /*tao*/, Vec x, Mat J, Mat Jpre, void * ctx)
  {
    LOG_SCOPE("inequality_constraints_jacobian()", "TaoOptimizationSolver");

    PetscErrorCode ierr = 0;

    libmesh_assert(x);
    libmesh_assert(J);
    libmesh_assert(Jpre);

    // ctx should be a pointer to the solver (it was passed in as void *)
    TaoOptimizationSolver<Number> * solver =
      static_cast<TaoOptimizationSolver<Number> *> (ctx);

    OptimizationSystem & sys = solver->system();

    // We'll use current_local_solution below, so let's ensure that it's consistent
    // with the vector x that was passed in.
    PetscVector<Number> & X_sys = *cast_ptr<PetscVector<Number> *>(sys.solution.get());
    PetscVector<Number> X(x, sys.comm());

    // Perform a swap so that sys.solution points to X
    X.swap(X_sys);
    // Impose constraints on X
    sys.get_dof_map().enforce_constraints_exactly(sys);
    // Update sys.current_local_solution based on X
    sys.update();
    // Swap back
    X.swap(X_sys);

    // Let's also wrap J and Jpre in PetscMatrix objects for convenience
    PetscMatrix<Number> J_petsc(J, sys.comm());
    PetscMatrix<Number> Jpre_petsc(Jpre, sys.comm());

    if (solver->inequality_constraints_jacobian_object != libmesh_nullptr)
      solver->inequality_constraints_jacobian_object->inequality_constraints_jacobian(*(sys.current_local_solution), J_petsc, sys);
    else
      libmesh_error_msg("Constraints function not defined in __libmesh_tao_inequality_constraints_jacobian");

    J_petsc.close();
    Jpre_petsc.close();

    return ierr;
  }

} // end extern "C"
//---------------------------------------------------------------------



//---------------------------------------------------------------------
// TaoOptimizationSolver<> methods
template <typename T>
TaoOptimizationSolver<T>::TaoOptimizationSolver (OptimizationSystem & system_in) :
  OptimizationSolver<T>(system_in),
  _reason(TAO_CONVERGED_USER) // Arbitrary initial value...
{
}



template <typename T>
TaoOptimizationSolver<T>::~TaoOptimizationSolver ()
{
  this->clear ();
}



template <typename T>
void TaoOptimizationSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      PetscErrorCode ierr=0;

      ierr = TaoDestroy(&_tao);
      LIBMESH_CHKERR(ierr);
    }
}



template <typename T>
void TaoOptimizationSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      PetscErrorCode ierr=0;

      ierr = TaoCreate(this->comm().get(),&_tao);
      LIBMESH_CHKERR(ierr);
    }
}

template <typename T>
void TaoOptimizationSolver<T>::solve ()
{
  LOG_SCOPE("solve()", "TaoOptimizationSolver");

  this->init ();

  this->system().solution->zero();

  PetscMatrix<T> * hessian  = cast_ptr<PetscMatrix<T> *>(this->system().matrix);
  // PetscVector<T> * gradient = cast_ptr<PetscVector<T> *>(this->system().rhs);
  PetscVector<T> * x         = cast_ptr<PetscVector<T> *>(this->system().solution.get());
  PetscVector<T> * ceq       = cast_ptr<PetscVector<T> *>(this->system().C_eq.get());
  PetscMatrix<T> * ceq_jac   = cast_ptr<PetscMatrix<T> *>(this->system().C_eq_jac.get());
  PetscVector<T> * cineq     = cast_ptr<PetscVector<T> *>(this->system().C_ineq.get());
  PetscMatrix<T> * cineq_jac = cast_ptr<PetscMatrix<T> *>(this->system().C_ineq_jac.get());
  PetscVector<T> * lb        = cast_ptr<PetscVector<T> *>(&this->system().get_vector("lower_bounds"));
  PetscVector<T> * ub        = cast_ptr<PetscVector<T> *>(&this->system().get_vector("upper_bounds"));

  // Set the starting guess to zero.
  x->zero();

  PetscErrorCode ierr = 0;

  // Workaround for bug where TaoSetFromOptions *reset*
  // programmatically set tolerance and max. function evaluation
  // values when "-tao_type ipm" was specified on the command line: we
  // call TaoSetFromOptions twice (both before and after setting
  // custom options programmatically)
  ierr = TaoSetFromOptions(_tao);
  LIBMESH_CHKERR(ierr);

  // Set convergence tolerances
  // f(X) - f(X*) (estimated)            <= fatol
  // |f(X) - f(X*)| (estimated) / |f(X)| <= frtol
  // ||g(X)||                            <= gatol
  // ||g(X)|| / |f(X)|                   <= grtol
  // ||g(X)|| / ||g(X0)||                <= gttol
  // Command line equivalents: -tao_fatol, -tao_frtol, -tao_gatol, -tao_grtol, -tao_gttol
  ierr = TaoSetTolerances(_tao,
#if PETSC_RELEASE_LESS_THAN(3,7,0)
                          // Releases up to 3.X.Y had fatol and frtol, after that they were removed.
                          // Hopefully we'll be able to know X and Y soon. Guessing at 3.7.0.
                          /*fatol=*/PETSC_DEFAULT,
                          /*frtol=*/PETSC_DEFAULT,
#endif
                          /*gatol=*/PETSC_DEFAULT,
                          /*grtol=*/this->objective_function_relative_tolerance,
                          /*gttol=*/PETSC_DEFAULT);
  LIBMESH_CHKERR(ierr);

  // Set the max-allowed number of objective function evaluations
  // Command line equivalent: -tao_max_funcs
  ierr = TaoSetMaximumFunctionEvaluations(_tao, this->max_objective_function_evaluations);
  LIBMESH_CHKERR(ierr);

  // Set the max-allowed number of optimization iterations.
  // Command line equivalent: -tao_max_it
  // Not implemented for now as it seems fairly similar to
  // ierr = TaoSetMaximumIterations(_tao, 4);
  // LIBMESH_CHKERR(ierr);

  // Set solution vec and an initial guess
  ierr = TaoSetInitialVector(_tao, x->vec());
  LIBMESH_CHKERR(ierr);

  // We have to have an objective function
  libmesh_assert( this->objective_object );

  // Set routines for objective, gradient, hessian evaluation
  ierr = TaoSetObjectiveRoutine(_tao, __libmesh_tao_objective, this);
  LIBMESH_CHKERR(ierr);

  if (this->gradient_object)
    {
      ierr = TaoSetGradientRoutine(_tao, __libmesh_tao_gradient, this);
      LIBMESH_CHKERR(ierr);
    }

  if (this->hessian_object)
    {
      ierr = TaoSetHessianRoutine(_tao, hessian->mat(), hessian->mat(), __libmesh_tao_hessian, this);
      LIBMESH_CHKERR(ierr);
    }

  if (this->lower_and_upper_bounds_object)
    {
      // Need to actually compute the bounds vectors first
      this->lower_and_upper_bounds_object->lower_and_upper_bounds(this->system());

      ierr = TaoSetVariableBounds(_tao,
                                  lb->vec(),
                                  ub->vec());
      LIBMESH_CHKERR(ierr);
    }

  if (this->equality_constraints_object)
    {
      ierr = TaoSetEqualityConstraintsRoutine(_tao, ceq->vec(), __libmesh_tao_equality_constraints, this);
      LIBMESH_CHKERR(ierr);
    }

  if (this->equality_constraints_jacobian_object)
    {
      ierr = TaoSetJacobianEqualityRoutine(_tao,
                                           ceq_jac->mat(),
                                           ceq_jac->mat(),
                                           __libmesh_tao_equality_constraints_jacobian,
                                           this);
      LIBMESH_CHKERR(ierr);
    }

  // Optionally set inequality constraints
  if (this->inequality_constraints_object)
    {
      ierr = TaoSetInequalityConstraintsRoutine(_tao, cineq->vec(), __libmesh_tao_inequality_constraints, this);
      LIBMESH_CHKERR(ierr);
    }

  // Optionally set inequality constraints Jacobian
  if (this->inequality_constraints_jacobian_object)
    {
      ierr = TaoSetJacobianInequalityRoutine(_tao,
                                             cineq_jac->mat(),
                                             cineq_jac->mat(),
                                             __libmesh_tao_inequality_constraints_jacobian,
                                             this);
      LIBMESH_CHKERR(ierr);
    }

  // Check for Tao command line options
  ierr = TaoSetFromOptions(_tao);
  LIBMESH_CHKERR(ierr);

  // Perform the optimization
  ierr = TaoSolve(_tao);
  LIBMESH_CHKERR(ierr);

  // Store the convergence/divergence reason
  ierr = TaoGetConvergedReason(_tao, &_reason);
  LIBMESH_CHKERR(ierr);
}


template <typename T>
void TaoOptimizationSolver<T>::get_dual_variables()
{
  LOG_SCOPE("get_dual_variables()", "TaoOptimizationSolver");

  PetscVector<T> * lambda_eq_petsc =
    cast_ptr<PetscVector<T> *>(this->system().lambda_eq.get());
  PetscVector<T> * lambda_ineq_petsc =
    cast_ptr<PetscVector<T> *>(this->system().lambda_ineq.get());

  Vec lambda_eq_petsc_vec = lambda_eq_petsc->vec();
  Vec lambda_ineq_petsc_vec = lambda_ineq_petsc->vec();

  PetscErrorCode ierr = 0;
  ierr = TaoGetDualVariables(_tao,
                             &lambda_eq_petsc_vec,
                             &lambda_ineq_petsc_vec);
  LIBMESH_CHKERR(ierr);
}


template <typename T>
void TaoOptimizationSolver<T>::print_converged_reason()
{
  libMesh::out << "Tao optimization solver convergence/divergence reason: "
               << TaoConvergedReasons[this->get_converged_reason()] << std::endl;
}



template <typename T>
int TaoOptimizationSolver<T>::get_converged_reason()
{
  PetscErrorCode ierr=0;

  if (this->initialized())
    {
      ierr = TaoGetConvergedReason(_tao, &_reason);
      LIBMESH_CHKERR(ierr);
    }

  return static_cast<int>(_reason);
}


//------------------------------------------------------------------
// Explicit instantiations
template class TaoOptimizationSolver<Number>;

} // namespace libMesh



#endif // #if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
