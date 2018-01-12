// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)


// C++ includes

// Local Includes
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nlopt_optimization_solver.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

double __libmesh_nlopt_objective(unsigned n,
                                 const double * x,
                                 double * gradient,
                                 void * data)
{
  LOG_SCOPE("objective()", "NloptOptimizationSolver");

  // ctx should be a pointer to the solver (it was passed in as void *)
  NloptOptimizationSolver<Number> * solver =
    static_cast<NloptOptimizationSolver<Number> *> (data);

  OptimizationSystem & sys = solver->system();

  // We'll use current_local_solution below, so let's ensure that it's consistent
  // with the vector x that was passed in.
  for (unsigned int i=sys.solution->first_local_index();
       i<sys.solution->last_local_index(); i++)
    sys.solution->set(i, x[i]);

  // Make sure the solution vector is parallel-consistent
  sys.solution->close();

  // Impose constraints on X
  sys.get_dof_map().enforce_constraints_exactly(sys);

  // Update sys.current_local_solution based on X
  sys.update();

  Real objective;
  if (solver->objective_object != libmesh_nullptr)
    {
      objective =
        solver->objective_object->objective(
                                            *(sys.current_local_solution), sys);
    }
  else
    {
      libmesh_error_msg("Objective function not defined in __libmesh_nlopt_objective");
    }

  // If the gradient has been requested, fill it in
  if (gradient)
    {
      if (solver->gradient_object != libmesh_nullptr)
        {
          solver->gradient_object->gradient(
                                            *(sys.current_local_solution), *(sys.rhs), sys);

          // we've filled up sys.rhs with the gradient data, now copy it
          // to the nlopt data structure
          libmesh_assert(sys.rhs->size() == n);

          std::vector<double> grad;
          sys.rhs->localize_to_one(grad);
          for (unsigned i=0; i<n; ++i)
            gradient[i] = grad[i];
        }
      else
        libmesh_error_msg("Gradient function not defined in __libmesh_nlopt_objective");
    }

  // Increment the iteration count.
  solver->get_iteration_count()++;

  // Possibly print the current value of the objective function
  if (solver->verbose)
    libMesh::out << objective << std::endl;

  return objective;
}



void __libmesh_nlopt_equality_constraints(unsigned m,
                                          double * result,
                                          unsigned n,
                                          const double * x,
                                          double * gradient,
                                          void * data)
{
  LOG_SCOPE("equality_constraints()", "NloptOptimizationSolver");

  libmesh_assert(data);

  // data should be a pointer to the solver (it was passed in as void *)
  NloptOptimizationSolver<Number> * solver =
    static_cast<NloptOptimizationSolver<Number> *> (data);

  OptimizationSystem & sys = solver->system();

  // We'll use current_local_solution below, so let's ensure that it's consistent
  // with the vector x that was passed in.
  if (sys.solution->size() != n)
    libmesh_error_msg("Error: Input vector x has different length than sys.solution!");

  for (unsigned int i=sys.solution->first_local_index(); i<sys.solution->last_local_index(); i++)
    sys.solution->set(i, x[i]);
  sys.solution->close();

  // Impose constraints on the solution vector
  sys.get_dof_map().enforce_constraints_exactly(sys);

  // Update sys.current_local_solution based on the solution vector
  sys.update();

  // Call the user's equality constraints function if there is one.
  OptimizationSystem::ComputeEqualityConstraints * eco = solver->equality_constraints_object;
  if (eco)
    {
      eco->equality_constraints(*sys.current_local_solution,
                                *sys.C_eq,
                                sys);

      sys.C_eq->close();

      // Copy the values out of eq_constraints into 'result'.
      // TODO: Even better would be if we could use 'result' directly
      // as the storage of eq_constraints.  Perhaps a serial-only
      // NumericVector variant which supports this option?
      for (unsigned i=0; i<m; ++i)
        result[i] = (*sys.C_eq)(i);

      // If gradient != NULL, then the Jacobian matrix of the equality
      // constraints has been requested.  The incoming 'gradient'
      // array is of length m*n and d(c_i)/d(x_j) = gradient[n*i+j].
      if (gradient)
        {
          OptimizationSystem::ComputeEqualityConstraintsJacobian * eco_jac =
            solver->equality_constraints_jacobian_object;

          if (eco_jac)
            {
              eco_jac->equality_constraints_jacobian(*sys.current_local_solution,
                                                     *sys.C_eq_jac,
                                                     sys);

              sys.C_eq_jac->close();

              // copy the Jacobian data to the gradient array
              for (numeric_index_type i=0; i<m; i++)
                {
                  std::set<numeric_index_type>::iterator it = sys.eq_constraint_jac_sparsity[i].begin();
                  std::set<numeric_index_type>::iterator it_end = sys.eq_constraint_jac_sparsity[i].end();
                  for ( ; it != it_end; ++it)
                    {
                      numeric_index_type dof_index = *it;
                      gradient[n*i+dof_index] = (*sys.C_eq_jac)(i,dof_index);
                    }
                }
            }
          else
            libmesh_error_msg("Jacobian function not defined in __libmesh_nlopt_equality_constraints");
        }

    }
  else
    libmesh_error_msg("Constraints function not defined in __libmesh_nlopt_equality_constraints");
}


void __libmesh_nlopt_inequality_constraints(unsigned m,
                                            double * result,
                                            unsigned n,
                                            const double * x,
                                            double * gradient,
                                            void * data)
{
  LOG_SCOPE("inequality_constraints()", "NloptOptimizationSolver");

  libmesh_assert(data);

  // data should be a pointer to the solver (it was passed in as void *)
  NloptOptimizationSolver<Number> * solver =
    static_cast<NloptOptimizationSolver<Number> *> (data);

  OptimizationSystem & sys = solver->system();

  // We'll use current_local_solution below, so let's ensure that it's consistent
  // with the vector x that was passed in.
  if (sys.solution->size() != n)
    libmesh_error_msg("Error: Input vector x has different length than sys.solution!");

  for (unsigned int i=sys.solution->first_local_index(); i<sys.solution->last_local_index(); i++)
    sys.solution->set(i, x[i]);
  sys.solution->close();

  // Impose constraints on the solution vector
  sys.get_dof_map().enforce_constraints_exactly(sys);

  // Update sys.current_local_solution based on the solution vector
  sys.update();

  // Call the user's inequality constraints function if there is one.
  OptimizationSystem::ComputeInequalityConstraints * ineco = solver->inequality_constraints_object;
  if (ineco)
    {
      ineco->inequality_constraints(*sys.current_local_solution,
                                    *sys.C_ineq,
                                    sys);

      sys.C_ineq->close();

      // Copy the values out of ineq_constraints into 'result'.
      // TODO: Even better would be if we could use 'result' directly
      // as the storage of ineq_constraints.  Perhaps a serial-only
      // NumericVector variant which supports this option?
      for (unsigned i=0; i<m; ++i)
        result[i] = (*sys.C_ineq)(i);

      // If gradient != NULL, then the Jacobian matrix of the equality
      // constraints has been requested.  The incoming 'gradient'
      // array is of length m*n and d(c_i)/d(x_j) = gradient[n*i+j].
      if (gradient)
        {
          OptimizationSystem::ComputeInequalityConstraintsJacobian * ineco_jac =
            solver->inequality_constraints_jacobian_object;

          if (ineco_jac)
            {
              ineco_jac->inequality_constraints_jacobian(*sys.current_local_solution,
                                                         *sys.C_ineq_jac,
                                                         sys);

              sys.C_ineq_jac->close();

              // copy the Jacobian data to the gradient array
              for (numeric_index_type i=0; i<m; i++)
                {
                  std::set<numeric_index_type>::iterator it = sys.ineq_constraint_jac_sparsity[i].begin();
                  std::set<numeric_index_type>::iterator it_end = sys.ineq_constraint_jac_sparsity[i].end();
                  for ( ; it != it_end; ++it)
                    {
                      numeric_index_type dof_index = *it;
                      gradient[n*i+dof_index] = (*sys.C_ineq_jac)(i,dof_index);
                    }
                }
            }
          else
            libmesh_error_msg("Jacobian function not defined in __libmesh_nlopt_inequality_constraints");
        }

    }
  else
    libmesh_error_msg("Constraints function not defined in __libmesh_nlopt_inequality_constraints");
}

//---------------------------------------------------------------------



//---------------------------------------------------------------------
// NloptOptimizationSolver<> methods
template <typename T>
NloptOptimizationSolver<T>::NloptOptimizationSolver (OptimizationSystem & system_in)
  :
  OptimizationSolver<T>(system_in),
  _opt(libmesh_nullptr),
  _result(NLOPT_SUCCESS),
  _iteration_count(0),
  _constraints_tolerance(1.e-8)
{
}



template <typename T>
NloptOptimizationSolver<T>::~NloptOptimizationSolver ()
{
  this->clear ();
}



template <typename T>
void NloptOptimizationSolver<T>::clear ()
{
  if (this->initialized())
    {
      this->_is_initialized = false;

      nlopt_destroy(_opt);
    }
}



template <typename T>
void NloptOptimizationSolver<T>::init ()
{
  // Initialize the data structures if not done so already.
  if (!this->initialized())
    {
      this->_is_initialized = true;

      // By default, use the LD_SLSQP solver
      std::string nlopt_algorithm_name = "LD_SLSQP";

      if (libMesh::on_command_line("--nlopt-algorithm"))
        nlopt_algorithm_name = libMesh::command_line_next ("--nlopt-algorithm",
                                                           nlopt_algorithm_name);

      // Convert string to an nlopt algorithm type
      std::map<std::string, nlopt_algorithm>::iterator it =
        _nlopt_algorithms.find(nlopt_algorithm_name);

      if (it == _nlopt_algorithms.end())
        libmesh_error_msg("Invalid nlopt algorithm requested on command line: " \
                          << nlopt_algorithm_name);

      _opt = nlopt_create(it->second, this->system().solution->size());
    }
}



template <typename T>
void NloptOptimizationSolver<T>::solve ()
{
  LOG_SCOPE("solve()", "NloptOptimizationSolver");

  this->init ();

  unsigned int nlopt_size = this->system().solution->size();

  // We have to have an objective function
  libmesh_assert( this->objective_object );

  // Set routine for objective and (optionally) gradient evaluation
  {
    nlopt_result ierr =
      nlopt_set_min_objective(_opt,
                              __libmesh_nlopt_objective,
                              this);
    if (ierr < 0)
      libmesh_error_msg("NLopt failed to set min objective: " << ierr);
  }

  if (this->lower_and_upper_bounds_object)
    {
      // Need to actually compute the bounds vectors first
      this->lower_and_upper_bounds_object->lower_and_upper_bounds(this->system());

      std::vector<Real> nlopt_lb(nlopt_size);
      std::vector<Real> nlopt_ub(nlopt_size);
      for (unsigned int i=0; i<nlopt_size; i++)
        {
          nlopt_lb[i] = this->system().get_vector("lower_bounds")(i);
          nlopt_ub[i] = this->system().get_vector("upper_bounds")(i);
        }

      nlopt_set_lower_bounds(_opt, &nlopt_lb[0]);
      nlopt_set_upper_bounds(_opt, &nlopt_ub[0]);
    }

  // If we have an equality constraints object, tell NLopt about it.
  if (this->equality_constraints_object)
    {
      // NLopt requires a vector to specify the tolerance for each constraint.
      // NLopt makes a copy of this vector internally, so it's safe for us to
      // let it go out of scope.
      std::vector<double> equality_constraints_tolerances(this->system().C_eq->size(),
                                                          _constraints_tolerance);

      // It would be nice to call the C interface directly, at least it should return an error
      // code we could parse... unfortunately, there does not seem to be a way to extract
      // the underlying nlopt_opt object from the nlopt::opt class!
      nlopt_result ierr =
        nlopt_add_equality_mconstraint(_opt,
                                       equality_constraints_tolerances.size(),
                                       __libmesh_nlopt_equality_constraints,
                                       this,
                                       &equality_constraints_tolerances[0]);

      if (ierr < 0)
        libmesh_error_msg("NLopt failed to add equality constraint: " << ierr);
    }

  // If we have an inequality constraints object, tell NLopt about it.
  if (this->inequality_constraints_object)
    {
      // NLopt requires a vector to specify the tolerance for each constraint
      std::vector<double> inequality_constraints_tolerances(this->system().C_ineq->size(),
                                                            _constraints_tolerance);

      nlopt_add_inequality_mconstraint(_opt,
                                       inequality_constraints_tolerances.size(),
                                       __libmesh_nlopt_inequality_constraints,
                                       this,
                                       &inequality_constraints_tolerances[0]);
    }

  // Set a relative tolerance on the optimization parameters
  nlopt_set_ftol_rel(_opt, this->objective_function_relative_tolerance);

  // Set the maximum number of allowed objective function evaluations
  nlopt_set_maxeval(_opt, this->max_objective_function_evaluations);

  // Reset internal iteration counter
  this->_iteration_count = 0;

  // Perform the optimization
  std::vector<Real> x(nlopt_size);
  Real min_val = 0.;
  _result = nlopt_optimize(_opt, &x[0], &min_val);

  if (_result < 0)
    libMesh::out << "NLopt failed!" << std::endl;
  else
    libMesh::out << "NLopt obtained optimal value: "
                 << min_val
                 << " in "
                 << this->get_iteration_count()
                 << " iterations."
                 << std::endl;
}


template <typename T>
void NloptOptimizationSolver<T>::print_converged_reason()
{
  libMesh::out << "NLopt optimization solver convergence/divergence reason: "
               << this->get_converged_reason() << std::endl;
}



template <typename T>
int NloptOptimizationSolver<T>::get_converged_reason()
{
  return static_cast<int>(_result);
}


//------------------------------------------------------------------
// Explicit instantiations
template class NloptOptimizationSolver<Number>;

} // namespace libMesh



#endif // #if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
