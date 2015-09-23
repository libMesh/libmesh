// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// C++ includes

// Local Includes
#include "libmesh/optimization_solver.h"
#include "libmesh/tao_optimization_solver.h"
#include "libmesh/nlopt_optimization_solver.h"

namespace libMesh
{

template <typename T>
inline
OptimizationSolver<T>::OptimizationSolver (sys_type& s) :
  ParallelObject(s),
  objective_object(NULL),
  gradient_object(NULL),
  hessian_object(NULL),
  equality_constraints_object(NULL),
  equality_constraints_jacobian_object(NULL),
  inequality_constraints_object(NULL),
  inequality_constraints_jacobian_object(NULL),
  lower_and_upper_bounds_object(NULL),
  max_objective_function_evaluations(500),
  objective_function_relative_tolerance(1.e-4),
  verbose(false),
  _system(s),
  _is_initialized (false)
{
}



template <typename T>
inline
OptimizationSolver<T>::~OptimizationSolver ()
{
  this->clear ();
}


template <typename T>
UniquePtr<OptimizationSolver<T> >
OptimizationSolver<T>::build(sys_type& s, const SolverPackage solver_package)
{
  // Prevent unused variables warnings when Tao is not available
  libmesh_ignore(s);

  // Build the appropriate solver
  switch (solver_package)
    {

#if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
    case PETSC_SOLVERS:
      return UniquePtr<OptimizationSolver<T> >(new TaoOptimizationSolver<T>(s));
#endif // #if defined(LIBMESH_HAVE_PETSC_TAO) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)

#if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
    case NLOPT_SOLVERS:
      return UniquePtr<OptimizationSolver<T> >(new NloptOptimizationSolver<T>(s));
#endif // #if defined(LIBMESH_HAVE_NLOPT) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<OptimizationSolver<T> >();
}


//------------------------------------------------------------------
// Explicit instantiations
template class OptimizationSolver<Number>;

} // namespace libMesh
