// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_SLEPC

// Local Includes
#include "libmesh/eigen_solver.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/enum_eigen_solver_type.h"

namespace libMesh
{


//------------------------------------------------------------------
// EigenSolver members
template <typename T>
EigenSolver<T>::EigenSolver (const Parallel::Communicator & comm_in) :
  ParallelObject(comm_in),
  _eigen_solver_type    (ARNOLDI),
  _eigen_problem_type   (NHEP),
  _position_of_spectrum (LARGEST_MAGNITUDE),
  _is_initialized       (false),
  _solver_configuration(nullptr),
  _close_matrix_before_solve(true)
{
}



template <typename T>
EigenSolver<T>::~EigenSolver ()
{
  this->clear ();
}



template <typename T>
std::unique_ptr<EigenSolver<T>>
EigenSolver<T>::build(const Parallel::Communicator & comm,
                      const SolverPackage solver_package)
{
  // Build the appropriate solver
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_SLEPC
    case SLEPC_SOLVERS:
      return libmesh_make_unique<SlepcEigenSolver<T>>(comm);
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized eigen solver package: " << solver_package);
    }

  return std::unique_ptr<EigenSolver<T>>();
}


template <typename T>
void EigenSolver<T>::set_solver_configuration(SolverConfiguration & solver_configuration)
{
  _solver_configuration = &solver_configuration;
}

template <typename T>
void EigenSolver<T>::set_position_of_spectrum (Real pos)
{
  if (pos >= 0)
    _position_of_spectrum = TARGET_MAGNITUDE;
  else
    _position_of_spectrum = TARGET_REAL;

  _target_val = pos;
}

template <typename T>
void EigenSolver<T>::set_position_of_spectrum (Real pos, PositionOfSpectrum target)
{
  _position_of_spectrum = target;
  _target_val = pos;
}



//------------------------------------------------------------------
// Explicit instantiations
template class EigenSolver<Number>;

} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC
