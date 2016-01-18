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



#ifndef LIBMESH_ENUM_CONVERGENCE_FLAGS
#define LIBMESH_ENUM_CONVERGENCE_FLAGS

// ------------------------------------------------------------
// enum ConvergenceFlags definition
namespace libMesh {

/**
 * Linear solver convergence flags (taken from the PETSc flags)
 */
enum LinearConvergenceReason {
  /* converged */
  CONVERGED_RTOL_NORMAL        =  1,
  CONVERGED_ATOL_NORMAL        =  9,
  CONVERGED_RTOL               =  2,
  CONVERGED_ATOL               =  3,
  CONVERGED_ITS                =  4,
  CONVERGED_CG_NEG_CURVE       =  5,
  CONVERGED_CG_CONSTRAINED     =  6,
  CONVERGED_STEP_LENGTH        =  7,
  CONVERGED_HAPPY_BREAKDOWN    =  8,
  /* diverged */
  DIVERGED_NULL                = -2,
  DIVERGED_ITS                 = -3,
  DIVERGED_DTOL                = -4,
  DIVERGED_BREAKDOWN           = -5,
  DIVERGED_BREAKDOWN_BICG      = -6,
  DIVERGED_NONSYMMETRIC        = -7,
  DIVERGED_INDEFINITE_PC       = -8,
  DIVERGED_NAN                 = -9,
  DIVERGED_INDEFINITE_MAT      = -10,

  CONVERGED_ITERATING          =  0,

  UNKNOWN_FLAG                 = -128};
}

#endif // LIBMESH_ENUM_CONVERGENCE_FLAGS
