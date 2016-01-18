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



#ifndef LIBMESH_ENUM_EIGENSOLVER_TYPE_H
#define LIBMESH_ENUM_EIGENSOLVER_TYPE_H

// ------------------------------------------------------------
// enum SolverType definition
namespace libMesh {

/**
 * Defines an \p enum for iterative eigenproblem solver types
 */
enum EigenSolverType {POWER=0,
                      LAPACK,
                      SUBSPACE,
                      ARNOLDI,
                      LANCZOS,
                      KRYLOVSCHUR,
                      // SLEPc optional packages
                      // EPSARPACK,
                      // EPSLAPACK,
                      // EPSBLZPACK,
                      // EPSPLANSO,
                      // EPSTRLAN,

                      INVALID_EIGENSOLVER};

/**
 * Defines an \p enum for  eigenproblem  types.
 * This can be Hermitian (HEP), generalized Hermitian (GHEP),
 * non-Hermitian (NHEP), generalized non-Hermitian (GNHEP), or
 * generalized indefinite Hermitian (GHIEP).
 */
enum EigenProblemType {NHEP=0,
                       HEP,
                       GNHEP,
                       GHEP,
                       GHIEP,

                       INVALID_EIGENPROBLEMTYPE};



/**
 * Defines an \p enum for the position of
 * the spectrum, i.e. the eigenvalues to be computed.
 */
enum PositionOfSpectrum {LARGEST_MAGNITUDE=0,
                         SMALLEST_MAGNITUDE,
                         LARGEST_REAL,
                         SMALLEST_REAL,
                         LARGEST_IMAGINARY,
                         SMALLEST_IMAGINARY,

                         INVALID_Postion_of_Spectrum};
}

#endif // LIBMESH_ENUM_EIGENSOLVER_TYPE_H
