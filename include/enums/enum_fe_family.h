// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ENUM_FE_FAMILY_H
#define LIBMESH_ENUM_FE_FAMILY_H

// ------------------------------------------------------------
// enum FEFamily definition
namespace libMesh {

/**
 * \enum libMesh::FEFamily defines an \p enum for finite
 * element families.
 */
// vanilla C0
enum FEFamily {LAGRANGE     = 0,
               HIERARCHIC   = 1,

               // discontinuous, in local coordinates
               MONOMIAL      = 2,
               L2_HIERARCHIC = 6,
               L2_LAGRANGE   = 7,

               // higher-order
               BERNSTEIN    = 3,
               SZABAB       = 4,

               // discontinuous, in global coordinates
               XYZ          = 5,

               // infinite element stuff
               INFINITE_MAP = 11,     //   for 1/r-map
               JACOBI_20_00 = 12,     //   i_max = 19
               JACOBI_30_00 = 13,     //   i_max = 19
               LEGENDRE     = 14,     //   i_max = 19

               // C1 elements
               CLOUGH       = 21,
               HERMITE      = 22,
               SUBDIV       = 23,

               // A scalar variable that couples to
               // all other DOFs in the system
               SCALAR       = 31,

               // Vector-valued elements
               LAGRANGE_VEC = 41,
               NEDELEC_ONE  = 42,

               INVALID_FE   = 99};

/**
 * \enum libMesh::FEContinuity defines an \p enum for finite element
 * types to libmesh_assert a certain level (or type? Hcurl?) of continuity.
 */
enum FEContinuity {DISCONTINUOUS,
                   C_ZERO,
                   C_ONE,
                   H_CURL};

/**
 * \enum libMesh::FEFieldType defines an \p enum for finite element
 * field types - i.e. is it a scalar element, vector, tensor, etc.
 */
enum FEFieldType {TYPE_SCALAR = 0,
                  TYPE_VECTOR};

}

#endif // LIBMESH_ENUM_FE_FAMILY_H
