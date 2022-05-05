// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ENUM_ERROR_ESTIMATOR_TYPE_H
#define LIBMESH_ENUM_ERROR_ESTIMATOR_TYPE_H

namespace libMesh {

/**
 * Defines an \p enum for the different types of error estimators which are available.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum ErrorEstimatorType : int;
 * reducing header file dependencies.
 */
enum ErrorEstimatorType : int {
                         INVALID                 = -1,
                         ADJOINT_REFINEMENT      =  0,
                         ADJOINT_RESIDUAL        =  1,
                         DISCONTINUITY_MEASURE   =  2,
                         EXACT                   =  3,
                         KELLY                   =  4,
                         LAPLACIAN               =  5,
                         PATCH_RECOVERY          =  6,
                         WEIGHTED_PATCH_RECOVERY =  7,
                         UNIFORM_REFINEMENT      =  8};

}

#endif
