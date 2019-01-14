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



#ifndef LIBMESH_ENUM_ELEM_QUALITY_H
#define LIBMESH_ENUM_ELEM_QUALITY_H

namespace libMesh
{

/**
 * Defines an \p enum for element quality metrics.
 *
 * The fixed type, i.e. ": int", enumeration syntax used here allows
 * this enum to be forward declared as
 * enum ElemQuality : int;
 * reducing header file dependencies.
 */
enum ElemQuality : int {
                  ASPECT_RATIO=0,
                  SKEW,
                  SHEAR,
                  SHAPE,
                  MAX_ANGLE,
                  MIN_ANGLE,
                  CONDITION,
                  DISTORTION,
                  TAPER,
                  WARP,
                  STRETCH,
                  DIAGONAL,
                  ASPECT_RATIO_BETA,
                  ASPECT_RATIO_GAMMA,
                  SIZE,
                  JACOBIAN};
}

#endif
