// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ENUM_NORM_TYPE_H
#define LIBMESH_ENUM_NORM_TYPE_H

// ------------------------------------------------------------
// enum FEMNormType definition
namespace libMeshEnums {

  /**
   * \enum libMeshEnums::FEMNormType defines an \p enum for norms
   * defined on vectors of finite element coefficients
   */
		    // Hilbert norms and seminorms in FE space
  enum FEMNormType {L2              = 0,
		    H1              = 1,
		    H2              = 2,
		    HCURL           = 3,
		    HDIV            = 4,

		    L1              = 5,
                    L_INF           = 6,

		    H1_SEMINORM     = 10,

		    H2_SEMINORM     = 11,

		    HCURL_SEMINORM  = 12,

		    HDIV_SEMINORM   = 13,

                    W1_INF_SEMINORM = 15,
                    W2_INF_SEMINORM = 16,

                    // discrete norms on coefficient vectors
                    DISCRETE_L1     = 20,
		    DISCRETE_L2     = 21,
		    DISCRETE_L_INF  = 22,

                    // Seminorms based on only individual gradient
                    // directional components
		    H1_X_SEMINORM    = 31,
		    H1_Y_SEMINORM    = 32,
		    H1_Z_SEMINORM    = 33,

		    INVALID_NORM    = 42};
}

using namespace libMeshEnums;

#endif // LIBMESH_ENUM_NORM_TYPE_H
