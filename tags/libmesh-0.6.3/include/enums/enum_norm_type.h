// $Id: enum_norm_type.h 2501 2007-11-20 02:33:29Z benkirk $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_norm_type_h__
#define __enum_norm_type_h__

// ------------------------------------------------------------
// enum NormType definition
namespace libMeshEnums {
  
  /**
   * \enum libMeshEnums::FEMNormType defines an \p enum for norms
   * defined on vectors of finite element coefficients
   */
		 // Hilbert norms and seminorms in FE space
  enum FEMNormType {L2             = 0,
		    H1             = 1,
		    H2             = 2,

		    H1_SEMINORM    = 10,
		    H2_SEMINORM    = 11,

                    // discrete vector norms
                    DISCRETE_L1    = 20,
		    DISCRETE_L2    = 21,
		    DISCRETE_L_INF = 22,

		    INVALID_NORM   = 42};
}

using namespace libMeshEnums;

#endif // #ifndef __norm_type_h__




