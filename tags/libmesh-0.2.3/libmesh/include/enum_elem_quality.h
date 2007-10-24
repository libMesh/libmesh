// $Id: enum_elem_quality.h,v 1.4 2003-01-24 17:24:38 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_elem_quality_h__
#define __enum_elem_quality_h__

// C++ includes

// Local includes




// ------------------------------------------------------------
// enum ElemType definition
namespace MeshEnums
{    
  /**
   * Defines an \p enum for element quality metrics.
   */
  enum ElemQuality {ASPECT_RATIO=0,
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
};


using namespace MeshEnums;




#endif
