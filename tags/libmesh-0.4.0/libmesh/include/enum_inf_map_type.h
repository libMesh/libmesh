// $Id: enum_inf_map_type.h,v 1.3 2003-05-15 23:34:33 benkirk Exp $

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



#ifndef __enum_inf_map_type_h__
#define __enum_inf_map_type_h__

// Local includes
#include "mesh_config.h"


#ifdef ENABLE_INFINITE_ELEMENTS

// ------------------------------------------------------------
// enum Order definition
namespace libMeshEnums {

  /**
   * \enum libMeshEnums::InfMapType defines an \p enum for the
   * types of coordinate mappings available in infinite elements.
   */
  enum InfMapType {CARTESIAN=0,
		   SPHERICAL,
		   ELLIPSOIDAL,
		   INVALID_INF_MAP};

}

using namespace libMeshEnums;

#endif //ifdef ENABLE_INFINITE_ELEMENTS

#endif
