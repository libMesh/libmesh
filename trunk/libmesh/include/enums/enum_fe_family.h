// $Id: enum_fe_family.h,v 1.2 2004-01-03 15:37:41 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __enum_fe_family_h__
#define __enum_fe_family_h__

// C++ includes

// Local includes
#include "libmesh_config.h"



// ------------------------------------------------------------
// enum FEFamily definition
namespace libMeshEnums {
  
  /**
   * \enum libMeshEnums::FEFamily defines an \p enum for finite element families.
   * Explicity assign numbers so that one can differentiate between enabled
   * and disabled infinite elements even when these are disabled.
   */
  enum FEFamily {LAGRANGE     = 0,
		 HIERARCHIC   = 1,
		 MONOMIAL     = 2,
		 
#ifdef ENABLE_INFINITE_ELEMENTS
		 INFINITE_MAP = 11,     //   for 1/r-map
		 JACOBI_20_00 = 12,     //   i_max = 19
                 JACOBI_30_00 = 13,     //   i_max = 19
		 LEGENDRE     = 14,     //   i_max = 19
#endif
		 
		 INVALID_FE   = 42};
}

using namespace libMeshEnums;




#endif // #ifndef __fe_type_h__




