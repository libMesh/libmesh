// $Id: enum_fe_family.h,v 1.4 2003-01-24 17:24:38 jwpeterson Exp $

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



#ifndef __enum_fe_family_h__
#define __enum_fe_family_h__

// C++ includes

// Local includes
#include "mesh_config.h"



// ------------------------------------------------------------
// enum FEFamily definition
namespace MeshEnums {
  
  /**
   * \enum MeshEnums::FEFamily defines an \p enum for finite element families.
   */
  enum FEFamily {LAGRANGE=0,
		 HIERARCHIC,
		 MONOMIAL,
		 
#ifdef ENABLE_INFINITE_ELEMENTS
		 INFINITE_MAP,     //   for 1/r-map
		 JACOBI_20_00,     //   i_max = 19
                 JACOBI_30_00,     //   i_max = 19
		 LEGENDRE,         //   i_max = 19
		 INF_LAGRANGE,     //   i_max = 14
#endif
		 
		 INVALID_FE};
};

using namespace MeshEnums;




#endif // #ifndef __fe_type_h__




