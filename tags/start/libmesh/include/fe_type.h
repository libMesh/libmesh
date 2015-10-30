// $Id: fe_type.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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



#ifndef __fe_type_h__
#define __fe_type_h__

// C++ includes

// Local includes
#include "mesh_config.h"
#include "order.h"



// ------------------------------------------------------------
// enum ElemType definition
namespace MeshEnums {
  
  /**
   * \enum MeshEnums::FEFamily defines an \p enum for finite element families.
   */
  enum FEFamily {LAGRANGE=0,
		 HIERARCHIC,
		 MONOMIAL,
		 
#ifdef ENABLE_INFINITE_ELEMENTS
		 
                 JACOBI_15_n05,    //   i_max = 17
                 JACOBI_15_00,     //   i_max = 19
                 JACOBI_15_05,     //   i_max = 17
                 JACOBI_20_n05,    //   i_max = 19
                 JACOBI_20_00,     //   i_max = 19
                 JACOBI_20_05,     //   i_max = 19
                 JACOBI_25_n05,    //   i_max = 17
                 JACOBI_25_00,     //   i_max = 17
                 JACOBI_25_05,     //   i_max = 17
                 JACOBI_30_n05,    //   i_max = 19
                 JACOBI_30_00,     //   i_max = 19
                 JACOBI_30_05,     //   i_max = 19
		 LEGENDRE,         //   i_max = 19
	         LAGRANGE_CONST,
	         LAGRANGE_FIRST,
	         LAGRANGE_SECOND,
	         LAGRANGE_THIRD,
	         LAGRANGE_FOURTH,
	         LAGRANGE_FIFTH,
	         LAGRANGE_SIXTH,
	         LAGRANGE_SEVENTH,
	         LAGRANGE_EIGHTH,
	         LAGRANGE_NINTH,
	         LAGRANGE_TENTH,
	         LAGRANGE_ELEVENTH,
	         LAGRANGE_TWELFTH,
	         LAGRANGE_THIRTEENTH,
	         LAGRANGE_FOURTEENTH,
		 
#endif
		 INVALID_FE};
};

using namespace MeshEnums;



/**
 * class FEType hides (possibly multiple) FEFamily and approximation
 * orders, thereby enabling specialized finite element families.
 */
class FEType
{
public:

#ifndef ENABLE_INFINITE_ELEMENTS
  
  /**
   * Constructor.  Optionally takes the approximation \p Order
   * and the finite element family \p FEFamily
   */
  FEType(const Order    o = FIRST,
	 const FEFamily f = LAGRANGE) :
    order(o),
    family(f)
  {};
    
#else
  
  /**
   * Constructor.  Optionally takes the approximation \p Order
   * and the finite element family \p FEFamily
   */
  FEType(const Order    o  = FIRST,
	 const FEFamily f  = LAGRANGE,
	 const FEFamily bf = INVALID_FE) :
    order(o),
    family(f),
    base_family(bf)
  {};
        
#endif


  /**
   * The approximation order of the element.  
   */
  Order order;
  
  /**
   * The type of finite element.  Valid types are \p LAGRANGE,
   * HIERARCHIC, etc...
   */
  FEFamily family;

#ifdef ENABLE_INFINITE_ELEMENTS
  /**
   * For InfFE, \p family contains the radial shape family, while 
   * \p base_family contains the base shape.
   */
  FEFamily base_family;
  
#endif

private:  
  
};


#endif // #ifndef __fe_type_h__




