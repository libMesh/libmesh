// $Id: fe_type.h,v 1.4 2003-01-21 19:24:34 benkirk Exp $

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
#include "enum_order.h"
#include "enum_fe_family.h"



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

  /**
   * The approximation order of the element.  
   */
  Order order;

  /**
   * The type of finite element.  Valid types are \p LAGRANGE,
   * \p HIERARCHIC, etc...
   */
  FEFamily family;
    
#else
  
  /**
   * Constructor.  Provides a uniform interface for non-infinite
   * elements.
   */
  FEType(const Order    o = FIRST,
	 const FEFamily f = LAGRANGE) :
    order(o),
    base_order(o),
    family(f),
    base_family(f)
  {};

  /**
   * Constructor.  Optionally takes the approximation \p Order
   * and the finite element family \p FEFamily
   */
  FEType(const Order    o,
	 const Order    bo,
	 const FEFamily f,
	 const FEFamily bf) :
    order(o),
    base_order(bo),
    family(f),
    base_family(bf)
  {};

  /**
   * The approximation order in radial direction of the infinite element.  
   */
  Order order;

  /**
   * The approximation order in the base of the infinite element.
   */
  Order base_order;

  /**
   * The type of radial approximation in radial direction.  Valid types are 
   * \p JACOBI_20_00, \p JACOBI_30_00, \p LEGENDRE, and \p LAGRANGE.
   */
  FEFamily family;

  /**
   * For InfFE, \p family contains the radial shape family, while
   * \p base_family contains the approximation type in circumferential
   * direction.  Valid types are \p LAGRANGE, \p HIERARCHIC, etc...
   */
  FEFamily base_family;
  
#endif

private:  
  
};


#endif // #ifndef __fe_type_h__




