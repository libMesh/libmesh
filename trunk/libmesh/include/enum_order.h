// $Id: enum_order.h,v 1.2 2003-01-20 17:06:08 jwpeterson Exp $

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



#ifndef __enum_order_h__
#define __enum_order_h__

// Local includes
#include "mesh_config.h"



// ------------------------------------------------------------
// enum Order definition
namespace MeshEnums {

  /**
   * \enum MeshEnums::Order defines an \p enum for polynomial orders.
   */
  enum Order {CONST=0,
	      FIRST,
	      SECOND,
	      THIRD,
	      FOURTH,
	      FIFTH,
	      SIXTH,
	      SEVENTH,
	      EIGHTH,
	      NINTH,
	      TENTH,
	      ELEVENTH,
	      TWELFTH,
	      THIRTEENTH,
	      FOURTEENTH,
	      FIFTEENTH,
	      SIXTEENTH,
	      SEVENTEENTH,
	      EIGHTTEENTH,
	      NINTEENTH,
	      TWENTIETH,
	      TWENTYFIRST,
	      TWENTYSECOND,
	      TWENTYTHIRD,
	      
	      INVALID_ORDER};


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * \enum MeshEnums::ShiftPattern for infinite elements
   */
  enum IfemShiftPattern {DD_SHIFT, ERROR};

#endif

};

using namespace MeshEnums;

#endif
