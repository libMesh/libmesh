// $Id: enum_order.h,v 1.6 2003-05-15 23:34:33 benkirk Exp $

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


// ------------------------------------------------------------
// enum Order definition
namespace libMeshEnums {

  /**
   * \enum libMeshEnums::Order defines an \p enum for polynomial orders.
   * Fixing each label to a specific int, since \p InfFE uses 
   * \p static_cast<unsigned int>.
   */
  enum Order {CONST        =  0,
	      FIRST        =  1,
	      SECOND       =  2,
	      THIRD        =  3,
	      FOURTH       =  4,
	      FIFTH        =  5,
	      SIXTH        =  6,
	      SEVENTH      =  7,
	      EIGHTH       =  8,
	      NINTH        =  9,
	      TENTH        = 10,
	      ELEVENTH     = 11,
	      TWELFTH      = 12,
	      THIRTEENTH   = 13,
	      FOURTEENTH   = 14,
	      FIFTEENTH    = 15,
	      SIXTEENTH    = 16,
	      SEVENTEENTH  = 17,
	      EIGHTTEENTH  = 18,
	      NINTEENTH    = 19,
	      TWENTIETH    = 20,
	      TWENTYFIRST  = 21,
	      TWENTYSECOND = 22,
	      TWENTYTHIRD  = 23,	      
	      INVALID_ORDER};

}

using namespace libMeshEnums;

#endif
