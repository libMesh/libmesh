// $Id$

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


// C++ inlcludes

// Local includes
#include "fe.h"
#include "elem.h"




template <>
Real FE<0,LAGRANGE>::shape(const ElemType,
			   const Order,
			   const unsigned int i,
			   const Point&)
{
  assert (i < 1);
  return 1.;
}



template <>
Real FE<0,LAGRANGE>::shape(const Elem*,
			   const Order,
			   const unsigned int i,
			   const Point&)
{
  assert (i < 1);
  return 1.;
}



template <>
Real FE<0,LAGRANGE>::shape_deriv(const ElemType,
				 const Order,
				 const unsigned int,
				 const unsigned int,
				 const Point&)
{
  // No spatial derivatives in 0D!
  error();
  return 0.;
}



template <>
Real FE<0,LAGRANGE>::shape_deriv(const Elem*,
				 const Order,
				 const unsigned int,
				 const unsigned int,
				 const Point&)
{
  // No spatial derivatives in 0D!
  error();
  return 0.;
}




template <>
Real FE<0,LAGRANGE>::shape_second_deriv(const ElemType,
					const Order,
					const unsigned int,
					const unsigned int,
					const Point&)
{
  // No spatial derivatives in 0D!
  error();
  return 0.;
}



template <>
Real FE<0,LAGRANGE>::shape_second_deriv(const Elem*,
				        const Order,
				        const unsigned int,
				        const unsigned int,
				        const Point&)
{
  // No spatial derivatives in 0D!
  error();
  return 0.;
}
