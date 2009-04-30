// $Id: fe_lagrange_shape_0D.C 3243 2009-01-30 23:04:41Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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


// C++ includes

// Local includes
#include "fe.h"
#include "elem.h"

template <>
Real FE<1,SCALAR>::shape(const ElemType,
			   const Order,
			   const unsigned int,
			   const Point&)
{
  libmesh_error();
  return 1.;
}

template <>
Real FE<1,SCALAR>::shape(const Elem*,
			   const Order,
			   const unsigned int,
			   const Point&)
{
  libmesh_error();
  return 1.;
}

template <>
Real FE<1,SCALAR>::shape_deriv(const ElemType,
				 const Order,
				 const unsigned int,
				 const unsigned int,
				 const Point&)
{
  libmesh_error();
  return 0.;
}

template <>
Real FE<1,SCALAR>::shape_deriv(const Elem*,
				 const Order,
				 const unsigned int,
				 const unsigned int,
				 const Point&)
{
  libmesh_error();
  return 0.;
}

template <>
Real FE<1,SCALAR>::shape_second_deriv(const ElemType,
					const Order,
					const unsigned int,
					const unsigned int,
					const Point&)
{
  libmesh_error();
  return 0.;
}

template <>
Real FE<1,SCALAR>::shape_second_deriv(const Elem*,
				        const Order,
				        const unsigned int,
				        const unsigned int,
				        const Point&)
{
  libmesh_error();
  return 0.;
}