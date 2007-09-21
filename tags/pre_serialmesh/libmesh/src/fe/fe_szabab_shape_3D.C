// $Id: fe_szabab_shape_3D.C,v 1.4 2005-06-12 18:36:40 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "libmesh_config.h"

#ifdef ENABLE_HIGHER_ORDER_SHAPES

#include "fe.h"


template <>
Real FE<3,SZABAB>::shape(const ElemType,
			 const Order,
			 const unsigned int,
			 const Point&)
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}



template <>
Real FE<3,SZABAB>::shape(const Elem*,
			 const Order,
			 const unsigned int,
			 const Point&)
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}


template <>
Real FE<3,SZABAB>::shape_deriv(const ElemType,
			       const Order,
			       const unsigned int,
			       const unsigned int,
			       const Point& )
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}



template <>
Real FE<3,SZABAB>::shape_deriv(const Elem*,
			       const Order,
			       const unsigned int,
			       const unsigned int,
			       const Point&)
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}



template <>
Real FE<3,SZABAB>::shape_second_deriv(const ElemType,
			              const Order,
			              const unsigned int,
			              const unsigned int,
			              const Point&)
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}



template <>
Real FE<3,SZABAB>::shape_second_deriv(const Elem*,
			              const Order,
			              const unsigned int,
			              const unsigned int,
			              const Point&)
{
  std::cerr << "Szabo-Babuska polynomials are not defined in 3D\n" << std::endl;
  error();
  return 0.;
}

#endif //ENABLE_HIGHER_ORDER_SHAPES

