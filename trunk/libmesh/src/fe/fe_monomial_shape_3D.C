// $Id: fe_monomial_shape_3D.C,v 1.5 2003-01-24 17:24:42 jwpeterson Exp $

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


// C++ inlcludes
#include <math.h>


// Local includes
#include "fe.h"
#include "elem.h"




template <>
real FE<3,MONOMIAL>::shape(const ElemType,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
#if DIM == 3
    
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONST:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 35);

	const real xi   = p(0);
	const real eta  = p(1);
	const real zeta = p(2);

	switch (i)
	  {
	    // constant
	  case 0:
	    return 1.;

	    // linears
	  case 1:
	    return xi;
	    
	  case 2:
	    return eta;
	    
	  case 3:
	    return zeta;

	    // quadratics
	  case 4:
	    return xi*xi;
	    
	  case 5:
	    return eta*eta;
	    
	  case 6:
	    return zeta*zeta;

	  case 7:
	    return xi*eta;

	  case 8:
	    return xi*zeta;

	  case 9:
	    return zeta*eta;

	    // cubics
	  case 10:
	    return xi*xi*xi;

	  case 11:
	    return eta*eta*eta;

	  case 12:
	    return zeta*zeta*zeta;

	  case 13:
	    return xi*eta*zeta;

	  case 14:
	    return eta*eta*zeta;

	  case 15:
	    return eta*zeta*zeta;

	  case 16:
	    return xi*zeta*zeta;

	  case 17:
	    return xi*xi*zeta;

	  case 18:
	    return xi*xi*eta;

	  case 19:
	    return xi*eta*eta;

	    // quartics
	  case 20:
	    return xi*xi*xi*xi;

	  case 21:
	    return eta*eta*eta*eta;

	  case 22:
	    return zeta*zeta*zeta*zeta;
	  case 23:
	    return xi*xi*xi*eta;

	  case 24:
	    return xi*xi*xi*zeta;

	  case 25:
	    return xi*xi*eta*eta;

	  case 26:
	    return xi*xi*eta*zeta;

	  case 27:
	    return xi*xi*zeta*zeta;

	  case 28:
	    return xi*eta*eta*eta;

	  case 29:
	    return xi*eta*eta*zeta;

	  case 30:
	    return xi*eta*zeta*zeta;

	  case 31:
	    return xi*zeta*zeta*zeta;

	  case 32:
	    return eta*eta*eta*zeta;

	  case 33:
	    return eta*eta*zeta*zeta;

	  case 34:
	    return eta*zeta*zeta*zeta;
	    	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  };
      };


            
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 3D FE order!: " << order
		  << std::endl;
	error();
      }
    };

#endif
  
  error();
  return 0.;
};



template <>
real FE<3,MONOMIAL>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  assert (elem != NULL);
      
  // call the orientation-independent shape functions
  return FE<3,MONOMIAL>::shape(elem->type(), order, i, p);
};



template <>
real FE<3,MONOMIAL>::shape_deriv(const ElemType,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
#if DIM == 3
  
  assert (j<3);
  
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONST:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 35);

	const real xi   = p(0);
	const real eta  = p(1);
	const real zeta = p(2);

	switch (j)
	  {
	    // d()/dxi
	  case 0:
	    {
	      switch (i)
		{
		  // constant
		case 0:
		  return 0.;
		  
		  // linear
		case 1:
		  return 1.;
		  
		case 2:
		  return 0.;
		  
		case 3:
		  return 0.;

		  // quadratic
		case 4:
		  return 2.*xi;
		  
		case 5:
		  return 0.;
		  
		case 6:
		  return 0.;
		  
		case 7:
		  return eta;
		  
		case 8:
		  return zeta;
		  
		case 9:
		  return 0.;

		  // cubic
		case 10:
		  return 3.*xi*xi;

		case 11:
		  return 0.;

		case 12:
		  return 0.;

		case 13:
		  return eta*zeta;

		case 14:
		  return 0.;

		case 15:
		  return 0.;

		case 16:
		  return zeta*zeta;

		case 17:
		  return 2.*xi*zeta;

		case 18:
		  return 2.*xi*eta;

		case 19:
		  return eta*eta;

		  // quartics
		case 20:
		  return 4.*xi*xi*xi;

		case 21:
		  return 0.;

		case 22:
		  return 0.;

		case 23:
		  return 3.*xi*xi*eta;

		case 24:
		  return 3.*xi*xi*zeta;

		case 25:
		  return 2.*xi*eta*eta;

		case 26:
		  return 2.*xi*eta*zeta;

		case 27:
		  return 2.*xi*zeta*zeta;

		case 28:
		  return eta*eta*eta;

		case 29:
		  return eta*eta*zeta;

		case 30:
		  return eta*zeta*zeta;

		case 31:
		  return zeta*zeta*zeta;

		case 32:
		  return 0.;

		case 33:
		  return 0.;

		case 34:
		  return 0.;
		  
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		};
	    };

	    
	    // d()/deta
	  case 1:
	    {
	      switch (i)
		{
		  // constant
		case 0:
		  return 0.;
		  
		  // linear
		case 1:
		  return 0.;
		  
		case 2:
		  return 1.;
		  
		case 3:
		  return 0.;

		  // quadratic
		case 4:
		  return 0.;
		  
		case 5:
		  return 2.*eta;
		  
		case 6:
		  return 0.;
		  
		case 7:
		  return xi;
		  
		case 8:
		  return 0.;
		  
		case 9:
		  return zeta;

		  // cubic
		case 10:
		  return 0.;

		case 11:
		  return 3.*eta*eta;

		case 12:
		  return 0.;

		case 13:
		  return xi*zeta;

		case 14:
		  return 2.*eta*zeta;

		case 15:
		  return zeta*zeta;

		case 16:
		  return 0.;

		case 17:
		  return 0.;

		case 18:
		  return xi*xi;

		case 19:
		  return 2.*xi*eta;

		  // quartics
		case 20:
		  return 0.;

		case 21:
		  return 4.*eta*eta*eta;

		case 22:
		  return 0.;

		case 23:
		  return xi*xi*xi;

		case 24:
		  return 0.;

		case 25:
		  return 2.*xi*xi*eta;

		case 26:
		  return xi*xi*zeta;

		case 27:
		  return 0.;

		case 28:
		  return 3.*xi*eta*eta;

		case 29:
		  return 2.*xi*eta*zeta;

		case 30:
		  return xi*zeta*zeta;

		case 31:
		  return 0.;

		case 32:
		  return 3.*eta*eta*zeta;

		case 33:
		  return 2.*eta*zeta*zeta;

		case 34:
		  return zeta*zeta*zeta;
		  
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		};
	    };

	    
	    // d()/dzeta
	  case 2:
	    {
	      switch (i)
		{
		  // constant
		case 0:
		  return 0.;
		  
		  // linear
		case 1:
		  return 0.;
		  
		case 2:
		  return 0.;
		  
		case 3:
		  return 1.;

		  // quadratic
		case 4:
		  return 0.;
		  
		case 5:
		  return 0.;
		  
		case 6:
		  return 2.*zeta;
		  
		case 7:
		  return 0.;
		  
		case 8:
		  return xi;
		  
		case 9:
		  return eta;

		  // cubic
		case 10:
		  return 0.;

		case 11:
		  return 0.;

		case 12:
		  return 3.*zeta*zeta;

		case 13:
		  return xi*eta;

		case 14:
		  return eta*eta;

		case 15:
		  return 2.*eta*zeta;

		case 16:
		  return 2.*xi*zeta;

		case 17:
		  return xi*xi;

		case 18:
		  return 0.;

		case 19:
		  return 0.;

		  // quartics
		case 20:
		  return 0.;

		case 21:
		  return 0.;

		case 22:
		  return 4.*zeta*zeta*zeta;

		case 23:
		  return 0.;

		case 24:
		  return xi*xi*xi;

		case 25:
		  return 0.;

		case 26:
		  return xi*xi*eta;

		case 27:
		  return 2.*xi*xi*zeta;

		case 28:
		  return 0.;

		case 29:
		  return xi*eta*eta;

		case 30:
		  return 2.*xi*eta*zeta;

		case 31:
		  return 3.*xi*zeta*zeta;

		case 32:
		  return eta*eta*eta;

		case 33:
		  return 2.*eta*eta*zeta;

		case 34:
		  return 3.*eta*zeta*zeta;
		  
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		};
	    };

	    
	  default:
	    error();
	  };
      };


            
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 3D FE order!: " << order
		  << std::endl;
	error();
      }
    };

#endif
  
  error();
  return 0.;  
};



template <>
real FE<3,MONOMIAL>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  assert (elem != NULL);
      
  // call the orientation-independent shape function derivatives
  return FE<3,MONOMIAL>::shape_deriv(elem->type(), order, i, j, p);
};
