// $Id: fe_monomial_shape_2D.C,v 1.3 2003-01-20 17:06:28 jwpeterson Exp $

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
real FE<2,MONOMIAL>::shape(const ElemType,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
#if DIM > 1
  
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONST:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 15);

	const real xi  = p(0);
	const real eta = p(1);

	switch (i)
	  {
	    // constant
	  case 0:
	    return 1.;

	    // linear
	  case 1:
	    return xi;
	    
	  case 2:
	    return eta;

	    // quadratics
	  case 3:
	    return xi*xi;
	    
	  case 4:
	    return xi*eta;
	    
	  case 5:
	    return eta*eta;

	    // cubics
	  case 6:
	    return xi*xi*xi;

	  case 7:
	    return eta*eta*eta;

	  case 8:
	    return xi*xi*eta;

	  case 9:
	    return xi*eta*eta;

	    // quartics
	  case 10:
	    return xi*xi*xi*xi;

	  case 11:
	    return eta*eta*eta*eta;

	  case 12:
	    return xi*xi*eta*eta;

	  case 13:
	    return xi*xi*xi*eta;

	  case 14:
	    return xi*eta*eta*eta;
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  };
      };

      
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      };
    };

  error();
  return 0.;

#endif
};



template <>
real FE<2,MONOMIAL>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  assert (elem != NULL);
  
  // by default call the orientation-independent shape functions
  return FE<2,MONOMIAL>::shape(elem->type(), order, i, p);
};



template <>
real FE<2,MONOMIAL>::shape_deriv(const ElemType,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
#if DIM > 1

  
  assert (j<2);

  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONST:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 15);

	const real xi  = p(0);
	const real eta = p(1);

	switch (j)
	  {
	    // d()/dxi
	  case 0:
	    {
	      switch (i)
		{
		  // constants
		case 0:
		  return 0.;
		    
		  // linears
		case 1:
		  return 1.;
		    
		case 2:
		  return 0.;

		  // quadratics
		case 3:
		  return 2.*xi;
		    
		case 4:
		  return eta;
		    
		case 5:
		  return 0.;

		  // cubics
		case 6:
		  return 3.*xi*xi;
		    
		case 7:
		  return 0.;
		    
		case 8:
		  return 2.*xi*eta;
		    
		case 9:
		  return eta*eta;
		    
		  // quartics
		case 10:
		  return 4.*xi*xi*xi;
		    
		case 11:
		  return 0.;
		    
		case 12:
		  return 2.*xi*eta*eta;
		    
		case 13:
		  return 3.*xi*xi*eta;
		    
		case 14:
		  return eta*eta*eta;
		    
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
		  // constants
		case 0:
		  return 0.;
		    
		  // linears
		case 1:
		  return 0.;
		    
		case 2:
		  return 1.;

		  // quadratics
		case 3:
		  return 0.;
		    
		case 4:
		  return xi;
		    
		case 5:
		  return 2.*eta;

		  // cubics
		case 6:
		  return 0.;
		    
		case 7:
		  return 3.*eta*eta;
		    
		case 8:
		  return xi*xi;
		    
		case 9:
		  return 2.*xi*eta;
		    
		  // quartics
		case 10:
		  return 0.;
		    
		case 11:
		  return 4.*eta*eta*eta;
		    
		case 12:
		  return 2.*xi*xi*eta;
		    
		case 13:
		  return xi*xi*xi;
		    
		case 14:
		  return 3.*xi*eta*eta;
		    
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
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      }
    };

  error();
  return 0.;

#endif
};



template <>
real FE<2,MONOMIAL>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  assert (elem != NULL);

  // by default call the orientation-independent shape functions
  return FE<2,MONOMIAL>::shape_deriv(elem->type(), order, i, j, p); 
};

