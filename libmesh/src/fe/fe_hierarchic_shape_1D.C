// $Id: fe_hierarchic_shape_1D.C,v 1.2 2003-01-20 16:31:32 jwpeterson Exp $

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
real FE<1,HIERARCHIC>::shape(const ElemType,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
  const real xi = p(0);

	
  switch (order)
    {
      // Hierarchics. since they are heirarchic we only need one case block.
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
      {
	assert (i < 6);
	
	switch (i)
	  {
	  case 0:
	    return .5*(1. - xi);

	  case 1:
	    return .5*(1.  + xi);

	    // All even-terms have the same form.
	    // (xi^p - 1.)/p!
	  case 2:
	    return (xi*xi - 1.)/2.;
	    
	  case 4:
	    return (pow(xi,4) - 1.)/24.;
	    
	  case 6:
	    return (pow(xi,6) - 1.)/720.;

	    // All odd-terms have the same form.
	    // (xi^p - xi)/p!
	  case 3:
	    return (xi*xi*xi - xi)/6.;

	  case 5:
	    return (pow(xi,5) - xi)/120.;

	  case 7:
	    return (pow(xi,7) - xi)/5040.;	    
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();	    
	  };
      };
      
    default:
      {
	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
	error();
      };
    };

  error();
  return 0.;
};



template <>
real FE<1,HIERARCHIC>::shape(const Elem* elem,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,HIERARCHIC>::shape(elem->type(), order, i, p);
};



template <>
real FE<1,HIERARCHIC>::shape_deriv(const ElemType,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
  // only d()/dxi in 1D!
  
  assert (j == 0);
	
  const real xi = p(0);

	
  switch (order)
    {      
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
      {
	assert (i < 6);
	
	switch (i)
	  {
	  case 0:
	    return -.5;
	    
	  case 1:
	    return  .5;

	    // All even-terms have the same form.
	    // xi^(p-1)/(p-1)!
	  case 2:
	    return xi;
	    
	  case 4:
	    return pow(xi,3)/6.;
	    
	  case 6:
	    return pow(xi,5)/120.;

	    // All odd-terms have the same form.
	    // (p*xi^(p-1) - 1.)/p!
	  case 3:
	    return (3*xi*xi - 1.)/6.;

	  case 5:
	    return (5.*pow(xi,4) - 1.)/120.;

	  case 7:
	    return (7.*pow(xi,6) - 1.)/5040.;	    
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();	    
	  };
      };


      
    default:
      {
	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
	error();
      };
    };

  error();
  return 0.;
};



template <>
real FE<1,HIERARCHIC>::shape_deriv(const Elem* elem,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,HIERARCHIC>::shape_deriv(elem->type(),
				       order, i, j, p);
};
