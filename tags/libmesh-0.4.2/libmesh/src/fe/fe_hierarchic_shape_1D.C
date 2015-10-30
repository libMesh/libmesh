// $Id: fe_hierarchic_shape_1D.C,v 1.12 2004-02-10 13:28:07 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "utility.h"




template <>
Real FE<1,HIERARCHIC>::shape(const ElemType,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
  const Real xi = p(0);
  
  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

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
	    return (pow<4>(xi) - 1.)/24.;
	    
	  case 6:
	    return (pow<6>(xi) - 1.)/720.;

	    // All odd-terms have the same form.
	    // (xi^p - xi)/p!
	  case 3:
	    return (xi*xi*xi - xi)/6.;

	  case 5:
	    return (pow<5>(xi) - xi)/120.;

	  case 7:
	    return (pow<7>(xi) - xi)/5040.;	    
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();	    
	  }
      }
      
    default:
      {
	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
	error();
      }
    }

  error();
  return 0.;
}



template <>
Real FE<1,HIERARCHIC>::shape(const Elem* elem,
			     const Order order,
			     const unsigned int i,
			     const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,HIERARCHIC>::shape(elem->type(), order, i, p);
}



template <>
Real FE<1,HIERARCHIC>::shape_deriv(const ElemType,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
  // only d()/dxi in 1D!
  
  assert (j == 0);

  // Declare that we are using our own special power function
  // from the Utility namespace.  This saves typing later.
  using Utility::pow;

  const Real xi = p(0);

	
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
	    return pow<3>(xi)/6.;
	    
	  case 6:
	    return pow<5>(xi)/120.;

	    // All odd-terms have the same form.
	    // (p*xi^(p-1) - 1.)/p!
	  case 3:
	    return (3*xi*xi - 1.)/6.;

	  case 5:
	    return (5.*pow<4>(xi) - 1.)/120.;

	  case 7:
	    return (7.*pow<6>(xi) - 1.)/5040.;	    
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();	    
	  }
      }


      
    default:
      {
	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
	error();
      }
    }

  error();
  return 0.;
}



template <>
Real FE<1,HIERARCHIC>::shape_deriv(const Elem* elem,
				   const Order order,
				   const unsigned int i,
				   const unsigned int j,
				   const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,HIERARCHIC>::shape_deriv(elem->type(),
				       order, i, j, p);
}
