// $Id: fe_lagrange_shape_1D.C,v 1.10 2004-01-03 15:37:42 benkirk Exp $

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




template <>
Real FE<1,LAGRANGE>::shape(const ElemType,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  const Real xi = p(0);

	
  switch (order)
    {
      // Lagrange linears
    case FIRST:
      {
	assert (i<2);
	
	switch (i)
	  {
	  case 0:
	    return .5*(1. - xi);
	  
	  case 1:
	    return .5*(1. + xi);
	   
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }

      
	      
      // Lagrange quadratics
    case SECOND:
      {
	assert (i<3);
	
	switch (i)
	  {
	  case 0:
	    return .5*xi*(xi - 1.);
	    
	  case 1:
	    return .5*xi*(xi + 1);
	   
	  case 2:
	    return (1. - xi*xi);
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }

      

      // Lagrange cubics
    case THIRD:
      {
	error();
	assert (i<4);
	
	switch (i)
	  {
	  case 0:
	    return 9./16.*(1./9.-xi*xi)*(xi-1.);
	   
	  case 1:
	    return 27./16.*(1.-xi*xi)*(1./3.-xi);
	   
	  case 2:
	    return 27./16.*(1.-xi*xi)*(1./3.+xi);
	   
	  case 3:
	    return -9./16.*(1./9.-xi*xi)*(xi+1.);
	   
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
Real FE<1,LAGRANGE>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,LAGRANGE>::shape(elem->type(), order, i, p);
}



template <>
Real FE<1,LAGRANGE>::shape_deriv(const ElemType,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  // only d()/dxi in 1D!
  
  assert (j == 0);
	
  const Real xi = p(0);

	
  switch (order)
    {      
      // Lagrange linear shape function derivatives
    case FIRST:
      {
	assert (i<2);
	
	switch (i)
	  {
	  case 0:
	    return -.5;
	   
	  case 1:
	    return .5;
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }
	
	      
      // Lagrange quadratic shape function derivatives
    case SECOND:
      {
	assert (i<3);
	
	switch (i)
	  {
	  case 0:
	    return xi-.5;
	    
	  case 1:
	    return xi+.5;
	    
	  case 2:
	    return -2.*xi;
	   
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }


      // Lagrange cubic shape function derivatives
    case THIRD:
      {
	error();
	assert (i<4);
	
	switch (i)
	  {
	  case 0:
	    return -9./16.*(3.*xi*xi-2.*xi-1./9.);
	    
	  case 1:
	    return 27./16.*(3.*xi*xi-2./3.*xi-1.);
	    
	  case 2:
	    return 27./16.*(-3.*xi*xi-2./3.*xi+1.);
	    
	  case 3:
	    return -9./16.*(-3.*xi*xi-2.*xi+1./9.);
	   
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
Real FE<1,LAGRANGE>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,LAGRANGE>::shape_deriv(elem->type(),
				     order, i, j, p);
}




template <>
Real FE<1,LAGRANGE>::shape_second_deriv(const ElemType,
					const Order order,
					const unsigned int i,
					const unsigned int j,
					const Point&)
{
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!

  assert (j == 0);
  
  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
	// All second derivatives of linears are zero....
	return 0.;
      }

      // quadratic Lagrange shape functions
    case SECOND:
      {
	switch (i)
	  {
	  case 0:
	    return 1.;
	    
	  case 1:
	    return 1.;
	    
	  case 2:
	    return -2.;
	    
	  default:
	    {
	      std::cerr << "Invalid shape function index requested!"
			<< std::endl;
	      error();
	    }
	  }
      } // end case SECOND

    default:
      {
	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
	error();
      }
    } // end switch (order)
  
  error();
  return 0.;
}
