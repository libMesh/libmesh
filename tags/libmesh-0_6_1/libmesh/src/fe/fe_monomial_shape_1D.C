// $Id: fe_monomial_shape_1D.C,v 1.15 2007-10-21 20:48:46 benkirk Exp $

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
Real FE<1,MONOMIAL>::shape(const ElemType,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  const Real xi = p(0);

	
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 5);

	switch (i)
	  {
	  case 0:
	    return 1.;

	  case 1:
	    return xi;
	    
	  case 2:
	    return xi*xi;
	    
	  case 3:
	    return xi*xi*xi;
	    
	  case 4:
	    return xi*xi*xi*xi;
	    
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
Real FE<1,MONOMIAL>::shape(const Elem* elem,
			   const Order order,
			   const unsigned int i,
			   const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,MONOMIAL>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
}



template <>
Real FE<1,MONOMIAL>::shape_deriv(const ElemType,
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
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 5);

	switch (i)
	  {
	  case 0:
	    return 0.;

	  case 1:
	    return 1.;
	    
	  case 2:
	    return 2.*xi;
	    
	  case 3:
	    return 3.*xi*xi;
	    
	  case 4:
	    return 4.*xi*xi*xi;
	    
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
Real FE<1,MONOMIAL>::shape_deriv(const Elem* elem,
				 const Order order,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,MONOMIAL>::shape_deriv(elem->type(),
				     static_cast<Order>(order + elem->p_level()), i, j, p);
}



template <>
Real FE<1,MONOMIAL>::shape_second_deriv(const ElemType,
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
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 5);

	switch (i)
	  {
	  case 0:
	  case 1:
	    return 0.;
	    
	  case 2:
	    return 2.;
	    
	  case 3:
	    return 6.*xi;
	    
	  case 4:
	    return 12.*xi*xi;
	    
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
Real FE<1,MONOMIAL>::shape_second_deriv(const Elem* elem,
				        const Order order,
				        const unsigned int i,
				        const unsigned int j,
				        const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,MONOMIAL>::shape_second_deriv(elem->type(),
				            static_cast<Order>(order + elem->p_level()), i, j, p);
}
