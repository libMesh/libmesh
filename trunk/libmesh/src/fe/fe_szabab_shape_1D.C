// $Id: fe_szabab_shape_1D.C,v 1.1 2004-01-09 19:33:22 spetersen Exp $

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



// C++ includes
#include <math.h>


// Local includes
#include "libmesh_config.h"

#ifdef ENABLE_HIGHER_ORDER_SHAPES

#include "fe.h"
#include "elem.h"


template <>
Real FE<1,SZABAB>::shape(const ElemType,
			 const Order order,
			 const unsigned int i,
			 const Point& p)
{
  const Real xi = p(0);
  	
  switch (order)
    {
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
    case SEVENTH:
      
      switch(i)
	{				
	  //nodal shape functions
	case 0: return (1.-xi)/2.;
	case 1:	return (1.+xi)/2.;
	case 2:	return 1./4.*sqrt(6.)*(pow(xi,2)-1.);
	case 3:	return 1./4.*sqrt(10.)*(pow(xi,3)-xi);
	case 4: return 1./16.*sqrt(14.)*(5.*pow(xi,4)-6.*pow(xi,2)+1.);
	case 5:	return 3./16.*sqrt(2.)*(7.*pow(xi,5)-10.*pow(xi,3)+3.*xi);
	case 6:	return 1./32.*sqrt(22.)*(21.*pow(xi,6)-35.*pow(xi,4)+15.*pow(xi,2)-1.);
	case 7: return 1./32.*sqrt(26.)*(33.*pow(xi,7)-63.*pow(xi,5)+35.*pow(xi,3)-5.*xi);
	case 8: return 1./256.*sqrt(30.)*(429.*pow(xi,8)-924.*pow(xi,6)+630.*pow(xi,4)-140.*pow(xi,2)+5.);
	  
	default:
	  std::cerr << "Invalid shape function index!" << std::endl;
	  error();	    
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
Real FE<1,SZABAB>::shape(const Elem* elem,
			 const Order order,
			 const unsigned int i,
			 const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,SZABAB>::shape(elem->type(), order, i, p);
}



template <>
Real FE<1,SZABAB>::shape_deriv(const ElemType,
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
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
    case SEVENTH:
      
      switch(i)
	{
	case 0:	return -0.5;
	case 1:	return  0.5;
	case 2:	return 1./ 2.*sqrt( 6.)*xi;
	case 3:	return 3./ 4.*sqrt(10.)*pow(xi,2)-1./4.*sqrt(10.);
	case 4:	return 1./ 4.*sqrt(14.)*(5.  *pow(xi,3)-3.*xi);
	case 5:	return 1./16.*sqrt( 2.)*(105.*pow(xi,4)-90.*pow(xi,2)+9.); 
	case 6:	return 1./16.*sqrt(22.)*(63. *pow(xi,5)-70.*pow(xi,3)+15.*xi);
	case 7: return 1./32.*sqrt(26.)*(231.*pow(xi,6)-315.*pow(xi,4)+105.*pow(xi,2)-5.);
	case 8: return 1./32.*sqrt(30.)*(429.*pow(xi,7)-693.*pow(xi,5)+315.*pow(xi,3)-35.*xi);
	  
	default:
	  std::cerr << "Invalid shape function index!" << std::endl;
	  error();
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
Real FE<1,SZABAB>::shape_deriv(const Elem* elem,
			       const Order order,
			       const unsigned int i,
			       const unsigned int j,
			       const Point& p)
{
  assert (elem != NULL);
  
  return FE<1,SZABAB>::shape_deriv(elem->type(),
				       order, i, j, p);
}


#endif //ENABLE_HIGHER_ORDER_SHAPES
