// $Id: fe_szabab_shape_1D.C,v 1.4 2004-02-10 13:28:07 benkirk Exp $

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
  const Real xi  = p(0);
  const Real xi2 = xi*xi;


  // Use this assert rather than a switch with a single entry...
  // It will go away in optimized mode, essentially has the same effect.
  assert (order <= SEVENTH);
  
//   switch (order)
//     {
//     case FIRST:
//     case SECOND:
//     case THIRD:
//     case FOURTH:
//     case FIFTH:
//     case SIXTH:
//     case SEVENTH:
      
      switch(i)
	{				
	  //nodal shape functions
	case 0: return 1./2.-1./2.*xi;
	case 1: return 1./2.+1./2.*xi;
	case 2: return 1./4.  *2.4494897427831780982*(xi2-1.);
	case 3: return 1./4.  *3.1622776601683793320*(xi2-1.)*xi;
	case 4: return 1./16. *3.7416573867739413856*((5.*xi2-6.)*xi2+1.);
	case 5: return 3./16. *1.4142135623730950488*(3.+(-10.+7.*xi2)*xi2)*xi;
	case 6: return 1./32. *4.6904157598234295546*(-1.+(15.+(-35.+21.*xi2)*xi2)*xi2);
	case 7: return 1./32. *5.0990195135927848300*(-5.+(35.+(-63.+33.*xi2)*xi2)*xi2)*xi;
	case 8: return 1./256.*5.4772255750516611346*(5.+(-140.+(630.+(-924.+429.*xi2)*xi2)*xi2)*xi2);
	  
	default:
	  std::cerr << "Invalid shape function index!" << std::endl;
	  error();	    
	}
      
//     default:
//       {
// 	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
// 	error();
//       }
//     }
  
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
	
  const Real xi  = p(0);
  const Real xi2 = xi*xi;
  
  // Use this assert rather than a switch with a single entry...
  // It will go away in optimized mode, essentially has the same effect.
  assert (order <= SEVENTH);
  
//   switch (order)
//     {      
//     case FIRST:
//     case SECOND:
//     case THIRD:
//     case FOURTH:
//     case FIFTH:
//     case SIXTH:
//     case SEVENTH:
      
      switch(i)
	{
	case 0:	return -1./2.;
	case 1:	return 1./2.;
	case 2:	return 1./2.*2.4494897427831780982*xi;
	case 3:	return -1./4.*3.1622776601683793320+3./4.*3.1622776601683793320*xi2;
	case 4:	return 1./16.*3.7416573867739413856*(-12.+20*xi2)*xi;
	case 5:	return 9./16.*1.4142135623730950488+(-45./8.*1.4142135623730950488+105./16.*1.4142135623730950488*xi2)*xi2;
	case 6:	return 1./32.*4.6904157598234295546*(30.+(-140.+126.*xi2)*xi2)*xi;
	case 7:	return -5./32.*5.0990195135927848300+(105./32.*5.0990195135927848300+(-315./32.*5.0990195135927848300+231./32.*5.0990195135927848300*xi2)*xi2)*xi2;
	case 8:	return 1./256.*5.4772255750516611346*(-280.+(2520.+(-5544.+3432.*xi2)*xi2)*xi2)*xi;
	  
	default:
	  std::cerr << "Invalid shape function index!" << std::endl;
	  error();
	}
 
//     default:
//       {
// 	std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
// 	error();
//       }
//     }

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
