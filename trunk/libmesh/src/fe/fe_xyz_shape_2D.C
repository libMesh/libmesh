// $Id: fe_xyz_shape_2D.C,v 1.4 2005-01-13 22:10:16 roystgnr Exp $

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

// Local includes
#include "fe.h"
#include "elem.h"




// Anonymous namespace for persistant variables.
// This allows us to determine when the centroid needs
// to be recalculated.
namespace
{
  static unsigned int old_elem_id = libMesh::invalid_uint;
  static Point centroid;
}



template <>
Real FE<2,XYZ>::shape(const ElemType,
		      const Order,
		      const unsigned int,
		      const Point&)
{
  std::cerr << "XYZ polynomials require the element\n"
            << "because the centroid is needed."
            << std::endl;

  error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
#if DIM > 1

  assert (elem != NULL);

  // Only recompute the centroid if the element
  // has changed from the last one we computed.
  // This avoids repeated centroid calculations
  // when called in succession with the same element.
  if (elem->id() != old_elem_id)
    {
      centroid = elem->centroid();
      old_elem_id = elem->id();
    }  
  
  const Real x  = p(0);
  const Real y  = p(1);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real dx = x - xc;
  const Real dy = y - yc;

  
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 15);

	switch (i)
	  {
	    // constant
	  case 0:
	    return 1.;

	    // linear
	  case 1:
	    return dx;
	    
	  case 2:
	    return dy;

	    // quadratics
	  case 3:
	    return dx*dx;
	    
	  case 4:
	    return dx*dy;
	    
	  case 5:
	    return dy*dy;

	    // cubics
	  case 6:
	    return dx*dx*dx;

	  case 7:
	    return dy*dy*dy;

	  case 8:
	    return dx*dx*dy;

	  case 9:
	    return dx*dy*dy;

	    // quartics
	  case 10:
	    return dx*dx*dx*dx;

	  case 11:
	    return dy*dy*dy*dy;

	  case 12:
	    return dx*dx*dy*dy;

	  case 13:
	    return dx*dx*dx*dy;

	  case 14:
	    return dx*dy*dy*dy;
	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }

      
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      }
    }

  error();
  return 0.;

#endif
}



template <>
Real FE<2,XYZ>::shape_deriv(const ElemType,
			    const Order,
			    const unsigned int,
			    const unsigned int,
			    const Point&)
{
  std::cerr << "XYZ polynomials require the element\n"
            << "because the centroid is needed."
            << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_deriv(const Elem* elem,
			    const Order order,
			    const unsigned int i,
			    const unsigned int j,
			    const Point& p)
{
#if DIM > 1

  
  assert (j<2);
  assert (elem != NULL);
  
  // Only recompute the centroid if the element
  // has changed from the last one we computed.
  // This avoids repeated centroid calculations
  // when called in succession with the same element.
  if (elem->id() != old_elem_id)
    {
      centroid = elem->centroid();
      old_elem_id = elem->id();
    }  
  
  const Real x  = p(0);
  const Real y  = p(1);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real dx = x - xc;
  const Real dy = y - yc;

  
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 15);

	switch (j)
	  {
	    // d()/dx
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
		  return 2.*dx;
		    
		case 4:
		  return dy;
		    
		case 5:
		  return 0.;

		  // cubics
		case 6:
		  return 3.*dx*dx;
		    
		case 7:
		  return 0.;
		    
		case 8:
		  return 2.*dx*dy;
		    
		case 9:
		  return dy*dy;
		    
		  // quartics
		case 10:
		  return 4.*dx*dx*dx;
		    
		case 11:
		  return 0.;
		    
		case 12:
		  return 2.*dx*dy*dy;
		    
		case 13:
		  return 3.*dx*dx*dy;
		    
		case 14:
		  return dy*dy*dy;
		    
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		}
	    }

	      
	    // d()/dy
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
		  return dx;
		    
		case 5:
		  return 2.*dy;

		  // cubics
		case 6:
		  return 0.;
		    
		case 7:
		  return 3.*dy*dy;
		    
		case 8:
		  return dx*dx;
		    
		case 9:
		  return 2.*dx*dy;
		    
		  // quartics
		case 10:
		  return 0.;
		    
		case 11:
		  return 4.*dy*dy*dy;
		    
		case 12:
		  return 2.*dx*dx*dy;
		    
		case 13:
		  return dx*dx*dx;
		    
		case 14:
		  return 3.*dx*dy*dy;
		    
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		}
	    }
	      
	      
	  default:
	    error();
	  }
      }

      
      
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 2D FE order!: " << order
		  << std::endl;
	error();
      }
    }

  error();
  return 0.;

#endif
}



template <>
Real FE<2,XYZ>::shape_second_deriv(const ElemType,
			           const Order,
			           const unsigned int,
			           const unsigned int,
			           const Point&)
{
  std::cerr << "XYZ polynomials require the element\n"
            << "because the centroid is needed."
            << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_second_deriv(const Elem*,
			           const Order,
			           const unsigned int,
			           const unsigned int,
			           const Point&)
{
  static bool warning_given = false;

  if (!warning_given)
  std::cerr << "Second derivatives for XYZ elements "
            << " are not yet implemented!"
            << std::endl;

  warning_given = true;
  return 0.;
}
