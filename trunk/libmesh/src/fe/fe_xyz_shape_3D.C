// $Id: fe_xyz_shape_3D.C,v 1.4 2005-01-13 22:10:16 roystgnr Exp $

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
Real FE<3,XYZ>::shape(const ElemType,
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
Real FE<3,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
#if DIM == 3
  
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
  const Real z  = p(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real dx = x - xc;
  const Real dy = y - yc;
  const Real dz = z - zc;
    
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 35);

	switch (i)
	  {
	    // constant
	  case 0:
	    return 1.;

	    // linears
	  case 1:
	    return dx;
	    
	  case 2:
	    return dy;
	    
	  case 3:
	    return dz;

	    // quadratics
	  case 4:
	    return dx*dx;
	    
	  case 5:
	    return dy*dy;
	    
	  case 6:
	    return dz*dz;

	  case 7:
	    return dx*dy;

	  case 8:
	    return dx*dz;

	  case 9:
	    return dz*dy;

	    // cubics
	  case 10:
	    return dx*dx*dx;

	  case 11:
	    return dy*dy*dy;

	  case 12:
	    return dz*dz*dz;

	  case 13:
	    return dx*dy*dz;

	  case 14:
	    return dy*dy*dz;

	  case 15:
	    return dy*dz*dz;

	  case 16:
	    return dx*dz*dz;

	  case 17:
	    return dx*dx*dz;

	  case 18:
	    return dx*dx*dy;

	  case 19:
	    return dx*dy*dy;

	    // quartics
	  case 20:
	    return dx*dx*dx*dx;

	  case 21:
	    return dy*dy*dy*dy;

	  case 22:
	    return dz*dz*dz*dz;
	  case 23:
	    return dx*dx*dx*dy;

	  case 24:
	    return dx*dx*dx*dz;

	  case 25:
	    return dx*dx*dy*dy;

	  case 26:
	    return dx*dx*dy*dz;

	  case 27:
	    return dx*dx*dz*dz;

	  case 28:
	    return dx*dy*dy*dy;

	  case 29:
	    return dx*dy*dy*dz;

	  case 30:
	    return dx*dy*dz*dz;

	  case 31:
	    return dx*dz*dz*dz;

	  case 32:
	    return dy*dy*dy*dz;

	  case 33:
	    return dy*dy*dz*dz;

	  case 34:
	    return dy*dz*dz*dz;
	    	    
	  default:
	    std::cerr << "Invalid shape function index!" << std::endl;
	    error();
	  }
      }


            
      // unsupported order
    default:
      {
	std::cerr << "ERROR: Unsupported 3D FE order!: " << order
		  << std::endl;
	error();
      }
    }

#endif
  
  error();
  return 0.;
}



template <>
Real FE<3,XYZ>::shape_deriv(const ElemType,
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
Real FE<3,XYZ>::shape_deriv(const Elem* elem,
			    const Order order,
			    const unsigned int i,
			    const unsigned int j,
			    const Point& p)
{
#if DIM == 3

  assert (elem != NULL);
  assert (j<3);
  
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
  const Real z  = p(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real dx = x - xc;
  const Real dy = y - yc;
  const Real dz = z - zc;
  
  switch (order)
    {
      // monomials. since they are heirarchic we only need one case block.
    case CONSTANT:
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
      {
	assert (i < 35);

	switch (j)
	  {
	    // d()/dx
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
		  return 2.*dx;
		  
		case 5:
		  return 0.;
		  
		case 6:
		  return 0.;
		  
		case 7:
		  return dy;
		  
		case 8:
		  return dz;
		  
		case 9:
		  return 0.;

		  // cubic
		case 10:
		  return 3.*dx*dx;

		case 11:
		  return 0.;

		case 12:
		  return 0.;

		case 13:
		  return dy*dz;

		case 14:
		  return 0.;

		case 15:
		  return 0.;

		case 16:
		  return dz*dz;

		case 17:
		  return 2.*dx*dz;

		case 18:
		  return 2.*dx*dy;

		case 19:
		  return dy*dy;

		  // quartics
		case 20:
		  return 4.*dx*dx*dx;

		case 21:
		  return 0.;

		case 22:
		  return 0.;

		case 23:
		  return 3.*dx*dx*dy;

		case 24:
		  return 3.*dx*dx*dz;

		case 25:
		  return 2.*dx*dy*dy;

		case 26:
		  return 2.*dx*dy*dz;

		case 27:
		  return 2.*dx*dz*dz;

		case 28:
		  return dy*dy*dy;

		case 29:
		  return dy*dy*dz;

		case 30:
		  return dy*dz*dz;

		case 31:
		  return dz*dz*dz;

		case 32:
		  return 0.;

		case 33:
		  return 0.;

		case 34:
		  return 0.;
		  
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
		  return 2.*dy;
		  
		case 6:
		  return 0.;
		  
		case 7:
		  return dx;
		  
		case 8:
		  return 0.;
		  
		case 9:
		  return dz;

		  // cubic
		case 10:
		  return 0.;

		case 11:
		  return 3.*dy*dy;

		case 12:
		  return 0.;

		case 13:
		  return dx*dz;

		case 14:
		  return 2.*dy*dz;

		case 15:
		  return dz*dz;

		case 16:
		  return 0.;

		case 17:
		  return 0.;

		case 18:
		  return dx*dx;

		case 19:
		  return 2.*dx*dy;

		  // quartics
		case 20:
		  return 0.;

		case 21:
		  return 4.*dy*dy*dy;

		case 22:
		  return 0.;

		case 23:
		  return dx*dx*dx;

		case 24:
		  return 0.;

		case 25:
		  return 2.*dx*dx*dy;

		case 26:
		  return dx*dx*dz;

		case 27:
		  return 0.;

		case 28:
		  return 3.*dx*dy*dy;

		case 29:
		  return 2.*dx*dy*dz;

		case 30:
		  return dx*dz*dz;

		case 31:
		  return 0.;

		case 32:
		  return 3.*dy*dy*dz;

		case 33:
		  return 2.*dy*dz*dz;

		case 34:
		  return dz*dz*dz;
		  
		default:
		  std::cerr << "Invalid shape function index!" << std::endl;
		  error();
		}
	    }

	    
	    // d()/dz
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
		  return 2.*dz;
		  
		case 7:
		  return 0.;
		  
		case 8:
		  return dx;
		  
		case 9:
		  return dy;

		  // cubic
		case 10:
		  return 0.;

		case 11:
		  return 0.;

		case 12:
		  return 3.*dz*dz;

		case 13:
		  return dx*dy;

		case 14:
		  return dy*dy;

		case 15:
		  return 2.*dy*dz;

		case 16:
		  return 2.*dx*dz;

		case 17:
		  return dx*dx;

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
		  return 4.*dz*dz*dz;

		case 23:
		  return 0.;

		case 24:
		  return dx*dx*dx;

		case 25:
		  return 0.;

		case 26:
		  return dx*dx*dy;

		case 27:
		  return 2.*dx*dx*dz;

		case 28:
		  return 0.;

		case 29:
		  return dx*dy*dy;

		case 30:
		  return 2.*dx*dy*dz;

		case 31:
		  return 3.*dx*dz*dz;

		case 32:
		  return dy*dy*dy;

		case 33:
		  return 2.*dy*dy*dz;

		case 34:
		  return 3.*dy*dz*dz;
		  
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
	std::cerr << "ERROR: Unsupported 3D FE order!: " << order
		  << std::endl;
	error();
      }
    }

#endif
  
  error();
  return 0.;  
}



template <>
Real FE<3,XYZ>::shape_second_deriv(const ElemType,
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
Real FE<3,XYZ>::shape_second_deriv(const Elem*,
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
