// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
#if LIBMESH_DIM > 1

  libmesh_assert (elem != NULL);

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

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert (i < (totalorder+1)*(totalorder+2)/2);

  
  // monomials. since they are hierarchic we only need one case block.
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
      return dx*dx*dy;

    case 8:
      return dx*dy*dy;

    case 9:
      return dy*dy*dy;

      // quartics
    case 10:
      return dx*dx*dx*dx;

    case 11:
      return dx*dx*dx*dy;

    case 12:
      return dx*dx*dy*dy;

    case 13:
      return dx*dy*dy*dy;

    case 14:
      return dy*dy*dy*dy;
    
    default:
      unsigned int o = 0;
      for (; i >= (o+1)*(o+2)/2; o++) { }
      unsigned int i2 = i - (o*(o+1)/2);
      Real val = 1.;
      for (unsigned int index=i2; index != o; index++)
        val *= dx;
      for (unsigned int index=0; index != i2; index++)
        val *= dy;
      return val;
    }

  libmesh_error();
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
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_deriv(const Elem* elem,
			    const Order order,
			    const unsigned int i,
			    const unsigned int j,
			    const Point& p)
{
#if LIBMESH_DIM > 1

  
  libmesh_assert (j<2);
  libmesh_assert (elem != NULL);
  
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

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert (i < (totalorder+1)*(totalorder+2)/2);
  
  // monomials. since they are hierarchic we only need one case block.

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
	  return 2.*dx*dy;
	    
	case 8:
	  return dy*dy;
	    
	case 9:
	  return 0.;
	    
	  // quartics
	case 10:
	  return 4.*dx*dx*dx;
	    
	case 11:
	  return 3.*dx*dx*dy;
	    
	case 12:
	  return 2.*dx*dy*dy;
	    
	case 13:
	  return dy*dy*dy;
	    
	case 14:
	  return 0.;
	    
        default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)/2; o++) { }
          unsigned int i2 = i - (o*(o+1)/2);
          Real val = o - i2;
          for (unsigned int index=i2+1; index < o; index++)
            val *= dx;
          for (unsigned int index=0; index != i2; index++)
            val *= dy;
          return val;
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
	  return dx*dx;
	    
	case 8:
	  return 2.*dx*dy;
	    
	case 9:
	  return 3.*dy*dy;
	    
	  // quartics
	case 10:
	  return 0.;
	    
	case 11:
	  return dx*dx*dx;
	    
	case 12:
	  return 2.*dx*dx*dy;
	    
	case 13:
	  return 3.*dx*dy*dy;
	    
	case 14:
	  return 4.*dy*dy*dy;
	    
        default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)/2; o++) { }
          unsigned int i2 = i - (o*(o+1)/2);
          Real val = i2;
          for (unsigned int index=i2; index != o; index++)
            val *= dx;
          for (unsigned int index=1; index <= i2; index++)
            val *= dy;
          return val;
	}
    }
	      
	      
  default:
    libmesh_error();
  }

  libmesh_error();
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
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<2,XYZ>::shape_second_deriv(const Elem* elem,
			           const Order order,
			           const unsigned int i,
			           const unsigned int j,
			           const Point& p)
{
#if LIBMESH_DIM > 1

  libmesh_assert (j<=2);
  libmesh_assert (elem != NULL);
  
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

#ifndef NDEBUG
  // totalorder is only used in the assertion below, so
  // we avoid declaring it when asserts are not active.
  const unsigned int totalorder = order + elem->p_level();
#endif
  libmesh_assert (i < (totalorder+1)*(totalorder+2)/2);
  
  // monomials. since they are hierarchic we only need one case block.

  switch (j)
    {
      // d^2()/dx^2
    case 0:
      {
        switch (i)
	  {
	    // constants
	  case 0:
	    // linears
	  case 1:
	  case 2:
	    return 0.;

	    // quadratics
	  case 3:
	    return 2.;
	    
	  case 4:
	  case 5:
	    return 0.;

	    // cubics
	  case 6:
	    return 6.*dx;
	    
	  case 7:
	    return 2.*dy;
	    
	  case 8:
	  case 9:
	    return 0.;
	    
	    // quartics
	  case 10:
	    return 12.*dx*dx;
	    
	  case 11:
	    return 6.*dx*dy;
	    
	  case 12:
	    return 2.*dy*dy;
	    
	  case 13:
	  case 14:
	    return 0.;
	    
          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = (o - i2) * (o - i2 - 1);
            for (unsigned int index=i2+2; index < o; index++)
              val *= dx;
            for (unsigned int index=0; index != i2; index++)
              val *= dy;
            return val;
	  }
      }

      // d^2()/dxdy
    case 1:
      {
        switch (i)
	  {
	    // constants
	  case 0:
	    
	    // linears
	  case 1:
	  case 2:
	    return 0.;

	    // quadratics
	  case 3:
	    return 0.;
	    
	  case 4:
	    return 1.;
	    
	  case 5:
	    return 0.;

	    // cubics
	  case 6:
	    return 0.;
	  case 7:
	    return 2.*dx;
	    
	  case 8:
	    return 2.*dy;
	    
	  case 9:
	    return 0.;
	    
	    // quartics
	  case 10:
	    return 0.;

	  case 11:
	    return 3.*dx*dx;
	    
	  case 12:
	    return 4.*dx*dy;
	    
	  case 13:
	    return 3.*dy*dy;
	    
	  case 14:
	    return 0.;
	    
          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = (o - i2) * i2;
            for (unsigned int index=i2+1; index < o; index++)
              val *= dx;
            for (unsigned int index=1; index < i2; index++)
              val *= dy;
            return val;
	  }
      }
	      
      // d^2()/dy^2
    case 2:
      {
        switch (i)
	  {
	    // constants
	  case 0:
	    
	    // linears
	  case 1:
	  case 2:
	    return 0.;

	    // quadratics
	  case 3:
	  case 4:
	    return 0.;
	    
	  case 5:
	    return 2.;

	    // cubics
	  case 6:
	    return 0.;
	    
	  case 7:
	    return 0.;
	    
	  case 8:
	    return 2.*dx;
	    
	  case 9:
	    return 6.*dy;
	    
	    // quartics
	  case 10:
	  case 11:
	    return 0.;
	    
	  case 12:
	    return 2.*dx*dx;
	    
	  case 13:
	    return 6.*dx*dy;
	    
	  case 14:
	    return 12.*dy*dy;
	    
          default:
            unsigned int o = 0;
            for (; i >= (o+1)*(o+2)/2; o++) { }
            unsigned int i2 = i - (o*(o+1)/2);
            Real val = i2 * (i2 - 1);
            for (unsigned int index=i2; index != o; index++)
              val *= dx;
            for (unsigned int index=2; index < i2; index++)
              val *= dy;
            return val;
	  }
      }
    }

  libmesh_error();
  return 0.;

#endif
}
