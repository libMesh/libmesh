// $Id: fe_xyz_shape_1D.C,v 1.4 2005-01-13 22:10:16 roystgnr Exp $

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
Real FE<1,XYZ>::shape(const ElemType,
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
Real FE<1,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
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
  const Real xc = centroid(0);
  const Real dx = x - xc;

	
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
	    return dx;
	    
	  case 2:
	    return dx*dx;
	    
	  case 3:
	    return dx*dx*dx;
	    
	  case 4:
	    return dx*dx*dx*dx;
	    
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
Real FE<1,XYZ>::shape_deriv(const ElemType,
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
Real FE<1,XYZ>::shape_deriv(const Elem* elem,
			    const Order order,
			    const unsigned int i,
			    const unsigned int j,
			    const Point& p)
{
  assert (elem != NULL);
  
  // only d()/dxi in 1D!
  
  assert (j == 0);
	
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
  const Real xc = centroid(0);
  const Real dx = x - xc;

	
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
	    return 2.*dx;
	    
	  case 3:
	    return 3.*dx*dx;
	    
	  case 4:
	    return 4.*dx*dx*dx;
	    
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
Real FE<1,XYZ>::shape_second_deriv(const ElemType,
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
Real FE<1,XYZ>::shape_second_deriv(const Elem*,
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
