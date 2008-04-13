// $Id$

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
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
  assert (elem != NULL);
  assert (i <= order + elem->p_level());

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

  // monomials. since they are hierarchic we only need one case block.
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
      Real val = 1.;
      for (unsigned int index = 0; index != i; ++index)
        val *= dx;
      return val;
    }
      
  libmesh_error();
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
  
  libmesh_error();
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
  assert (i <= order + elem->p_level());
  
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

  // monomials. since they are hierarchic we only need one case block.
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
      Real val = i;
      for (unsigned int index = 1; index != i; ++index)
        val *= dx;
      return val;
    }

  libmesh_error();
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
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<1,XYZ>::shape_second_deriv(const Elem* elem,
			           const Order order,
			           const unsigned int i,
			           const unsigned int j,
			           const Point& p)
{
  assert (elem != NULL);
  assert (i <= order + elem->p_level());
  
  // only d2()/dxi2 in 1D!
  
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

  // monomials. since they are hierarchic we only need one case block.
  switch (i)
    {
    case 0:
    case 1:
      return 0.;
    
    case 2:
      return 2.;
    
    case 3:
      return 6.*dx;
    
    case 4:
      return 12.*dx*dx;
    
    default:
      Real val = 2.;
      for (unsigned int index = 2; index != i; ++index)
        val *= (index+1) * dx;
      return val;
    }

  libmesh_error();
  return 0.;
}
