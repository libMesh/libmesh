// $Id: fe_hermite_shape_2D.C,v 1.1 2005-08-25 18:31:37 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
// This allows us to cache the global-to-local mapping transformation
// FIXME: This should also screw up multithreading royally
namespace
{
  static unsigned int old_elem_id = libMesh::invalid_uint;
  // Coefficient naming: d1xd1x is the derivative of x with respect to
  // xi at xi=0, d2xd2x is at xi=1, similarly for y/eta
  static Real d1xd1x, d2xd2x, d1yd1y, d2yd2y;


// Compute the static coefficients for an element
void hermite_compute_coefs(const Elem* elem)
{
  // Coefficients are cached from old elements
  if (elem->id() == old_elem_id)
    return;

  old_elem_id = elem->id();

  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());
  const int n_mapping_shape_functions =
    FE<2,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  std::vector<Point> dofpt;
  dofpt.push_back(Point(1,0));
  dofpt.push_back(Point(0,1));

  // Mapping functions - derivatives at each dofpt
  std::vector<Real> dxdxi(2), dydeta(2);

  for (int p = 0; p != 2; ++p)
    {
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = FE<2,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 0, dofpt[p]);
          const Real ddeta = FE<2,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 1, dofpt[p]);

// dxdeta and dydxi should be 0!
          dxdxi[p] += elem->point(i)(0) * ddxi;
          dydeta[p] += elem->point(i)(1) * ddeta;
        }
    }

  d1xd1x = dxdxi[0];
  d2xd2x = dxdxi[1];
  d1yd1y = dydeta[0];
  d2yd2y = dydeta[1];
}


} // end anonymous namespace



template <>
Real FE<2,HERMITE>::shape(const ElemType,
			  const Order,
			  const unsigned int,
			  const Point&)
{
  std::cerr << "Hermite elements require the real element\n"
	    << "to construct gradient-based degrees of freedom."
	    << std::endl;
  
  error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 3rd-order bicubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      assert (i<16);

	      switch (i)
		{
		case 0:
		  return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                         FEHermite<1>::hermite_raw_shape(0,p(1));
		case 1:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                         FEHermite<1>::hermite_raw_shape(0,p(1));
		case 2:
		  return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                         d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
		case 3:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                         d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
		case 4:
		  return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                         FEHermite<1>::hermite_raw_shape(0,p(1));
		case 5:
		  return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                         FEHermite<1>::hermite_raw_shape(0,p(1));
		case 6:
		  return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                         d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
		case 7:
		  return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                         d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
		case 8:
		  return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                         FEHermite<1>::hermite_raw_shape(1,p(1));
		case 9:
		  return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                         FEHermite<1>::hermite_raw_shape(1,p(1));
		case 10:
		  return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                         d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
		case 11:
		  return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                         d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
		case 12:
		  return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                         FEHermite<1>::hermite_raw_shape(1,p(1));
		case 13:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                         FEHermite<1>::hermite_raw_shape(1,p(1));
		case 14:
		  return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                         d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
		case 15:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                         d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
		  
		default:
		  error();
		}
	    }
	  default:
            std::cerr << "ERROR: Unsupported element type!" << std::endl;
	    error();
	  }
      }
      // by default throw an error
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();
    }
  
  error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_deriv(const ElemType,
				const Order,			    
				const unsigned int,
				const unsigned int,
				const Point&)
{
  std::cerr << "Hermite elements require the real element\n"
	    << "to construct gradient-based degrees of freedom."
	    << std::endl;

  error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_deriv(const Elem* elem,
				const Order order,
				const unsigned int i,
				const unsigned int j,
				const Point& p)
{
  assert (elem != NULL);
  assert (j == 0 || j == 1);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 3rd-order bicubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      assert (i<16);

	      switch (i)
		{
		case 0:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                    }
		case 1:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                    }
		case 2:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                    }
		case 3:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                    }
		case 4:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                    }
		case 5:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                    }
		case 6:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                    }
		case 7:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                    }
		case 8:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                    }
		case 9:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                    }
		case 10:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                    }
		case 11:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                    }
		case 12:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                    }
		case 13:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                    }
		case 14:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                    }
		case 15:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                    }
		default:
		      error();
		}
	    }
	  default:
            std::cerr << "ERROR: Unsupported element type!" << std::endl;
	    error();
	  }
      }
      // by default throw an error
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();
    }
  
  error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_second_deriv(const Elem* elem,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point& p)
{
  assert (elem != NULL);
  assert (j == 0 || j == 1 || j == 2);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 3rd-order bicubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case QUAD4:
	  case QUAD8:
	  case QUAD9:
	    {
	      assert (i<16);

	      switch (i)
		{
		case 0:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(0,p(1));
                    }
		case 1:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                  case 2:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(0,p(1));
                    }
		case 2:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(0,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(1));
                    }
		case 3:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                  case 2:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(1));
                    }
		case 4:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(0,p(1));
                    }
		case 5:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape(0,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(0,p(1));
                  case 2:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(0,p(1));
                    }
		case 6:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(1,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(1));
                    }
		case 7:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape(2,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_deriv(2,p(1));
                  case 2:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             d1yd1y * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(1));
                    }
		case 8:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(1,p(1));
                    }
		case 9:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                  case 2:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(1,p(1));
                    }
		case 10:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(1,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(1,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(1,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(1));
                    }
		case 11:
	          switch (j)
		    {
                  case 0:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                  case 2:
		      return d2xd2x * FEHermite<1>::hermite_raw_shape(3,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(1));
                    }
		case 12:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(1,p(1));
                    }
		case 13:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape(1,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_deriv(1,p(1));
                  case 2:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             FEHermite<1>::hermite_raw_shape_second_deriv(1,p(1));
                    }
		case 14:
	          switch (j)
		    {
                  case 0:
		      return FEHermite<1>::hermite_raw_shape_second_deriv(0,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return FEHermite<1>::hermite_raw_shape_deriv(0,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                  case 2:
		      return FEHermite<1>::hermite_raw_shape(0,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(1));
                    }
		case 15:
	          switch (j)
		    {
                  case 0:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape(3,p(1));
                  case 1:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_deriv(3,p(1));
                  case 2:
		      return d1xd1x * FEHermite<1>::hermite_raw_shape(2,p(0)) *
                             d2yd2y * FEHermite<1>::hermite_raw_shape_second_deriv(3,p(1));
                    }
		default:
		      error();
		}
	    }
	  default:
            std::cerr << "ERROR: Unsupported element type!" << std::endl;
	    error();
	  }
      }
      // by default throw an error
    default:
      std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
      error();
    }
  
  error();
  return 0.;
}
