// $Id: fe_hermite_shape_1D.C,v 1.4 2005-09-02 18:48:55 roystgnr Exp $

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
// This should also screw up multithreading royally
namespace
{
  static unsigned int old_elem_id = libMesh::invalid_uint;
  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static Real d1xd1x, d2xd2x;

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
    FE<1,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  // Degrees of freedom are at vertices and edge midpoints
  std::vector<Point> dofpt;
  dofpt.push_back(Point(0));
  dofpt.push_back(Point(1));

  // Mapping functions - first derivatives at each dofpt
  std::vector<Real> dxdxi(2);
  std::vector<Real> dxidx(2);

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[p] = 0;
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = FE<1,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 0, dofpt[p]);
          dxdxi[p] += dofpt[p](0) * ddxi;
        }
    }

  // Calculate derivative scaling factors
  
  d1xd1x = dxdxi[0];
  d2xd2x = dxdxi[1];
}

} // end anonymous namespace



template<>
Real FEHermite<1>::hermite_raw_shape_second_deriv
 (const unsigned int basis_num, const Real xi)
{
  switch (basis_num)
    {
      case 0:
        return 1.5 * xi;
      case 1:
        return -1.5 * xi;
      case 2:
        return 0.5 * (-1. + 3.*xi);
      case 3:
        return 0.5 * (1. + 3.*xi);
    }

  error();
  return 0.;
}



template<>
Real FEHermite<1>::hermite_raw_shape_deriv
 (const unsigned int basis_num, const Real xi)
{
  switch (basis_num)
    {
      case 0:
        return 0.75 * (-1. + xi*xi);
      case 1:
        return 0.75 * (1. - xi*xi);
      case 2:
        return 0.25 * (-1. - 2.*xi + 3.*xi*xi);
      case 3:
        return 0.25 * (-1. + 2.*xi + 3.*xi*xi);
    }

  error();
  return 0.;
}

template<>
Real FEHermite<1>::hermite_raw_shape
 (const unsigned int basis_num, const Real xi)
{
  switch (basis_num)
    {
      case 0:
        return 0.25 * (2. - 3.*xi + xi*xi*xi);
      case 1:
        return 0.25 * (2. + 3.*xi - xi*xi*xi);
      case 2:
        return 0.25 * (1. - xi - xi*xi + xi*xi*xi);
      case 3:
        return 0.25 * (-1. - xi + xi*xi + xi*xi*xi);
    }

  error();
  return 0.;
}


template<>
Real FEHermite<1>::hermite_bases
 (std::vector<unsigned int> &bases1D,
  const std::vector<std::vector<Real> > &dxdxi,
  unsigned int index,
  unsigned int dim)
{
  unsigned int pow2[4] = {1, 2, 4, 8};
  assert (dim < 4);

  bases1D.clear();
  bases1D.resize(dim,0);

  switch (index/pow2[dim]) // Node number
    {
    case 0:
      break;
    case 1:
      bases1D[0] = 1; break;
    case 2:
      bases1D[0] = 1; bases1D[1] = 1; break;
    case 3:
      bases1D[1] = 1; break;
    case 4:
      bases1D[2] = 1; break;
    case 5:
      bases1D[0] = 1; bases1D[2] = 1; break;
    case 6:
      bases1D[0] = 1; bases1D[1] = 1; bases1D[2] = 1; break;
    case 7:
      bases1D[1] = 1; bases1D[2] = 1; break;
    default:
      error();
    }

  Real coef;
  switch (index%pow2[dim]) // DoF type
    {
    case 0: // DoF = value at node
      coef = 1.0;
      break;
    case 1: // DoF = x derivative at node
      coef = dxdxi[0][bases1D[0]];
      bases1D[0] += 2; break;
    case 2: // DoF = y derivative at node
      coef = dxdxi[1][bases1D[1]];
      bases1D[1] += 2; break;
    case 3: // DoF = xy derivative at node
      coef = dxdxi[0][bases1D[0]] * dxdxi[1][bases1D[1]];
      bases1D[0] += 2; bases1D[1] += 2; break;
    case 4: // DoF = z derivative at node
      coef = dxdxi[2][bases1D[2]];
      bases1D[2] += 2; break;
    case 5: // DoF = xz derivative at node
      coef = dxdxi[0][bases1D[0]] * dxdxi[2][bases1D[2]];
      bases1D[0] += 2; bases1D[2] += 2; break;
    case 6: // DoF = yz derivative at node
      coef = dxdxi[1][bases1D[1]] * dxdxi[2][bases1D[2]];
      bases1D[1] += 2; bases1D[2] += 2; break;
    case 7: // DoF = xyz derivative at node
      coef = dxdxi[0][bases1D[0]] * dxdxi[1][bases1D[1]] * dxdxi[2][bases1D[2]];
      bases1D[0] += 2; bases1D[1] += 2; bases1D[2] += 2; break;
    }

  // No singular elements
  assert(coef);
  return coef;
}

  

template <>
Real FE<1,HERMITE>::shape(const ElemType,
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
Real FE<1,HERMITE>::shape(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // Hermite cubic shape functions
    case THIRD:
      {
	switch (type)
	  {
	    // C1 functions on the C1 cubic edge
	  case EDGE2:
	  case EDGE3:
	    {
	      assert (i<4);

	      switch (i)
		{
		case 0:
		  return FEHermite<1>::hermite_raw_shape(0, p(0));
		case 1:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape(2, p(0));
		case 2:
		  return FEHermite<1>::hermite_raw_shape(1, p(0));
		case 3:
                  return d2xd2x * FEHermite<1>::hermite_raw_shape(3, p(0));
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
Real FE<1,HERMITE>::shape_deriv(const ElemType,
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
Real FE<1,HERMITE>::shape_deriv(const Elem* elem,
				const Order order,
				const unsigned int i,
				const unsigned int,
				const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // Hermite cubic shape functions
    case THIRD:
      {
	switch (type)
	  {
	    // C1 functions on the C1 cubic edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (i)
		{
		case 0:
		  return FEHermite<1>::hermite_raw_shape_deriv(0, p(0));
		case 1:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape_deriv(2, p(0));
		case 2:
		  return FEHermite<1>::hermite_raw_shape_deriv(1, p(0));
		case 3:
                  return d2xd2x * FEHermite<1>::hermite_raw_shape_deriv(3, p(0));
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
Real FE<1,HERMITE>::shape_second_deriv(const Elem* elem,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int,
                                       const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // Hermite cubic shape functions
    case THIRD:
      {
	switch (type)
	  {
	    // C1 functions on the C1 cubic edge
	  case EDGE2:
	  case EDGE3:
	    {
	      switch (i)
		{
		case 0:
		  return FEHermite<1>::hermite_raw_shape_second_deriv(0, p(0));
		case 1:
		  return d1xd1x * FEHermite<1>::hermite_raw_shape_second_deriv(2, p(0));
		case 2:
		  return FEHermite<1>::hermite_raw_shape_second_deriv(1, p(0));
		case 3:
                  return d2xd2x * FEHermite<1>::hermite_raw_shape_second_deriv(3, p(0));
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
