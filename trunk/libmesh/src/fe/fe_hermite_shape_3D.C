// $Id: fe_hermite_shape_3D.C,v 1.2 2005-08-26 06:26:29 roystgnr Exp $

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
  static std::vector<Real> dxdxi(2), dydeta(2), dzdzeta(2);


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
    FE<3,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  std::vector<Point> dofpt;
  dofpt.push_back(Point(1,0,0));
  dofpt.push_back(Point(0,1,0));
  dofpt.push_back(Point(0,0,1));

  // Mapping functions - derivatives at each dofpt
  std::vector<Real> dxdxi(2), dydeta(2), dzdzeta(2);

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[0] = 0;
      dydeta[0] = 0;
      dzdzeta[0] = 0;
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = FE<3,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 0, dofpt[p]);
          const Real ddeta = FE<3,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 1, dofpt[p]);
          const Real ddzeta = FE<3,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 2, dofpt[p]);

	  // dxdeta, dxdzeta, dydxi, dydzeta, dzdxi, dzdeta should all
          // be 0!
          dxdxi[p] += elem->point(i)(0) * ddxi;
          dydeta[p] += elem->point(i)(1) * ddeta;
          dzdzeta[p] += elem->point(i)(2) * ddzeta;
        }
    }
}


} // end anonymous namespace



template <>
Real FE<3,HERMITE>::shape(const ElemType,
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
Real FE<3,HERMITE>::shape(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case HEX8:
	  case HEX20:
	  case HEX27:
	    {
	      assert (i<64);

              // Is node at maximum x/y/z plane?
              unsigned int xmax, ymax, zmax, xbasis, ybasis, zbasis;
              switch (i/8) // Node number
                {
                case 0:
                  xmax = 0; ymax = 0; zmax = 0; break;
                case 1:
                  xmax = 1; ymax = 0; zmax = 0; break;
                case 2:
                  xmax = 1; ymax = 1; zmax = 0; break;
                case 3:
                  xmax = 0; ymax = 1; zmax = 0; break;
                case 4:
                  xmax = 0; ymax = 0; zmax = 1; break;
                case 5:
                  xmax = 1; ymax = 0; zmax = 1; break;
                case 6:
                  xmax = 1; ymax = 1; zmax = 1; break;
                case 7:
                  xmax = 0; ymax = 1; zmax = 1; break;
		default:
		  error();
                }

              Real coef;
	      switch (i%8) // DoF type
		{
		case 0:
                  xbasis = xmax; ybasis = ymax; zbasis = zmax;
                  coef = 1.0;
                  break;
		case 1:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = zmax;
                  coef = dxdxi[xmax];
                  break;
		case 2:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dydeta[ymax];
                  break;
		case 3:
                  xbasis = xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dzdzeta[zmax];
                  break;
		case 4:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dxdxi[xmax] * dydeta[ymax];
                  break;
		case 5:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dxdxi[xmax] * dzdzeta[zmax];
                  break;
		case 6:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = 2+zmax;
                  coef = dydeta[xmax] * dzdzeta[zmax];
                  break;
		case 7:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = 2+zmax;
		  coef = dxdxi[xmax] * dydeta[xmax] * dzdzeta[zmax];
                  break;
		default:
		  error();
		}

	      return coef * FEHermite<1>::hermite_raw_shape(xbasis,p(0)) *
                     FEHermite<1>::hermite_raw_shape(ybasis,p(1)) *
                     FEHermite<1>::hermite_raw_shape(zbasis,p(2));
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
Real FE<3,HERMITE>::shape_deriv(const ElemType,
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
Real FE<3,HERMITE>::shape_deriv(const Elem* elem,
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
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case HEX8:
	  case HEX20:
	  case HEX27:
	    {
	      assert (i<64);

              // Is node at maximum x/y/z plane?
              unsigned int xmax, ymax, zmax, xbasis, ybasis, zbasis;
              switch (i/8) // Node number
                {
                case 0:
                  xmax = 0; ymax = 0; zmax = 0; break;
                case 1:
                  xmax = 1; ymax = 0; zmax = 0; break;
                case 2:
                  xmax = 1; ymax = 1; zmax = 0; break;
                case 3:
                  xmax = 0; ymax = 1; zmax = 0; break;
                case 4:
                  xmax = 0; ymax = 0; zmax = 1; break;
                case 5:
                  xmax = 1; ymax = 0; zmax = 1; break;
                case 6:
                  xmax = 1; ymax = 1; zmax = 1; break;
                case 7:
                  xmax = 0; ymax = 1; zmax = 1; break;
		default:
		  error();
                }

              Real coef;
	      switch (i%8) // DoF type
		{
		case 0:
                  xbasis = xmax; ybasis = ymax; zbasis = zmax;
                  coef = 1.0;
                  break;
		case 1:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = zmax;
                  coef = dxdxi[xmax];
                  break;
		case 2:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dydeta[ymax];
                  break;
		case 3:
                  xbasis = xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dzdzeta[zmax];
                  break;
		case 4:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dxdxi[xmax] * dydeta[ymax];
                  break;
		case 5:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dxdxi[xmax] * dzdzeta[zmax];
                  break;
		case 6:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = 2+zmax;
                  coef = dydeta[xmax] * dzdzeta[zmax];
                  break;
		case 7:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = 2+zmax;
		  coef = dxdxi[xmax] * dydeta[xmax] * dzdzeta[zmax];
                  break;
		default:
		  error();
		}

              switch (j) // Derivative type
		{
		case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape(zbasis,p(2));
                  break;
		case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape(zbasis,p(2));
                  break;
		case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(zbasis,p(2));
                  break;
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
Real FE<3,HERMITE>::shape_second_deriv(const Elem* elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point& p)
{
  assert (elem != NULL);

  hermite_compute_coefs(elem);

  const ElemType type = elem->type();
  
  switch (order)
    {      
      // 3rd-order tricubic Hermite functions
    case THIRD:
      {
	switch (type)
	  {
	  case HEX8:
	  case HEX20:
	  case HEX27:
	    {
	      assert (i<64);

              // Is node at maximum x/y/z plane?
              unsigned int xmax, ymax, zmax, xbasis, ybasis, zbasis;
              switch (i/8) // Node number
                {
                case 0:
                  xmax = 0; ymax = 0; zmax = 0; break;
                case 1:
                  xmax = 1; ymax = 0; zmax = 0; break;
                case 2:
                  xmax = 1; ymax = 1; zmax = 0; break;
                case 3:
                  xmax = 0; ymax = 1; zmax = 0; break;
                case 4:
                  xmax = 0; ymax = 0; zmax = 1; break;
                case 5:
                  xmax = 1; ymax = 0; zmax = 1; break;
                case 6:
                  xmax = 1; ymax = 1; zmax = 1; break;
                case 7:
                  xmax = 0; ymax = 1; zmax = 1; break;
		default:
		  error();
                }

              Real coef;
	      switch (i%8) // DoF type
		{
		case 0:
                  xbasis = xmax; ybasis = ymax; zbasis = zmax;
                  coef = 1.0;
                  break;
		case 1:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = zmax;
                  coef = dxdxi[xmax];
                  break;
		case 2:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dydeta[ymax];
                  break;
		case 3:
                  xbasis = xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dzdzeta[zmax];
                  break;
		case 4:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = zmax;
                  coef = dxdxi[xmax] * dydeta[ymax];
                  break;
		case 5:
                  xbasis = 2+xmax; ybasis = ymax; zbasis = 2+zmax;
                  coef = dxdxi[xmax] * dzdzeta[zmax];
                  break;
		case 6:
                  xbasis = xmax; ybasis = 2+ymax; zbasis = 2+zmax;
                  coef = dydeta[xmax] * dzdzeta[zmax];
                  break;
		case 7:
                  xbasis = 2+xmax; ybasis = 2+ymax; zbasis = 2+zmax;
		  coef = dxdxi[xmax] * dydeta[xmax] * dzdzeta[zmax];
                  break;
		default:
		  error();
		}

              switch (j) // Derivative type
		{
		case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_second_deriv(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape(zbasis,p(2));
                  break;
		case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape(zbasis,p(2));
                  break;
		case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_second_deriv(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape(zbasis,p(2));
                  break;
		case 3:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(zbasis,p(2));
                  break;
		case 4:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(zbasis,p(2));
                  break;
		case 5:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(xbasis,p(0)) * 
                    FEHermite<1>::hermite_raw_shape(ybasis,p(1)) * 
                    FEHermite<1>::hermite_raw_shape_second_deriv(zbasis,p(2));
                  break;
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
