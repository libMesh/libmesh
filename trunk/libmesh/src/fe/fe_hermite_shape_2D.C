// $Id: fe_hermite_shape_2D.C,v 1.4 2006-03-29 18:47:23 roystgnr Exp $

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
  // Mapping functions - derivatives at each dofpt
  std::vector<std::vector<Real> > dxdxi(2, std::vector<Real>(2, 0));
#ifdef DEBUG
  std::vector<Real> dxdeta(2), dydxi(2);
#endif



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
  dofpt.push_back(Point(-1,-1));
  dofpt.push_back(Point(1,1));

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[0][p] = 0;
      dxdxi[1][p] = 0;
#ifdef DEBUG
      dxdeta[p] = 0;
      dydxi[p] = 0;
#endif
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = FE<2,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 0, dofpt[p]);
          const Real ddeta = FE<2,LAGRANGE>::shape_deriv 
            (mapping_elem_type, mapping_order, i, 1, dofpt[p]);

          dxdxi[0][p] += elem->point(i)(0) * ddxi;
          dxdxi[1][p] += elem->point(i)(1) * ddeta;
// dxdeta and dydxi should be 0!
#ifdef DEBUG
          dxdeta[p] += elem->point(i)(0) * ddeta;
          dydxi[p] += elem->point(i)(1) * ddxi;
#endif
        }
      // No singular elements!
      assert(dxdxi[0][p]);
      assert(dxdxi[1][p]);
      // No non-rectilinear or non-axis-aligned elements!
#ifdef DEBUG
      assert(std::abs(dxdeta[p]) < 1e-9);
      assert(std::abs(dydxi[p]) < 1e-9);
#endif
    }
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = FEHermite<1>::hermite_bases(bases1D, dxdxi, i, 2);

              return coef * FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                     FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = FEHermite<1>::hermite_bases(bases1D, dxdxi, i, 2);

              switch (j)
                {
                case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
                case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1));
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = FEHermite<1>::hermite_bases(bases1D, dxdxi, i, 2);

              switch (j)
                {
                case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
                case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1));
                case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[1],p(1));
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
