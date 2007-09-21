// $Id: fe_hermite_shape_3D.C,v 1.7 2006-04-25 22:31:34 roystgnr Exp $

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
#include "number_lookups.h"


// Anonymous namespace for persistant variables.
// This allows us to cache the global-to-local mapping transformation
// FIXME: This should also screw up multithreading royally
namespace
{
  static unsigned int old_elem_id = libMesh::invalid_uint;
  static std::vector<std::vector<Real> > dxdxi(3, std::vector<Real>(2, 0));
#ifdef DEBUG
  static std::vector<Real> dydxi(2), dzdeta(2), dxdzeta(2);
  static std::vector<Real> dzdxi(2), dxdeta(2), dydzeta(2);
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
    FE<3,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  std::vector<Point> dofpt;
  dofpt.push_back(Point(-1,-1,-1));
  dofpt.push_back(Point(1,1,1));

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[0][p] = 0;
      dxdxi[1][p] = 0;
      dxdxi[2][p] = 0;
#ifdef DEBUG
      dydxi[p] = 0;
      dzdeta[p] = 0;
      dxdzeta[p] = 0;
      dzdxi[p] = 0;
      dxdeta[p] = 0;
      dydzeta[p] = 0;
#endif
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
          dxdxi[0][p] += elem->point(i)(0) * ddxi;
          dxdxi[1][p] += elem->point(i)(1) * ddeta;
          dxdxi[2][p] += elem->point(i)(2) * ddzeta;
#ifdef DEBUG
          dydxi[p] += elem->point(i)(1) * ddxi;
          dzdeta[p] += elem->point(i)(2) * ddeta;
          dxdzeta[p] += elem->point(i)(0) * ddzeta;
          dzdxi[p] += elem->point(i)(2) * ddxi;
          dxdeta[p] += elem->point(i)(0) * ddeta;
          dydzeta[p] += elem->point(i)(1) * ddzeta;
#endif
        }

      // No singular elements!
      assert(dxdxi[0][p]);
      assert(dxdxi[1][p]);
      assert(dxdxi[2][p]);
      // No non-rectilinear or non-axis-aligned elements!
#ifdef DEBUG
      assert(std::abs(dydxi[p]) < TOLERANCE);
      assert(std::abs(dzdeta[p]) < TOLERANCE);
      assert(std::abs(dxdzeta[p]) < TOLERANCE);
      assert(std::abs(dzdxi[p]) < TOLERANCE);
      assert(std::abs(dxdeta[p]) < TOLERANCE);
      assert(std::abs(dydzeta[p]) < TOLERANCE);
#endif
    }
}



Real hermite_bases_3D
 (std::vector<unsigned int> &bases1D,
  const std::vector<std::vector<Real> > &dxdxi,
  const Order &o,
  unsigned int i)
{
  bases1D.clear();
  bases1D.resize(3,0);
  Real coef = 1.0;

  unsigned int e = o-2;

  // Nodes
  if (i < 64)
    {
      switch (i / 8)
        {
        case 0:
          break;
        case 1:
          bases1D[0] = 1;
          break;
        case 2:
          bases1D[0] = 1;
          bases1D[1] = 1;
          break;
        case 3:
          bases1D[1] = 1;
          break;
        case 4:
          bases1D[2] = 1;
          break;
        case 5:
          bases1D[0] = 1;
          bases1D[2] = 1;
          break;
        case 6:
          bases1D[0] = 1;
          bases1D[1] = 1;
          bases1D[2] = 1;
          break;
        case 7:
          bases1D[1] = 1;
          bases1D[2] = 1;
          break;
        }

      unsigned int basisnum = i%8;
      switch (basisnum) // DoF type
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
    }
  // Edges
  else if (i < 64 + 12*4*e)
    {
      unsigned int basisnum = (i - 64) % (4*e);
      switch ((i - 64) / (4*e))
        {
        case 0:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][0];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 1:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 2:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2 + 1;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][1];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 3:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 4:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 5:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 6:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 7:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = basisnum / 4 + 4;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 8:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][0];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 9:
          bases1D[0] = basisnum % 4 / 2 * 2 + 1;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][1];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 10:
          bases1D[0] = basisnum / 4 + 4;
          bases1D[1] = basisnum % 4 / 2 * 2 + 1;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[1][1];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        case 11:
          bases1D[0] = basisnum % 4 / 2 * 2;
          bases1D[1] = basisnum / 4 + 4;
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 4 / 2)
            coef *= dxdxi[0][0];
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        }
    }
  // Faces
  else if (i < 64 + 12*4*e + 6*2*e*e)
    {
      unsigned int basisnum = (i - 64 - 12*4*e) % (2*e*e);
      switch ((i - 64 - 12*4*e) / (2*e*e))
        {
        case 0:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = square_number_row[basisnum / 2];
          bases1D[2] = basisnum % 2 * 2;
          if (basisnum % 2)
            coef *= dxdxi[2][0];
          break;
        case 1:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = basisnum % 2 * 2;
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[1][0];
          break;
        case 2:
          bases1D[0] = basisnum % 2 * 2 + 1;
          bases1D[1] = square_number_column[basisnum / 2];
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[0][1];
          break;
        case 3:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = basisnum % 2 * 2 + 1;
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[1][1];
          break;
        case 4:
          bases1D[0] = basisnum % 2 * 2;
          bases1D[1] = square_number_column[basisnum / 2];
          bases1D[2] = square_number_row[basisnum / 2];
          if (basisnum % 2)
            coef *= dxdxi[0][0];
          break;
        case 5:
          bases1D[0] = square_number_column[basisnum / 2];
          bases1D[1] = square_number_row[basisnum / 2];
          bases1D[2] = basisnum % 2 * 2 + 1;
          if (basisnum % 2)
            coef *= dxdxi[2][1];
          break;
        }
    }
  // Interior
  else
    {
      unsigned int basisnum = i - 64 - 12*4*e;
      bases1D[0] = cube_number_column[basisnum] + 4;
      bases1D[1] = cube_number_row[basisnum] + 4;
      bases1D[2] = cube_number_page[basisnum] + 4;
    }

  // No singular elements
  assert(coef);
  return coef;
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

	      return coef *
                     FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                     FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) *
                     FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

              switch (j) // Derivative type
		{
		case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
		case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
		case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
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
  
  const Order totalorder = static_cast<Order>(order + elem->p_level());
  
  switch (totalorder)
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

              std::vector<unsigned int> bases1D;

              Real coef = hermite_bases_3D(bases1D, dxdxi, totalorder, i);

              switch (j) // Derivative type
		{
		case 0:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
		case 1:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
		case 2:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[2],p(2));
                  break;
		case 3:
                  return coef *
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
                  break;
		case 4:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape_deriv(bases1D[2],p(2));
                  break;
		case 5:
                  return coef *
                    FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) * 
                    FEHermite<1>::hermite_raw_shape(bases1D[1],p(1)) * 
                    FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[2],p(2));
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
