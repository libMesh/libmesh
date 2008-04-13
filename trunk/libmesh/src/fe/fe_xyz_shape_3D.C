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
Real FE<3,XYZ>::shape(const ElemType,
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
Real FE<3,XYZ>::shape(const Elem* elem,
		      const Order order,
		      const unsigned int i,
		      const Point& p)
{
#if DIM == 3
  
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
  const Real z  = p(2);
  const Real xc = centroid(0);
  const Real yc = centroid(1);
  const Real zc = centroid(2);
  const Real dx = x - xc;
  const Real dy = y - yc;
  const Real dz = z - zc;

  const unsigned int totalorder = order + elem->p_level();
  libmesh_assert (i < (static_cast<unsigned int>(totalorder)+1)*
              (static_cast<unsigned int>(totalorder)+2)*
              (static_cast<unsigned int>(totalorder)+2)/6);
    
  // monomials. since they are hierarchic we only need one case block.
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
      return dx*dy;
      
    case 6:
      return dy*dy;

    case 7:
      return dx*dz;

    case 8:
      return dz*dy;

    case 9:
      return dz*dz;

      // cubics
    case 10:
      return dx*dx*dx;

    case 11:
      return dx*dx*dy;

    case 12:
      return dx*dy*dy;

    case 13:
      return dy*dy*dy;

    case 14:
      return dx*dx*dz;

    case 15:
      return dx*dy*dz;

    case 16:
      return dy*dy*dz;

    case 17:
      return dx*dz*dz;

    case 18:
      return dy*dz*dz;

    case 19:
      return dz*dz*dz;

      // quartics
    case 20:
      return dx*dx*dx*dx;

    case 21:
      return dx*dx*dx*dy;

    case 22:
      return dx*dx*dy*dy;

    case 23:
      return dx*dy*dy*dy;

    case 24:
      return dy*dy*dy*dy;

    case 25:
      return dx*dx*dx*dz;

    case 26:
      return dx*dx*dy*dz;

    case 27:
      return dx*dy*dy*dz;

    case 28:
      return dy*dy*dy*dz;

    case 29:
      return dx*dx*dz*dz;

    case 30:
      return dx*dy*dz*dz;

    case 31:
      return dy*dy*dz*dz;

    case 32:
      return dx*dz*dz*dz;

    case 33:
      return dy*dz*dz*dz;

    case 34:
      return dz*dz*dz*dz;
      	    
    default:
      unsigned int o = 0;
      for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
      unsigned int i2 = i - (o*(o+1)*(o+2)/6);
      unsigned int block=o, nz = 0;
      for (; block < i2; block += (o-nz+1)) { nz++; }
      const unsigned int nx = block - i2;
      const unsigned int ny = o - nx - nz;
      Real val = 1.;
      for (unsigned int index=0; index != nx; index++)
        val *= dx;
      for (unsigned int index=0; index != ny; index++)
        val *= dy;
      for (unsigned int index=0; index != nz; index++)
        val *= dz;
      return val;
    }

#endif
  
  libmesh_error();
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
  
  libmesh_error();
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

  libmesh_assert (elem != NULL);
  libmesh_assert (j<3);
  
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

  const unsigned int totalorder = static_cast<Order>(order + elem->p_level());
  libmesh_assert (i < (static_cast<unsigned int>(totalorder)+1)*
              (static_cast<unsigned int>(totalorder)+2)*
              (static_cast<unsigned int>(totalorder)+2)/6);

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
  	  return dy;
  	  
  	case 6:
  	  return 0.;
  	  
  	case 7:
  	  return dz;
  	  
  	case 8:
  	  return 0.;
  	  
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	  return 3.*dx*dx;

  	case 11:
  	  return 2.*dx*dy;

  	case 12:
  	  return dy*dy;

  	case 13:
  	  return 0.;

  	case 14:
  	  return 2.*dx*dz;

  	case 15:
  	  return dy*dz;

  	case 16:
  	  return 0.;

  	case 17:
  	  return dz*dz;

  	case 18:
  	  return 0.;

  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	  return 4.*dx*dx*dx;

  	case 21:
  	  return 3.*dx*dx*dy;

  	case 22:
  	  return 2.*dx*dy*dy;

  	case 23:
  	  return dy*dy*dy;

  	case 24:
  	  return 0.;

  	case 25:
  	  return 3.*dx*dx*dz;

  	case 26:
  	  return 2.*dx*dy*dz;

  	case 27:
  	  return dy*dy*dz;

  	case 28:
  	  return 0.;

  	case 29:
  	  return 2.*dx*dz*dz;

  	case 30:
  	  return dy*dz*dz;

  	case 31:
  	  return 0.;

  	case 32:
  	  return dz*dz*dz;

  	case 33:
  	  return 0.;

  	case 34:
  	  return 0.;
  	  
        default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nx;
          for (unsigned int index=1; index < nx; index++)
            val *= dx;
          for (unsigned int index=0; index != ny; index++)
            val *= dy;
          for (unsigned int index=0; index != nz; index++)
            val *= dz;
          return val;
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
  	  return dx;
  	  
  	case 6:
  	  return 2.*dy;
  	  
  	case 7:
  	  return 0.;
  	  
  	case 8:
  	  return dz;
  	  
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	  return 0.;

  	case 11:
  	  return dx*dx;

  	case 12:
  	  return 2.*dx*dy;

  	case 13:
  	  return 3.*dy*dy;

  	case 14:
  	  return 0.;

  	case 15:
  	  return dx*dz;

  	case 16:
  	  return 2.*dy*dz;

  	case 17:
  	  return 0.;

  	case 18:
  	  return dz*dz;

  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	  return 0.;

  	case 21:
  	  return dx*dx*dx;

  	case 22:
  	  return 2.*dx*dx*dy;

  	case 23:
  	  return 3.*dx*dy*dy;

  	case 24:
  	  return 4.*dy*dy*dy;

  	case 25:
  	  return 0.;

  	case 26:
  	  return dx*dx*dz;

  	case 27:
  	  return 2.*dx*dy*dz;

  	case 28:
  	  return 3.*dy*dy*dz;

  	case 29:
  	  return 0.;

  	case 30:
  	  return dx*dz*dz;

  	case 31:
  	  return 2.*dy*dz*dz;

  	case 32:
  	  return 0.;

  	case 33:
  	  return dz*dz*dz;

  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = ny;
          for (unsigned int index=0; index != nx; index++)
            val *= dx;
          for (unsigned int index=1; index < ny; index++)
            val *= dy;
          for (unsigned int index=0; index != nz; index++)
            val *= dz;
          return val;
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
  	  return 0.;
  	  
  	case 7:
  	  return dx;
  	  
  	case 8:
  	  return dy;
  	  
  	case 9:
  	  return 2.*dz;

  	  // cubic
  	case 10:
  	  return 0.;

  	case 11:
  	  return 0.;

  	case 12:
  	  return 0.;

  	case 13:
  	  return 0.;

  	case 14:
  	  return dx*dx;

  	case 15:
  	  return dx*dy;

  	case 16:
  	  return dy*dy;

  	case 17:
  	  return 2.*dx*dz;

  	case 18:
  	  return 2.*dy*dz;

  	case 19:
  	  return 3.*dz*dz;

  	  // quartics
  	case 20:
  	  return 0.;

  	case 21:
  	  return 0.;

  	case 22:
  	  return 0.;

  	case 23:
  	  return 0.;

  	case 24:
  	  return 0.;

  	case 25:
  	  return dx*dx*dx;

  	case 26:
  	  return dx*dx*dy;

  	case 27:
  	  return dx*dy*dy;

  	case 28:
  	  return dy*dy*dy;

  	case 29:
  	  return 2.*dx*dx*dz;

  	case 30:
  	  return 2.*dx*dy*dz;

  	case 31:
  	  return 2.*dy*dy*dz;

  	case 32:
  	  return 3.*dx*dz*dz;

  	case 33:
  	  return 3.*dy*dz*dz;

  	case 34:
  	  return 4.*dz*dz*dz;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nz;
          for (unsigned int index=0; index != nx; index++)
            val *= dx;
          for (unsigned int index=0; index != ny; index++)
            val *= dy;
          for (unsigned int index=1; index < nz; index++)
            val *= dz;
          return val;
  	}
      }

      
    default:
      libmesh_error();
    }

#endif
  
  libmesh_error();
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
  
  libmesh_error();
  return 0.;
}



template <>
Real FE<3,XYZ>::shape_second_deriv(const Elem* elem,
			           const Order order,
			           const unsigned int i,
			           const unsigned int j,
			           const Point& p)
{
#if DIM == 3

  libmesh_assert (elem != NULL);
  libmesh_assert (j<3);
  
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

  const unsigned int totalorder = static_cast<Order>(order + elem->p_level());
  libmesh_assert (i < (static_cast<unsigned int>(totalorder)+1)*
              (static_cast<unsigned int>(totalorder)+2)*
              (static_cast<unsigned int>(totalorder)+2)/6);

    // monomials. since they are hierarchic we only need one case block.
  switch (j)
    {
      // d^2()/dx^2
    case 0:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	  return 2.;
  	  
  	case 5:
  	case 6:
  	case 7:
  	case 8:
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	  return 6.*dx;

  	case 11:
  	  return 2.*dy;

  	case 12:
  	case 13:
  	  return 0.;

  	case 14:
  	  return 2.*dz;

  	case 15:
  	case 16:
  	case 17:
  	case 18:
  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	  return 12.*dx*dx;

  	case 21:
  	  return 6.*dx*dy;

  	case 22:
  	  return 2.*dy*dy;

  	case 23:
  	case 24:
  	  return 0.;

  	case 25:
  	  return 6.*dx*dz;

  	case 26:
  	  return 2.*dy*dz;

  	case 27:
  	case 28:
  	  return 0.;

  	case 29:
  	  return 2.*dz*dz;

  	case 30:
  	case 31:
  	case 32:
  	case 33:
  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nx * (nx - 1);
          for (unsigned int index=2; index < nx; index++)
            val *= dx;
          for (unsigned int index=0; index != ny; index++)
            val *= dy;
          for (unsigned int index=0; index != nz; index++)
            val *= dz;
          return val;
  	}
      }


      // d^2()/dxdy
    case 1:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	  return 0.;

  	case 5:
  	  return 1.;

  	case 6:
  	case 7:
  	case 8:
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	  return 0.;

  	case 11:
  	  return 2.*dx;

  	case 12:
  	  return 2.*dy;

  	case 13:
  	case 14:
  	  return 0.;

  	case 15:
  	  return dz;

  	case 16:
  	case 17:
  	case 18:
  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	  return 0.;

  	case 21:
  	  return 3.*dx*dx;

  	case 22:
  	  return 4.*dx*dy;

  	case 23:
  	  return 3.*dy*dy;

  	case 24:
  	case 25:
  	  return 0.;

  	case 26:
  	  return 2.*dx*dz;

  	case 27:
  	  return 2.*dy*dz;

  	case 28:
  	case 29:
  	  return 0.;

  	case 30:
  	  return dz*dz;

  	case 31:
  	case 32:
  	case 33:
  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nx * ny;
          for (unsigned int index=1; index < nx; index++)
            val *= dx;
          for (unsigned int index=1; index < ny; index++)
            val *= dy;
          for (unsigned int index=0; index != nz; index++)
            val *= dz;
          return val;
  	}
      }

      
      // d^2()/dy^2
    case 2:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	case 5:
  	  return 0.;
  	  
  	case 6:
  	  return 2.;

  	case 7:
  	case 8:
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	case 11:
  	  return 0.;

  	case 12:
  	  return 2.*dx;
  	case 13:
  	  return 6.*dy;

  	case 14:
  	case 15:
  	  return 0.;

  	case 16:
  	  return 2.*dz;

  	case 17:
  	case 18:
  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	case 21:
  	  return 0.;

  	case 22:
  	  return 2.*dx*dx;

  	case 23:
  	  return 6.*dx*dy;

  	case 24:
  	  return 12.*dy*dy;

  	case 25:
  	case 26:
  	  return 0.;

  	case 27:
  	  return 2.*dx*dz;

  	case 28:
  	  return 6.*dy*dz;

  	case 29:
  	case 30:
  	  return 0.;

  	case 31:
  	  return 2.*dz*dz;

  	case 32:
  	case 33:
  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = ny * (ny - 1);
          for (unsigned int index=0; index != nx; index++)
            val *= dx;
          for (unsigned int index=2; index < ny; index++)
            val *= dy;
          for (unsigned int index=0; index != nz; index++)
            val *= dz;
          return val;
  	}
      }

      
      // d^2()/dxdz
    case 3:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	case 5:
  	case 6:
  	  return 0.;

  	case 7:
  	  return 1.;
  	  
  	case 8:
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	case 11:
  	case 12:
  	case 13:
  	  return 0.;

  	case 14:
  	  return 2.*dx;

  	case 15:
  	  return dy;

  	case 16:
  	  return 0.;

  	case 17:
  	  return 2.*dz;

  	case 18:
  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	case 21:
  	case 22:
  	case 23:
  	case 24:
  	  return 0.;

  	case 25:
  	  return 3.*dx*dx;

  	case 26:
  	  return 2.*dx*dy;

  	case 27:
  	  return dy*dy;

  	case 28:
  	  return 0.;

  	case 29:
  	  return 4.*dx*dz;

  	case 30:
  	  return 2.*dy*dz;

  	case 31:
  	  return 0.;

  	case 32:
  	  return 3.*dz*dz;

  	case 33:
  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nx * nz;
          for (unsigned int index=1; index < nx; index++)
            val *= dx;
          for (unsigned int index=0; index != ny; index++)
            val *= dy;
          for (unsigned int index=1; index < nz; index++)
            val *= dz;
          return val;
  	}
      }

      // d^2()/dydz
    case 4:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	case 5:
  	case 6:
  	case 7:
  	  return 0.;

  	case 8:
  	  return 1.;
  	  
  	case 9:
  	  return 0.;

  	  // cubic
  	case 10:
  	case 11:
  	case 12:
  	case 13:
  	case 14:
  	  return 0.;

  	case 15:
  	  return dx;

  	case 16:
  	  return 2.*dy;

  	case 17:
  	  return 0.;

  	case 18:
  	  return 2.*dz;

  	case 19:
  	  return 0.;

  	  // quartics
  	case 20:
  	case 21:
  	case 22:
  	case 23:
  	case 24:
  	case 25:
  	  return 0.;

  	case 26:
  	  return dx*dx;

  	case 27:
  	  return 2.*dx*dy;

  	case 28:
  	  return 3.*dy*dy;

  	case 29:
  	  return 0.;

  	case 30:
  	  return 2.*dx*dz;

  	case 31:
  	  return 4.*dy*dz;

  	case 32:
  	  return 0.;

  	case 33:
  	  return 3.*dz*dz;

  	case 34:
  	  return 0.;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = ny * nz;
          for (unsigned int index=0; index != nx; index++)
            val *= dx;
          for (unsigned int index=1; index < ny; index++)
            val *= dy;
          for (unsigned int index=1; index < nz; index++)
            val *= dz;
          return val;
  	}
      }


      // d^2()/dz^2
    case 5:
      {
        switch (i)
  	{
  	  // constant
  	case 0:
  	  
  	  // linear
  	case 1:
  	case 2:
  	case 3:
  	  return 0.;

  	  // quadratic
  	case 4:
  	case 5:
  	case 6:
  	case 7:
  	case 8:
  	  return 0.;

  	case 9:
  	  return 2.;

  	  // cubic
  	case 10:
  	case 11:
  	case 12:
  	case 13:
  	case 14:
  	case 15:
  	case 16:
  	  return 0.;

  	case 17:
  	  return 2.*dx;

  	case 18:
  	  return 2.*dy;

  	case 19:
  	  return 6.*dz;

  	  // quartics
  	case 20:
  	case 21:
  	case 22:
  	case 23:
  	case 24:
  	case 25:
  	case 26:
  	case 27:
  	case 28:
  	  return 0.;

  	case 29:
  	  return 2.*dx*dx;

  	case 30:
  	  return 2.*dx*dy;

  	case 31:
  	  return 2.*dy*dy;

  	case 32:
  	  return 6.*dx*dz;

  	case 33:
  	  return 6.*dy*dz;

  	case 34:
  	  return 12.*dz*dz;
  	  
  	default:
          unsigned int o = 0;
          for (; i >= (o+1)*(o+2)*(o+3)/6; o++) { }
          unsigned int i2 = i - (o*(o+1)*(o+2)/6);
          unsigned int block=o, nz = 0;
          for (; block < i2; block += (o-nz+1)) { nz++; }
          const unsigned int nx = block - i2;
          const unsigned int ny = o - nx - nz;
          Real val = nz * (nz - 1);
          for (unsigned int index=0; index != nx; index++)
            val *= dx;
          for (unsigned int index=0; index != ny; index++)
            val *= dy;
          for (unsigned int index=2; index < nz; index++)
            val *= dz;
          return val;
  	}
      }

      
    default:
      libmesh_error();
    }

#endif
  
  libmesh_error();
  return 0.;  
}
