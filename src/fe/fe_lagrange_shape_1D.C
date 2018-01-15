// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{




template <>
Real FE<1,LAGRANGE>::shape(const ElemType,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  const Real xi = p(0);


  switch (order)
    {
      // Lagrange linears
    case FIRST:
      {
        libmesh_assert_less (i, 2);

        switch (i)
          {
          case 0:
            return .5*(1. - xi);

          case 1:
            return .5*(1. + xi);

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }



      // Lagrange quadratics
    case SECOND:
      {
        libmesh_assert_less (i, 3);

        switch (i)
          {
          case 0:
            return .5*xi*(xi - 1.);

          case 1:
            return .5*xi*(xi + 1);

          case 2:
            return (1. - xi*xi);

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }



      // Lagrange cubics
    case THIRD:
      {
        libmesh_assert_less (i, 4);

        switch (i)
          {
          case 0:
            return 9./16.*(1./9.-xi*xi)*(xi-1.);

          case 1:
            return -9./16.*(1./9.-xi*xi)*(xi+1.);

          case 2:
            return 27./16.*(1.-xi*xi)*(1./3.-xi);

          case 3:
            return 27./16.*(1.-xi*xi)*(1./3.+xi);

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }

    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}



template <>
Real FE<1,LAGRANGE>::shape(const Elem * elem,
                           const Order order,
                           const unsigned int i,
                           const Point & p)
{
  libmesh_assert(elem);

  return FE<1,LAGRANGE>::shape(elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
}



template <>
Real FE<1,LAGRANGE>::shape_deriv(const ElemType,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int libmesh_dbg_var(j),
                                 const Point & p)
{
  // only d()/dxi in 1D!

  libmesh_assert_equal_to (j, 0);

  const Real xi = p(0);


  switch (order)
    {
      // Lagrange linear shape function derivatives
    case FIRST:
      {
        libmesh_assert_less (i, 2);

        switch (i)
          {
          case 0:
            return -.5;

          case 1:
            return .5;

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }


      // Lagrange quadratic shape function derivatives
    case SECOND:
      {
        libmesh_assert_less (i, 3);

        switch (i)
          {
          case 0:
            return xi-.5;

          case 1:
            return xi+.5;

          case 2:
            return -2.*xi;

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }


      // Lagrange cubic shape function derivatives
    case THIRD:
      {
        libmesh_assert_less (i, 4);

        switch (i)
          {
          case 0:
            return -9./16.*(3.*xi*xi-2.*xi-1./9.);

          case 1:
            return -9./16.*(-3.*xi*xi-2.*xi+1./9.);

          case 2:
            return 27./16.*(3.*xi*xi-2./3.*xi-1.);

          case 3:
            return 27./16.*(-3.*xi*xi-2./3.*xi+1.);

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      }


    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}



template <>
Real FE<1,LAGRANGE>::shape_deriv(const Elem * elem,
                                 const Order order,
                                 const unsigned int i,
                                 const unsigned int j,
                                 const Point & p)
{
  libmesh_assert(elem);

  return FE<1,LAGRANGE>::shape_deriv(elem->type(),
                                     static_cast<Order>(order + elem->p_level()), i, j, p);
}




template <>
Real FE<1,LAGRANGE>::shape_second_deriv(const ElemType,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int libmesh_dbg_var(j),
                                        const Point & p)
{
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!

  const Real xi = p(0);
  libmesh_assert_equal_to (j, 0);

  switch (order)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        // All second derivatives of linears are zero....
        return 0.;
      }

      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (i)
          {
          case 0:
            return 1.;

          case 1:
            return 1.;

          case 2:
            return -2.;

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      } // end case SECOND

    case THIRD:
      {
        switch (i)
          {
          case 0:
            return -9./16.*(6.*xi-2);

          case 1:
            return -9./16.*(-6*xi-2.);

          case 2:
            return 27./16.*(6*xi-2./3.);

          case 3:
            return 27./16.*(-6*xi-2./3.);

          default:
            libmesh_error_msg("Invalid shape function index i = " << i);
          }
      } // end case THIRD


    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    } // end switch (order)
}



template <>
Real FE<1,LAGRANGE>::shape_second_deriv(const Elem * elem,
                                        const Order order,
                                        const unsigned int i,
                                        const unsigned int j,
                                        const Point & p)
{
  libmesh_assert(elem);

  return FE<1,LAGRANGE>::shape_second_deriv(elem->type(),
                                            static_cast<Order>(order + elem->p_level()), i, j, p);
}

} // namespace libMesh
