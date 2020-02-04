// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FE_LAGRANGE_SHAPE_1D_H
#define LIBMESH_FE_LAGRANGE_SHAPE_1D_H

// Local includes
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/point.h"

// Inline functions useful to inline on tensor elements.

namespace libMesh
{

inline
Real fe_lagrange_1D_linear_shape(const unsigned int i,
                                 const Real xi)
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



inline
Real fe_lagrange_1D_quadratic_shape(const unsigned int i,
                                    const Real xi)
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



inline
Real fe_lagrange_1D_cubic_shape(const unsigned int i,
                                const Real xi)
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



inline
Real fe_lagrange_1D_shape(const Order order,
                          const unsigned int i,
                          const Real xi)
{
  switch (order)
    {
      // Lagrange linears
    case FIRST:
      return fe_lagrange_1D_linear_shape(i, xi);

      // Lagrange quadratics
    case SECOND:
      return fe_lagrange_1D_quadratic_shape(i, xi);

      // Lagrange cubics
    case THIRD:
      return fe_lagrange_1D_cubic_shape(i, xi);

    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}



inline
Real fe_lagrange_1D_linear_shape_deriv(const unsigned int i,
                                       const unsigned int libmesh_dbg_var(j),
                                       const Real)
{
  // only d()/dxi in 1D!
  libmesh_assert_equal_to (j, 0);

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


inline
Real fe_lagrange_1D_quadratic_shape_deriv(const unsigned int i,
                                          const unsigned int libmesh_dbg_var(j),
                                          const Real xi)
{
  // only d()/dxi in 1D!
  libmesh_assert_equal_to (j, 0);

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


inline
Real fe_lagrange_1D_cubic_shape_deriv(const unsigned int i,
                                      const unsigned int libmesh_dbg_var(j),
                                      const Real xi)
{
  // only d()/dxi in 1D!
  libmesh_assert_equal_to (j, 0);

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



inline
Real fe_lagrange_1D_shape_deriv(const Order order,
                                const unsigned int i,
                                const unsigned int j,
                                const Real xi)
{
  switch (order)
    {
    case FIRST:
      return fe_lagrange_1D_linear_shape_deriv(i, j, xi);

    case SECOND:
      return fe_lagrange_1D_quadratic_shape_deriv(i, j, xi);

    case THIRD:
      return fe_lagrange_1D_cubic_shape_deriv(i, j, xi);

    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

// fe_lagrange_1D_linear_shape_second_deriv is 0


inline
Real fe_lagrange_1D_quadratic_shape_second_deriv(const unsigned int i,
                                                 const unsigned int libmesh_dbg_var(j),
                                                 const Real)
{
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!
  libmesh_assert_equal_to (j, 0);

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
}


inline
Real fe_lagrange_1D_cubic_shape_second_deriv(const unsigned int i,
                                             const unsigned int libmesh_dbg_var(j),
                                             const Real xi)
{
  // Don't need to switch on j.  1D shape functions
  // depend on xi only!
  libmesh_assert_equal_to (j, 0);

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
}



inline
Real fe_lagrange_1D_shape_second_deriv(const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Real xi)
{
  switch (order)
    {
    // All second derivatives of linears are zero....
    case FIRST:
      return 0.;

    case SECOND:
      return fe_lagrange_1D_quadratic_shape_second_deriv(i, j, xi);

    case THIRD:
      return fe_lagrange_1D_cubic_shape_second_deriv(i, j, xi);

    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    } // end switch (order)
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

}

#endif // LIBMESH_FE_LAGRANGE_SHAPE_1D_H
