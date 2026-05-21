// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

#ifndef LIBMESH_FE_TENSOR_PRODUCT_LAGRANGE_H
#define LIBMESH_FE_TENSOR_PRODUCT_LAGRANGE_H

#include "libmesh/fe_lagrange_shape_1D.h"

namespace libMesh
{
namespace detail
{

LIBMESH_DEVICE_INLINE
unsigned int quad4_i0(const unsigned int i)
{
  libmesh_assert_less(i, 4);
  return (i == 0 || i == 3) ? 0u : 1u;
}

LIBMESH_DEVICE_INLINE
unsigned int quad4_i1(const unsigned int i)
{
  libmesh_assert_less(i, 4);
  return i < 2 ? 0u : 1u;
}

LIBMESH_DEVICE_INLINE
unsigned int quad9_i0(const unsigned int i)
{
  libmesh_assert_less(i, 9);

  switch (i)
    {
    case 0:
    case 3:
    case 7:
      return 0;
    case 1:
    case 2:
    case 5:
      return 1;
    default:
      return 2;
    }
}

LIBMESH_DEVICE_INLINE
unsigned int quad9_i1(const unsigned int i)
{
  libmesh_assert_less(i, 9);

  switch (i)
    {
    case 0:
    case 1:
    case 4:
      return 0;
    case 2:
    case 3:
    case 6:
      return 1;
    default:
      return 2;
    }
}

LIBMESH_DEVICE_INLINE
unsigned int hex8_i0(const unsigned int i)
{
  libmesh_assert_less(i, 8);
  return (i == 0 || i == 3 || i == 4 || i == 7) ? 0u : 1u;
}

LIBMESH_DEVICE_INLINE
unsigned int hex8_i1(const unsigned int i)
{
  libmesh_assert_less(i, 8);
  return (i == 0 || i == 1 || i == 4 || i == 5) ? 0u : 1u;
}

LIBMESH_DEVICE_INLINE
unsigned int hex8_i2(const unsigned int i)
{
  libmesh_assert_less(i, 8);
  return i < 4 ? 0u : 1u;
}

LIBMESH_DEVICE_INLINE
unsigned int hex27_i0(const unsigned int i)
{
  libmesh_assert_less(i, 27);

  switch (i)
    {
    case 0:
    case 3:
    case 4:
    case 7:
    case 11:
    case 12:
    case 15:
    case 19:
    case 24:
      return 0;
    case 1:
    case 2:
    case 5:
    case 6:
    case 9:
    case 13:
    case 14:
    case 17:
    case 22:
      return 1;
    default:
      return 2;
    }
}

LIBMESH_DEVICE_INLINE
unsigned int hex27_i1(const unsigned int i)
{
  libmesh_assert_less(i, 27);

  switch (i)
    {
    case 0:
    case 1:
    case 4:
    case 5:
    case 8:
    case 12:
    case 13:
    case 16:
    case 21:
      return 0;
    case 2:
    case 3:
    case 6:
    case 7:
    case 10:
    case 14:
    case 15:
    case 18:
    case 23:
      return 1;
    default:
      return 2;
    }
}

LIBMESH_DEVICE_INLINE
unsigned int hex27_i2(const unsigned int i)
{
  libmesh_assert_less(i, 27);

  switch (i)
    {
    case 0:
    case 1:
    case 2:
    case 3:
    case 8:
    case 9:
    case 10:
    case 11:
    case 20:
      return 0;
    case 4:
    case 5:
    case 6:
    case 7:
    case 16:
    case 17:
    case 18:
    case 19:
    case 25:
      return 1;
    default:
      return 2;
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad4_shape(const unsigned int i,
                             const Real xi,
                             const Real eta)
{
  libmesh_assert_less(i, 4);

  return fe_lagrange_1D_linear_shape(quad4_i0(i), xi) *
         fe_lagrange_1D_linear_shape(quad4_i1(i), eta);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad4_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 2);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_linear_shape_deriv(quad4_i0(i), 0, xi) *
             fe_lagrange_1D_linear_shape(quad4_i1(i), eta);

    default:
      return fe_lagrange_1D_linear_shape(quad4_i0(i), xi) *
             fe_lagrange_1D_linear_shape_deriv(quad4_i1(i), 0, eta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad9_shape(const unsigned int i,
                             const Real xi,
                             const Real eta)
{
  libmesh_assert_less(i, 9);

  return fe_lagrange_1D_quadratic_shape(quad9_i0(i), xi) *
         fe_lagrange_1D_quadratic_shape(quad9_i1(i), eta);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad9_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta)
{
  libmesh_assert_less(i, 9);
  libmesh_assert_less(j, 2);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_quadratic_shape_deriv(quad9_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape(quad9_i1(i), eta);

    default:
      return fe_lagrange_1D_quadratic_shape(quad9_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape_deriv(quad9_i1(i), 0, eta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex8_shape(const unsigned int i,
                            const Real xi,
                            const Real eta,
                            const Real zeta)
{
  libmesh_assert_less(i, 8);

  return fe_lagrange_1D_linear_shape(hex8_i0(i), xi) *
         fe_lagrange_1D_linear_shape(hex8_i1(i), eta) *
         fe_lagrange_1D_linear_shape(hex8_i2(i), zeta);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex8_shape_deriv(const unsigned int i,
                                  const unsigned int j,
                                  const Real xi,
                                  const Real eta,
                                  const Real zeta)
{
  libmesh_assert_less(i, 8);
  libmesh_assert_less(j, 3);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_linear_shape_deriv(hex8_i0(i), 0, xi) *
             fe_lagrange_1D_linear_shape(hex8_i1(i), eta) *
             fe_lagrange_1D_linear_shape(hex8_i2(i), zeta);

    case 1:
      return fe_lagrange_1D_linear_shape(hex8_i0(i), xi) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i1(i), 0, eta) *
             fe_lagrange_1D_linear_shape(hex8_i2(i), zeta);

    default:
      return fe_lagrange_1D_linear_shape(hex8_i0(i), xi) *
             fe_lagrange_1D_linear_shape(hex8_i1(i), eta) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i2(i), 0, zeta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex27_shape(const unsigned int i,
                             const Real xi,
                             const Real eta,
                             const Real zeta)
{
  libmesh_assert_less(i, 27);

  return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
         fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
         fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex27_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta,
                                   const Real zeta)
{
  libmesh_assert_less(i, 27);
  libmesh_assert_less(j, 3);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_quadratic_shape_deriv(hex27_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
             fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);

    case 1:
      return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i1(i), 0, eta) *
             fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);

    default:
      return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i2(i), 0, zeta);
    }
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad4_shape_second_deriv(const unsigned int i,
                                          const unsigned int j,
                                          const Real xi,
                                          const Real eta)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 3);

  switch (j)
    {
    case 0:
    case 2:
      return 0.;

    default:
      return fe_lagrange_1D_linear_shape_deriv(quad4_i0(i), 0, xi) *
             fe_lagrange_1D_linear_shape_deriv(quad4_i1(i), 0, eta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_quad9_shape_second_deriv(const unsigned int i,
                                          const unsigned int j,
                                          const Real xi,
                                          const Real eta)
{
  libmesh_assert_less(i, 9);
  libmesh_assert_less(j, 3);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_quadratic_shape_second_deriv(quad9_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape(quad9_i1(i), eta);

    case 1:
      return fe_lagrange_1D_quadratic_shape_deriv(quad9_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape_deriv(quad9_i1(i), 0, eta);

    default:
      return fe_lagrange_1D_quadratic_shape(quad9_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape_second_deriv(quad9_i1(i), 0, eta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex8_shape_second_deriv(const unsigned int i,
                                         const unsigned int j,
                                         const Real xi,
                                         const Real eta,
                                         const Real zeta)
{
  libmesh_assert_less(i, 8);
  libmesh_assert_less(j, 6);

  switch (j)
    {
    case 0:
    case 2:
    case 5:
      return 0.;

    case 1:
      return fe_lagrange_1D_linear_shape_deriv(hex8_i0(i), 0, xi) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i1(i), 0, eta) *
             fe_lagrange_1D_linear_shape(hex8_i2(i), zeta);

    case 3:
      return fe_lagrange_1D_linear_shape_deriv(hex8_i0(i), 0, xi) *
             fe_lagrange_1D_linear_shape(hex8_i1(i), eta) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i2(i), 0, zeta);

    default:
      return fe_lagrange_1D_linear_shape(hex8_i0(i), xi) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i1(i), 0, eta) *
             fe_lagrange_1D_linear_shape_deriv(hex8_i2(i), 0, zeta);
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_hex27_shape_second_deriv(const unsigned int i,
                                          const unsigned int j,
                                          const Real xi,
                                          const Real eta,
                                          const Real zeta)
{
  libmesh_assert_less(i, 27);
  libmesh_assert_less(j, 6);

  switch (j)
    {
    case 0:
      return fe_lagrange_1D_quadratic_shape_second_deriv(hex27_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
             fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);

    case 1:
      return fe_lagrange_1D_quadratic_shape_deriv(hex27_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i1(i), 0, eta) *
             fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);

    case 2:
      return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape_second_deriv(hex27_i1(i), 0, eta) *
             fe_lagrange_1D_quadratic_shape(hex27_i2(i), zeta);

    case 3:
      return fe_lagrange_1D_quadratic_shape_deriv(hex27_i0(i), 0, xi) *
             fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i2(i), 0, zeta);

    case 4:
      return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i1(i), 0, eta) *
             fe_lagrange_1D_quadratic_shape_deriv(hex27_i2(i), 0, zeta);

    default:
      return fe_lagrange_1D_quadratic_shape(hex27_i0(i), xi) *
             fe_lagrange_1D_quadratic_shape(hex27_i1(i), eta) *
             fe_lagrange_1D_quadratic_shape_second_deriv(hex27_i2(i), 0, zeta);
    }
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_FE_TENSOR_PRODUCT_LAGRANGE_H
