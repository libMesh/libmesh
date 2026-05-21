// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

#ifndef LIBMESH_FE_SIMPLEX_LAGRANGE_H
#define LIBMESH_FE_SIMPLEX_LAGRANGE_H

#include "libmesh/point.h"

namespace libMesh
{
namespace detail
{

LIBMESH_DEVICE_INLINE
Real tri_dzeta(const unsigned int i,
               const unsigned int j)
{
  libmesh_assert_less(i, 3);
  libmesh_assert_less(j, 2);

  switch (i)
    {
    case 0:
      return -1.;
    case 1:
      return j == 0 ? 1. : 0.;
    default:
      return j == 0 ? 0. : 1.;
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tri6_zeta_index(const unsigned int i,
                               const unsigned int j)
{
  libmesh_assert_less(i, 6);
  libmesh_assert_less(j, 2);

  switch (i)
    {
    case 0: return 0;
    case 1: return 1;
    case 2: return 2;
    case 3: return j == 0 ? 0 : 1;
    case 4: return j == 0 ? 1 : 2;
    default: return j == 0 ? 2 : 0;
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tri7_bubble_zeta_index(const unsigned int j)
{
  libmesh_assert_less(j, 3);
  return static_cast<unsigned short>(j);
}

LIBMESH_DEVICE_INLINE
Real tet_dzeta(const unsigned int i,
               const unsigned int j)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 3);

  switch (i)
    {
    case 0:
      return -1.;
    case 1:
      return j == 0 ? 1. : 0.;
    case 2:
      return j == 1 ? 1. : 0.;
    default:
      return j == 2 ? 1. : 0.;
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tet10_zeta_index(const unsigned int i,
                                const unsigned int j)
{
  libmesh_assert_less(i, 10);
  libmesh_assert_less(j, 2);

  switch (i)
    {
    case 0: return 0;
    case 1: return 1;
    case 2: return 2;
    case 3: return 3;
    case 4: return j == 0 ? 0 : 1;
    case 5: return j == 0 ? 1 : 2;
    case 6: return j == 0 ? 2 : 0;
    case 7: return j == 0 ? 0 : 3;
    case 8: return j == 0 ? 1 : 3;
    default: return j == 0 ? 2 : 3;
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tet14_bubble_zeta_index(const unsigned int i,
                                       const unsigned int j)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 3);

  switch (i)
    {
    case 0:
      return static_cast<unsigned short>(j);
    case 1:
      return j < 2 ? static_cast<unsigned short>(j) : 3;
    case 2:
      return j == 0 ? 1 : static_cast<unsigned short>(j + 1);
    default:
      return j == 0 ? 0 : static_cast<unsigned short>(j + 1);
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tet14_vertex_bubble_index(const unsigned int i,
                                         const unsigned int j)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 3);

  switch (i)
    {
    case 0:
      return j < 2 ? static_cast<unsigned short>(j) : 3;
    case 1:
      return static_cast<unsigned short>(j);
    case 2:
      return j == 0 ? 0 : static_cast<unsigned short>(j + 1);
    default:
      return static_cast<unsigned short>(j + 1);
    }
}

LIBMESH_DEVICE_INLINE
unsigned short tet14_edge_bubble_index(const unsigned int i,
                                       const unsigned int j)
{
  libmesh_assert_less(i, 6);
  libmesh_assert_less(j, 2);

  switch (i)
    {
    case 0: return static_cast<unsigned short>(j);
    case 1: return j == 0 ? 0 : 2;
    case 2: return j == 0 ? 0 : 3;
    case 3: return j == 0 ? 1 : 3;
    case 4: return j == 0 ? 1 : 2;
    default: return j == 0 ? 3 : 2;
    }
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
LIBMESH_DEVICE_INLINE
unsigned short tet_second_deriv_index(const unsigned int i,
                                      const unsigned int j)
{
  libmesh_assert_less(i, 6);
  libmesh_assert_less(j, 2);

  switch (i)
    {
    case 0: return 0;
    case 1: return static_cast<unsigned short>(j);
    case 2: return 1;
    case 3: return j == 0 ? 0 : 2;
    case 4: return j == 0 ? 1 : 2;
    default: return 2;
    }
}
#endif

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri3_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 3);

  switch (i)
    {
    case 0: return 1. - xi - eta;
    case 1: return xi;
    default: return eta;
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri3_shape_deriv(const unsigned int i,
                                  const unsigned int j)
{
  libmesh_assert_less(i, 3);
  libmesh_assert_less(j, 2);

  return tri_dzeta(i, j);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri6_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 6);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const unsigned short m = tri6_zeta_index(i, 0);
  const unsigned short n = tri6_zeta_index(i, 1);

  if (i < 3)
    return bary[m] * (2. * bary[m] - 1.);

  return 4. * bary[m] * bary[n];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri6_shape_deriv(const unsigned int i,
                                  const unsigned int j,
                                  const Real xi,
                                  const Real eta)
{
  libmesh_assert_less(i, 6);
  libmesh_assert_less(j, 2);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const unsigned short m = tri6_zeta_index(i, 0);
  const unsigned short n = tri6_zeta_index(i, 1);

  if (i < 3)
    return (4. * bary[m] - 1.) * tri_dzeta(m, j);

  return 4. * bary[n] * tri_dzeta(m, j) + 4. * bary[m] * tri_dzeta(n, j);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri7_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 7);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const Real bubble =
    bary[tri7_bubble_zeta_index(0)] *
    bary[tri7_bubble_zeta_index(1)] *
    bary[tri7_bubble_zeta_index(2)];

  if (i < 3)
    return fe_lagrange_tri6_shape(i, xi, eta) + 3. * bubble;

  if (i < 6)
    return fe_lagrange_tri6_shape(i, xi, eta) - 12. * bubble;

  return 27. * bubble;
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri7_shape_deriv(const unsigned int i,
                                  const unsigned int j,
                                  const Real xi,
                                  const Real eta)
{
  libmesh_assert_less(i, 7);
  libmesh_assert_less(j, 2);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const Real bubble_deriv =
    tri_dzeta(tri7_bubble_zeta_index(0), j) * bary[tri7_bubble_zeta_index(1)] * bary[tri7_bubble_zeta_index(2)] +
    bary[tri7_bubble_zeta_index(0)] * tri_dzeta(tri7_bubble_zeta_index(1), j) * bary[tri7_bubble_zeta_index(2)] +
    bary[tri7_bubble_zeta_index(0)] * bary[tri7_bubble_zeta_index(1)] * tri_dzeta(tri7_bubble_zeta_index(2), j);

  if (i < 3)
    return fe_lagrange_tri6_shape_deriv(i, j, xi, eta) + 3. * bubble_deriv;

  if (i < 6)
    return fe_lagrange_tri6_shape_deriv(i, j, xi, eta) - 12. * bubble_deriv;

  return 27. * bubble_deriv;
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri6_shape_second_deriv(const unsigned int i,
                                         const unsigned int j)
{
  libmesh_assert_less(i, 6);
  libmesh_assert_less(j, 3);

  const unsigned short my_j = j == 2 ? 1 : 0;
  const unsigned short my_k = j == 0 ? 0 : 1;

  if (i < 3)
    return 4. * tri_dzeta(i, my_j) * tri_dzeta(i, my_k);

  const unsigned short m = tri6_zeta_index(i, 0);
  const unsigned short n = tri6_zeta_index(i, 1);

  return 4. * (tri_dzeta(n, my_j) * tri_dzeta(m, my_k) +
               tri_dzeta(m, my_j) * tri_dzeta(n, my_k));
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri7_shape_second_deriv(const unsigned int i,
                                         const unsigned int j,
                                         const Real xi,
                                         const Real eta)
{
  libmesh_assert_less(i, 7);
  libmesh_assert_less(j, 3);

  const unsigned short my_j = j == 2 ? 1 : 0;
  const unsigned short my_k = j == 0 ? 0 : 1;
  const Real bary[3] = {1. - xi - eta, xi, eta};
  const Real bubble_second_deriv =
    tri_dzeta(tri7_bubble_zeta_index(0), my_j) * tri_dzeta(tri7_bubble_zeta_index(1), my_k) * bary[tri7_bubble_zeta_index(2)] +
    tri_dzeta(tri7_bubble_zeta_index(0), my_j) * bary[tri7_bubble_zeta_index(1)] * tri_dzeta(tri7_bubble_zeta_index(2), my_k) +
    bary[tri7_bubble_zeta_index(0)] * tri_dzeta(tri7_bubble_zeta_index(1), my_j) * tri_dzeta(tri7_bubble_zeta_index(2), my_k) +
    tri_dzeta(tri7_bubble_zeta_index(0), my_k) * tri_dzeta(tri7_bubble_zeta_index(1), my_j) * bary[tri7_bubble_zeta_index(2)] +
    tri_dzeta(tri7_bubble_zeta_index(0), my_k) * bary[tri7_bubble_zeta_index(1)] * tri_dzeta(tri7_bubble_zeta_index(2), my_j) +
    bary[tri7_bubble_zeta_index(0)] * tri_dzeta(tri7_bubble_zeta_index(1), my_k) * tri_dzeta(tri7_bubble_zeta_index(2), my_j);

  if (i < 3)
    return fe_lagrange_tri6_shape_second_deriv(i, j) + 3. * bubble_second_deriv;

  if (i < 6)
    return fe_lagrange_tri6_shape_second_deriv(i, j) - 12. * bubble_second_deriv;

  return 27. * bubble_second_deriv;
}
#endif

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet4_shape(const unsigned int i,
                            const Real xi,
                            const Real eta,
                            const Real zeta)
{
  libmesh_assert_less(i, 4);

  switch (i)
    {
    case 0: return 1. - xi - eta - zeta;
    case 1: return xi;
    case 2: return eta;
    default: return zeta;
    }
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet4_shape_deriv(const unsigned int i,
                                  const unsigned int j)
{
  libmesh_assert_less(i, 4);
  libmesh_assert_less(j, 3);

  return tet_dzeta(i, j);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet10_shape(const unsigned int i,
                             const Real xi,
                             const Real eta,
                             const Real zeta)
{
  libmesh_assert_less(i, 10);

  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  const unsigned short m = tet10_zeta_index(i, 0);
  const unsigned short n = tet10_zeta_index(i, 1);

  if (i < 4)
    return bary[m] * (2. * bary[m] - 1.);

  return 4. * bary[m] * bary[n];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet10_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta,
                                   const Real zeta)
{
  libmesh_assert_less(i, 10);
  libmesh_assert_less(j, 3);

  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  const unsigned short m = tet10_zeta_index(i, 0);
  const unsigned short n = tet10_zeta_index(i, 1);

  if (i < 4)
    return (4. * bary[m] - 1.) * tet_dzeta(m, j);

  return 4. * bary[n] * tet_dzeta(m, j) + 4. * bary[m] * tet_dzeta(n, j);
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet14_shape(const unsigned int i,
                             const Real xi,
                             const Real eta,
                             const Real zeta)
{
  libmesh_assert_less(i, 14);

  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  Real bubble[4];

  for (unsigned short b = 0; b != 4; ++b)
    bubble[b] =
      bary[tet14_bubble_zeta_index(b, 0)] *
      bary[tet14_bubble_zeta_index(b, 1)] *
      bary[tet14_bubble_zeta_index(b, 2)];

  if (i < 4)
    {
      return fe_lagrange_tet10_shape(i, xi, eta, zeta) +
             3. * (bubble[tet14_vertex_bubble_index(i, 0)] +
                   bubble[tet14_vertex_bubble_index(i, 1)] +
                   bubble[tet14_vertex_bubble_index(i, 2)]);
    }

  if (i < 10)
    {
      return fe_lagrange_tet10_shape(i, xi, eta, zeta) -
             12. * (bubble[tet14_edge_bubble_index(i - 4, 0)] +
                    bubble[tet14_edge_bubble_index(i - 4, 1)]);
    }

  return 27. * bubble[i - 10];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet14_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta,
                                   const Real zeta)
{
  libmesh_assert_less(i, 14);
  libmesh_assert_less(j, 3);

  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  Real bubble_deriv[4];

  for (unsigned short b = 0; b != 4; ++b)
    {
      bubble_deriv[b] =
        tet_dzeta(tet14_bubble_zeta_index(b, 0), j) * bary[tet14_bubble_zeta_index(b, 1)] * bary[tet14_bubble_zeta_index(b, 2)] +
        bary[tet14_bubble_zeta_index(b, 0)] * tet_dzeta(tet14_bubble_zeta_index(b, 1), j) * bary[tet14_bubble_zeta_index(b, 2)] +
        bary[tet14_bubble_zeta_index(b, 0)] * bary[tet14_bubble_zeta_index(b, 1)] * tet_dzeta(tet14_bubble_zeta_index(b, 2), j);
    }

  if (i < 4)
    {
      return fe_lagrange_tet10_shape_deriv(i, j, xi, eta, zeta) +
             3. * (bubble_deriv[tet14_vertex_bubble_index(i, 0)] +
                   bubble_deriv[tet14_vertex_bubble_index(i, 1)] +
                   bubble_deriv[tet14_vertex_bubble_index(i, 2)]);
    }

  if (i < 10)
    {
      return fe_lagrange_tet10_shape_deriv(i, j, xi, eta, zeta) -
             12. * (bubble_deriv[tet14_edge_bubble_index(i - 4, 0)] +
                    bubble_deriv[tet14_edge_bubble_index(i - 4, 1)]);
    }

  return 27. * bubble_deriv[i - 10];
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet10_shape_second_deriv(const unsigned int i,
                                          const unsigned int j)
{
  libmesh_assert_less(i, 10);
  libmesh_assert_less(j, 6);

  const unsigned short my_j = tet_second_deriv_index(j, 0);
  const unsigned short my_k = tet_second_deriv_index(j, 1);

  if (i < 4)
    return 4. * tet_dzeta(i, my_j) * tet_dzeta(i, my_k);

  const unsigned short m = tet10_zeta_index(i, 0);
  const unsigned short n = tet10_zeta_index(i, 1);

  return 4. * (tet_dzeta(n, my_j) * tet_dzeta(m, my_k) +
               tet_dzeta(m, my_j) * tet_dzeta(n, my_k));
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet14_shape_second_deriv(const unsigned int i,
                                          const unsigned int j,
                                          const Real xi,
                                          const Real eta,
                                          const Real zeta)
{
  libmesh_assert_less(i, 14);
  libmesh_assert_less(j, 6);

  const unsigned short my_j = tet_second_deriv_index(j, 0);
  const unsigned short my_k = tet_second_deriv_index(j, 1);
  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  Real bubble_second_deriv[4];

  for (unsigned short b = 0; b != 4; ++b)
    {
      bubble_second_deriv[b] =
        tet_dzeta(tet14_bubble_zeta_index(b, 0), my_j) * tet_dzeta(tet14_bubble_zeta_index(b, 1), my_k) * bary[tet14_bubble_zeta_index(b, 2)] +
        tet_dzeta(tet14_bubble_zeta_index(b, 0), my_j) * bary[tet14_bubble_zeta_index(b, 1)] * tet_dzeta(tet14_bubble_zeta_index(b, 2), my_k) +
        bary[tet14_bubble_zeta_index(b, 0)] * tet_dzeta(tet14_bubble_zeta_index(b, 1), my_j) * tet_dzeta(tet14_bubble_zeta_index(b, 2), my_k) +
        tet_dzeta(tet14_bubble_zeta_index(b, 0), my_k) * tet_dzeta(tet14_bubble_zeta_index(b, 1), my_j) * bary[tet14_bubble_zeta_index(b, 2)] +
        tet_dzeta(tet14_bubble_zeta_index(b, 0), my_k) * bary[tet14_bubble_zeta_index(b, 1)] * tet_dzeta(tet14_bubble_zeta_index(b, 2), my_j) +
        bary[tet14_bubble_zeta_index(b, 0)] * tet_dzeta(tet14_bubble_zeta_index(b, 1), my_k) * tet_dzeta(tet14_bubble_zeta_index(b, 2), my_j);
    }

  if (i < 4)
    {
      return fe_lagrange_tet10_shape_second_deriv(i, j) +
             3. * (bubble_second_deriv[tet14_vertex_bubble_index(i, 0)] +
                   bubble_second_deriv[tet14_vertex_bubble_index(i, 1)] +
                   bubble_second_deriv[tet14_vertex_bubble_index(i, 2)]);
    }

  if (i < 10)
    {
      return fe_lagrange_tet10_shape_second_deriv(i, j) -
             12. * (bubble_second_deriv[tet14_edge_bubble_index(i - 4, 0)] +
                    bubble_second_deriv[tet14_edge_bubble_index(i - 4, 1)]);
    }

  return 27. * bubble_second_deriv[i - 10];
}
#endif

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_FE_SIMPLEX_LAGRANGE_H
