// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

#ifndef LIBMESH_FE_SIMPLEX_LAGRANGE_H
#define LIBMESH_FE_SIMPLEX_LAGRANGE_H

#include "libmesh/point.h"

namespace libMesh
{
namespace detail
{

constexpr Real tri_dzeta[3][2] =
  {
    {-1., -1.},
    { 1.,  0.},
    { 0.,  1.}
  };

constexpr unsigned short tri6_zeta_indices[6][2] =
  {
    {0, 0},
    {1, 1},
    {2, 2},
    {0, 1},
    {1, 2},
    {2, 0}
  };

constexpr unsigned short tri7_bubble_zeta_indices[1][3] =
  {
    {0, 1, 2}
  };

constexpr Real tet_dzeta[4][3] =
  {
    {-1., -1., -1.},
    { 1.,  0.,  0.},
    { 0.,  1.,  0.},
    { 0.,  0.,  1.}
  };

constexpr unsigned short tet10_zeta_indices[10][2] =
  {
    {0, 0},
    {1, 1},
    {2, 2},
    {3, 3},
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3}
  };

constexpr unsigned short tet14_bubble_zeta_indices[4][3] =
  {
    {0, 1, 2},
    {0, 1, 3},
    {1, 2, 3},
    {0, 2, 3}
  };

constexpr unsigned short tet14_vertex_bubble_indices[4][3] =
  {
    {0, 1, 3},
    {0, 1, 2},
    {0, 2, 3},
    {1, 2, 3}
  };

constexpr unsigned short tet14_edge_bubble_indices[6][2] =
  {
    {0, 1},
    {0, 2},
    {0, 3},
    {1, 3},
    {1, 2},
    {3, 2}
  };

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
constexpr unsigned short tet_second_deriv_indices[6][2] =
  {
    {0, 0},
    {0, 1},
    {1, 1},
    {0, 2},
    {1, 2},
    {2, 2}
  };
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

  return tri_dzeta[i][j];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri6_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 6);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const unsigned short m = tri6_zeta_indices[i][0];
  const unsigned short n = tri6_zeta_indices[i][1];

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
  const unsigned short m = tri6_zeta_indices[i][0];
  const unsigned short n = tri6_zeta_indices[i][1];

  if (i < 3)
    return (4. * bary[m] - 1.) * tri_dzeta[m][j];

  return 4. * bary[n] * tri_dzeta[m][j] + 4. * bary[m] * tri_dzeta[n][j];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tri7_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 7);

  const Real bary[3] = {1. - xi - eta, xi, eta};
  const auto & bubble_indices = tri7_bubble_zeta_indices[0];
  const Real bubble =
    bary[bubble_indices[0]] * bary[bubble_indices[1]] * bary[bubble_indices[2]];

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
  const auto & bubble_indices = tri7_bubble_zeta_indices[0];
  const Real bubble_deriv =
    tri_dzeta[bubble_indices[0]][j] * bary[bubble_indices[1]] * bary[bubble_indices[2]] +
    bary[bubble_indices[0]] * tri_dzeta[bubble_indices[1]][j] * bary[bubble_indices[2]] +
    bary[bubble_indices[0]] * bary[bubble_indices[1]] * tri_dzeta[bubble_indices[2]][j];

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
    return 4. * tri_dzeta[i][my_j] * tri_dzeta[i][my_k];

  const unsigned short m = tri6_zeta_indices[i][0];
  const unsigned short n = tri6_zeta_indices[i][1];

  return 4. * (tri_dzeta[n][my_j] * tri_dzeta[m][my_k] +
               tri_dzeta[m][my_j] * tri_dzeta[n][my_k]);
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
  const auto & bubble_indices = tri7_bubble_zeta_indices[0];
  const Real bubble_second_deriv =
    tri_dzeta[bubble_indices[0]][my_j] * tri_dzeta[bubble_indices[1]][my_k] * bary[bubble_indices[2]] +
    tri_dzeta[bubble_indices[0]][my_j] * bary[bubble_indices[1]] * tri_dzeta[bubble_indices[2]][my_k] +
    bary[bubble_indices[0]] * tri_dzeta[bubble_indices[1]][my_j] * tri_dzeta[bubble_indices[2]][my_k] +
    tri_dzeta[bubble_indices[0]][my_k] * tri_dzeta[bubble_indices[1]][my_j] * bary[bubble_indices[2]] +
    tri_dzeta[bubble_indices[0]][my_k] * bary[bubble_indices[1]] * tri_dzeta[bubble_indices[2]][my_j] +
    bary[bubble_indices[0]] * tri_dzeta[bubble_indices[1]][my_k] * tri_dzeta[bubble_indices[2]][my_j];

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

  return tet_dzeta[i][j];
}

LIBMESH_DEVICE_INLINE
Real fe_lagrange_tet10_shape(const unsigned int i,
                             const Real xi,
                             const Real eta,
                             const Real zeta)
{
  libmesh_assert_less(i, 10);

  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  const unsigned short m = tet10_zeta_indices[i][0];
  const unsigned short n = tet10_zeta_indices[i][1];

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
  const unsigned short m = tet10_zeta_indices[i][0];
  const unsigned short n = tet10_zeta_indices[i][1];

  if (i < 4)
    return (4. * bary[m] - 1.) * tet_dzeta[m][j];

  return 4. * bary[n] * tet_dzeta[m][j] + 4. * bary[m] * tet_dzeta[n][j];
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
      bary[tet14_bubble_zeta_indices[b][0]] *
      bary[tet14_bubble_zeta_indices[b][1]] *
      bary[tet14_bubble_zeta_indices[b][2]];

  if (i < 4)
    {
      const auto & bubble_ids = tet14_vertex_bubble_indices[i];
      return fe_lagrange_tet10_shape(i, xi, eta, zeta) +
             3. * (bubble[bubble_ids[0]] + bubble[bubble_ids[1]] + bubble[bubble_ids[2]]);
    }

  if (i < 10)
    {
      const auto & bubble_ids = tet14_edge_bubble_indices[i - 4];
      return fe_lagrange_tet10_shape(i, xi, eta, zeta) -
             12. * (bubble[bubble_ids[0]] + bubble[bubble_ids[1]]);
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
      const auto & bubble_ids = tet14_bubble_zeta_indices[b];
      bubble_deriv[b] =
        tet_dzeta[bubble_ids[0]][j] * bary[bubble_ids[1]] * bary[bubble_ids[2]] +
        bary[bubble_ids[0]] * tet_dzeta[bubble_ids[1]][j] * bary[bubble_ids[2]] +
        bary[bubble_ids[0]] * bary[bubble_ids[1]] * tet_dzeta[bubble_ids[2]][j];
    }

  if (i < 4)
    {
      const auto & bubble_ids = tet14_vertex_bubble_indices[i];
      return fe_lagrange_tet10_shape_deriv(i, j, xi, eta, zeta) +
             3. * (bubble_deriv[bubble_ids[0]] + bubble_deriv[bubble_ids[1]] + bubble_deriv[bubble_ids[2]]);
    }

  if (i < 10)
    {
      const auto & bubble_ids = tet14_edge_bubble_indices[i - 4];
      return fe_lagrange_tet10_shape_deriv(i, j, xi, eta, zeta) -
             12. * (bubble_deriv[bubble_ids[0]] + bubble_deriv[bubble_ids[1]]);
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

  const unsigned short my_j = tet_second_deriv_indices[j][0];
  const unsigned short my_k = tet_second_deriv_indices[j][1];

  if (i < 4)
    return 4. * tet_dzeta[i][my_j] * tet_dzeta[i][my_k];

  const unsigned short m = tet10_zeta_indices[i][0];
  const unsigned short n = tet10_zeta_indices[i][1];

  return 4. * (tet_dzeta[n][my_j] * tet_dzeta[m][my_k] +
               tet_dzeta[m][my_j] * tet_dzeta[n][my_k]);
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

  const unsigned short my_j = tet_second_deriv_indices[j][0];
  const unsigned short my_k = tet_second_deriv_indices[j][1];
  const Real bary[4] = {1. - xi - eta - zeta, xi, eta, zeta};
  Real bubble_second_deriv[4];

  for (unsigned short b = 0; b != 4; ++b)
    {
      const auto & bubble_ids = tet14_bubble_zeta_indices[b];
      bubble_second_deriv[b] =
        tet_dzeta[bubble_ids[0]][my_j] * tet_dzeta[bubble_ids[1]][my_k] * bary[bubble_ids[2]] +
        tet_dzeta[bubble_ids[0]][my_j] * bary[bubble_ids[1]] * tet_dzeta[bubble_ids[2]][my_k] +
        bary[bubble_ids[0]] * tet_dzeta[bubble_ids[1]][my_j] * tet_dzeta[bubble_ids[2]][my_k] +
        tet_dzeta[bubble_ids[0]][my_k] * tet_dzeta[bubble_ids[1]][my_j] * bary[bubble_ids[2]] +
        tet_dzeta[bubble_ids[0]][my_k] * bary[bubble_ids[1]] * tet_dzeta[bubble_ids[2]][my_j] +
        bary[bubble_ids[0]] * tet_dzeta[bubble_ids[1]][my_k] * tet_dzeta[bubble_ids[2]][my_j];
    }

  if (i < 4)
    {
      const auto & bubble_ids = tet14_vertex_bubble_indices[i];
      return fe_lagrange_tet10_shape_second_deriv(i, j) +
             3. * (bubble_second_deriv[bubble_ids[0]] + bubble_second_deriv[bubble_ids[1]] + bubble_second_deriv[bubble_ids[2]]);
    }

  if (i < 10)
    {
      const auto & bubble_ids = tet14_edge_bubble_indices[i - 4];
      return fe_lagrange_tet10_shape_second_deriv(i, j) -
             12. * (bubble_second_deriv[bubble_ids[0]] + bubble_second_deriv[bubble_ids[1]]);
    }

  return 27. * bubble_second_deriv[i - 10];
}
#endif

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_FE_SIMPLEX_LAGRANGE_H
