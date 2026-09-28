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

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_FE_SIMPLEX_LAGRANGE_H
#define LIBMESH_FE_SIMPLEX_LAGRANGE_H

// Lagrange shape functions and their derivatives for the simplex
// topologies, in one place rather than spelled out per element order in
// fe_lagrange_shape_2D.C and fe_lagrange_shape_3D.C.
//
// They are inline and call nothing but each other and the 1D shapes,
// which is what the Kokkos work needs in order to evaluate the same
// expressions on a device.  The assertions are the exception: those are
// ours, host-only, and a device build will have to say what it wants
// done with them.
//
// Everything here is expressed in the barycentric coordinates of the
// reference simplex -- zeta_0 = 1 - xi - eta (- zeta) and the remaining
// zeta_i the reference coordinates themselves -- because that is the
// form in which the higher-order members stay readable: each is the
// order below it plus a correction, rather than a fresh set of
// polynomials.

#include "libmesh/libmesh_device.h"
#include "libmesh/point.h"

namespace libMesh
{
namespace detail
{

/**
 * \returns The derivative of barycentric coordinate \p i with respect
 * to reference coordinate \p j on the triangle; a constant, since the
 * coordinates are affine in xi and eta.
 */
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

/**
 * \returns Which barycentric coordinate the \p j th factor of Tri6
 * shape \p i uses: a vertex shape uses its own coordinate twice, and a
 * mid-edge shape the two coordinates of the vertices it lies between.
 */
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
/**
 * The linear triangle's shape functions are its barycentric
 * coordinates.
 */
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
/**
 * Quadratic on the triangle: zeta*(2*zeta - 1) at a vertex, and
 * 4*zeta_a*zeta_b at the node between vertices a and b, which
 * tri6_zeta_index() names.
 */
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
/**
 * Tri6 plus the interior bubble 27*zeta_0*zeta_1*zeta_2, with the vertex
 * and edge shapes corrected by three and minus twelve bubbles so that
 * each still vanishes at every node but its own.  The Tri6 halves are
 * written out rather than obtained from fe_lagrange_tri6_shape() so that
 * the arithmetic, and so the rounding, is the same as the .C file this
 * came from.
 */
Real fe_lagrange_tri7_shape(const unsigned int i,
                            const Real xi,
                            const Real eta)
{
  libmesh_assert_less(i, 7);

  const Real zeta1 = xi;
  const Real zeta2 = eta;
  const Real zeta0 = 1. - zeta1 - zeta2;
  const Real bubble_27th = zeta0*zeta1*zeta2;

  switch (i)
    {
    case 0: return 2.*zeta0*(zeta0-0.5) + 3.*bubble_27th;
    case 1: return 2.*zeta1*(zeta1-0.5) + 3.*bubble_27th;
    case 2: return 2.*zeta2*(zeta2-0.5) + 3.*bubble_27th;
    case 3: return 4.*zeta0*zeta1 - 12.*bubble_27th;
    case 4: return 4.*zeta1*zeta2 - 12.*bubble_27th;
    case 5: return 4.*zeta2*zeta0 - 12.*bubble_27th;
    default: return 27.*bubble_27th;
    }
}

LIBMESH_DEVICE_INLINE
/**
 * Differentiating zeta_0*zeta_1*zeta_2 and substituting zeta_0 = 1 - xi -
 * eta, zeta_1 = xi, zeta_2 = eta collapses the product rule to one term
 * per direction.
 */
Real fe_lagrange_tri7_shape_deriv(const unsigned int i,
                                  const unsigned int j,
                                  const Real xi,
                                  const Real eta)
{
  libmesh_assert_less(i, 7);
  libmesh_assert_less(j, 2);

  const Real zeta1 = xi;
  const Real zeta2 = eta;
  const Real zeta0 = 1. - zeta1 - zeta2;

  const Real dzeta0dxi  = -1.;
  const Real dzeta1dxi  = 1.;
  const Real dzeta2dxi  = 0.;
  const Real dbubbledxi = zeta2 * (1. - 2.*zeta1 - zeta2);

  const Real dzeta0deta = -1.;
  const Real dzeta1deta = 0.;
  const Real dzeta2deta = 1.;
  const Real dbubbledeta= zeta1 * (1. - zeta1 - 2.*zeta2);

  if (j == 0)
    switch (i)
      {
      case 0: return (4.*zeta0-1.)*dzeta0dxi + 3.*dbubbledxi;
      case 1: return (4.*zeta1-1.)*dzeta1dxi + 3.*dbubbledxi;
      case 2: return (4.*zeta2-1.)*dzeta2dxi + 3.*dbubbledxi;
      case 3: return 4.*zeta1*dzeta0dxi + 4.*zeta0*dzeta1dxi - 12.*dbubbledxi;
      case 4: return 4.*zeta2*dzeta1dxi + 4.*zeta1*dzeta2dxi - 12.*dbubbledxi;
      case 5: return 4.*zeta2*dzeta0dxi + 4*zeta0*dzeta2dxi - 12.*dbubbledxi;
      default: return 27.*dbubbledxi;
      }

  switch (i)
    {
    case 0: return (4.*zeta0-1.)*dzeta0deta + 3.*dbubbledeta;
    case 1: return (4.*zeta1-1.)*dzeta1deta + 3.*dbubbledeta;
    case 2: return (4.*zeta2-1.)*dzeta2deta + 3.*dbubbledeta;
    case 3: return 4.*zeta1*dzeta0deta + 4.*zeta0*dzeta1deta - 12.*dbubbledeta;
    case 4: return 4.*zeta2*dzeta1deta + 4.*zeta1*dzeta2deta - 12.*dbubbledeta;
    case 5: return 4.*zeta2*dzeta0deta + 4*zeta0*dzeta2deta - 12.*dbubbledeta;
    default: return 27.*dbubbledeta;
    }
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
/**
 * The bubble is quadratic in each direction at most, so its second
 * derivatives are linear: -2*eta in xi, -2*xi in eta, and 1 - 2*xi -
 * 2*eta mixed.  The Tri6 halves are constants, since a quadratic's
 * second derivative does not depend on the point.
 */
Real fe_lagrange_tri7_shape_second_deriv(const unsigned int i,
                                         const unsigned int j,
                                         const Real xi,
                                         const Real eta)
{
  libmesh_assert_less(i, 7);
  libmesh_assert_less(j, 3);

  const Real zeta1 = xi;
  const Real zeta2 = eta;

  const Real dzeta0dxi  = -1.;
  const Real dzeta1dxi  = 1.;
  const Real dzeta2dxi  = 0.;
  const Real d2bubbledxi2 = -2. * zeta2;

  const Real dzeta0deta = -1.;
  const Real dzeta1deta = 0.;
  const Real dzeta2deta = 1.;
  const Real d2bubbledeta2 = -2. * zeta1;

  const Real d2bubbledxideta = (1. - 2.*zeta1 - 2.*zeta2);

  if (j == 0)
    switch (i)
      {
      case 0: return 4.*dzeta0dxi*dzeta0dxi + 3.*d2bubbledxi2;
      case 1: return 4.*dzeta1dxi*dzeta1dxi + 3.*d2bubbledxi2;
      case 2: return 4.*dzeta2dxi*dzeta2dxi + 3.*d2bubbledxi2;
      case 3: return 8.*dzeta0dxi*dzeta1dxi - 12.*d2bubbledxi2;
      case 4: return 8.*dzeta1dxi*dzeta2dxi - 12.*d2bubbledxi2;
      case 5: return 8.*dzeta0dxi*dzeta2dxi - 12.*d2bubbledxi2;
      default: return 27.*d2bubbledxi2;
      }

  if (j == 1)
    switch (i)
      {
      case 0: return 4.*dzeta0dxi*dzeta0deta + 3.*d2bubbledxideta;
      case 1: return 4.*dzeta1dxi*dzeta1deta + 3.*d2bubbledxideta;
      case 2: return 4.*dzeta2dxi*dzeta2deta + 3.*d2bubbledxideta;
      case 3: return 4.*dzeta1deta*dzeta0dxi + 4.*dzeta0deta*dzeta1dxi - 12.*d2bubbledxideta;
      case 4: return 4.*dzeta2deta*dzeta1dxi + 4.*dzeta1deta*dzeta2dxi - 12.*d2bubbledxideta;
      case 5: return 4.*dzeta2deta*dzeta0dxi + 4.*dzeta0deta*dzeta2dxi - 12.*d2bubbledxideta;
      default: return 27.*d2bubbledxideta;
      }

  switch (i)
    {
    case 0: return 4.*dzeta0deta*dzeta0deta + 3.*d2bubbledeta2;
    case 1: return 4.*dzeta1deta*dzeta1deta + 3.*d2bubbledeta2;
    case 2: return 4.*dzeta2deta*dzeta2deta + 3.*d2bubbledeta2;
    case 3: return 8.*dzeta0deta*dzeta1deta - 12.*d2bubbledeta2;
    case 4: return 8.*dzeta1deta*dzeta2deta - 12.*d2bubbledeta2;
    case 5: return 8.*dzeta0deta*dzeta2deta - 12.*d2bubbledeta2;
    default: return 27.*d2bubbledeta2;
    }
}

#endif

LIBMESH_DEVICE_INLINE
/**
 * The linear tet's shape functions are its barycentric coordinates,
 * as for the triangle.
 */
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
/**
 * Quadratic on the tet, in the same two forms as Tri6:
 * zeta*(2*zeta - 1) at a vertex and 4*zeta_a*zeta_b at a mid-edge node.
 */
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
/**
 * Tet10 plus one bubble per face, zeta_a*zeta_b*zeta_c over the face's
 * three vertices.  A vertex shape picks up the three bubbles of the faces
 * it touches, an edge shape loses twelve of each bubble on the two faces
 * sharing that edge, and each face node is 27 times its own bubble.  The
 * Tet10 halves are written out rather than obtained from
 * fe_lagrange_tet10_shape() so that the arithmetic, and so the rounding,
 * is the same as the .C file this came from.
 */
Real fe_lagrange_tet14_shape(const unsigned int i,
                             const Real xi,
                             const Real eta,
                             const Real zeta)
{
  libmesh_assert_less(i, 14);

  // Area coordinates, pg. 205, Vol. I, Carey, Oden, Becker FEM
  const Real zeta1 = xi;
  const Real zeta2 = eta;
  const Real zeta3 = zeta;
  const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

  // Bubble functions (not yet scaled) on side nodes
  const Real bubble_012 = zeta0*zeta1*zeta2;
  const Real bubble_013 = zeta0*zeta1*zeta3;
  const Real bubble_123 = zeta1*zeta2*zeta3;
  const Real bubble_023 = zeta0*zeta2*zeta3;

  switch (i)
    {
    case 0: return zeta0*(2.*zeta0 - 1.) + 3.*(bubble_012+bubble_013+bubble_023);
    case 1: return zeta1*(2.*zeta1 - 1.) + 3.*(bubble_012+bubble_013+bubble_123);
    case 2: return zeta2*(2.*zeta2 - 1.) + 3.*(bubble_012+bubble_023+bubble_123);
    case 3: return zeta3*(2.*zeta3 - 1.) + 3.*(bubble_013+bubble_023+bubble_123);
    case 4: return 4.*zeta0*zeta1 - 12.*(bubble_012+bubble_013);
    case 5: return 4.*zeta1*zeta2 - 12.*(bubble_012+bubble_123);
    case 6: return 4.*zeta2*zeta0 - 12.*(bubble_012+bubble_023);
    case 7: return 4.*zeta0*zeta3 - 12.*(bubble_013+bubble_023);
    case 8: return 4.*zeta1*zeta3 - 12.*(bubble_013+bubble_123);
    case 9: return 4.*zeta2*zeta3 - 12.*(bubble_023+bubble_123);
    case 10: return 27.*bubble_012;
    case 11: return 27.*bubble_013;
    case 12: return 27.*bubble_123;
    default: return 27.*bubble_023;
    }
}

LIBMESH_DEVICE_INLINE
/**
 * Tet10 plus one cubic bubble per face, each bubble differentiated in
 * place: substituting zeta_0 = 1 - xi - eta - zeta leaves a single
 * product per face and direction, so the face opposite the direction's
 * vertex loses a factor outright while the faces containing it pick up
 * the difference of the two coordinates that vary.
 */
Real fe_lagrange_tet14_shape_deriv(const unsigned int i,
                                   const unsigned int j,
                                   const Real xi,
                                   const Real eta,
                                   const Real zeta)
{
  libmesh_assert_less(i, 14);
  libmesh_assert_less(j, 3);

  const Real zeta1 = xi;
  const Real zeta2 = eta;
  const Real zeta3 = zeta;
  const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

  const Real dzeta0dxi = -1.;
  const Real dzeta1dxi =  1.;
  const Real dzeta2dxi =  0.;
  const Real dzeta3dxi =  0.;
  const Real dbubble012dxi = (zeta0-zeta1)*zeta2;
  const Real dbubble013dxi = (zeta0-zeta1)*zeta3;
  const Real dbubble123dxi = zeta2*zeta3;
  const Real dbubble023dxi = -zeta2*zeta3;

  const Real dzeta0deta = -1.;
  const Real dzeta1deta =  0.;
  const Real dzeta2deta =  1.;
  const Real dzeta3deta =  0.;
  const Real dbubble012deta = (zeta0-zeta2)*zeta1;
  const Real dbubble013deta = -zeta1*zeta3;
  const Real dbubble123deta = zeta1*zeta3;
  const Real dbubble023deta = (zeta0-zeta2)*zeta3;

  const Real dzeta0dzeta = -1.;
  const Real dzeta1dzeta =  0.;
  const Real dzeta2dzeta =  0.;
  const Real dzeta3dzeta =  1.;
  const Real dbubble012dzeta = -zeta1*zeta2;
  const Real dbubble013dzeta = (zeta0-zeta3)*zeta1;
  const Real dbubble123dzeta = zeta1*zeta2;
  const Real dbubble023dzeta = (zeta0-zeta3)*zeta2;

  if (j == 0)
    switch (i)
      {
        case 0: return (4.*zeta0 - 1.)*dzeta0dxi + 3.*(dbubble012dxi+dbubble013dxi+dbubble023dxi);
        case 1: return (4.*zeta1 - 1.)*dzeta1dxi + 3.*(dbubble012dxi+dbubble013dxi+dbubble123dxi);
        case 2: return (4.*zeta2 - 1.)*dzeta2dxi + 3.*(dbubble012dxi+dbubble023dxi+dbubble123dxi);
        case 3: return (4.*zeta3 - 1.)*dzeta3dxi + 3.*(dbubble013dxi+dbubble023dxi+dbubble123dxi);
        case 4: return 4.*(zeta0*dzeta1dxi + dzeta0dxi*zeta1) - 12.*(dbubble012dxi+dbubble013dxi);
        case 5: return 4.*(zeta1*dzeta2dxi + dzeta1dxi*zeta2) - 12.*(dbubble012dxi+dbubble123dxi);
        case 6: return 4.*(zeta0*dzeta2dxi + dzeta0dxi*zeta2) - 12.*(dbubble012dxi+dbubble023dxi);
        case 7: return 4.*(zeta0*dzeta3dxi + dzeta0dxi*zeta3) - 12.*(dbubble013dxi+dbubble023dxi);
        case 8: return 4.*(zeta1*dzeta3dxi + dzeta1dxi*zeta3) - 12.*(dbubble013dxi+dbubble123dxi);
        case 9: return 4.*(zeta2*dzeta3dxi + dzeta2dxi*zeta3) - 12.*(dbubble023dxi+dbubble123dxi);
        case 10: return 27.*dbubble012dxi;
        case 11: return 27.*dbubble013dxi;
        case 12: return 27.*dbubble123dxi;
        default: return 27.*dbubble023dxi;
      }

  if (j == 1)
    switch (i)
      {
        case 0: return (4.*zeta0 - 1.)*dzeta0deta + 3.*(dbubble012deta+dbubble013deta+dbubble023deta);
        case 1: return (4.*zeta1 - 1.)*dzeta1deta + 3.*(dbubble012deta+dbubble013deta+dbubble123deta);
        case 2: return (4.*zeta2 - 1.)*dzeta2deta + 3.*(dbubble012deta+dbubble023deta+dbubble123deta);
        case 3: return (4.*zeta3 - 1.)*dzeta3deta + 3.*(dbubble013deta+dbubble023deta+dbubble123deta);
        case 4: return 4.*(zeta0*dzeta1deta + dzeta0deta*zeta1) - 12.*(dbubble012deta+dbubble013deta);
        case 5: return 4.*(zeta1*dzeta2deta + dzeta1deta*zeta2) - 12.*(dbubble012deta+dbubble123deta);
        case 6: return 4.*(zeta0*dzeta2deta + dzeta0deta*zeta2) - 12.*(dbubble012deta+dbubble023deta);
        case 7: return 4.*(zeta0*dzeta3deta + dzeta0deta*zeta3) - 12.*(dbubble013deta+dbubble023deta);
        case 8: return 4.*(zeta1*dzeta3deta + dzeta1deta*zeta3) - 12.*(dbubble013deta+dbubble123deta);
        case 9: return 4.*(zeta2*dzeta3deta + dzeta2deta*zeta3) - 12.*(dbubble023deta+dbubble123deta);
        case 10: return 27.*dbubble012deta;
        case 11: return 27.*dbubble013deta;
        case 12: return 27.*dbubble123deta;
        default: return 27.*dbubble023deta;
      }

  switch (i)
    {
      case 0: return (4.*zeta0 - 1.)*dzeta0dzeta + 3.*(dbubble012dzeta+dbubble013dzeta+dbubble023dzeta);
      case 1: return (4.*zeta1 - 1.)*dzeta1dzeta + 3.*(dbubble012dzeta+dbubble013dzeta+dbubble123dzeta);
      case 2: return (4.*zeta2 - 1.)*dzeta2dzeta + 3.*(dbubble012dzeta+dbubble023dzeta+dbubble123dzeta);
      case 3: return (4.*zeta3 - 1.)*dzeta3dzeta + 3.*(dbubble013dzeta+dbubble023dzeta+dbubble123dzeta);
      case 4: return 4.*(zeta0*dzeta1dzeta + dzeta0dzeta*zeta1) - 12.*(dbubble012dzeta+dbubble013dzeta);
      case 5: return 4.*(zeta1*dzeta2dzeta + dzeta1dzeta*zeta2) - 12.*(dbubble012dzeta+dbubble123dzeta);
      case 6: return 4.*(zeta0*dzeta2dzeta + dzeta0dzeta*zeta2) - 12.*(dbubble012dzeta+dbubble023dzeta);
      case 7: return 4.*(zeta0*dzeta3dzeta + dzeta0dzeta*zeta3) - 12.*(dbubble013dzeta+dbubble023dzeta);
      case 8: return 4.*(zeta1*dzeta3dzeta + dzeta1dzeta*zeta3) - 12.*(dbubble013dzeta+dbubble123dzeta);
      case 9: return 4.*(zeta2*dzeta3dzeta + dzeta2dzeta*zeta3) - 12.*(dbubble023dzeta+dbubble123dzeta);
      case 10: return 27.*dbubble012dzeta;
      case 11: return 27.*dbubble013dzeta;
      case 12: return 27.*dbubble123dzeta;
      default: return 27.*dbubble023dzeta;
    }
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
/**
 * The Tet10 half of each second derivative is a constant, so it comes
 * from fe_lagrange_tet10_shape_second_deriv(); the face nodes have no
 * Tet10 half at all.  Each bubble is trilinear in the barycentric
 * coordinates, so its second derivatives are linear in the point.
 */
Real fe_lagrange_tet14_shape_second_deriv(const unsigned int i,
                                          const unsigned int j,
                                          const Real xi,
                                          const Real eta,
                                          const Real zeta)
{
  libmesh_assert_less(i, 14);
  libmesh_assert_less(j, 6);

  const Real returnval = (i < 10) ? fe_lagrange_tet10_shape_second_deriv(i, j) : 0.;

  const Real zeta1 = xi;
  const Real zeta2 = eta;
  const Real zeta3 = zeta;
  const Real zeta0 = 1. - zeta1 - zeta2 - zeta3;

  // Fill these with whichever derivative we're concerned with
  Real d2bubble012, d2bubble013, d2bubble023, d2bubble123;
  switch (j)
    {
      // d^2()/dxi^2
    case 0:
      d2bubble012 = -2.*zeta2;
      d2bubble013 = -2.*zeta3;
      d2bubble023 = 0.;
      d2bubble123 = 0.;
      break;

      // d^2()/dxideta
    case 1:
      d2bubble012 = (zeta0-zeta1)-zeta2;
      d2bubble013 = -zeta3;
      d2bubble123 = zeta3;
      d2bubble023 = -zeta3;
      break;

      // d^2()/deta^2
    case 2:
      d2bubble012 = -2.*zeta1;
      d2bubble013 = 0.;
      d2bubble123 = 0.;
      d2bubble023 = -2.*zeta3;
      break;

      // d^2()/dxi dzeta
    case 3:
      d2bubble012 = -zeta2;
      d2bubble013 = (zeta0-zeta3)-zeta1;
      d2bubble123 = zeta2;
      d2bubble023 = -zeta2;
      break;

      // d^2()/deta dzeta
    case 4:
      d2bubble012 = -zeta1;
      d2bubble013 = -zeta1;
      d2bubble123 = zeta1;
      d2bubble023 = (zeta0-zeta3)-zeta2;
      break;

      // d^2()/dzeta^2
    default:
      d2bubble012 = 0.;
      d2bubble013 = -2.*zeta1;
      d2bubble123 = 0.;
      d2bubble023 = -2.*zeta2;
      break;
    }

  switch (i)
    {
    case 0: return returnval + 3.*(d2bubble012+d2bubble013+d2bubble023);
    case 1: return returnval + 3.*(d2bubble012+d2bubble013+d2bubble123);
    case 2: return returnval + 3.*(d2bubble012+d2bubble023+d2bubble123);
    case 3: return returnval + 3.*(d2bubble013+d2bubble023+d2bubble123);
    case 4: return returnval - 12.*(d2bubble012+d2bubble013);
    case 5: return returnval - 12.*(d2bubble012+d2bubble123);
    case 6: return returnval - 12.*(d2bubble012+d2bubble023);
    case 7: return returnval - 12.*(d2bubble013+d2bubble023);
    case 8: return returnval - 12.*(d2bubble013+d2bubble123);
    case 9: return returnval - 12.*(d2bubble023+d2bubble123);
    case 10: return 27.*d2bubble012;
    case 11: return 27.*d2bubble013;
    case 12: return 27.*d2bubble123;
    default: return 27.*d2bubble023;
    }
}

#endif

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_FE_SIMPLEX_LAGRANGE_H
