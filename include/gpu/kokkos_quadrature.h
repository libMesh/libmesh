// Kokkos device-compatible Gauss quadrature rules.
//
// All evaluation functions are LIBMESH_DEVICE_INLINE — callable from both
// host and GPU device code.
//
// GaussLegendre1D:  1-D Gauss-Legendre on [-1,1], 1-7 point rules.
// GaussQuadrature:  Full quadrature dispatcher for all supported topologies.
//   - n_points(topo, order): number of quadrature points
//   - point(topo, order, qp): reference coordinate of qp-th point
//   - weight(topo, order, qp): weight of qp-th point
//
// Values match the libMesh QGauss implementation.

#ifndef LIBMESH_KOKKOS_QUADRATURE_H
#define LIBMESH_KOKKOS_QUADRATURE_H

#include "gpu/kokkos_scalar_types.h"
#include "libmesh/enum_elem_type.h"
#include <cmath>
#include <vector>

namespace libMesh::Kokkos
{

// ---------------------------------------------------------------------------
//  1-D Gauss-Legendre quadrature on [-1, 1]
// ---------------------------------------------------------------------------

struct GaussLegendre1D
{
  LIBMESH_DEVICE_INLINE static unsigned int n_points(unsigned int alg_order)
  {
    const unsigned int n = (alg_order + 2u) / 2u;
    return (n < 1u) ? 1u : (n > 7u ? 7u : n);
  }

  LIBMESH_DEVICE_INLINE static Real point(unsigned int n, unsigned int i)
  {
    switch (n)
    {
      case 1: return 0.0;
      case 2:
        switch (i)
        {
          case 0: return -5.7735026918962576450914878050196e-01;
          case 1: return  5.7735026918962576450914878050196e-01;
          default: return 0.0;
        }
      case 3:
        switch (i)
        {
          case 0: return -7.7459666924148337703585307995648e-01;
          case 1: return  0.0;
          case 2: return  7.7459666924148337703585307995648e-01;
          default: return 0.0;
        }
      case 4:
        switch (i)
        {
          case 0: return -8.6113631159405257522394648889281e-01;
          case 1: return -3.3998104358485626480266575910324e-01;
          case 2: return  3.3998104358485626480266575910324e-01;
          case 3: return  8.6113631159405257522394648889281e-01;
          default: return 0.0;
        }
      case 5:
        switch (i)
        {
          case 0: return -9.0617984593866399279762687829939e-01;
          case 1: return -5.3846931010568309103631442070021e-01;
          case 2: return  0.0;
          case 3: return  5.3846931010568309103631442070021e-01;
          case 4: return  9.0617984593866399279762687829939e-01;
          default: return 0.0;
        }
      case 6:
        switch (i)
        {
          case 0: return -9.3246951420315202781230155449399e-01;
          case 1: return -6.6120938646626451366139959501991e-01;
          case 2: return -2.3861918608319690863050172168071e-01;
          case 3: return  2.3861918608319690863050172168071e-01;
          case 4: return  6.6120938646626451366139959501991e-01;
          case 5: return  9.3246951420315202781230155449399e-01;
          default: return 0.0;
        }
      case 7:
        switch (i)
        {
          case 0: return -9.4910791234275852452618968404785e-01;
          case 1: return -7.4153118559939443986386477328079e-01;
          case 2: return -4.0584515137739716690660641207696e-01;
          case 3: return  0.0;
          case 4: return  4.0584515137739716690660641207696e-01;
          case 5: return  7.4153118559939443986386477328079e-01;
          case 6: return  9.4910791234275852452618968404785e-01;
          default: return 0.0;
        }
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static Real weight(unsigned int n, unsigned int i)
  {
    switch (n)
    {
      case 1: return 2.0;
      case 2: return 1.0;
      case 3:
        switch (i)
        {
          case 0: case 2: return 5.5555555555555555555555555555556e-01;
          case 1:         return 8.8888888888888888888888888888889e-01;
          default: return 0.0;
        }
      case 4:
        switch (i)
        {
          case 0: case 3: return 3.4785484513745385737306394922200e-01;
          case 1: case 2: return 6.5214515486254614262693605077800e-01;
          default: return 0.0;
        }
      case 5:
        switch (i)
        {
          case 0: case 4: return 2.3692688505618908751426404071992e-01;
          case 1: case 3: return 4.7862867049936646804129151483564e-01;
          case 2:         return 5.6888888888888888888888888888889e-01;
          default: return 0.0;
        }
      case 6:
        switch (i)
        {
          case 0: case 5: return 1.7132449237917034504029614217273e-01;
          case 1: case 4: return 3.6076157304813860756983351383772e-01;
          case 2: case 3: return 4.6791393457269104738987034398955e-01;
          default: return 0.0;
        }
      case 7:
        switch (i)
        {
          case 0: case 6: return 1.2948496616886969327061143267908e-01;
          case 1: case 5: return 2.7970539148927666790146777142378e-01;
          case 2: case 4: return 3.8183005050511894495036977548898e-01;
          case 3:         return 4.1795918367346938775510204081633e-01;
          default: return 0.0;
        }
      default: return 0.0;
    }
  }
};

// ---------------------------------------------------------------------------
//  GaussQuadrature — device-callable quadrature for all supported topologies
//
//  Coordinate conventions (same as libMesh):
//    EDGE:  xi in [-1,1]
//    QUAD:  (xi,eta) in [-1,1]^2, tensor product
//    HEX:   (xi,eta,zeta) in [-1,1]^3, tensor product
//    TRI:   (x,y) on unit triangle {(0,0),(1,0),(0,1)}
//    TET:   (x,y,z) on unit tet {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}
// ---------------------------------------------------------------------------

struct GaussQuadrature
{
  /// Number of quadrature points for a given topology and polynomial order.
  LIBMESH_DEVICE_INLINE static unsigned int
  n_points(libMesh::ElemType topo, unsigned int order)
  {
    switch (topo)
    {
      case libMesh::EDGE2: case libMesh::EDGE3:
        return GaussLegendre1D::n_points(order);

      case libMesh::QUAD4: case libMesh::QUAD8: case libMesh::QUAD9:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return n * n;
      }

      case libMesh::HEX8: case libMesh::HEX20: case libMesh::HEX27:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return n * n * n;
      }

      case libMesh::TRI3: case libMesh::TRI6:
        switch (order)
        {
          case 0: case 1: return 1;
          case 2: return 3;
          case 3: return 4;
          case 4: return 6;
          case 5: return 7;
          default: return 12;
        }

      case libMesh::TET4: case libMesh::TET10:
        switch (order)
        {
          case 0: case 1: return 1;
          case 2: return 4;
          case 3: return 5;
          case 4: return 11;
          case 5: return 14;
          default: return 24;
        }

      default: return 0;
    }
  }

  /// Reference coordinate of the qp-th quadrature point.
  LIBMESH_DEVICE_INLINE static RealVector
  point(libMesh::ElemType topo, unsigned int order, unsigned int qp)
  {
    switch (topo)
    {
      case libMesh::EDGE2: case libMesh::EDGE3:
        return make_vector(GaussLegendre1D::point(GaussLegendre1D::n_points(order), qp), 0.0, 0.0);

      case libMesh::QUAD4: case libMesh::QUAD8: case libMesh::QUAD9:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        const unsigned int i = qp % n;
        const unsigned int j = qp / n;
        return make_vector(GaussLegendre1D::point(n, i),
                           GaussLegendre1D::point(n, j), 0.0);
      }

      case libMesh::HEX8: case libMesh::HEX20: case libMesh::HEX27:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        const unsigned int i = qp % n;
        const unsigned int j = (qp / n) % n;
        const unsigned int k = qp / (n * n);
        return make_vector(GaussLegendre1D::point(n, i),
                           GaussLegendre1D::point(n, j),
                           GaussLegendre1D::point(n, k));
      }

      case libMesh::TRI3: case libMesh::TRI6:
        return tri_point(order, qp);

      case libMesh::TET4: case libMesh::TET10:
        return tet_point(order, qp);

      default: return zero_vector();
    }
  }

  /// Weight of the qp-th quadrature point.
  LIBMESH_DEVICE_INLINE static Real
  weight(libMesh::ElemType topo, unsigned int order, unsigned int qp)
  {
    switch (topo)
    {
      case libMesh::EDGE2: case libMesh::EDGE3:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return GaussLegendre1D::weight(n, qp);
      }

      case libMesh::QUAD4: case libMesh::QUAD8: case libMesh::QUAD9:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return GaussLegendre1D::weight(n, qp % n) *
               GaussLegendre1D::weight(n, qp / n);
      }

      case libMesh::HEX8: case libMesh::HEX20: case libMesh::HEX27:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return GaussLegendre1D::weight(n, qp % n) *
               GaussLegendre1D::weight(n, (qp / n) % n) *
               GaussLegendre1D::weight(n, qp / (n * n));
      }

      case libMesh::TRI3: case libMesh::TRI6:
        return tri_weight(order, qp);

      case libMesh::TET4: case libMesh::TET10:
        return tet_weight(order, qp);

      default: return 0.0;
    }
  }

private:
  // ── Triangle rules ────────────────────────────────────────────────────────

  LIBMESH_DEVICE_INLINE static RealVector
  tri_point(unsigned int order, unsigned int qp)
  {
    switch (order)
    {
      case 0: case 1:
        return make_vector(1.0 / 3.0, 1.0 / 3.0, 0.0);

      case 2:
        switch (qp)
        {
          case 0: return make_vector(2.0 / 3.0, 1.0 / 6.0, 0.0);
          case 1: return make_vector(1.0 / 6.0, 2.0 / 3.0, 0.0);
          case 2: return make_vector(1.0 / 6.0, 1.0 / 6.0, 0.0);
          default: return zero_vector();
        }

      case 3:
        switch (qp)
        {
          case 0: return make_vector(1.5505102572168219018e-01, 1.7855872826361642312e-01, 0.0);
          case 1: return make_vector(6.4494897427831780982e-01, 7.5031110222608118177e-02, 0.0);
          case 2: return make_vector(1.5505102572168219018e-01, 6.6639024601470138670e-01, 0.0);
          case 3: return make_vector(6.4494897427831780982e-01, 2.8001991549907407200e-01, 0.0);
          default: return zero_vector();
        }

      case 4:
      {
        constexpr Real a1 = 4.4594849091596488632e-01, b1 = 1.0 - 2.0 * a1;
        constexpr Real a2 = 9.1576213509770743460e-02, b2 = 1.0 - 2.0 * a2;
        switch (qp)
        {
          case 0: return make_vector(a1, a1, 0.0);
          case 1: return make_vector(a1, b1, 0.0);
          case 2: return make_vector(b1, a1, 0.0);
          case 3: return make_vector(a2, a2, 0.0);
          case 4: return make_vector(a2, b2, 0.0);
          case 5: return make_vector(b2, a2, 0.0);
          default: return zero_vector();
        }
      }

      case 5:
      {
        const Real sq15 = 3.872983346207417; // sqrt(15)
        const Real a1 = 2.0 / 7.0 + sq15 / 21.0;
        const Real a2 = 2.0 / 7.0 - sq15 / 21.0;
        const Real b1 = 1.0 - 2.0 * a1, b2 = 1.0 - 2.0 * a2;
        switch (qp)
        {
          case 0: return make_vector(1.0 / 3.0, 1.0 / 3.0, 0.0);
          case 1: return make_vector(a1, a1, 0.0);
          case 2: return make_vector(a1, b1, 0.0);
          case 3: return make_vector(b1, a1, 0.0);
          case 4: return make_vector(a2, a2, 0.0);
          case 5: return make_vector(a2, b2, 0.0);
          case 6: return make_vector(b2, a2, 0.0);
          default: return zero_vector();
        }
      }

      case 6:
      {
        constexpr Real a1 = 2.4928674517091042129163855310701908e-01;
        constexpr Real a2 = 6.3089014491502228340331602870819157e-02;
        constexpr Real a3 = 3.1035245103378440541660773395655215e-01;
        constexpr Real b1 = 1.0 - 2.0 * a1;
        constexpr Real b2 = 1.0 - 2.0 * a2;
        constexpr Real b3 = 6.3650249912139864723014259441204970e-01;
        constexpr Real c3 = 1.0 - a3 - b3;
        switch (qp)
        {
          case 0: return make_vector(a1, a1, 0.0);
          case 1: return make_vector(a1, b1, 0.0);
          case 2: return make_vector(b1, a1, 0.0);
          case 3: return make_vector(a2, a2, 0.0);
          case 4: return make_vector(a2, b2, 0.0);
          case 5: return make_vector(b2, a2, 0.0);
          case 6: return make_vector(a3, b3, 0.0);
          case 7: return make_vector(b3, a3, 0.0);
          case 8: return make_vector(a3, c3, 0.0);
          case 9: return make_vector(c3, a3, 0.0);
          case 10: return make_vector(b3, c3, 0.0);
          case 11: return make_vector(c3, b3, 0.0);
          default: return zero_vector();
        }
      }

      default: // order >= 7: 12-point Ro3-invariant rule
      {
        constexpr Real rd[4][2] = {
          {6.2382265094402118174e-02, 6.7517867073916085443e-02},
          {5.5225456656926611737e-02, 3.2150249385198182267e-01},
          {3.4324302945097146470e-02, 6.6094919618673565761e-01},
          {5.1584233435359177926e-01, 2.7771616697639178257e-01}
        };
        const unsigned int row = qp / 3;
        const unsigned int sub = qp % 3;
        if (row >= 4)
          return zero_vector();
        const Real z1 = rd[row][0], z2 = rd[row][1], z3 = 1.0 - z1 - z2;
        switch (sub)
        {
          case 0: return make_vector(z1, z2, 0.0);
          case 1: return make_vector(z3, z1, 0.0);
          case 2: return make_vector(z2, z3, 0.0);
          default: return zero_vector();
        }
      }
    }
  }

  LIBMESH_DEVICE_INLINE static Real
  tri_weight(unsigned int order, unsigned int qp)
  {
    switch (order)
    {
      case 0: case 1: return 0.5;
      case 2: return 1.0 / 6.0;
      case 3: return (qp % 2 == 0) ? 1.5902069087198858470e-01 : 9.0979309128011415303e-02;
      case 4: return (qp < 3) ? 1.1169079483900573285e-01 : 5.4975871827660933819e-02;
      case 5:
      {
        if (qp == 0)
          return 9.0 / 80.0;
        const Real sq15 = 3.872983346207417;
        return (qp <= 3) ? (31.0 / 480.0 + sq15 / 2400.0) : (31.0 / 480.0 - sq15 / 2400.0);
      }
      case 6:
      {
        if (qp <= 2)
          return 5.8393137863189683012644805692789721e-02;
        if (qp <= 5)
          return 2.5422453185103408460468404553434492e-02;
        return 4.1425537809186787596776728210221227e-02;
      }
      default:
      {
        constexpr Real wts[4] = {
          2.6517028157436251429e-02, 4.3881408714446055037e-02,
          2.8775042784981585738e-02, 6.7493187009802774463e-02
        };
        return (qp / 3 < 4) ? wts[qp / 3] : 0.0;
      }
    }
  }

  // ── Tetrahedral rules ─────────────────────────────────────────────────────

  LIBMESH_DEVICE_INLINE static RealVector
  tet_point(unsigned int order, unsigned int qp)
  {
    switch (order)
    {
      case 0: case 1:
        return make_vector(0.25, 0.25, 0.25);

      case 2:
      {
        const Real b = 0.25 * (1.0 - 1.0 / 2.2360679774997896964); // 1/sqrt(5)
        const Real a = 1.0 - 3.0 * b;
        switch (qp)
        {
          case 0: return make_vector(a, b, b);
          case 1: return make_vector(b, a, b);
          case 2: return make_vector(b, b, a);
          case 3: return make_vector(b, b, b);
          default: return zero_vector();
        }
      }

      case 3:
        switch (qp)
        {
          case 0: return make_vector(0.25, 0.25, 0.25);
          case 1: return make_vector(0.5, 1.0 / 6.0, 1.0 / 6.0);
          case 2: return make_vector(1.0 / 6.0, 0.5, 1.0 / 6.0);
          case 3: return make_vector(1.0 / 6.0, 1.0 / 6.0, 0.5);
          case 4: return make_vector(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0);
          default: return zero_vector();
        }

      case 4:
      {
        constexpr Real a1 = 2.5e-01;
        constexpr Real a2 = 7.85714285714285714e-01, b2 = 7.14285714285714285e-02;
        constexpr Real a3 = 3.99403576166799219e-01, b3 = 1.00596423833200785e-01;
        switch (qp)
        {
          case 0:  return make_vector(a1, a1, a1);
          case 1:  return make_vector(a2, b2, b2);
          case 2:  return make_vector(b2, a2, b2);
          case 3:  return make_vector(b2, b2, a2);
          case 4:  return make_vector(b2, b2, b2);
          case 5:  return make_vector(a3, a3, b3);
          case 6:  return make_vector(a3, b3, b3);
          case 7:  return make_vector(b3, b3, a3);
          case 8:  return make_vector(b3, a3, b3);
          case 9:  return make_vector(b3, a3, a3);
          case 10: return make_vector(a3, b3, a3);
          default: return zero_vector();
        }
      }

      case 5:
      {
        constexpr Real af[3] = {3.1088591926330060980e-01,
                                9.2735250310891226402e-02,
                                4.5503704125649649492e-02};
        if (qp < 8)
        {
          const unsigned int g = qp / 4;
          const unsigned int sub = qp % 4;
          const Real ag = af[g], bg = 1.0 - 3.0 * ag;
          switch (sub)
          {
            case 0: return make_vector(ag, ag, ag);
            case 1: return make_vector(ag, bg, ag);
            case 2: return make_vector(bg, ag, ag);
            case 3: return make_vector(ag, ag, bg);
            default: return zero_vector();
          }
        }
        else
        {
          const Real a2 = af[2], b2 = 0.5 * (1.0 - 2.0 * a2);
          switch (qp - 8)
          {
            case 0: return make_vector(b2, b2, a2);
            case 1: return make_vector(b2, a2, a2);
            case 2: return make_vector(a2, a2, b2);
            case 3: return make_vector(a2, b2, a2);
            case 4: return make_vector(b2, a2, b2);
            case 5: return make_vector(a2, b2, b2);
            default: return zero_vector();
          }
        }
      }

      default: // order >= 6: 24-point Keast rule
      {
        constexpr Real data[4][3] = {
          {3.56191386222544953e-01, 2.14602871259151684e-01, 0.0},
          {8.77978124396165982e-01, 4.06739585346113397e-02, 0.0},
          {3.29863295731730594e-02, 3.22337890142275646e-01, 0.0},
          {0.0, 0.0, 0.0} // 12-perm group handled separately
        };

        if (qp < 12)
        {
          // Three 4-permutation groups
          const unsigned int grp = qp / 4;
          const unsigned int sub = qp % 4;
          const Real a = data[grp][0], b = data[grp][1];
          switch (sub)
          {
            case 0: return make_vector(a, b, b);
            case 1: return make_vector(b, a, b);
            case 2: return make_vector(b, b, a);
            case 3: return make_vector(b, b, b);
            default: return zero_vector();
          }
        }
        else
        {
          // 12-permutation group
          constexpr Real a4 = 6.36610018750175299e-02;
          constexpr Real b4 = 2.69672331458315867e-01;
          constexpr Real c4 = 6.03005664791649076e-01;
          switch (qp - 12)
          {
            case 0:  return make_vector(a4, a4, b4);
            case 1:  return make_vector(a4, a4, c4);
            case 2:  return make_vector(b4, a4, a4);
            case 3:  return make_vector(c4, a4, a4);
            case 4:  return make_vector(a4, b4, a4);
            case 5:  return make_vector(a4, c4, a4);
            case 6:  return make_vector(a4, b4, c4);
            case 7:  return make_vector(a4, c4, b4);
            case 8:  return make_vector(b4, a4, c4);
            case 9:  return make_vector(b4, c4, a4);
            case 10: return make_vector(c4, a4, b4);
            case 11: return make_vector(c4, b4, a4);
            default: return zero_vector();
          }
        }
      }
    }
  }

  LIBMESH_DEVICE_INLINE static Real
  tet_weight(unsigned int order, unsigned int qp)
  {
    switch (order)
    {
      case 0: case 1: return 1.0 / 6.0;
      case 2: return 1.0 / 24.0;
      case 3: return (qp == 0) ? -2.0 / 15.0 : 0.075;
      case 4:
      {
        if (qp == 0)
          return -1.31555555555555556e-02;
        if (qp <= 4)
          return 7.62222222222222222e-03;
        return 2.48888888888888889e-02;
      }
      case 5:
      {
        constexpr Real wf[3] = {1.8781320953002641800e-02,
                                1.2248840519393658257e-02,
                                7.0910034628469110730e-03};
        if (qp < 4)
          return wf[0];
        if (qp < 8)
          return wf[1];
        return wf[2];
      }
      default:
      {
        constexpr Real wts[4] = {6.65379170969464506e-03,
                                 1.67953517588677620e-03,
                                 9.22619692394239843e-03,
                                 8.03571428571428248e-03};
        if (qp < 4)
          return wts[0];
        if (qp < 8)
          return wts[1];
        if (qp < 12)
          return wts[2];
        return wts[3];
      }
    }
  }
};

// ---------------------------------------------------------------------------
//  fill_quadrature — host-side convenience wrapper
//
//  Fills std::vectors using the device-callable GaussQuadrature functions.
// ---------------------------------------------------------------------------

inline void
fill_quadrature(libMesh::ElemType topo,
                unsigned int order,
                std::vector<RealVector> & qpts,
                std::vector<Real> & weights)
{
  const unsigned int nqp = GaussQuadrature::n_points(topo, order);
  qpts.resize(nqp);
  weights.resize(nqp);
  for (unsigned int q = 0; q < nqp; ++q)
  {
    qpts[q] = GaussQuadrature::point(topo, order, q);
    weights[q] = GaussQuadrature::weight(topo, order, q);
  }
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_QUADRATURE_H
