// Kokkos FEEvaluator specializations for MONOMIAL elements.
//
// MONOMIAL uses the complete total-degree polynomial space P_p.  Following
// libMesh's FE<Dim, MONOMIAL>, the basis is parameterised by spatial dimension,
// not element class — TRI and QUAD share MonomialImpl2D; TET/HEX/PRISM/PYRAMID
// share MonomialImpl3D.  This gives 3 x 6 = 18 impl specializations (dims 1/2/3,
// orders 0-5), then per-topology FEEvaluator delegating specializations wire each
// libMesh::ElemType to the matching impl.
//
// Basis ordering: graded-lex (total degree first, then lexicographic by
// decreasing xi exponent).  Matches libMesh::FE<Dim, MONOMIAL>::shape ordering.

#ifndef LIBMESH_KOKKOS_FE_MONOMIAL_H
#define LIBMESH_KOKKOS_FE_MONOMIAL_H

#include "kokkos_fe_base.h"
#include "libmesh/enum_elem_type.h"

namespace libMesh::Kokkos
{

// ═══════════════════════════════════════════════════════════════════════════
// MonomialImpl1D<N> — 1-D MONOMIAL basis, order N
// n_dofs = N + 1
// Basis: {1, xi, xi², xi³, ...}
// ═══════════════════════════════════════════════════════════════════════════

template <unsigned int N>
struct MonomialImpl1D;

template <>
struct MonomialImpl1D<0>
{
  static constexpr unsigned int n_dofs() { return 1; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return 1.0;
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return zero_vector();
  }
};

template <>
struct MonomialImpl1D<1>
{
  static constexpr unsigned int n_dofs() { return 2; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl1D<2>
{
  static constexpr unsigned int n_dofs() { return 3; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return xi * xi;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(2.0 * xi, 0.0, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl1D<3>
{
  static constexpr unsigned int n_dofs() { return 4; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return xi * xi;
      case 3: return xi * xi * xi;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(2.0 * xi, 0.0, 0.0);
      case 3: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl1D<4>
{
  static constexpr unsigned int n_dofs() { return 5; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return xi * xi;
      case 3: return xi * xi * xi;
      case 4: return xi * xi * xi * xi;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(2.0 * xi, 0.0, 0.0);
      case 3: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 4: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl1D<5>
{
  static constexpr unsigned int n_dofs() { return 6; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return xi * xi;
      case 3: return xi * xi * xi;
      case 4: return xi * xi * xi * xi;
      case 5: return xi * xi * xi * xi * xi;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(2.0 * xi, 0.0, 0.0);
      case 3: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 4: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      case 5: return make_vector(5.0 * xi * xi * xi * xi, 0.0, 0.0);
      default: return zero_vector();
    }
  }
};

// ═══════════════════════════════════════════════════════════════════════════
// MonomialImpl2D<N> — 2-D MONOMIAL basis, order N
// n_dofs = (N+1)(N+2)/2
// Graded-lex basis: {1, xi, eta, xi², xi·eta, eta², ...}
// Shared by TRI and QUAD element classes.
// ═══════════════════════════════════════════════════════════════════════════

template <unsigned int N>
struct MonomialImpl2D;

template <>
struct MonomialImpl2D<0>
{
  static constexpr unsigned int n_dofs() { return 1; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return 1.0;
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return zero_vector();
  }
};

template <>
struct MonomialImpl2D<1>
{
  static constexpr unsigned int n_dofs() { return 3; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return eta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 1.0, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl2D<2>
{
  static constexpr unsigned int n_dofs() { return 6; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return eta;
      case 3: return xi * xi;
      case 4: return xi * eta;
      case 5: return eta * eta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 1.0, 0.0);
      case 3: return make_vector(2.0 * xi, 0.0, 0.0);
      case 4: return make_vector(eta, xi, 0.0);
      case 5: return make_vector(0.0, 2.0 * eta, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl2D<3>
{
  static constexpr unsigned int n_dofs() { return 10; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return eta;
      case 3: return xi * xi;
      case 4: return xi * eta;
      case 5: return eta * eta;
      case 6: return xi * xi * xi;
      case 7: return xi * xi * eta;
      case 8: return xi * eta * eta;
      case 9: return eta * eta * eta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 1.0, 0.0);
      case 3: return make_vector(2.0 * xi, 0.0, 0.0);
      case 4: return make_vector(eta, xi, 0.0);
      case 5: return make_vector(0.0, 2.0 * eta, 0.0);
      case 6: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 7: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case 8: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case 9: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl2D<4>
{
  static constexpr unsigned int n_dofs() { return 15; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case  0: return 1.0;
      case  1: return xi;
      case  2: return eta;
      case  3: return xi * xi;
      case  4: return xi * eta;
      case  5: return eta * eta;
      case  6: return xi * xi * xi;
      case  7: return xi * xi * eta;
      case  8: return xi * eta * eta;
      case  9: return eta * eta * eta;
      case 10: return xi * xi * xi * xi;
      case 11: return xi * xi * xi * eta;
      case 12: return xi * xi * eta * eta;
      case 13: return xi * eta * eta * eta;
      case 14: return eta * eta * eta * eta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case  0: return zero_vector();
      case  1: return make_vector(1.0, 0.0, 0.0);
      case  2: return make_vector(0.0, 1.0, 0.0);
      case  3: return make_vector(2.0 * xi, 0.0, 0.0);
      case  4: return make_vector(eta, xi, 0.0);
      case  5: return make_vector(0.0, 2.0 * eta, 0.0);
      case  6: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case  7: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case  8: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case  9: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      case 10: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      case 11: return make_vector(3.0 * xi * xi * eta, xi * xi * xi, 0.0);
      case 12: return make_vector(2.0 * xi * eta * eta, 2.0 * xi * xi * eta, 0.0);
      case 13: return make_vector(eta * eta * eta, 3.0 * xi * eta * eta, 0.0);
      case 14: return make_vector(0.0, 4.0 * eta * eta * eta, 0.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl2D<5>
{
  static constexpr unsigned int n_dofs() { return 21; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case  0: return 1.0;
      case  1: return xi;
      case  2: return eta;
      case  3: return xi * xi;
      case  4: return xi * eta;
      case  5: return eta * eta;
      case  6: return xi * xi * xi;
      case  7: return xi * xi * eta;
      case  8: return xi * eta * eta;
      case  9: return eta * eta * eta;
      case 10: return xi * xi * xi * xi;
      case 11: return xi * xi * xi * eta;
      case 12: return xi * xi * eta * eta;
      case 13: return xi * eta * eta * eta;
      case 14: return eta * eta * eta * eta;
      case 15: return xi * xi * xi * xi * xi;
      case 16: return xi * xi * xi * xi * eta;
      case 17: return xi * xi * xi * eta * eta;
      case 18: return xi * xi * eta * eta * eta;
      case 19: return xi * eta * eta * eta * eta;
      case 20: return eta * eta * eta * eta * eta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case  0: return zero_vector();
      case  1: return make_vector(1.0, 0.0, 0.0);
      case  2: return make_vector(0.0, 1.0, 0.0);
      case  3: return make_vector(2.0 * xi, 0.0, 0.0);
      case  4: return make_vector(eta, xi, 0.0);
      case  5: return make_vector(0.0, 2.0 * eta, 0.0);
      case  6: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case  7: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case  8: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case  9: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      case 10: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      case 11: return make_vector(3.0 * xi * xi * eta, xi * xi * xi, 0.0);
      case 12: return make_vector(2.0 * xi * eta * eta, 2.0 * xi * xi * eta, 0.0);
      case 13: return make_vector(eta * eta * eta, 3.0 * xi * eta * eta, 0.0);
      case 14: return make_vector(0.0, 4.0 * eta * eta * eta, 0.0);
      case 15: return make_vector(5.0 * xi * xi * xi * xi, 0.0, 0.0);
      case 16: return make_vector(4.0 * xi * xi * xi * eta, xi * xi * xi * xi, 0.0);
      case 17: return make_vector(3.0 * xi * xi * eta * eta, 2.0 * xi * xi * xi * eta, 0.0);
      case 18: return make_vector(2.0 * xi * eta * eta * eta, 3.0 * xi * xi * eta * eta, 0.0);
      case 19: return make_vector(eta * eta * eta * eta, 4.0 * xi * eta * eta * eta, 0.0);
      case 20: return make_vector(0.0, 5.0 * eta * eta * eta * eta, 0.0);
      default: return zero_vector();
    }
  }
};

// ═══════════════════════════════════════════════════════════════════════════
// MonomialImpl3D<N> — 3-D MONOMIAL basis, order N
// n_dofs = (N+1)(N+2)(N+3)/6
// Basis ordering: graded-lex; for each total degree d, iterate c (zeta
// exponent) from 0 to d, then a (xi exponent) from d-c down to 0 (b=d-c-a).
// Shared by TET, HEX, PRISM, and PYRAMID element classes.
// ═══════════════════════════════════════════════════════════════════════════

template <unsigned int N>
struct MonomialImpl3D;

template <>
struct MonomialImpl3D<0>
{
  static constexpr unsigned int n_dofs() { return 1; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return 1.0;
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int /*i*/, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return zero_vector();
  }
};

template <>
struct MonomialImpl3D<1>
{
  static constexpr unsigned int n_dofs() { return 4; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return eta;
      case 3: return zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 1.0, 0.0);
      case 3: return make_vector(0.0, 0.0, 1.0);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl3D<2>
{
  static constexpr unsigned int n_dofs() { return 10; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0: return 1.0;
      case 1: return xi;
      case 2: return eta;
      case 3: return zeta;
      case 4: return xi * xi;
      case 5: return xi * eta;
      case 6: return eta * eta;
      case 7: return xi * zeta;
      case 8: return eta * zeta;
      case 9: return zeta * zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0: return zero_vector();
      case 1: return make_vector(1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 1.0, 0.0);
      case 3: return make_vector(0.0, 0.0, 1.0);
      case 4: return make_vector(2.0 * xi, 0.0, 0.0);
      case 5: return make_vector(eta, xi, 0.0);
      case 6: return make_vector(0.0, 2.0 * eta, 0.0);
      case 7: return make_vector(zeta, 0.0, xi);
      case 8: return make_vector(0.0, zeta, eta);
      case 9: return make_vector(0.0, 0.0, 2.0 * zeta);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl3D<3>
{
  static constexpr unsigned int n_dofs() { return 20; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return 1.0;
      case  1: return xi;
      case  2: return eta;
      case  3: return zeta;
      case  4: return xi * xi;
      case  5: return xi * eta;
      case  6: return eta * eta;
      case  7: return xi * zeta;
      case  8: return eta * zeta;
      case  9: return zeta * zeta;
      case 10: return xi * xi * xi;
      case 11: return xi * xi * eta;
      case 12: return xi * eta * eta;
      case 13: return eta * eta * eta;
      case 14: return xi * xi * zeta;
      case 15: return xi * eta * zeta;
      case 16: return eta * eta * zeta;
      case 17: return xi * zeta * zeta;
      case 18: return eta * zeta * zeta;
      case 19: return zeta * zeta * zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return zero_vector();
      case  1: return make_vector(1.0, 0.0, 0.0);
      case  2: return make_vector(0.0, 1.0, 0.0);
      case  3: return make_vector(0.0, 0.0, 1.0);
      case  4: return make_vector(2.0 * xi, 0.0, 0.0);
      case  5: return make_vector(eta, xi, 0.0);
      case  6: return make_vector(0.0, 2.0 * eta, 0.0);
      case  7: return make_vector(zeta, 0.0, xi);
      case  8: return make_vector(0.0, zeta, eta);
      case  9: return make_vector(0.0, 0.0, 2.0 * zeta);
      case 10: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 11: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case 12: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case 13: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      case 14: return make_vector(2.0 * xi * zeta, 0.0, xi * xi);
      case 15: return make_vector(eta * zeta, xi * zeta, xi * eta);
      case 16: return make_vector(0.0, 2.0 * eta * zeta, eta * eta);
      case 17: return make_vector(zeta * zeta, 0.0, 2.0 * xi * zeta);
      case 18: return make_vector(0.0, zeta * zeta, 2.0 * eta * zeta);
      case 19: return make_vector(0.0, 0.0, 3.0 * zeta * zeta);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl3D<4>
{
  static constexpr unsigned int n_dofs() { return 35; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return 1.0;
      case  1: return xi;
      case  2: return eta;
      case  3: return zeta;
      case  4: return xi * xi;
      case  5: return xi * eta;
      case  6: return eta * eta;
      case  7: return xi * zeta;
      case  8: return eta * zeta;
      case  9: return zeta * zeta;
      case 10: return xi * xi * xi;
      case 11: return xi * xi * eta;
      case 12: return xi * eta * eta;
      case 13: return eta * eta * eta;
      case 14: return xi * xi * zeta;
      case 15: return xi * eta * zeta;
      case 16: return eta * eta * zeta;
      case 17: return xi * zeta * zeta;
      case 18: return eta * zeta * zeta;
      case 19: return zeta * zeta * zeta;
      case 20: return xi * xi * xi * xi;
      case 21: return xi * xi * xi * eta;
      case 22: return xi * xi * eta * eta;
      case 23: return xi * eta * eta * eta;
      case 24: return eta * eta * eta * eta;
      case 25: return xi * xi * xi * zeta;
      case 26: return xi * xi * eta * zeta;
      case 27: return xi * eta * eta * zeta;
      case 28: return eta * eta * eta * zeta;
      case 29: return xi * xi * zeta * zeta;
      case 30: return xi * eta * zeta * zeta;
      case 31: return eta * eta * zeta * zeta;
      case 32: return xi * zeta * zeta * zeta;
      case 33: return eta * zeta * zeta * zeta;
      case 34: return zeta * zeta * zeta * zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return zero_vector();
      case  1: return make_vector(1.0, 0.0, 0.0);
      case  2: return make_vector(0.0, 1.0, 0.0);
      case  3: return make_vector(0.0, 0.0, 1.0);
      case  4: return make_vector(2.0 * xi, 0.0, 0.0);
      case  5: return make_vector(eta, xi, 0.0);
      case  6: return make_vector(0.0, 2.0 * eta, 0.0);
      case  7: return make_vector(zeta, 0.0, xi);
      case  8: return make_vector(0.0, zeta, eta);
      case  9: return make_vector(0.0, 0.0, 2.0 * zeta);
      case 10: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 11: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case 12: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case 13: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      case 14: return make_vector(2.0 * xi * zeta, 0.0, xi * xi);
      case 15: return make_vector(eta * zeta, xi * zeta, xi * eta);
      case 16: return make_vector(0.0, 2.0 * eta * zeta, eta * eta);
      case 17: return make_vector(zeta * zeta, 0.0, 2.0 * xi * zeta);
      case 18: return make_vector(0.0, zeta * zeta, 2.0 * eta * zeta);
      case 19: return make_vector(0.0, 0.0, 3.0 * zeta * zeta);
      case 20: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      case 21: return make_vector(3.0 * xi * xi * eta, xi * xi * xi, 0.0);
      case 22: return make_vector(2.0 * xi * eta * eta, 2.0 * xi * xi * eta, 0.0);
      case 23: return make_vector(eta * eta * eta, 3.0 * xi * eta * eta, 0.0);
      case 24: return make_vector(0.0, 4.0 * eta * eta * eta, 0.0);
      case 25: return make_vector(3.0 * xi * xi * zeta, 0.0, xi * xi * xi);
      case 26: return make_vector(2.0 * xi * eta * zeta, xi * xi * zeta, xi * xi * eta);
      case 27: return make_vector(eta * eta * zeta, 2.0 * xi * eta * zeta, xi * eta * eta);
      case 28: return make_vector(0.0, 3.0 * eta * eta * zeta, eta * eta * eta);
      case 29: return make_vector(2.0 * xi * zeta * zeta, 0.0, 2.0 * xi * xi * zeta);
      case 30: return make_vector(eta * zeta * zeta, xi * zeta * zeta, 2.0 * xi * eta * zeta);
      case 31: return make_vector(0.0, 2.0 * eta * zeta * zeta, 2.0 * eta * eta * zeta);
      case 32: return make_vector(zeta * zeta * zeta, 0.0, 3.0 * xi * zeta * zeta);
      case 33: return make_vector(0.0, zeta * zeta * zeta, 3.0 * eta * zeta * zeta);
      case 34: return make_vector(0.0, 0.0, 4.0 * zeta * zeta * zeta);
      default: return zero_vector();
    }
  }
};

template <>
struct MonomialImpl3D<5>
{
  static constexpr unsigned int n_dofs() { return 56; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return 1.0;
      case  1: return xi;
      case  2: return eta;
      case  3: return zeta;
      case  4: return xi*xi;
      case  5: return xi*eta;
      case  6: return eta*eta;
      case  7: return xi*zeta;
      case  8: return eta*zeta;
      case  9: return zeta*zeta;
      case 10: return xi*xi*xi;
      case 11: return xi*xi*eta;
      case 12: return xi*eta*eta;
      case 13: return eta*eta*eta;
      case 14: return xi*xi*zeta;
      case 15: return xi*eta*zeta;
      case 16: return eta*eta*zeta;
      case 17: return xi*zeta*zeta;
      case 18: return eta*zeta*zeta;
      case 19: return zeta*zeta*zeta;
      case 20: return xi*xi*xi*xi;
      case 21: return xi*xi*xi*eta;
      case 22: return xi*xi*eta*eta;
      case 23: return xi*eta*eta*eta;
      case 24: return eta*eta*eta*eta;
      case 25: return xi*xi*xi*zeta;
      case 26: return xi*xi*eta*zeta;
      case 27: return xi*eta*eta*zeta;
      case 28: return eta*eta*eta*zeta;
      case 29: return xi*xi*zeta*zeta;
      case 30: return xi*eta*zeta*zeta;
      case 31: return eta*eta*zeta*zeta;
      case 32: return xi*zeta*zeta*zeta;
      case 33: return eta*zeta*zeta*zeta;
      case 34: return zeta*zeta*zeta*zeta;
      case 35: return xi*xi*xi*xi*xi;
      case 36: return xi*xi*xi*xi*eta;
      case 37: return xi*xi*xi*eta*eta;
      case 38: return xi*xi*eta*eta*eta;
      case 39: return xi*eta*eta*eta*eta;
      case 40: return eta*eta*eta*eta*eta;
      case 41: return xi*xi*xi*xi*zeta;
      case 42: return xi*xi*xi*eta*zeta;
      case 43: return xi*xi*eta*eta*zeta;
      case 44: return xi*eta*eta*eta*zeta;
      case 45: return eta*eta*eta*eta*zeta;
      case 46: return xi*xi*xi*zeta*zeta;
      case 47: return xi*xi*eta*zeta*zeta;
      case 48: return xi*eta*eta*zeta*zeta;
      case 49: return eta*eta*eta*zeta*zeta;
      case 50: return xi*xi*zeta*zeta*zeta;
      case 51: return xi*eta*zeta*zeta*zeta;
      case 52: return eta*eta*zeta*zeta*zeta;
      case 53: return xi*zeta*zeta*zeta*zeta;
      case 54: return eta*zeta*zeta*zeta*zeta;
      case 55: return zeta*zeta*zeta*zeta*zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case  0: return zero_vector();
      case  1: return make_vector(1.0, 0.0, 0.0);
      case  2: return make_vector(0.0, 1.0, 0.0);
      case  3: return make_vector(0.0, 0.0, 1.0);
      case  4: return make_vector(2.0 * xi, 0.0, 0.0);
      case  5: return make_vector(eta, xi, 0.0);
      case  6: return make_vector(0.0, 2.0 * eta, 0.0);
      case  7: return make_vector(zeta, 0.0, xi);
      case  8: return make_vector(0.0, zeta, eta);
      case  9: return make_vector(0.0, 0.0, 2.0 * zeta);
      case 10: return make_vector(3.0 * xi * xi, 0.0, 0.0);
      case 11: return make_vector(2.0 * xi * eta, xi * xi, 0.0);
      case 12: return make_vector(eta * eta, 2.0 * xi * eta, 0.0);
      case 13: return make_vector(0.0, 3.0 * eta * eta, 0.0);
      case 14: return make_vector(2.0 * xi * zeta, 0.0, xi * xi);
      case 15: return make_vector(eta * zeta, xi * zeta, xi * eta);
      case 16: return make_vector(0.0, 2.0 * eta * zeta, eta * eta);
      case 17: return make_vector(zeta * zeta, 0.0, 2.0 * xi * zeta);
      case 18: return make_vector(0.0, zeta * zeta, 2.0 * eta * zeta);
      case 19: return make_vector(0.0, 0.0, 3.0 * zeta * zeta);
      case 20: return make_vector(4.0 * xi * xi * xi, 0.0, 0.0);
      case 21: return make_vector(3.0 * xi * xi * eta, xi * xi * xi, 0.0);
      case 22: return make_vector(2.0 * xi * eta * eta, 2.0 * xi * xi * eta, 0.0);
      case 23: return make_vector(eta * eta * eta, 3.0 * xi * eta * eta, 0.0);
      case 24: return make_vector(0.0, 4.0 * eta * eta * eta, 0.0);
      case 25: return make_vector(3.0 * xi * xi * zeta, 0.0, xi * xi * xi);
      case 26: return make_vector(2.0 * xi * eta * zeta, xi * xi * zeta, xi * xi * eta);
      case 27: return make_vector(eta * eta * zeta, 2.0 * xi * eta * zeta, xi * eta * eta);
      case 28: return make_vector(0.0, 3.0 * eta * eta * zeta, eta * eta * eta);
      case 29: return make_vector(2.0 * xi * zeta * zeta, 0.0, 2.0 * xi * xi * zeta);
      case 30: return make_vector(eta * zeta * zeta, xi * zeta * zeta, 2.0 * xi * eta * zeta);
      case 31: return make_vector(0.0, 2.0 * eta * zeta * zeta, 2.0 * eta * eta * zeta);
      case 32: return make_vector(zeta * zeta * zeta, 0.0, 3.0 * xi * zeta * zeta);
      case 33: return make_vector(0.0, zeta * zeta * zeta, 3.0 * eta * zeta * zeta);
      case 34: return make_vector(0.0, 0.0, 4.0 * zeta * zeta * zeta);
      case 35: return make_vector(5.0 * xi * xi * xi * xi, 0.0, 0.0);
      case 36: return make_vector(4.0 * xi * xi * xi * eta, xi * xi * xi * xi, 0.0);
      case 37: return make_vector(3.0 * xi * xi * eta * eta, 2.0 * xi * xi * xi * eta, 0.0);
      case 38: return make_vector(2.0 * xi * eta * eta * eta, 3.0 * xi * xi * eta * eta, 0.0);
      case 39: return make_vector(eta * eta * eta * eta, 4.0 * xi * eta * eta * eta, 0.0);
      case 40: return make_vector(0.0, 5.0 * eta * eta * eta * eta, 0.0);
      case 41: return make_vector(4.0 * xi * xi * xi * zeta, 0.0, xi * xi * xi * xi);
      case 42: return make_vector(3.0 * xi * xi * eta * zeta, xi * xi * xi * zeta, xi * xi * xi * eta);
      case 43: return make_vector(2.0 * xi * eta * eta * zeta, 2.0 * xi * xi * eta * zeta, xi * xi * eta * eta);
      case 44: return make_vector(eta * eta * eta * zeta, 3.0 * xi * eta * eta * zeta, xi * eta * eta * eta);
      case 45: return make_vector(0.0, 4.0 * eta * eta * eta * zeta, eta * eta * eta * eta);
      case 46: return make_vector(3.0 * xi * xi * zeta * zeta, 0.0, 2.0 * xi * xi * xi * zeta);
      case 47: return make_vector(2.0 * xi * eta * zeta * zeta, xi * xi * zeta * zeta, 2.0 * xi * xi * eta * zeta);
      case 48: return make_vector(eta * eta * zeta * zeta, 2.0 * xi * eta * zeta * zeta, 2.0 * xi * eta * eta * zeta);
      case 49: return make_vector(0.0, 3.0 * eta * eta * zeta * zeta, 2.0 * eta * eta * eta * zeta);
      case 50: return make_vector(2.0 * xi * zeta * zeta * zeta, 0.0, 3.0 * xi * xi * zeta * zeta);
      case 51: return make_vector(eta * zeta * zeta * zeta, xi * zeta * zeta * zeta, 3.0 * xi * eta * zeta * zeta);
      case 52: return make_vector(0.0, 2.0 * eta * zeta * zeta * zeta, 3.0 * eta * eta * zeta * zeta);
      case 53: return make_vector(zeta * zeta * zeta * zeta, 0.0, 4.0 * xi * zeta * zeta * zeta);
      case 54: return make_vector(0.0, zeta * zeta * zeta * zeta, 4.0 * eta * zeta * zeta * zeta);
      case 55: return make_vector(0.0, 0.0, 5.0 * zeta * zeta * zeta * zeta);
      default: return zero_vector();
    }
  }
};

// ═══════════════════════════════════════════════════════════════════════════
// Per-topology FEEvaluator delegating specializations
//
// Each partial specialization fixes family=MONOMIAL and elem_type, leaving the
// polynomial Order as a template parameter, then inherits the matching impl.
// ═══════════════════════════════════════════════════════════════════════════

// ── 1-D ──────────────────────────────────────────────────────────────────────

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::EDGE2, N> : MonomialImpl1D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::EDGE3, N> : MonomialImpl1D<N> {};

// ── 2-D ──────────────────────────────────────────────────────────────────────

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI3, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI6, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI7, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD4, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD8, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD9, N> : MonomialImpl2D<N> {};

// ── 3-D ──────────────────────────────────────────────────────────────────────

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET4, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET10, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET14, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX8, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX20, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX27, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM6, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM15, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM18, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID5, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID13, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID14, N> : MonomialImpl3D<N> {};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_MONOMIAL_H
