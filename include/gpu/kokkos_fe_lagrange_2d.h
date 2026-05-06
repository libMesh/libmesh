// Kokkos FEEvaluator specializations for 2-D Lagrange elements.
//
// Covers TRI3, TRI6, QUAD4, QUAD8, QUAD9.
// Reference-element coordinate conventions (libMesh-compatible):
//   Tri:   xi >= 0, eta >= 0, xi+eta <= 1  (unit triangle)
//   Quad:  (xi, eta) in [-1,1]²

#ifndef LIBMESH_KOKKOS_FE_LAGRANGE_2D_H
#define LIBMESH_KOKKOS_FE_LAGRANGE_2D_H

#include "kokkos_fe_base.h"

namespace libMesh::Kokkos
{

// ── TRI3 (linear triangle, 3 nodes) ──────────────────────────────────────────
// Barycentric: zeta0 = 1-xi-eta,  zeta1 = xi,  zeta2 = eta

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>
{
  static constexpr unsigned int n_dofs() { return 3; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 1.0 - xi - eta;
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
      case 0: return make_vector(-1.0, -1.0, 0.0);
      case 1: return make_vector( 1.0,  0.0, 0.0);
      case 2: return make_vector( 0.0,  1.0, 0.0);
      default: return zero_vector();
    }
  }
#endif
};

// ── TRI6 (quadratic triangle, 6 nodes) ───────────────────────────────────────
// Barycentric: z0=1-xi-eta, z1=xi, z2=eta
//   phi_0 = z0*(2*z0-1) = (1-xi-eta)*(1-2*xi-2*eta)
//   phi_1 = z1*(2*z1-1) = xi*(2*xi-1)
//   phi_2 = z2*(2*z2-1) = eta*(2*eta-1)
//   phi_3 = 4*z0*z1     = 4*(1-xi-eta)*xi
//   phi_4 = 4*z1*z2     = 4*xi*eta
//   phi_5 = 4*z2*z0     = 4*eta*(1-xi-eta)

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>
{
  static constexpr unsigned int n_dofs() { return 6; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    const Real z0 = 1.0 - xi - eta;
    switch (i)
    {
      case 0: return z0 * (2.0 * z0 - 1.0);
      case 1: return xi * (2.0 * xi - 1.0);
      case 2: return eta * (2.0 * eta - 1.0);
      case 3: return 4.0 * z0 * xi;
      case 4: return 4.0 * xi * eta;
      case 5: return 4.0 * eta * z0;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return make_vector(4.0*xi + 4.0*eta - 3.0, 4.0*xi + 4.0*eta - 3.0, 0.0);
      case 1: return make_vector(4.0*xi - 1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 4.0*eta - 1.0, 0.0);
      case 3: return make_vector(4.0*(1.0 - 2.0*xi - eta), -4.0*xi, 0.0);
      case 4: return make_vector(4.0*eta, 4.0*xi, 0.0);
      case 5: return make_vector(-4.0*eta, 4.0*(1.0 - xi - 2.0*eta), 0.0);
      default: return zero_vector();
    }
  }
#endif
};

// ── QUAD4 (bilinear quadrilateral, 4 nodes) ───────────────────────────────────
// Tensor product of two EDGE2 bases. libMesh node ordering:
//   node 0: (-1,-1)   node 1: (+1,-1)
//   node 2: (+1,+1)   node 3: (-1,+1)

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>
{
  static constexpr unsigned int n_dofs() { return 4; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 0.25 * (1.0 - xi) * (1.0 - eta);
      case 1: return 0.25 * (1.0 + xi) * (1.0 - eta);
      case 2: return 0.25 * (1.0 + xi) * (1.0 + eta);
      case 3: return 0.25 * (1.0 - xi) * (1.0 + eta);
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return make_vector(-0.25*(1.0-eta), -0.25*(1.0-xi), 0.0);
      case 1: return make_vector( 0.25*(1.0-eta), -0.25*(1.0+xi), 0.0);
      case 2: return make_vector( 0.25*(1.0+eta),  0.25*(1.0+xi), 0.0);
      case 3: return make_vector(-0.25*(1.0+eta),  0.25*(1.0-xi), 0.0);
      default: return zero_vector();
    }
  }
#endif
};

// ── QUAD8 (serendipity quadrilateral, 8 nodes) ────────────────────────────────
// Node ordering:
//   0: (-1,-1)   1: (+1,-1)   2: (+1,+1)   3: (-1,+1)
//   4: ( 0,-1)   5: (+1, 0)   6: ( 0,+1)   7: (-1, 0)

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>
{
  static constexpr unsigned int n_dofs() { return 8; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 0.25 * (1.0-xi) * (1.0-eta) * (-1.0-xi-eta);
      case 1: return 0.25 * (1.0+xi) * (1.0-eta) * (-1.0+xi-eta);
      case 2: return 0.25 * (1.0+xi) * (1.0+eta) * (-1.0+xi+eta);
      case 3: return 0.25 * (1.0-xi) * (1.0+eta) * (-1.0-xi+eta);
      case 4: return 0.5  * (1.0-xi*xi) * (1.0-eta);
      case 5: return 0.5  * (1.0+xi)    * (1.0-eta*eta);
      case 6: return 0.5  * (1.0-xi*xi) * (1.0+eta);
      case 7: return 0.5  * (1.0-xi)    * (1.0-eta*eta);
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return make_vector(0.25*(1.0-eta)*(2.0*xi+eta),
                                 0.25*(1.0-xi)*(xi+2.0*eta),
                                 0.0);
      case 1: return make_vector(0.25*(1.0-eta)*(2.0*xi-eta),
                                 0.25*(1.0+xi)*(2.0*eta-xi),
                                 0.0);
      case 2: return make_vector(0.25*(1.0+eta)*(2.0*xi+eta),
                                 0.25*(1.0+xi)*(xi+2.0*eta),
                                 0.0);
      case 3: return make_vector(0.25*(1.0+eta)*(2.0*xi-eta),
                                 0.25*(1.0-xi)*(2.0*eta-xi),
                                 0.0);
      case 4: return make_vector(-xi*(1.0-eta), -0.5*(1.0-xi*xi), 0.0);
      case 5: return make_vector(0.5*(1.0-eta*eta), -eta*(1.0+xi), 0.0);
      case 6: return make_vector(-xi*(1.0+eta), 0.5*(1.0-xi*xi), 0.0);
      case 7: return make_vector(-0.5*(1.0-eta*eta), -eta*(1.0-xi), 0.0);
      default: return zero_vector();
    }
  }
#endif
};

// ── QUAD9 (biquadratic quadrilateral, 9 nodes) ────────────────────────────────
// Tensor product of two EDGE3 bases. libMesh node ordering:
//   i0[] = {0,1,1,0, 2,1,2,0, 2}
//   i1[] = {0,0,1,1, 0,2,1,2, 2}
//
// 1D basis (libMesh non-sequential ordering):
//   L_0(t) = 0.5*t*(t-1)   dL_0/dt = t - 0.5
//   L_1(t) = 0.5*t*(t+1)   dL_1/dt = t + 0.5
//   L_2(t) = 1 - t²        dL_2/dt = -2*t

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>
{
  static constexpr unsigned int n_dofs() { return 9; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real L(unsigned int k, Real t)
  {
    switch (k)
    {
      case 0: return 0.5 * t * (t - 1.0);
      case 1: return 0.5 * t * (t + 1.0);
      case 2: return 1.0 - t * t;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static Real dL(unsigned int k, Real t)
  {
    switch (k)
    {
      case 0: return t - 0.5;
      case 1: return t + 0.5;
      case 2: return -2.0 * t;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
    static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
    return L(i0[i], xi) * L(i1[i], eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    static const unsigned int i0[] = {0, 1, 1, 0, 2, 1, 2, 0, 2};
    static const unsigned int i1[] = {0, 0, 1, 1, 0, 2, 1, 2, 2};
    const Real dxi  = dL(i0[i], xi)  * L(i1[i], eta);
    const Real deta = L(i0[i], xi)   * dL(i1[i], eta);
    return make_vector(dxi, deta, 0.0);
  }
#endif
};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_LAGRANGE_2D_H
