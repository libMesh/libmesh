// Kokkos FEEvaluator specializations for 2-D Lagrange elements.
//
// Covers TRI3, TRI6, QUAD4, QUAD8, QUAD9.
// Reference-element coordinate conventions (libMesh-compatible):
//   Tri:   xi >= 0, eta >= 0, xi+eta <= 1  (unit triangle)
//   Quad:  (xi, eta) in [-1,1]²

#ifndef LIBMESH_KOKKOS_FE_LAGRANGE_2D_H
#define LIBMESH_KOKKOS_FE_LAGRANGE_2D_H

#include "kokkos_fe_base.h"
#include "libmesh/fe_serendipity_lagrange.h"
#include "libmesh/fe_simplex_lagrange.h"
#include "libmesh/fe_tensor_product_lagrange.h"

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
    return libMesh::detail::fe_lagrange_tri3_shape(i, xi, eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_tri3_shape_deriv(i, 0),
                       libMesh::detail::fe_lagrange_tri3_shape_deriv(i, 1),
                       0.0);
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
    return libMesh::detail::fe_lagrange_tri6_shape(i, xi, eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_tri6_shape_deriv(i, 0, xi, eta),
                       libMesh::detail::fe_lagrange_tri6_shape_deriv(i, 1, xi, eta),
                       0.0);
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
    return libMesh::detail::fe_lagrange_quad4_shape(i, xi, eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_quad4_shape_deriv(i, 0, xi, eta),
                       libMesh::detail::fe_lagrange_quad4_shape_deriv(i, 1, xi, eta),
                       0.0);
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
    return libMesh::detail::fe_lagrange_quad8_shape(i, xi, eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_quad8_shape_deriv(i, 0, xi, eta),
                       libMesh::detail::fe_lagrange_quad8_shape_deriv(i, 1, xi, eta),
                       0.0);
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
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    return libMesh::detail::fe_lagrange_quad9_shape(i, xi, eta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_quad9_shape_deriv(i, 0, xi, eta),
                       libMesh::detail::fe_lagrange_quad9_shape_deriv(i, 1, xi, eta),
                       0.0);
  }
#endif
};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_LAGRANGE_2D_H
