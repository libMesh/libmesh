// Kokkos FEEvaluator specializations for 3-D Lagrange elements.
//
// Covers TET4, TET10, HEX8, HEX20, HEX27.
// Reference-element coordinate conventions (libMesh-compatible):
//   Tet: xi >= 0, eta >= 0, zeta >= 0, xi+eta+zeta <= 1  (unit tetrahedron)
//   Hex: (xi, eta, zeta) in [-1,1]³

#ifndef LIBMESH_KOKKOS_FE_LAGRANGE_3D_H
#define LIBMESH_KOKKOS_FE_LAGRANGE_3D_H

#include "kokkos_fe_base.h"
#include "libmesh/fe_serendipity_lagrange.h"
#include "libmesh/fe_simplex_lagrange.h"
#include "libmesh/fe_tensor_product_lagrange.h"

namespace libMesh::Kokkos
{

// ── TET4 (linear tetrahedron, 4 nodes) ───────────────────────────────────────
// Barycentric: z0=1-xi-eta-zeta, z1=xi, z2=eta, z3=zeta

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>
{
  static constexpr unsigned int n_dofs() { return 4; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return libMesh::detail::fe_lagrange_tet4_shape(i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    return make_vector(libMesh::detail::fe_lagrange_tet4_shape_deriv(i, 0),
                       libMesh::detail::fe_lagrange_tet4_shape_deriv(i, 1),
                       libMesh::detail::fe_lagrange_tet4_shape_deriv(i, 2));
  }
#endif
};

// ── TET10 (quadratic tetrahedron, 10 nodes) ───────────────────────────────────
// Barycentric: z0=1-xi-eta-zeta, z1=xi, z2=eta, z3=zeta
//   phi_0  = z0*(2*z0-1)
//   phi_1  = z1*(2*z1-1) = xi*(2*xi-1)
//   phi_2  = z2*(2*z2-1) = eta*(2*eta-1)
//   phi_3  = z3*(2*z3-1) = zeta*(2*zeta-1)
//   phi_4  = 4*z0*z1     = 4*(1-xi-eta-zeta)*xi
//   phi_5  = 4*z1*z2     = 4*xi*eta
//   phi_6  = 4*z2*z0     = 4*eta*(1-xi-eta-zeta)
//   phi_7  = 4*z0*z3     = 4*(1-xi-eta-zeta)*zeta
//   phi_8  = 4*z1*z3     = 4*xi*zeta
//   phi_9  = 4*z2*z3     = 4*eta*zeta

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>
{
  static constexpr unsigned int n_dofs() { return 10; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return libMesh::detail::fe_lagrange_tet10_shape(i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return make_vector(libMesh::detail::fe_lagrange_tet10_shape_deriv(i, 0, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_tet10_shape_deriv(i, 1, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_tet10_shape_deriv(i, 2, xi, eta, zeta));
  }
#endif
};

// ── HEX8 (trilinear hexahedron, 8 nodes) ─────────────────────────────────────
// Tensor product of three EDGE2 bases.
// Node ordering (same as libMesh):
//   0:(-1,-1,-1)  1:(+1,-1,-1)  2:(+1,+1,-1)  3:(-1,+1,-1)
//   4:(-1,-1,+1)  5:(+1,-1,+1)  6:(+1,+1,+1)  7:(-1,+1,+1)

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>
{
  static constexpr unsigned int n_dofs() { return 8; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return libMesh::detail::fe_lagrange_hex8_shape(i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return make_vector(libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 0, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 1, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex8_shape_deriv(i, 2, xi, eta, zeta));
  }
#endif
};

// ── HEX20 (serendipity hexahedron, 20 nodes) ─────────────────────────────────
// Corner nodes: phi = 0.125*(1+sx*xi)*(1+sy*eta)*(1+sz*zeta)*(sx*xi+sy*eta+sz*zeta-2)
// Node ordering follows libMesh (nodes 0-7 corners, 8-19 midside).

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>
{
  static constexpr unsigned int n_dofs() { return 20; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return libMesh::detail::fe_lagrange_hex20_shape(i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return make_vector(libMesh::detail::fe_lagrange_hex20_shape_deriv(i, 0, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex20_shape_deriv(i, 1, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex20_shape_deriv(i, 2, xi, eta, zeta));
  }
#endif
};

// ── HEX27 (triquadratic hexahedron, 27 nodes) ─────────────────────────────────
// Tensor product of three EDGE3 bases.
// Index tables (libMesh fe_lagrange_shape_3D.C):
//   i0[] = {0,1,1,0, 0,1,1,0, 2,1,2,0, 0,1,1,0, 2,1,2,0, 2,2,1,2,0,2,2}
//   i1[] = {0,0,1,1, 0,0,1,1, 0,2,1,2, 0,0,1,1, 0,2,1,2, 2,0,2,1,2,2,2}
//   i2[] = {0,0,0,0, 1,1,1,1, 0,0,0,0, 2,2,2,2, 1,1,1,1, 0,2,2,2,2,1,2}

template <>
struct FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>
{
  static constexpr unsigned int n_dofs() { return 27; }

#ifdef LIBMESH_HAVE_KOKKOS
  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return libMesh::detail::fe_lagrange_hex27_shape(i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    return make_vector(libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 0, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 1, xi, eta, zeta),
                       libMesh::detail::fe_lagrange_hex27_shape_deriv(i, 2, xi, eta, zeta));
  }
#endif
};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_LAGRANGE_3D_H
