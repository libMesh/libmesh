// Kokkos FEEvaluator specializations for 3-D Lagrange elements.
//
// Covers TET4, TET10, HEX8, HEX20, HEX27.
// Reference-element coordinate conventions (libMesh-compatible):
//   Tet: xi >= 0, eta >= 0, zeta >= 0, xi+eta+zeta <= 1  (unit tetrahedron)
//   Hex: (xi, eta, zeta) in [-1,1]³

#ifndef LIBMESH_KOKKOS_FE_LAGRANGE_3D_H
#define LIBMESH_KOKKOS_FE_LAGRANGE_3D_H

#include "gpu/kokkos_fe_base.h"

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
    switch (i)
    {
      case 0: return 1.0 - xi - eta - zeta;
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
      case 0: return make_vector(-1.0, -1.0, -1.0);
      case 1: return make_vector( 1.0,  0.0,  0.0);
      case 2: return make_vector( 0.0,  1.0,  0.0);
      case 3: return make_vector( 0.0,  0.0,  1.0);
      default: return zero_vector();
    }
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
    const Real z0 = 1.0 - xi - eta - zeta;
    switch (i)
    {
      case 0: return z0  * (2.0*z0   - 1.0);
      case 1: return xi  * (2.0*xi   - 1.0);
      case 2: return eta * (2.0*eta  - 1.0);
      case 3: return zeta* (2.0*zeta - 1.0);
      case 4: return 4.0 * z0 * xi;
      case 5: return 4.0 * xi * eta;
      case 6: return 4.0 * eta * z0;
      case 7: return 4.0 * z0 * zeta;
      case 8: return 4.0 * xi * zeta;
      case 9: return 4.0 * eta * zeta;
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0:
      {
        const Real v = 4.0*(xi + eta + zeta) - 3.0;
        return make_vector(v, v, v);
      }
      case 1: return make_vector(4.0*xi - 1.0, 0.0, 0.0);
      case 2: return make_vector(0.0, 4.0*eta - 1.0, 0.0);
      case 3: return make_vector(0.0, 0.0, 4.0*zeta - 1.0);
      case 4: return make_vector( 4.0*(1.0-2.0*xi-eta-zeta), -4.0*xi, -4.0*xi);
      case 5: return make_vector( 4.0*eta, 4.0*xi, 0.0);
      case 6: return make_vector(-4.0*eta, 4.0*(1.0-xi-2.0*eta-zeta), -4.0*eta);
      case 7: return make_vector(-4.0*zeta, -4.0*zeta, 4.0*(1.0-xi-eta-2.0*zeta));
      case 8: return make_vector(4.0*zeta, 0.0, 4.0*xi);
      case 9: return make_vector(0.0, 4.0*zeta, 4.0*eta);
      default: return zero_vector();
    }
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
    switch (i)
    {
      case 0: return 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta);
      case 1: return 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta);
      case 2: return 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta);
      case 3: return 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta);
      case 4: return 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta);
      case 5: return 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta);
      case 6: return 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta);
      case 7: return 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta);
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0: return make_vector(-0.125*(1.0-eta)*(1.0-zeta),
                                 -0.125*(1.0-xi) *(1.0-zeta),
                                 -0.125*(1.0-xi) *(1.0-eta));
      case 1: return make_vector( 0.125*(1.0-eta)*(1.0-zeta),
                                 -0.125*(1.0+xi) *(1.0-zeta),
                                 -0.125*(1.0+xi) *(1.0-eta));
      case 2: return make_vector( 0.125*(1.0+eta)*(1.0-zeta),
                                  0.125*(1.0+xi) *(1.0-zeta),
                                 -0.125*(1.0+xi) *(1.0+eta));
      case 3: return make_vector(-0.125*(1.0+eta)*(1.0-zeta),
                                  0.125*(1.0-xi) *(1.0-zeta),
                                 -0.125*(1.0-xi) *(1.0+eta));
      case 4: return make_vector(-0.125*(1.0-eta)*(1.0+zeta),
                                 -0.125*(1.0-xi) *(1.0+zeta),
                                  0.125*(1.0-xi) *(1.0-eta));
      case 5: return make_vector( 0.125*(1.0-eta)*(1.0+zeta),
                                 -0.125*(1.0+xi) *(1.0+zeta),
                                  0.125*(1.0+xi) *(1.0-eta));
      case 6: return make_vector( 0.125*(1.0+eta)*(1.0+zeta),
                                  0.125*(1.0+xi) *(1.0+zeta),
                                  0.125*(1.0+xi) *(1.0+eta));
      case 7: return make_vector(-0.125*(1.0+eta)*(1.0+zeta),
                                  0.125*(1.0-xi) *(1.0+zeta),
                                  0.125*(1.0-xi) *(1.0+eta));
      default: return zero_vector();
    }
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
    switch (i)
    {
      case 0:  return 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta)*(-xi-eta-zeta-2.0);
      case 1:  return 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta)*( xi-eta-zeta-2.0);
      case 2:  return 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta)*( xi+eta-zeta-2.0);
      case 3:  return 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta)*(-xi+eta-zeta-2.0);
      case 4:  return 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta)*(-xi-eta+zeta-2.0);
      case 5:  return 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta)*( xi-eta+zeta-2.0);
      case 6:  return 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta)*( xi+eta+zeta-2.0);
      case 7:  return 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta)*(-xi+eta+zeta-2.0);
      case 8:  return 0.25*(1.0-xi*xi)*(1.0-eta)*(1.0-zeta);
      case 10: return 0.25*(1.0-xi*xi)*(1.0+eta)*(1.0-zeta);
      case 16: return 0.25*(1.0-xi*xi)*(1.0-eta)*(1.0+zeta);
      case 18: return 0.25*(1.0-xi*xi)*(1.0+eta)*(1.0+zeta);
      case 9:  return 0.25*(1.0+xi)*(1.0-eta*eta)*(1.0-zeta);
      case 11: return 0.25*(1.0-xi)*(1.0-eta*eta)*(1.0-zeta);
      case 17: return 0.25*(1.0+xi)*(1.0-eta*eta)*(1.0+zeta);
      case 19: return 0.25*(1.0-xi)*(1.0-eta*eta)*(1.0+zeta);
      case 12: return 0.25*(1.0-xi)*(1.0-eta)*(1.0-zeta*zeta);
      case 13: return 0.25*(1.0+xi)*(1.0-eta)*(1.0-zeta*zeta);
      case 14: return 0.25*(1.0+xi)*(1.0+eta)*(1.0-zeta*zeta);
      case 15: return 0.25*(1.0-xi)*(1.0+eta)*(1.0-zeta*zeta);
      default: return 0.0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    switch (i)
    {
      case 0: return make_vector(
        -0.125*(1.0-eta)*(1.0-zeta)*(-2.0*xi-eta-zeta-1.0),
        -0.125*(1.0-xi) *(1.0-zeta)*(-xi-2.0*eta-zeta-1.0),
        -0.125*(1.0-xi) *(1.0-eta) *(-xi-eta-2.0*zeta-1.0));
      case 1: return make_vector(
         0.125*(1.0-eta)*(1.0-zeta)*(2.0*xi-eta-zeta-1.0),
        -0.125*(1.0+xi) *(1.0-zeta)*(xi-2.0*eta-zeta-1.0),
        -0.125*(1.0+xi) *(1.0-eta) *(xi-eta-2.0*zeta-1.0));
      case 2: return make_vector(
         0.125*(1.0+eta)*(1.0-zeta)*(2.0*xi+eta-zeta-1.0),
         0.125*(1.0+xi) *(1.0-zeta)*(xi+2.0*eta-zeta-1.0),
        -0.125*(1.0+xi) *(1.0+eta) *(xi+eta-2.0*zeta-1.0));
      case 3: return make_vector(
        -0.125*(1.0+eta)*(1.0-zeta)*(-2.0*xi+eta-zeta-1.0),
         0.125*(1.0-xi) *(1.0-zeta)*(-xi+2.0*eta-zeta-1.0),
        -0.125*(1.0-xi) *(1.0+eta) *(-xi+eta-2.0*zeta-1.0));
      case 4: return make_vector(
        -0.125*(1.0-eta)*(1.0+zeta)*(-2.0*xi-eta+zeta-1.0),
        -0.125*(1.0-xi) *(1.0+zeta)*(-xi-2.0*eta+zeta-1.0),
         0.125*(1.0-xi) *(1.0-eta) *(-xi-eta+2.0*zeta-1.0));
      case 5: return make_vector(
         0.125*(1.0-eta)*(1.0+zeta)*(2.0*xi-eta+zeta-1.0),
        -0.125*(1.0+xi) *(1.0+zeta)*(xi-2.0*eta+zeta-1.0),
         0.125*(1.0+xi) *(1.0-eta) *(xi-eta+2.0*zeta-1.0));
      case 6: return make_vector(
         0.125*(1.0+eta)*(1.0+zeta)*(2.0*xi+eta+zeta-1.0),
         0.125*(1.0+xi) *(1.0+zeta)*(xi+2.0*eta+zeta-1.0),
         0.125*(1.0+xi) *(1.0+eta) *(xi+eta+2.0*zeta-1.0));
      case 7: return make_vector(
        -0.125*(1.0+eta)*(1.0+zeta)*(-2.0*xi+eta+zeta-1.0),
         0.125*(1.0-xi) *(1.0+zeta)*(-xi+2.0*eta+zeta-1.0),
         0.125*(1.0-xi) *(1.0+eta) *(-xi+eta+2.0*zeta-1.0));
      case 8:  return make_vector(-0.5*xi*(1.0-eta)*(1.0-zeta),
                            -0.25*(1.0-xi*xi)*(1.0-zeta),
                            -0.25*(1.0-xi*xi)*(1.0-eta));
      case 10: return make_vector(-0.5*xi*(1.0+eta)*(1.0-zeta),
                             0.25*(1.0-xi*xi)*(1.0-zeta),
                            -0.25*(1.0-xi*xi)*(1.0+eta));
      case 16: return make_vector(-0.5*xi*(1.0-eta)*(1.0+zeta),
                            -0.25*(1.0-xi*xi)*(1.0+zeta),
                             0.25*(1.0-xi*xi)*(1.0-eta));
      case 18: return make_vector(-0.5*xi*(1.0+eta)*(1.0+zeta),
                             0.25*(1.0-xi*xi)*(1.0+zeta),
                             0.25*(1.0-xi*xi)*(1.0+eta));
      case 9:  return make_vector( 0.25*(1.0-eta*eta)*(1.0-zeta),
                            -0.5*eta*(1.0+xi)*(1.0-zeta),
                            -0.25*(1.0+xi)*(1.0-eta*eta));
      case 11: return make_vector(-0.25*(1.0-eta*eta)*(1.0-zeta),
                            -0.5*eta*(1.0-xi)*(1.0-zeta),
                            -0.25*(1.0-xi)*(1.0-eta*eta));
      case 17: return make_vector( 0.25*(1.0-eta*eta)*(1.0+zeta),
                            -0.5*eta*(1.0+xi)*(1.0+zeta),
                             0.25*(1.0+xi)*(1.0-eta*eta));
      case 19: return make_vector(-0.25*(1.0-eta*eta)*(1.0+zeta),
                            -0.5*eta*(1.0-xi)*(1.0+zeta),
                             0.25*(1.0-xi)*(1.0-eta*eta));
      case 12: return make_vector(-0.25*(1.0-eta)*(1.0-zeta*zeta),
                            -0.25*(1.0-xi)*(1.0-zeta*zeta),
                            -0.5*zeta*(1.0-xi)*(1.0-eta));
      case 13: return make_vector( 0.25*(1.0-eta)*(1.0-zeta*zeta),
                            -0.25*(1.0+xi)*(1.0-zeta*zeta),
                            -0.5*zeta*(1.0+xi)*(1.0-eta));
      case 14: return make_vector( 0.25*(1.0+eta)*(1.0-zeta*zeta),
                             0.25*(1.0+xi)*(1.0-zeta*zeta),
                            -0.5*zeta*(1.0+xi)*(1.0+eta));
      case 15: return make_vector(-0.25*(1.0+eta)*(1.0-zeta*zeta),
                             0.25*(1.0-xi)*(1.0-zeta*zeta),
                            -0.5*zeta*(1.0-xi)*(1.0+eta));
      default: return zero_vector();
    }
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
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    static const unsigned int i0[] =
      {0,1,1,0, 0,1,1,0, 2,1,2,0, 0,1,1,0, 2,1,2,0, 2,2,1,2,0,2,2};
    static const unsigned int i1[] =
      {0,0,1,1, 0,0,1,1, 0,2,1,2, 0,0,1,1, 0,2,1,2, 2,0,2,1,2,2,2};
    static const unsigned int i2[] =
      {0,0,0,0, 1,1,1,1, 0,0,0,0, 2,2,2,2, 1,1,1,1, 0,2,2,2,2,1,2};
    return L(i0[i], xi) * L(i1[i], eta) * L(i2[i], zeta);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    static const unsigned int i0[] =
      {0,1,1,0, 0,1,1,0, 2,1,2,0, 0,1,1,0, 2,1,2,0, 2,2,1,2,0,2,2};
    static const unsigned int i1[] =
      {0,0,1,1, 0,0,1,1, 0,2,1,2, 0,0,1,1, 0,2,1,2, 2,0,2,1,2,2,2};
    static const unsigned int i2[] =
      {0,0,0,0, 1,1,1,1, 0,0,0,0, 2,2,2,2, 1,1,1,1, 0,2,2,2,2,1,2};
    const Real lxi   = L(i0[i], xi);
    const Real leta  = L(i1[i], eta);
    const Real lzeta = L(i2[i], zeta);
    return make_vector(dL(i0[i], xi)  * leta  * lzeta,
                       lxi * dL(i1[i], eta)   * lzeta,
                       lxi * leta  * dL(i2[i], zeta));
  }
#endif
};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_LAGRANGE_3D_H
