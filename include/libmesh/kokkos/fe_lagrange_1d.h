// Kokkos FEEvaluator specializations for 1-D Lagrange elements.
//
// Covers EDGE2 (linear) and EDGE3 (quadratic).
// Reference-element coordinate convention (libMesh-compatible):
//   EDGE2/EDGE3: xi in [-1, 1]
//
// EDGE3 node ordering (libMesh non-sequential):
//   index 0 -> xi = -1   (left node)
//   index 1 -> xi = +1   (right node)
//   index 2 -> xi =  0   (midpoint)

#pragma once

#include "libmesh/kokkos/fe_base.h"

namespace libMesh::Kokkos
{

// ── EDGE2 (linear edge, 2 nodes) ─────────────────────────────────────────────

template <>
struct FEEvaluator<LagrangeTag, Edge2Tag>
{
  static constexpr unsigned int n_dofs() { return 2; }

#ifdef MOOSE_KOKKOS_SCOPE
  KOKKOS_INLINE_FUNCTION static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 0.5 * (1.0 - xi);
      case 1: return 0.5 * (1.0 + xi);
      default: return 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION static Real3
  grad_shape(unsigned int i, Real /*xi*/, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return Real3(-0.5, 0.0, 0.0);
      case 1: return Real3( 0.5, 0.0, 0.0);
      default: return Real3(0.0, 0.0, 0.0);
    }
  }
#endif
};

// ── EDGE3 (quadratic edge, 3 nodes) ──────────────────────────────────────────
// Node ordering matches libMesh: 0->left(-1), 1->right(+1), 2->mid(0)
//   L_0(xi) = 0.5*xi*(xi-1)   dL_0/dxi = xi - 0.5
//   L_1(xi) = 0.5*xi*(xi+1)   dL_1/dxi = xi + 0.5
//   L_2(xi) = 1 - xi²         dL_2/dxi = -2*xi

template <>
struct FEEvaluator<LagrangeTag, Edge3Tag>
{
  static constexpr unsigned int n_dofs() { return 3; }

#ifdef MOOSE_KOKKOS_SCOPE
  KOKKOS_INLINE_FUNCTION static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return 0.5 * xi * (xi - 1.0);
      case 1: return 0.5 * xi * (xi + 1.0);
      case 2: return 1.0 - xi * xi;
      default: return 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION static Real3
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    switch (i)
    {
      case 0: return Real3(xi - 0.5,  0.0, 0.0);
      case 1: return Real3(xi + 0.5,  0.0, 0.0);
      case 2: return Real3(-2.0 * xi, 0.0, 0.0);
      default: return Real3(0.0, 0.0, 0.0);
    }
  }
#endif
};

} // namespace libMesh::Kokkos
