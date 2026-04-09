// Kokkos on-device FE shape function dispatch (fe_evaluator.h).
//
// Provides:
//   nativeMapShape    — isoparametric Lagrange shape (topology-based)
//   nativeGradMapShape — isoparametric Lagrange gradient (topology-based)
//   nativeShape        — physics FE shape (FEShapeKey-based)
//   nativeGradShape    — physics FE gradient (FEShapeKey-based)
//
// All functions are KOKKOS_INLINE_FUNCTION and dispatch via switch statements
// that compile to fast GPU branch logic.
//
// Compiled only when LIBMESH_HAVE_KOKKOS is defined (i.e. from .K translation
// units compiled by the Kokkos device compiler).

#pragma once

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos/fe_base.h"
#include "kokkos/fe_types.h"
#include "kokkos/fe_lagrange_1d.h"
#include "kokkos/fe_lagrange_2d.h"
#include "kokkos/fe_lagrange_3d.h"
#include "kokkos/fe_monomial.h"

namespace libMesh::Kokkos
{

// ── On-device helpers: element class -> spatial dimension ─────────────────────

KOKKOS_INLINE_FUNCTION unsigned int
dimFromClass(FEElemClass cls)
{
  switch (cls)
  {
    case FEElemClass::EDGE:
      return 1;
    case FEElemClass::TRI:
    case FEElemClass::QUAD:
      return 2;
    case FEElemClass::TET:
    case FEElemClass::HEX:
    case FEElemClass::PRISM:
    case FEElemClass::PYRAMID:
      return 3;
    default:
      return 0;
  }
}

// ── On-device helper: (class, order) -> canonical Lagrange topology ────────────

KOKKOS_INLINE_FUNCTION FEElemTopology
lagrangeTopologyForClassAndOrder(FEElemClass cls, unsigned int order)
{
  switch (cls)
  {
    case FEElemClass::EDGE:
      switch (order)
      {
        case 1: return FEElemTopology::EDGE2;
        case 2: return FEElemTopology::EDGE3;
        default: return FEElemTopology::EDGE2;
      }
    case FEElemClass::TRI:
      switch (order)
      {
        case 1: return FEElemTopology::TRI3;
        case 2: return FEElemTopology::TRI6;
        default: return FEElemTopology::TRI3;
      }
    case FEElemClass::QUAD:
      switch (order)
      {
        case 1: return FEElemTopology::QUAD4;
        case 2: return FEElemTopology::QUAD9;
        default: return FEElemTopology::QUAD4;
      }
    case FEElemClass::TET:
      switch (order)
      {
        case 1: return FEElemTopology::TET4;
        case 2: return FEElemTopology::TET10;
        default: return FEElemTopology::TET4;
      }
    case FEElemClass::HEX:
      switch (order)
      {
        case 1: return FEElemTopology::HEX8;
        case 2: return FEElemTopology::HEX27;
        default: return FEElemTopology::HEX8;
      }
    default:
      return FEElemTopology::EDGE2;
  }
}

// ── Geometry-only shape dispatch (topology-based) ─────────────────────────────
//
// Used by mapFaceQpToParent() for the isoparametric mapping from face reference
// coordinates to parent reference coordinates.

/// Evaluate the i-th isoparametric Lagrange shape function at (xi, eta, zeta).
KOKKOS_INLINE_FUNCTION Real
nativeMapShape(FEElemTopology topo, unsigned int i, Real xi, Real eta, Real zeta)
{
  switch (topo)
  {
    case FEElemTopology::EDGE2:
      return FEEvaluator<LagrangeTag, Edge2Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::EDGE3:
      return FEEvaluator<LagrangeTag, Edge3Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::TRI3:
      return FEEvaluator<LagrangeTag, Tri3Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::TRI6:
      return FEEvaluator<LagrangeTag, Tri6Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD4:
      return FEEvaluator<LagrangeTag, Quad4Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD8:
      return FEEvaluator<LagrangeTag, Quad8Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD9:
      return FEEvaluator<LagrangeTag, Quad9Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::TET4:
      return FEEvaluator<LagrangeTag, Tet4Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::TET10:
      return FEEvaluator<LagrangeTag, Tet10Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::HEX8:
      return FEEvaluator<LagrangeTag, Hex8Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::HEX20:
      return FEEvaluator<LagrangeTag, Hex20Tag>::shape(i, xi, eta, zeta);
    case FEElemTopology::HEX27:
      return FEEvaluator<LagrangeTag, Hex27Tag>::shape(i, xi, eta, zeta);
    default:
      return Real(0);
  }
}

/// Evaluate the reference-space gradient of the i-th isoparametric Lagrange shape function.
KOKKOS_INLINE_FUNCTION Real3
nativeGradMapShape(FEElemTopology topo, unsigned int i, Real xi, Real eta, Real zeta)
{
  switch (topo)
  {
    case FEElemTopology::EDGE2:
      return FEEvaluator<LagrangeTag, Edge2Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::EDGE3:
      return FEEvaluator<LagrangeTag, Edge3Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::TRI3:
      return FEEvaluator<LagrangeTag, Tri3Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::TRI6:
      return FEEvaluator<LagrangeTag, Tri6Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD4:
      return FEEvaluator<LagrangeTag, Quad4Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD8:
      return FEEvaluator<LagrangeTag, Quad8Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::QUAD9:
      return FEEvaluator<LagrangeTag, Quad9Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::TET4:
      return FEEvaluator<LagrangeTag, Tet4Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::TET10:
      return FEEvaluator<LagrangeTag, Tet10Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::HEX8:
      return FEEvaluator<LagrangeTag, Hex8Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::HEX20:
      return FEEvaluator<LagrangeTag, Hex20Tag>::grad_shape(i, xi, eta, zeta);
    case FEElemTopology::HEX27:
      return FEEvaluator<LagrangeTag, Hex27Tag>::grad_shape(i, xi, eta, zeta);
    default:
      return Real3(0, 0, 0);
  }
}

// ── Physics shape dispatch (FEShapeKey-based) ─────────────────────────────────

/// Evaluate the i-th physics shape function at (xi, eta, zeta).
KOKKOS_INLINE_FUNCTION Real
nativeShape(FEShapeKey key, unsigned int i, Real xi, Real eta, Real zeta)
{
  switch (key.family)
  {
    case FEFamily::LAGRANGE:
    {
      switch (lagrangeTopologyForClassAndOrder(key.cls, key.order))
      {
        case FEElemTopology::EDGE2:
          return FEEvaluator<LagrangeTag, Edge2Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::EDGE3:
          return FEEvaluator<LagrangeTag, Edge3Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::TRI3:
          return FEEvaluator<LagrangeTag, Tri3Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::TRI6:
          return FEEvaluator<LagrangeTag, Tri6Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD4:
          return FEEvaluator<LagrangeTag, Quad4Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD8:
          return FEEvaluator<LagrangeTag, Quad8Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD9:
          return FEEvaluator<LagrangeTag, Quad9Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::TET4:
          return FEEvaluator<LagrangeTag, Tet4Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::TET10:
          return FEEvaluator<LagrangeTag, Tet10Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::HEX8:
          return FEEvaluator<LagrangeTag, Hex8Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::HEX20:
          return FEEvaluator<LagrangeTag, Hex20Tag>::shape(i, xi, eta, zeta);
        case FEElemTopology::HEX27:
          return FEEvaluator<LagrangeTag, Hex27Tag>::shape(i, xi, eta, zeta);
        default:
          return Real(0);
      }
    }

    case FEFamily::MONOMIAL:
    {
      switch (dimFromClass(key.cls))
      {
        case 1:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim1Tag, 0>::shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim1Tag, 1>::shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim1Tag, 2>::shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim1Tag, 3>::shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim1Tag, 4>::shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim1Tag, 5>::shape(i, xi, eta, zeta);
            default: return Real(0);
          }
        case 2:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim2Tag, 0>::shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim2Tag, 1>::shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim2Tag, 2>::shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim2Tag, 3>::shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim2Tag, 4>::shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim2Tag, 5>::shape(i, xi, eta, zeta);
            default: return Real(0);
          }
        case 3:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim3Tag, 0>::shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim3Tag, 1>::shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim3Tag, 2>::shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim3Tag, 3>::shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim3Tag, 4>::shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim3Tag, 5>::shape(i, xi, eta, zeta);
            default: return Real(0);
          }
        default:
          return Real(0);
      }
    }

    default:
      return Real(0);
  }
}

/// Evaluate the reference-space gradient of the i-th physics shape function.
/// The physical gradient is obtained by: grad_u_physical = J * nativeGradShape(key, ...)
KOKKOS_INLINE_FUNCTION Real3
nativeGradShape(FEShapeKey key, unsigned int i, Real xi, Real eta, Real zeta)
{
  switch (key.family)
  {
    case FEFamily::LAGRANGE:
    {
      switch (lagrangeTopologyForClassAndOrder(key.cls, key.order))
      {
        case FEElemTopology::EDGE2:
          return FEEvaluator<LagrangeTag, Edge2Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::EDGE3:
          return FEEvaluator<LagrangeTag, Edge3Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::TRI3:
          return FEEvaluator<LagrangeTag, Tri3Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::TRI6:
          return FEEvaluator<LagrangeTag, Tri6Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD4:
          return FEEvaluator<LagrangeTag, Quad4Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD8:
          return FEEvaluator<LagrangeTag, Quad8Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::QUAD9:
          return FEEvaluator<LagrangeTag, Quad9Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::TET4:
          return FEEvaluator<LagrangeTag, Tet4Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::TET10:
          return FEEvaluator<LagrangeTag, Tet10Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::HEX8:
          return FEEvaluator<LagrangeTag, Hex8Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::HEX20:
          return FEEvaluator<LagrangeTag, Hex20Tag>::grad_shape(i, xi, eta, zeta);
        case FEElemTopology::HEX27:
          return FEEvaluator<LagrangeTag, Hex27Tag>::grad_shape(i, xi, eta, zeta);
        default:
          return Real3(0, 0, 0);
      }
    }

    case FEFamily::MONOMIAL:
    {
      switch (dimFromClass(key.cls))
      {
        case 1:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim1Tag, 0>::grad_shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim1Tag, 1>::grad_shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim1Tag, 2>::grad_shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim1Tag, 3>::grad_shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim1Tag, 4>::grad_shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim1Tag, 5>::grad_shape(i, xi, eta, zeta);
            default: return Real3(0, 0, 0);
          }
        case 2:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim2Tag, 0>::grad_shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim2Tag, 1>::grad_shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim2Tag, 2>::grad_shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim2Tag, 3>::grad_shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim2Tag, 4>::grad_shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim2Tag, 5>::grad_shape(i, xi, eta, zeta);
            default: return Real3(0, 0, 0);
          }
        case 3:
          switch (key.order)
          {
            case 0: return FEEvaluator<MonomialTag, Dim3Tag, 0>::grad_shape(i, xi, eta, zeta);
            case 1: return FEEvaluator<MonomialTag, Dim3Tag, 1>::grad_shape(i, xi, eta, zeta);
            case 2: return FEEvaluator<MonomialTag, Dim3Tag, 2>::grad_shape(i, xi, eta, zeta);
            case 3: return FEEvaluator<MonomialTag, Dim3Tag, 3>::grad_shape(i, xi, eta, zeta);
            case 4: return FEEvaluator<MonomialTag, Dim3Tag, 4>::grad_shape(i, xi, eta, zeta);
            case 5: return FEEvaluator<MonomialTag, Dim3Tag, 5>::grad_shape(i, xi, eta, zeta);
            default: return Real3(0, 0, 0);
          }
        default:
          return Real3(0, 0, 0);
      }
    }

    default:
      return Real3(0, 0, 0);
  }
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_HAVE_KOKKOS
