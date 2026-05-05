// Kokkos on-device FE shape function dispatch (fe_evaluator.h).
//
// Provides:
//   map_shape           — isoparametric Lagrange shape (topology-based)
//   grad_map_shape      — isoparametric Lagrange gradient (topology-based)
//   shape               — physics FE shape (FEShapeKey-based)
//   grad_shape          — physics FE gradient (FEShapeKey-based)
//
// All functions are LIBMESH_DEVICE_INLINE and dispatch via switch statements
// that compile to fast GPU branch logic.
//
// These helpers are intended for Kokkos-enabled code paths. Device execution
// happens from .K translation units, but the header is also parsed by host code.

#ifndef LIBMESH_KOKKOS_FE_EVALUATOR_H
#define LIBMESH_KOKKOS_FE_EVALUATOR_H

#include "gpu/kokkos_fe_base.h"
#include "gpu/kokkos_fe_types.h"
#include "gpu/kokkos_fe_lagrange_1d.h"
#include "gpu/kokkos_fe_lagrange_2d.h"
#include "gpu/kokkos_fe_lagrange_3d.h"
#include "gpu/kokkos_fe_monomial.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"

namespace libMesh::Kokkos
{

// ── On-device helpers: element class -> spatial dimension ─────────────────────

LIBMESH_DEVICE_INLINE unsigned int
dim_from_class(FEElemClass cls)
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
      detail::abort_unsupported("dim_from_class(): unsupported element class");
      return 0;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
dim_from_topology(libMesh::ElemType topo)
{
  return dim_from_class(class_from_topology(topo));
}

// ── On-device helper: exact libMesh Lagrange key -> evaluator topology ─────────

LIBMESH_DEVICE_INLINE libMesh::ElemType
lagrange_shape_topology_for_key(FEShapeKey key)
{
  const libMesh::ElemType topo = lagrange_shape_topology_or_invalid(key);

  if (topo == libMesh::INVALID_ELEM)
  {
    detail::abort_unsupported("lagrange_shape_topology_for_key(): unsupported LAGRANGE key for current Kokkos evaluator support boundary");
    return libMesh::INVALID_ELEM;
  }

  return topo;
}

LIBMESH_DEVICE_INLINE Real
eval_lagrange_shape(libMesh::ElemType topo,
                    unsigned int i,
                    Real xi,
                    Real eta,
                    Real zeta)
{
  switch (topo)
  {
    case libMesh::EDGE2:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::shape(i, xi, eta, zeta);
    case libMesh::EDGE3:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::shape(i, xi, eta, zeta);
    case libMesh::TRI3:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::shape(i, xi, eta, zeta);
    case libMesh::TRI6:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::shape(i, xi, eta, zeta);
    case libMesh::QUAD4:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::shape(i, xi, eta, zeta);
    case libMesh::QUAD8:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::shape(i, xi, eta, zeta);
    case libMesh::QUAD9:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::shape(i, xi, eta, zeta);
    case libMesh::TET4:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::shape(i, xi, eta, zeta);
    case libMesh::TET10:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::shape(i, xi, eta, zeta);
    case libMesh::HEX8:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::shape(i, xi, eta, zeta);
    case libMesh::HEX20:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::shape(i, xi, eta, zeta);
    case libMesh::HEX27:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::shape(i, xi, eta, zeta);
    default:
      detail::abort_unsupported("eval_lagrange_shape(): unsupported evaluator topology");
      return Real(0);
  }
}

LIBMESH_DEVICE_INLINE RealVector
eval_lagrange_grad_shape(libMesh::ElemType topo,
                         unsigned int i,
                         Real xi,
                         Real eta,
                         Real zeta)
{
  switch (topo)
  {
    case libMesh::EDGE2:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE2>::grad_shape(i, xi, eta, zeta);
    case libMesh::EDGE3:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::EDGE3>::grad_shape(i, xi, eta, zeta);
    case libMesh::TRI3:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TRI3>::grad_shape(i, xi, eta, zeta);
    case libMesh::TRI6:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TRI6>::grad_shape(i, xi, eta, zeta);
    case libMesh::QUAD4:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD4>::grad_shape(i, xi, eta, zeta);
    case libMesh::QUAD8:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD8>::grad_shape(i, xi, eta, zeta);
    case libMesh::QUAD9:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::QUAD9>::grad_shape(i, xi, eta, zeta);
    case libMesh::TET4:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TET4>::grad_shape(i, xi, eta, zeta);
    case libMesh::TET10:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::TET10>::grad_shape(i, xi, eta, zeta);
    case libMesh::HEX8:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX8>::grad_shape(i, xi, eta, zeta);
    case libMesh::HEX20:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX20>::grad_shape(i, xi, eta, zeta);
    case libMesh::HEX27:
      return FEEvaluator<libMesh::LAGRANGE, libMesh::HEX27>::grad_shape(i, xi, eta, zeta);
    default:
      detail::abort_unsupported("eval_lagrange_grad_shape(): unsupported evaluator topology");
      return zero_vector();
  }
}

// ── Geometry-only shape dispatch (mapping-type + topology) ────────────────────
//
// Used by map_face_qp_to_parent() for the isoparametric mapping from face reference
// coordinates to parent reference coordinates.
//
// The mapping_type parameter selects the geometric map family.  Currently only
// LAGRANGE_MAP is supported; RATIONAL_BERNSTEIN_MAP requires additional
// rational-weight data that is not yet threaded through the device path.

// ── Compile-time topology versions (preferred for GPU) ───────────────────
// Template on FEFamily and ElemType so gpu compiler only instantiates the specific
// FEEvaluator specialization. No topology switch means no stack pressure.

/// Compile-time map shape evaluation.
template <libMesh::FEFamily family, libMesh::ElemType topo>
LIBMESH_DEVICE_INLINE Real
map_shape(unsigned int i, Real xi, Real eta, Real zeta)
{
  return FEEvaluator<family, topo>::shape(i, xi, eta, zeta);
}

/// Compile-time map gradient evaluation.
template <libMesh::FEFamily family, libMesh::ElemType topo>
LIBMESH_DEVICE_INLINE RealVector
grad_map_shape(unsigned int i, Real xi, Real eta, Real zeta)
{
  return FEEvaluator<family, topo>::grad_shape(i, xi, eta, zeta);
}

// ── Runtime topology versions (larger GPU stack usage) ───────────────────

/// Evaluate the i-th geometric map shape function at (xi, eta, zeta).
LIBMESH_DEVICE_INLINE Real
map_shape(libMesh::ElemMappingType mapping_type,
          libMesh::ElemType topo,
          unsigned int i,
          Real xi,
          Real eta,
          Real zeta)
{
  switch (mapping_type)
  {
    case libMesh::LAGRANGE_MAP:
      return eval_lagrange_shape(topo, i, xi, eta, zeta);
    default:
      detail::abort_unsupported("map_shape(): only LAGRANGE_MAP is implemented");
      return Real(0);
  }
}

/// Evaluate the reference-space gradient of the i-th geometric map shape function.
LIBMESH_DEVICE_INLINE RealVector
grad_map_shape(libMesh::ElemMappingType mapping_type,
               libMesh::ElemType topo,
               unsigned int i,
               Real xi,
               Real eta,
               Real zeta)
{
  switch (mapping_type)
  {
    case libMesh::LAGRANGE_MAP:
      return eval_lagrange_grad_shape(topo, i, xi, eta, zeta);
    default:
      detail::abort_unsupported("grad_map_shape(): only LAGRANGE_MAP is implemented");
      return zero_vector();
  }
}

// ── Physics shape dispatch (FEShapeKey-based) ─────────────────────────────────

/// Evaluate the i-th physics shape function at (xi, eta, zeta).
LIBMESH_DEVICE_INLINE Real
shape(FEShapeKey key, unsigned int i, Real xi, Real eta, Real zeta)
{
  if (!supports_shape(key))
  {
    detail::abort_unsupported("shape(): unsupported FE key for current Kokkos evaluator support boundary");
    return Real(0);
  }

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return eval_lagrange_shape(lagrange_shape_topology_for_key(key), i, xi, eta, zeta);

    case libMesh::MONOMIAL:
    {
      switch (monomial_evaluator_dim_or_zero(key.elem_type))
      {
        case 1:
          switch (key.order)
          {
            case 0: return MonomialImpl1D<0>::shape(i, xi, eta, zeta);
            case 1: return MonomialImpl1D<1>::shape(i, xi, eta, zeta);
            case 2: return MonomialImpl1D<2>::shape(i, xi, eta, zeta);
            case 3: return MonomialImpl1D<3>::shape(i, xi, eta, zeta);
            case 4: return MonomialImpl1D<4>::shape(i, xi, eta, zeta);
            case 5: return MonomialImpl1D<5>::shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("shape(): unsupported 1D MONOMIAL order");
              return Real(0);
          }
        case 2:
          switch (key.order)
          {
            case 0: return MonomialImpl2D<0>::shape(i, xi, eta, zeta);
            case 1: return MonomialImpl2D<1>::shape(i, xi, eta, zeta);
            case 2: return MonomialImpl2D<2>::shape(i, xi, eta, zeta);
            case 3: return MonomialImpl2D<3>::shape(i, xi, eta, zeta);
            case 4: return MonomialImpl2D<4>::shape(i, xi, eta, zeta);
            case 5: return MonomialImpl2D<5>::shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("shape(): unsupported 2D MONOMIAL order");
              return Real(0);
          }
        case 3:
          switch (key.order)
          {
            case 0: return MonomialImpl3D<0>::shape(i, xi, eta, zeta);
            case 1: return MonomialImpl3D<1>::shape(i, xi, eta, zeta);
            case 2: return MonomialImpl3D<2>::shape(i, xi, eta, zeta);
            case 3: return MonomialImpl3D<3>::shape(i, xi, eta, zeta);
            case 4: return MonomialImpl3D<4>::shape(i, xi, eta, zeta);
            case 5: return MonomialImpl3D<5>::shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("shape(): unsupported 3D MONOMIAL order");
              return Real(0);
          }
        default:
          detail::abort_unsupported("shape(): unsupported MONOMIAL element topology");
          return Real(0);
      }
    }

    default:
      detail::abort_unsupported("shape(): unsupported FE family");
      return Real(0);
  }
}

/// Evaluate the reference-space gradient of the i-th physics shape function.
/// With J from jacobian(), rows are reference derivatives, so
/// grad_ref = J * grad_phys and grad_phys = J.inverse(dim) * grad_ref.
LIBMESH_DEVICE_INLINE RealVector
grad_shape(FEShapeKey key, unsigned int i, Real xi, Real eta, Real zeta)
{
  if (!supports_grad_shape(key))
  {
    detail::abort_unsupported("grad_shape(): unsupported FE key for current Kokkos evaluator support boundary");
    return zero_vector();
  }

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return eval_lagrange_grad_shape(lagrange_shape_topology_for_key(key), i, xi, eta, zeta);

    case libMesh::MONOMIAL:
    {
      switch (monomial_evaluator_dim_or_zero(key.elem_type))
      {
        case 1:
          switch (key.order)
          {
            case 0: return MonomialImpl1D<0>::grad_shape(i, xi, eta, zeta);
            case 1: return MonomialImpl1D<1>::grad_shape(i, xi, eta, zeta);
            case 2: return MonomialImpl1D<2>::grad_shape(i, xi, eta, zeta);
            case 3: return MonomialImpl1D<3>::grad_shape(i, xi, eta, zeta);
            case 4: return MonomialImpl1D<4>::grad_shape(i, xi, eta, zeta);
            case 5: return MonomialImpl1D<5>::grad_shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("grad_shape(): unsupported 1D MONOMIAL order");
              return zero_vector();
          }
        case 2:
          switch (key.order)
          {
            case 0: return MonomialImpl2D<0>::grad_shape(i, xi, eta, zeta);
            case 1: return MonomialImpl2D<1>::grad_shape(i, xi, eta, zeta);
            case 2: return MonomialImpl2D<2>::grad_shape(i, xi, eta, zeta);
            case 3: return MonomialImpl2D<3>::grad_shape(i, xi, eta, zeta);
            case 4: return MonomialImpl2D<4>::grad_shape(i, xi, eta, zeta);
            case 5: return MonomialImpl2D<5>::grad_shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("grad_shape(): unsupported 2D MONOMIAL order");
              return zero_vector();
          }
        case 3:
          switch (key.order)
          {
            case 0: return MonomialImpl3D<0>::grad_shape(i, xi, eta, zeta);
            case 1: return MonomialImpl3D<1>::grad_shape(i, xi, eta, zeta);
            case 2: return MonomialImpl3D<2>::grad_shape(i, xi, eta, zeta);
            case 3: return MonomialImpl3D<3>::grad_shape(i, xi, eta, zeta);
            case 4: return MonomialImpl3D<4>::grad_shape(i, xi, eta, zeta);
            case 5: return MonomialImpl3D<5>::grad_shape(i, xi, eta, zeta);
            default:
              detail::abort_unsupported("grad_shape(): unsupported 3D MONOMIAL order");
              return zero_vector();
          }
        default:
          detail::abort_unsupported("grad_shape(): unsupported MONOMIAL element topology");
          return zero_vector();
      }
    }

    default:
      detail::abort_unsupported("grad_shape(): unsupported FE family");
      return zero_vector();
  }
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_EVALUATOR_H
