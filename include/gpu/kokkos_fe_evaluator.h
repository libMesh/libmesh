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

#include "kokkos_fe_base.h"
#include "kokkos_fe_types.h"
#include "kokkos_fe_lagrange_1d.h"
#include "kokkos_fe_lagrange_2d.h"
#include "kokkos_fe_lagrange_3d.h"
#include "kokkos_fe_monomial.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"

namespace libMesh::Kokkos
{

LIBMESH_DEVICE_INLINE libMesh::ElemType
lagrange_shape_topology_for_key(FEShapeKey key);

LIBMESH_DEVICE_INLINE Real
eval_lagrange_shape(libMesh::ElemType topo,
                    unsigned int i,
                    Real xi,
                    Real eta,
                    Real zeta);

LIBMESH_DEVICE_INLINE RealVector
eval_lagrange_grad_shape(libMesh::ElemType topo,
                         unsigned int i,
                         Real xi,
                         Real eta,
                         Real zeta);

namespace detail
{

template <typename Op>
LIBMESH_DEVICE_INLINE auto
dispatch_lagrange_topology(libMesh::ElemType topo, const Op & op)
  -> decltype(op.template operator()<libMesh::EDGE2>())
{
  return libMesh::dispatch_lagrange_map_topology_or(
    topo,
    op,
    [&](libMesh::ElemType) -> decltype(op.template operator()<libMesh::EDGE2>())
    {
      detail::abort_unsupported("dispatch_lagrange_topology(): unsupported evaluator topology");
      return op.template operator()<libMesh::EDGE2>();
    });
}

template <unsigned int Dim, typename Op>
LIBMESH_DEVICE_INLINE auto
dispatch_monomial_order(libMesh::Order order, const Op & op)
  -> decltype(op.template operator()<Dim, 0>())
{
  switch (order)
  {
    case libMesh::CONSTANT: return op.template operator()<Dim, 0>();
    case libMesh::FIRST: return op.template operator()<Dim, 1>();
    case libMesh::SECOND: return op.template operator()<Dim, 2>();
    case libMesh::THIRD: return op.template operator()<Dim, 3>();
    case libMesh::FOURTH: return op.template operator()<Dim, 4>();
    case libMesh::FIFTH: return op.template operator()<Dim, 5>();
    default:
      detail::abort_unsupported("dispatch_monomial_order(): unsupported MONOMIAL order");
      return op.template operator()<Dim, 0>();
  }
}

template <typename Op>
LIBMESH_DEVICE_INLINE auto
dispatch_monomial(libMesh::ElemType elem_type, libMesh::Order order, const Op & op)
  -> decltype(op.template operator()<1, 0>())
{
  switch (monomial_evaluator_dim_or_zero(elem_type))
  {
    case 1: return dispatch_monomial_order<1>(order, op);
    case 2: return dispatch_monomial_order<2>(order, op);
    case 3: return dispatch_monomial_order<3>(order, op);
    default:
      detail::abort_unsupported("dispatch_monomial(): unsupported MONOMIAL element topology");
      return op.template operator()<1, 0>();
  }
}

struct LagrangeShapeOp
{
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  template <libMesh::ElemType topo>
  LIBMESH_DEVICE_INLINE Real operator()() const
  {
    return FEEvaluator<libMesh::LAGRANGE, topo>::shape(i, xi, eta, zeta);
  }
};

struct LagrangeGradShapeOp
{
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  template <libMesh::ElemType topo>
  LIBMESH_DEVICE_INLINE RealVector operator()() const
  {
    return FEEvaluator<libMesh::LAGRANGE, topo>::grad_shape(i, xi, eta, zeta);
  }
};

struct MonomialShapeOp
{
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  template <unsigned int Dim, unsigned int Order>
  LIBMESH_DEVICE_INLINE Real operator()() const
  {
    if constexpr (Dim == 1)
      return MonomialImpl1D<Order>::shape(i, xi, eta, zeta);
    else if constexpr (Dim == 2)
      return MonomialImpl2D<Order>::shape(i, xi, eta, zeta);
    else
      return MonomialImpl3D<Order>::shape(i, xi, eta, zeta);
  }
};

struct MonomialGradShapeOp
{
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  template <unsigned int Dim, unsigned int Order>
  LIBMESH_DEVICE_INLINE RealVector operator()() const
  {
    if constexpr (Dim == 1)
      return MonomialImpl1D<Order>::grad_shape(i, xi, eta, zeta);
    else if constexpr (Dim == 2)
      return MonomialImpl2D<Order>::grad_shape(i, xi, eta, zeta);
    else
      return MonomialImpl3D<Order>::grad_shape(i, xi, eta, zeta);
  }
};

template <typename LagrangeOp, typename MonomialOp>
LIBMESH_DEVICE_INLINE auto
dispatch_shape_family(libMesh::FEShapeKey key,
                      const LagrangeOp & lagrange_op,
                      const MonomialOp & monomial_op,
                      const char * unsupported_message)
  -> decltype(lagrange_op())
{
  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return lagrange_op();

    case libMesh::MONOMIAL:
      return monomial_op();

    default:
      detail::abort_unsupported(unsupported_message);
      return lagrange_op();
  }
}

struct KeyedLagrangeShapeOp
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  LIBMESH_DEVICE_INLINE Real operator()() const
  {
    return eval_lagrange_shape(lagrange_shape_topology_for_key(key), i, xi, eta, zeta);
  }
};

struct KeyedLagrangeGradShapeOp
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  LIBMESH_DEVICE_INLINE RealVector operator()() const
  {
    return eval_lagrange_grad_shape(lagrange_shape_topology_for_key(key), i, xi, eta, zeta);
  }
};

struct KeyedMonomialShapeOp
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  LIBMESH_DEVICE_INLINE Real operator()() const
  {
    return detail::dispatch_monomial(key.elem_type,
                                     key.order,
                                     detail::MonomialShapeOp{i, xi, eta, zeta});
  }
};

struct KeyedMonomialGradShapeOp
{
  libMesh::Kokkos::FEShapeKey key;
  unsigned int i;
  Real xi;
  Real eta;
  Real zeta;

  LIBMESH_DEVICE_INLINE RealVector operator()() const
  {
    return detail::dispatch_monomial(key.elem_type,
                                     key.order,
                                     detail::MonomialGradShapeOp{i, xi, eta, zeta});
  }
};

} // namespace detail

// ── On-device helpers: element class -> spatial dimension ─────────────────────

LIBMESH_DEVICE_INLINE unsigned int
dim_from_class(FEElemClass cls)
{
  const unsigned int dim = libMesh::elem_class_dim_or_zero(cls);

  if (!dim)
  {
    detail::abort_unsupported("dim_from_class(): unsupported element class");
    return 0;
  }

  return dim;
}

LIBMESH_DEVICE_INLINE unsigned int
dim_from_topology(libMesh::ElemType topo)
{
  const unsigned int dim = libMesh::topology_dim_or_zero(topo);

  if (!dim)
  {
    detail::abort_unsupported("dim_from_topology(): unsupported element type");
    return 0;
  }

  return dim;
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
  return detail::dispatch_lagrange_topology(topo, detail::LagrangeShapeOp{i, xi, eta, zeta});
}

LIBMESH_DEVICE_INLINE RealVector
eval_lagrange_grad_shape(libMesh::ElemType topo,
                         unsigned int i,
                         Real xi,
                         Real eta,
                         Real zeta)
{
  return detail::dispatch_lagrange_topology(topo, detail::LagrangeGradShapeOp{i, xi, eta, zeta});
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

  return detail::dispatch_shape_family(
    key,
    detail::KeyedLagrangeShapeOp{key, i, xi, eta, zeta},
    detail::KeyedMonomialShapeOp{key, i, xi, eta, zeta},
    "shape(): unsupported FE family");
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

  return detail::dispatch_shape_family(
    key,
    detail::KeyedLagrangeGradShapeOp{key, i, xi, eta, zeta},
    detail::KeyedMonomialGradShapeOp{key, i, xi, eta, zeta},
    "grad_shape(): unsupported FE family");
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_EVALUATOR_H
