// Shared Kokkos FE shape dispatch helpers.
//
// These helpers capture the supported Kokkos FE evaluator boundary in one
// place so production code and oracle tests can dispatch exact FE keys without
// duplicating the support matrix.

#ifndef LIBMESH_KOKKOS_FE_SHAPE_DISPATCH_H
#define LIBMESH_KOKKOS_FE_SHAPE_DISPATCH_H

#include "libmesh/fe_shape_traits.h"
#include "libmesh/kokkos_fe_evaluator.h"

namespace libMesh::Kokkos
{

template <unsigned int Dim, libMesh::Order ExactOrder>
struct monomial_order_evaluator;

#define LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(dim_value, exact_order, impl_suffix, impl_order) \
  template <>                                                                                \
  struct monomial_order_evaluator<dim_value, libMesh::exact_order>                           \
  {                                                                                          \
    LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,                         \
                                                      libMesh::Real xi,                      \
                                                      libMesh::Real eta,                     \
                                                      libMesh::Real zeta)                    \
    {                                                                                        \
      return libMesh::Kokkos::impl_suffix<impl_order>::shape(i, xi, eta, zeta);             \
    }                                                                                        \
                                                                                             \
    LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(                     \
      unsigned int i, libMesh::Real xi, libMesh::Real eta, libMesh::Real zeta)              \
    {                                                                                        \
      return libMesh::Kokkos::impl_suffix<impl_order>::grad_shape(i, xi, eta, zeta);        \
    }                                                                                        \
  }

LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, CONSTANT, MonomialImpl1D, 0);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, FIRST, MonomialImpl1D, 1);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, SECOND, MonomialImpl1D, 2);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, THIRD, MonomialImpl1D, 3);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, FOURTH, MonomialImpl1D, 4);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(1, FIFTH, MonomialImpl1D, 5);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, CONSTANT, MonomialImpl2D, 0);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, FIRST, MonomialImpl2D, 1);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, SECOND, MonomialImpl2D, 2);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, THIRD, MonomialImpl2D, 3);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, FOURTH, MonomialImpl2D, 4);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(2, FIFTH, MonomialImpl2D, 5);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, CONSTANT, MonomialImpl3D, 0);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, FIRST, MonomialImpl3D, 1);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, SECOND, MonomialImpl3D, 2);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, THIRD, MonomialImpl3D, 3);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, FOURTH, MonomialImpl3D, 4);
LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE(3, FIFTH, MonomialImpl3D, 5);

#undef LIBMESH_KOKKOS_MONOMIAL_ORDER_CASE

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator;

template <libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator<libMesh::LAGRANGE, ExactTopo, ExactOrder>
{
  static constexpr libMesh::FEShapeKey exact_key{ libMesh::LAGRANGE, ExactTopo, ExactOrder };
  static constexpr libMesh::ElemType evaluator_topology =
    libMesh::lagrange_shape_topology_or_invalid(exact_key);

  LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,
                                                   libMesh::Real xi,
                                                   libMesh::Real eta,
                                                   libMesh::Real zeta)
  {
    return map_shape<libMesh::LAGRANGE, evaluator_topology>(
      i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(unsigned int i,
                                                                      libMesh::Real xi,
                                                                      libMesh::Real eta,
                                                                      libMesh::Real zeta)
  {
    return grad_map_shape<libMesh::LAGRANGE, evaluator_topology>(
      i, xi, eta, zeta);
  }
};

template <libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator<libMesh::MONOMIAL, ExactTopo, ExactOrder>
{
  static constexpr unsigned int evaluator_dim =
    libMesh::monomial_evaluator_dim_or_zero(ExactTopo);

  LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,
                                                   libMesh::Real xi,
                                                   libMesh::Real eta,
                                                   libMesh::Real zeta)
  {
    return monomial_order_evaluator<evaluator_dim, ExactOrder>::shape(
      i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(unsigned int i,
                                                                      libMesh::Real xi,
                                                                      libMesh::Real eta,
                                                                      libMesh::Real zeta)
  {
    return monomial_order_evaluator<evaluator_dim, ExactOrder>::grad_shape(
      i, xi, eta, zeta);
  }
};

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
LIBMESH_DEVICE_INLINE libMesh::Real
shape_for_key(unsigned int i, libMesh::Real xi, libMesh::Real eta, libMesh::Real zeta)
{
  return exact_shape_evaluator<Family, ExactTopo, ExactOrder>::shape(i, xi, eta, zeta);
}

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
LIBMESH_DEVICE_INLINE libMesh::Kokkos::RealVector
grad_shape_for_key(unsigned int i, libMesh::Real xi, libMesh::Real eta, libMesh::Real zeta)
{
  return exact_shape_evaluator<Family, ExactTopo, ExactOrder>::grad_shape(i, xi, eta, zeta);
}

template <libMesh::ElemType ExactTopo, typename Dispatcher>
inline int
dispatch_supported_monomial_order(libMesh::Order order, const Dispatcher & dispatcher)
{
  switch (order)
  {
    case libMesh::CONSTANT:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::CONSTANT>();
    case libMesh::FIRST:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::FIRST>();
    case libMesh::SECOND:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::SECOND>();
    case libMesh::THIRD:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::THIRD>();
    case libMesh::FOURTH:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::FOURTH>();
    case libMesh::FIFTH:
      return dispatcher.template operator()<libMesh::MONOMIAL, ExactTopo, libMesh::FIFTH>();
    default:
      return dispatcher.unsupported_key(libMesh::FEShapeKey{ libMesh::MONOMIAL, ExactTopo, order });
  }
}

template <typename Dispatcher>
inline int
dispatch_exact_lagrange_shape_key(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  switch (key.elem_type)
  {
    case libMesh::EDGE2:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE2, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::EDGE3:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE3, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE3, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::EDGE4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE4, libMesh::FIRST>();
        case libMesh::THIRD:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE4, libMesh::THIRD>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TRI3:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI3, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TRI6:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI6, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI6, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TRI7:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI7, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI7, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD4, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD8:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD8, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD8, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD9:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD9, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD9, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TET4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET4, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TET10:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET10, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET10, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TET14:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET14, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET14, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX8:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX8, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX20:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX20, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX20, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX27:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX27, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX27, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    default:
      return dispatcher.unsupported_key(key);
  }
}

template <typename Dispatcher>
inline int
dispatch_exact_lagrange_shape_key_with_map(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  switch (key.elem_type)
  {
    case libMesh::EDGE2:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE2, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::EDGE3:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE3, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE3, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::EDGE4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE4, libMesh::FIRST>();
        case libMesh::THIRD:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::EDGE4, libMesh::THIRD>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TRI3:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI3, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TRI6:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI6, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TRI6, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD4, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD8:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD8, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD8, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::QUAD9:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD9, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::QUAD9, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TET4:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET4, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::TET10:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET10, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::TET10, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX8:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX8, libMesh::FIRST>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX20:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX20, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX20, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    case libMesh::HEX27:
      switch (key.order)
      {
        case libMesh::FIRST:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX27, libMesh::FIRST>();
        case libMesh::SECOND:
          return dispatcher.template operator()<libMesh::LAGRANGE, libMesh::HEX27, libMesh::SECOND>();
        default:
          return dispatcher.unsupported_key(key);
      }
    default:
      return dispatcher.unsupported_key(key);
  }
}

template <typename Dispatcher>
inline int
dispatch_exact_monomial_shape_key(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  switch (key.elem_type)
  {
    case libMesh::EDGE2:
      return dispatch_supported_monomial_order<libMesh::EDGE2>(key.order, dispatcher);
    case libMesh::EDGE3:
      return dispatch_supported_monomial_order<libMesh::EDGE3>(key.order, dispatcher);
    case libMesh::EDGE4:
      return dispatch_supported_monomial_order<libMesh::EDGE4>(key.order, dispatcher);
    case libMesh::TRI3:
      return dispatch_supported_monomial_order<libMesh::TRI3>(key.order, dispatcher);
    case libMesh::TRI6:
      return dispatch_supported_monomial_order<libMesh::TRI6>(key.order, dispatcher);
    case libMesh::TRI7:
      return dispatch_supported_monomial_order<libMesh::TRI7>(key.order, dispatcher);
    case libMesh::QUAD4:
      return dispatch_supported_monomial_order<libMesh::QUAD4>(key.order, dispatcher);
    case libMesh::QUAD8:
      return dispatch_supported_monomial_order<libMesh::QUAD8>(key.order, dispatcher);
    case libMesh::QUAD9:
      return dispatch_supported_monomial_order<libMesh::QUAD9>(key.order, dispatcher);
    case libMesh::TET4:
      return dispatch_supported_monomial_order<libMesh::TET4>(key.order, dispatcher);
    case libMesh::TET10:
      return dispatch_supported_monomial_order<libMesh::TET10>(key.order, dispatcher);
    case libMesh::TET14:
      return dispatch_supported_monomial_order<libMesh::TET14>(key.order, dispatcher);
    case libMesh::HEX8:
      return dispatch_supported_monomial_order<libMesh::HEX8>(key.order, dispatcher);
    case libMesh::HEX20:
      return dispatch_supported_monomial_order<libMesh::HEX20>(key.order, dispatcher);
    case libMesh::HEX27:
      return dispatch_supported_monomial_order<libMesh::HEX27>(key.order, dispatcher);
    case libMesh::PRISM6:
      return dispatch_supported_monomial_order<libMesh::PRISM6>(key.order, dispatcher);
    case libMesh::PRISM15:
      return dispatch_supported_monomial_order<libMesh::PRISM15>(key.order, dispatcher);
    case libMesh::PRISM18:
      return dispatch_supported_monomial_order<libMesh::PRISM18>(key.order, dispatcher);
    case libMesh::PRISM20:
      return dispatch_supported_monomial_order<libMesh::PRISM20>(key.order, dispatcher);
    case libMesh::PRISM21:
      return dispatch_supported_monomial_order<libMesh::PRISM21>(key.order, dispatcher);
    case libMesh::PYRAMID5:
      return dispatch_supported_monomial_order<libMesh::PYRAMID5>(key.order, dispatcher);
    case libMesh::PYRAMID13:
      return dispatch_supported_monomial_order<libMesh::PYRAMID13>(key.order, dispatcher);
    case libMesh::PYRAMID14:
      return dispatch_supported_monomial_order<libMesh::PYRAMID14>(key.order, dispatcher);
    case libMesh::PYRAMID18:
      return dispatch_supported_monomial_order<libMesh::PYRAMID18>(key.order, dispatcher);
    default:
      return dispatcher.unsupported_key(key);
  }
}

template <typename Dispatcher>
inline int
dispatch_exact_monomial_shape_key_with_map(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  switch (key.elem_type)
  {
    case libMesh::EDGE2:
      return dispatch_supported_monomial_order<libMesh::EDGE2>(key.order, dispatcher);
    case libMesh::EDGE3:
      return dispatch_supported_monomial_order<libMesh::EDGE3>(key.order, dispatcher);
    case libMesh::EDGE4:
      return dispatch_supported_monomial_order<libMesh::EDGE4>(key.order, dispatcher);
    case libMesh::TRI3:
      return dispatch_supported_monomial_order<libMesh::TRI3>(key.order, dispatcher);
    case libMesh::TRI6:
      return dispatch_supported_monomial_order<libMesh::TRI6>(key.order, dispatcher);
    case libMesh::QUAD4:
      return dispatch_supported_monomial_order<libMesh::QUAD4>(key.order, dispatcher);
    case libMesh::QUAD8:
      return dispatch_supported_monomial_order<libMesh::QUAD8>(key.order, dispatcher);
    case libMesh::QUAD9:
      return dispatch_supported_monomial_order<libMesh::QUAD9>(key.order, dispatcher);
    case libMesh::TET4:
      return dispatch_supported_monomial_order<libMesh::TET4>(key.order, dispatcher);
    case libMesh::TET10:
      return dispatch_supported_monomial_order<libMesh::TET10>(key.order, dispatcher);
    case libMesh::HEX8:
      return dispatch_supported_monomial_order<libMesh::HEX8>(key.order, dispatcher);
    case libMesh::HEX20:
      return dispatch_supported_monomial_order<libMesh::HEX20>(key.order, dispatcher);
    case libMesh::HEX27:
      return dispatch_supported_monomial_order<libMesh::HEX27>(key.order, dispatcher);
    default:
      return dispatcher.unsupported_key(key);
  }
}

template <typename Dispatcher>
inline int
dispatch_exact_shape_key(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return dispatch_exact_lagrange_shape_key(key, dispatcher);

    case libMesh::MONOMIAL:
      return dispatch_exact_monomial_shape_key(key, dispatcher);

    default:
      return dispatcher.unsupported_key(key);
  }
}

template <typename Dispatcher>
inline int
dispatch_supported_shape_key(libMesh::FEShapeKey key, const Dispatcher & dispatcher)
{
  if (!libMesh::supports_shape(key))
    return dispatcher.unsupported_key(key);

  return dispatch_exact_shape_key(key, dispatcher);
}

inline bool
is_supported_lagrange_map_topology(libMesh::ElemType topo)
{
  return libMesh::supports_lagrange_map_topology(topo);
}

inline bool
supports_shape_key_with_lagrange_map(libMesh::FEShapeKey key)
{
  return libMesh::supports_shape_with_lagrange_map(key);
}

template <typename Dispatcher>
inline int
dispatch_supported_shape_key_with_lagrange_map(libMesh::FEShapeKey key,
                                               const Dispatcher & dispatcher)
{
  if (!supports_shape_key_with_lagrange_map(key))
    return dispatcher.unsupported_key(key);

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return dispatch_exact_lagrange_shape_key_with_map(key, dispatcher);

    case libMesh::MONOMIAL:
      return dispatch_exact_monomial_shape_key_with_map(key, dispatcher);

    default:
      return dispatcher.unsupported_key(key);
  }
}

inline bool
is_supported_lagrange_face_map_topology(libMesh::ElemType topo)
{
  return libMesh::supports_lagrange_face_map_topology(topo);
}

template <typename Dispatcher>
inline int
dispatch_supported_lagrange_map_topology(libMesh::ElemType topo,
                                         const Dispatcher & dispatcher)
{
  return libMesh::dispatch_lagrange_map_topology_or(
    topo,
    dispatcher,
    [&](libMesh::ElemType unsupported) { return dispatcher.unsupported_topology(unsupported); });
}

template <typename Dispatcher>
inline int
dispatch_supported_lagrange_face_map_topology(libMesh::ElemType topo,
                                              const Dispatcher & dispatcher)
{
  if (!is_supported_lagrange_face_map_topology(topo))
    return dispatcher.unsupported_topology(topo);

  return dispatch_supported_lagrange_map_topology(topo, dispatcher);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_SHAPE_DISPATCH_H
