#ifndef LIBMESH_TESTS_FE_KOKKOS_FE_ORACLE_TEST_UTILS_H
#define LIBMESH_TESTS_FE_KOKKOS_FE_ORACLE_TEST_UTILS_H

#include "gpu/kokkos_fe_evaluator.h"
#include "gpu/kokkos_fe_face_map.h"
#include "gpu/kokkos_fe_map.h"
#include "gpu/kokkos_fe_types.h"

#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_map.h"
#include "libmesh/node.h"
#include "libmesh/quadrature_gauss.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

// This header is intended for the standalone Kokkos test executables, which
// include Kokkos before pulling in these helpers.

namespace kokkos_test_utils
{

LIBMESH_DEVICE_INLINE libMesh::Real
vector_component(const libMesh::Kokkos::RealVector & v, unsigned int component)
{
  switch (component)
  {
    case 0:
      return v(0);
    case 1:
#if LIBMESH_DIM > 1
      return v(1);
#else
      return 0.0;
#endif
    case 2:
#if LIBMESH_DIM > 2
      return v(2);
#else
      return 0.0;
#endif
    default:
      return 0.0;
  }
}

LIBMESH_DEVICE_INLINE libMesh::Real
tensor_component(const libMesh::Kokkos::RealTensor & T, unsigned int i, unsigned int j)
{
#if LIBMESH_DIM > 2
  return T(i, j);
#elif LIBMESH_DIM > 1
  if (i < 2 && j < 2)
    return T(i, j);
  return 0.0;
#else
  if (i == 0 && j == 0)
    return T(0, 0);
  return 0.0;
#endif
}

struct element_fixture
{
  std::unique_ptr<libMesh::Elem>              elem;
  std::vector<std::unique_ptr<libMesh::Node>> nodes;
};

struct map_helper_context
{
  std::vector<libMesh::Real> ref_values;
  Kokkos::View<double *>     d_coords;
  Kokkos::View<double *>     d_xi;
  Kokkos::View<double *>     d_eta;
  Kokkos::View<double *>     d_zeta;
  Kokkos::View<double *>     d_w;
  unsigned int               nqp;
  unsigned int               dim;
  unsigned int               n_nodes;
};

struct face_helper_context
{
  std::vector<libMesh::Real> ref_values;
  Kokkos::View<double *>     d_face_coords;
  Kokkos::View<double *>     d_parent_coords;
  Kokkos::View<double *>     d_xi;
  Kokkos::View<double *>     d_eta;
  Kokkos::View<double *>     d_zeta;
  Kokkos::View<double *>     d_w;
  Kokkos::View<double *>     d_parent_xi;
  Kokkos::View<double *>     d_parent_eta;
  Kokkos::View<double *>     d_parent_zeta;
  unsigned int               nqp;
  unsigned int               parent_dim;
  unsigned int               n_parent_nodes;
  unsigned int               n_face_nodes;
};

template <libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct lagrange_evaluator_topology
{
  static const libMesh::ElemType value = libMesh::INVALID_ELEM;
};

#define KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(exact_topo, exact_order, evaluator_topo)      \
  template <>                                                                                    \
  struct lagrange_evaluator_topology<libMesh::exact_topo, libMesh::exact_order>                 \
  {                                                                                              \
    static const libMesh::ElemType value = libMesh::evaluator_topo;                             \
  }

KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(EDGE2, FIRST, EDGE2);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(EDGE3, FIRST, EDGE2);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(EDGE3, SECOND, EDGE3);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(EDGE4, FIRST, EDGE2);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TRI3, FIRST, TRI3);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TRI6, FIRST, TRI3);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TRI6, SECOND, TRI6);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TRI7, FIRST, TRI3);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TRI7, SECOND, TRI6);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(QUAD4, FIRST, QUAD4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(QUAD8, FIRST, QUAD4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(QUAD8, SECOND, QUAD8);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(QUAD9, FIRST, QUAD4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(QUAD9, SECOND, QUAD9);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TET4, FIRST, TET4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TET10, FIRST, TET4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TET10, SECOND, TET10);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TET14, FIRST, TET4);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(TET14, SECOND, TET10);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(HEX8, FIRST, HEX8);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(HEX20, FIRST, HEX8);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(HEX20, SECOND, HEX20);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(HEX27, FIRST, HEX8);
KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE(HEX27, SECOND, HEX27);

#undef KOKKOS_TEST_UTILS_LAGRANGE_TOPOLOGY_CASE

template <libMesh::ElemType ExactTopo>
struct monomial_evaluator_dim
{
  static const unsigned int value = 0;
};

#define KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(exact_topo, dim_value)                              \
  template <>                                                                                    \
  struct monomial_evaluator_dim<libMesh::exact_topo>                                            \
  {                                                                                              \
    static const unsigned int value = dim_value;                                                \
  }

KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(EDGE2, 1);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(EDGE3, 1);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(EDGE4, 1);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TRI3, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TRI6, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TRI7, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(QUAD4, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(QUAD8, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(QUAD9, 2);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TET4, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TET10, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(TET14, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(HEX8, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(HEX20, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(HEX27, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PRISM6, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PRISM15, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PRISM18, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PRISM20, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PRISM21, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PYRAMID5, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PYRAMID13, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PYRAMID14, 3);
KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE(PYRAMID18, 3);

#undef KOKKOS_TEST_UTILS_MONOMIAL_DIM_CASE

template <unsigned int Dim, libMesh::Order ExactOrder>
struct monomial_order_evaluator;

#define KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(dim_value, exact_order, impl_suffix, impl_order)  \
  template <>                                                                                    \
  struct monomial_order_evaluator<dim_value, libMesh::exact_order>                              \
  {                                                                                              \
    LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,                             \
                                                      libMesh::Real xi,                          \
                                                      libMesh::Real eta,                         \
                                                      libMesh::Real zeta)                        \
    {                                                                                            \
      return libMesh::Kokkos::impl_suffix<impl_order>::shape(i, xi, eta, zeta);                 \
    }                                                                                            \
                                                                                                 \
    LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(unsigned int i,         \
                                                                        libMesh::Real xi,        \
                                                                        libMesh::Real eta,       \
                                                                        libMesh::Real zeta)      \
    {                                                                                            \
      return libMesh::Kokkos::impl_suffix<impl_order>::grad_shape(i, xi, eta, zeta);            \
    }                                                                                            \
  }

KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, CONSTANT, MonomialImpl1D, 0);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, FIRST, MonomialImpl1D, 1);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, SECOND, MonomialImpl1D, 2);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, THIRD, MonomialImpl1D, 3);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, FOURTH, MonomialImpl1D, 4);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(1, FIFTH, MonomialImpl1D, 5);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, CONSTANT, MonomialImpl2D, 0);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, FIRST, MonomialImpl2D, 1);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, SECOND, MonomialImpl2D, 2);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, THIRD, MonomialImpl2D, 3);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, FOURTH, MonomialImpl2D, 4);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(2, FIFTH, MonomialImpl2D, 5);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, CONSTANT, MonomialImpl3D, 0);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, FIRST, MonomialImpl3D, 1);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, SECOND, MonomialImpl3D, 2);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, THIRD, MonomialImpl3D, 3);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, FOURTH, MonomialImpl3D, 4);
KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE(3, FIFTH, MonomialImpl3D, 5);

#undef KOKKOS_TEST_UTILS_MONOMIAL_ORDER_CASE

template <libMesh::FEFamily Family, libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator;

template <libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator<libMesh::LAGRANGE, ExactTopo, ExactOrder>
{
  LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,
                                                   libMesh::Real xi,
                                                   libMesh::Real eta,
                                                   libMesh::Real zeta)
  {
    return libMesh::Kokkos::map_shape<libMesh::LAGRANGE,
                                      lagrange_evaluator_topology<ExactTopo, ExactOrder>::value>(
      i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(unsigned int i,
                                                                      libMesh::Real xi,
                                                                      libMesh::Real eta,
                                                                      libMesh::Real zeta)
  {
    return libMesh::Kokkos::grad_map_shape<libMesh::LAGRANGE,
                                           lagrange_evaluator_topology<ExactTopo, ExactOrder>::value>(
      i, xi, eta, zeta);
  }
};

template <libMesh::ElemType ExactTopo, libMesh::Order ExactOrder>
struct exact_shape_evaluator<libMesh::MONOMIAL, ExactTopo, ExactOrder>
{
  LIBMESH_DEVICE_INLINE static libMesh::Real shape(unsigned int i,
                                                   libMesh::Real xi,
                                                   libMesh::Real eta,
                                                   libMesh::Real zeta)
  {
    return monomial_order_evaluator<monomial_evaluator_dim<ExactTopo>::value, ExactOrder>::shape(
      i, xi, eta, zeta);
  }

  LIBMESH_DEVICE_INLINE static libMesh::Kokkos::RealVector grad_shape(unsigned int i,
                                                                      libMesh::Real xi,
                                                                      libMesh::Real eta,
                                                                      libMesh::Real zeta)
  {
    return monomial_order_evaluator<monomial_evaluator_dim<ExactTopo>::value, ExactOrder>::grad_shape(
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
      return dispatcher.unsupported_key(libMesh::Kokkos::FEShapeKey{ libMesh::MONOMIAL, ExactTopo, order });
  }
}

inline bool
is_supported_lagrange_map_topology(libMesh::ElemType topo);

template <typename Dispatcher>
inline int
dispatch_supported_lagrange_shape_key(libMesh::Kokkos::FEShapeKey key, const Dispatcher & dispatcher)
{
  if (key.family != libMesh::LAGRANGE || !libMesh::Kokkos::supports_shape(key))
    return dispatcher.unsupported_key(key);

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
dispatch_supported_lagrange_shape_key_with_map(libMesh::Kokkos::FEShapeKey key,
                                               const Dispatcher & dispatcher)
{
  if (key.family != libMesh::LAGRANGE ||
      !libMesh::Kokkos::supports_shape(key) ||
      !is_supported_lagrange_map_topology(key.elem_type))
    return dispatcher.unsupported_key(key);

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
dispatch_supported_shape_key(libMesh::Kokkos::FEShapeKey key, const Dispatcher & dispatcher)
{
  if (!libMesh::Kokkos::supports_shape(key))
    return dispatcher.unsupported_key(key);

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return dispatch_supported_lagrange_shape_key(key, dispatcher);

    case libMesh::MONOMIAL:
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

    default:
      return dispatcher.unsupported_key(key);
  }
}

inline bool
is_supported_lagrange_map_topology(libMesh::ElemType topo)
{
  switch (topo)
  {
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::TRI3:
    case libMesh::TRI6:
    case libMesh::QUAD4:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
    case libMesh::TET4:
    case libMesh::TET10:
    case libMesh::HEX8:
    case libMesh::HEX20:
    case libMesh::HEX27:
      return true;

    default:
      return false;
  }
}

inline bool
supports_shape_key_with_lagrange_map(libMesh::Kokkos::FEShapeKey key)
{
  return libMesh::Kokkos::supports_shape(key) &&
         is_supported_lagrange_map_topology(key.elem_type);
}

template <typename Dispatcher>
inline int
dispatch_supported_shape_key_with_lagrange_map(libMesh::Kokkos::FEShapeKey key,
                                               const Dispatcher & dispatcher)
{
  if (!supports_shape_key_with_lagrange_map(key))
    return dispatcher.unsupported_key(key);

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return dispatch_supported_lagrange_shape_key_with_map(key, dispatcher);

    case libMesh::MONOMIAL:
      switch (key.elem_type)
      {
        case libMesh::EDGE2:
          return dispatch_supported_monomial_order<libMesh::EDGE2>(key.order, dispatcher);
        case libMesh::EDGE3:
          return dispatch_supported_monomial_order<libMesh::EDGE3>(key.order, dispatcher);
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

    default:
      return dispatcher.unsupported_key(key);
  }
}

inline bool
is_supported_lagrange_face_map_topology(libMesh::ElemType topo)
{
  switch (topo)
  {
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::TRI3:
    case libMesh::TRI6:
    case libMesh::QUAD4:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
      return true;

    default:
      return false;
  }
}

template <typename Dispatcher>
inline int
dispatch_supported_lagrange_map_topology(libMesh::ElemType topo,
                                         const Dispatcher & dispatcher)
{
  switch (topo)
  {
    case libMesh::EDGE2:
      return dispatcher.template operator()<libMesh::EDGE2>();
    case libMesh::EDGE3:
      return dispatcher.template operator()<libMesh::EDGE3>();
    case libMesh::TRI3:
      return dispatcher.template operator()<libMesh::TRI3>();
    case libMesh::TRI6:
      return dispatcher.template operator()<libMesh::TRI6>();
    case libMesh::QUAD4:
      return dispatcher.template operator()<libMesh::QUAD4>();
    case libMesh::QUAD8:
      return dispatcher.template operator()<libMesh::QUAD8>();
    case libMesh::QUAD9:
      return dispatcher.template operator()<libMesh::QUAD9>();
    case libMesh::TET4:
      return dispatcher.template operator()<libMesh::TET4>();
    case libMesh::TET10:
      return dispatcher.template operator()<libMesh::TET10>();
    case libMesh::HEX8:
      return dispatcher.template operator()<libMesh::HEX8>();
    case libMesh::HEX20:
      return dispatcher.template operator()<libMesh::HEX20>();
    case libMesh::HEX27:
      return dispatcher.template operator()<libMesh::HEX27>();

    default:
      return dispatcher.unsupported_topology(topo);
  }
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

inline int
compare_device_values(const Kokkos::View<double *> & d_values,
                      const std::vector<libMesh::Real> & ref_values,
                      double tol = 1.0e-13)
{
  auto h_values = Kokkos::create_mirror_view(d_values);
  Kokkos::deep_copy(h_values, d_values);

  int fail = 0;
  for (std::size_t i = 0; i < ref_values.size(); ++i)
    if (std::fabs(h_values(i) - ref_values[i]) > tol)
      ++fail;

  return fail;
}

inline std::unique_ptr<libMesh::Elem>
build_reference_elem(libMesh::ElemType elem_type)
{
  auto elem = libMesh::Elem::build(elem_type);
  elem->set_mapping_type(libMesh::LAGRANGE_MAP);
  return elem;
}

inline unsigned int
build_qps(libMesh::ElemType elem_type,
          unsigned int dim,
          unsigned int quadrature_order,
          std::vector<libMesh::Real> & xi_h,
          std::vector<libMesh::Real> & eta_h,
          std::vector<libMesh::Real> & zeta_h)
{
  libMesh::QGauss qr(dim, static_cast<libMesh::Order>(quadrature_order));
  qr.allow_rules_with_negative_weights = true;
  qr.init(elem_type);

  const unsigned int nqp = qr.n_points();
  xi_h.resize(nqp);
  eta_h.resize(nqp);
  zeta_h.resize(nqp);

  for (unsigned int q = 0; q < nqp; ++q)
  {
    xi_h[q] = qr.qp(q)(0);
    eta_h[q] = (dim >= 2) ? qr.qp(q)(1) : libMesh::Real(0);
    zeta_h[q] = (dim >= 3) ? qr.qp(q)(2) : libMesh::Real(0);
  }

  return nqp;
}

inline unsigned int
build_qps(libMesh::ElemType elem_type,
          unsigned int dim,
          std::vector<libMesh::Real> & xi_h,
          std::vector<libMesh::Real> & eta_h,
          std::vector<libMesh::Real> & zeta_h)
{
  return build_qps(elem_type, dim, /*quadrature_order=*/4, xi_h, eta_h, zeta_h);
}

inline unsigned int
build_host_qgauss(libMesh::ElemType topo,
                  unsigned int dim,
                  unsigned int order,
                  std::vector<libMesh::Real> & x_ref,
                  std::vector<libMesh::Real> & y_ref,
                  std::vector<libMesh::Real> & z_ref,
                  std::vector<libMesh::Real> & w_ref)
{
  libMesh::QGauss qr(dim, static_cast<libMesh::Order>(order));
  qr.allow_rules_with_negative_weights = true;
  qr.init(topo);

  const unsigned int nqp = qr.n_points();
  x_ref.resize(nqp);
  y_ref.resize(nqp);
  z_ref.resize(nqp);
  w_ref.resize(nqp);

  for (unsigned int q = 0; q < nqp; ++q)
  {
    x_ref[q] = qr.qp(q)(0);
    y_ref[q] = (dim >= 2) ? qr.qp(q)(1) : libMesh::Real(0);
    z_ref[q] = (dim >= 3) ? qr.qp(q)(2) : libMesh::Real(0);
    w_ref[q] = qr.w(q);
  }

  return nqp;
}

inline Kokkos::View<double *>
upload_real(const std::vector<libMesh::Real> & values, const char * label)
{
  Kokkos::View<double *> d(std::string(label), values.size());
  auto h = Kokkos::create_mirror_view(d);
  for (std::size_t i = 0; i < values.size(); ++i)
    h(i) = values[i];
  Kokkos::deep_copy(d, h);
  return d;
}

inline Kokkos::View<double *>
upload_point_coordinates(const libMesh::Elem & elem, const char * label)
{
  Kokkos::View<double *> d(std::string(label), 3 * elem.n_nodes());
  auto h = Kokkos::create_mirror_view(d);
  for (unsigned int i = 0; i < elem.n_nodes(); ++i)
  {
    h(3 * i + 0) = elem.point(i)(0);
    h(3 * i + 1) = elem.point(i)(1);
    h(3 * i + 2) = elem.point(i)(2);
  }
  Kokkos::deep_copy(d, h);
  return d;
}

inline std::string
make_label(const char * prefix, const char * suffix)
{
  return std::string(prefix) + suffix;
}

inline element_fixture
build_reference_fixture(libMesh::ElemType elem_type)
{
  element_fixture fixture;
  fixture.elem = build_reference_elem(elem_type);
  fixture.nodes.reserve(fixture.elem->n_nodes());

  const unsigned int dim = fixture.elem->dim();

  for (unsigned int i = 0; i < fixture.elem->n_nodes(); ++i)
  {
    const libMesh::Point master = fixture.elem->master_point(i);
    const libMesh::Real xi = master(0);
    const libMesh::Real eta = master(1);
    const libMesh::Real zeta = master(2);

    libMesh::Point xyz;
    switch (dim)
    {
      case 1:
        xyz = libMesh::Point(
          0.7 + 0.8 * xi + 0.06 * xi * xi,
          -0.3 + 0.25 * xi + 0.04 * xi * xi,
          0.2 + 0.1 * xi - 0.03 * xi * xi);
        break;

      case 2:
        xyz = libMesh::Point(
          0.4 + 0.9 * xi + 0.15 * eta + 0.04 * xi * eta + 0.03 * eta * eta,
          -0.2 + 0.2 * xi + 0.85 * eta + 0.05 * xi * xi + 0.03 * xi * eta,
          0.1 + 0.12 * xi - 0.08 * eta + 0.02 * xi * eta);
        break;

      case 3:
        xyz = libMesh::Point(
          0.3 + 0.9 * xi + 0.12 * eta + 0.08 * zeta + 0.03 * xi * eta + 0.02 * zeta * zeta,
          -0.1 + 0.18 * xi + 0.8 * eta + 0.11 * zeta + 0.02 * eta * zeta,
          0.2 + 0.10 * xi + 0.14 * eta + 0.85 * zeta + 0.02 * xi * zeta + 0.01 * xi * eta);
        break;

      default:
        xyz = libMesh::Point();
        break;
    }

    fixture.nodes.push_back(libMesh::Node::build(xyz(0), xyz(1), xyz(2), i));
    fixture.elem->set_node(i, fixture.nodes.back().get());
  }

  return fixture;
}

inline element_fixture
build_flat_reference_fixture(libMesh::ElemType elem_type)
{
  element_fixture fixture;
  fixture.elem = build_reference_elem(elem_type);
  fixture.nodes.reserve(fixture.elem->n_nodes());

  const unsigned int dim = fixture.elem->dim();

  for (unsigned int i = 0; i < fixture.elem->n_nodes(); ++i)
  {
    const libMesh::Point master = fixture.elem->master_point(i);
    const libMesh::Real xi = master(0);
    const libMesh::Real eta = master(1);
    const libMesh::Real zeta = master(2);

    libMesh::Point xyz;
    switch (dim)
    {
      case 1:
        xyz = libMesh::Point(0.7 + 0.8 * xi + 0.06 * xi * xi,
                             0.0,
                             0.0);
        break;

      case 2:
        xyz = libMesh::Point(0.4 + 0.9 * xi + 0.15 * eta + 0.04 * xi * eta + 0.03 * eta * eta,
                             -0.2 + 0.2 * xi + 0.85 * eta + 0.05 * xi * xi + 0.03 * xi * eta,
                             0.0);
        break;

      case 3:
        xyz = libMesh::Point(
          0.3 + 0.9 * xi + 0.12 * eta + 0.08 * zeta + 0.03 * xi * eta + 0.02 * zeta * zeta,
          -0.1 + 0.18 * xi + 0.8 * eta + 0.11 * zeta + 0.02 * eta * zeta,
          0.2 + 0.10 * xi + 0.14 * eta + 0.85 * zeta + 0.02 * xi * zeta + 0.01 * xi * eta);
        break;

      default:
        xyz = libMesh::Point();
        break;
    }

    fixture.nodes.push_back(libMesh::Node::build(xyz(0), xyz(1), xyz(2), i));
    fixture.elem->set_node(i, fixture.nodes.back().get());
  }

  return fixture;
}

inline element_fixture
build_permuted_reference_fixture(libMesh::ElemType elem_type,
                                 unsigned int perm_num)
{
  element_fixture fixture = build_reference_fixture(elem_type);
  fixture.elem->permute(perm_num);
  return fixture;
}

inline map_helper_context
build_map_helper_context(const element_fixture & fixture,
                         libMesh::ElemType topo,
                         const char * label_prefix)
{
  map_helper_context context;

  const unsigned int dim = fixture.elem->dim();
  const unsigned int n_nodes = fixture.elem->n_nodes();
  const libMesh::FEType fe_type(fixture.elem->default_order(),
                                libMesh::FEMap::map_fe_type(*fixture.elem));
  auto fe = libMesh::FEBase::build(dim, fe_type);

  libMesh::QGauss qr(dim, libMesh::FOURTH);
  qr.allow_rules_with_negative_weights = true;
  qr.init(topo);

  fe->attach_quadrature_rule(&qr);
  fe->get_xyz();
  fe->get_dxyzdxi();
  if (dim >= 2)
    fe->get_dxyzdeta();
  if (dim >= 3)
    fe->get_dxyzdzeta();
  fe->get_JxW();
  fe->reinit(fixture.elem.get());

  const unsigned int nqp = qr.n_points();
  const auto & xyz = fe->get_xyz();
  const auto & dxyzdxi = fe->get_dxyzdxi();
  const auto & jxw = fe->get_JxW();

  context.ref_values.resize(13 * nqp);
  std::vector<libMesh::Real> xi_h(nqp), eta_h(nqp), zeta_h(nqp), w_h(nqp);
  for (unsigned int q = 0; q < nqp; ++q)
  {
    libMesh::RealGradient dxyzdeta(0.0);
    libMesh::RealGradient dxyzdzeta(0.0);
    if (dim >= 2)
      dxyzdeta = fe->get_dxyzdeta()[q];
    if (dim >= 3)
      dxyzdzeta = fe->get_dxyzdzeta()[q];

    const unsigned int base = 13 * q;
    context.ref_values[base + 0] = xyz[q](0);
    context.ref_values[base + 1] = xyz[q](1);
    context.ref_values[base + 2] = xyz[q](2);
    context.ref_values[base + 3] = dxyzdxi[q](0);
    context.ref_values[base + 4] = dxyzdxi[q](1);
    context.ref_values[base + 5] = dxyzdxi[q](2);
    context.ref_values[base + 6] = dxyzdeta(0);
    context.ref_values[base + 7] = dxyzdeta(1);
    context.ref_values[base + 8] = dxyzdeta(2);
    context.ref_values[base + 9] = dxyzdzeta(0);
    context.ref_values[base + 10] = dxyzdzeta(1);
    context.ref_values[base + 11] = dxyzdzeta(2);
    context.ref_values[base + 12] = jxw[q];

    xi_h[q] = qr.qp(q)(0);
    eta_h[q] = (dim >= 2) ? qr.qp(q)(1) : libMesh::Real(0);
    zeta_h[q] = (dim >= 3) ? qr.qp(q)(2) : libMesh::Real(0);
    w_h[q] = qr.w(q);
  }

  context.d_coords = upload_point_coordinates(*fixture.elem, make_label(label_prefix, "_coords").c_str());
  context.d_xi = upload_real(xi_h, make_label(label_prefix, "_xi").c_str());
  context.d_eta = upload_real(eta_h, make_label(label_prefix, "_eta").c_str());
  context.d_zeta = upload_real(zeta_h, make_label(label_prefix, "_zeta").c_str());
  context.d_w = upload_real(w_h, make_label(label_prefix, "_w").c_str());
  context.nqp = nqp;
  context.dim = dim;
  context.n_nodes = n_nodes;

  return context;
}

template <libMesh::ElemType Topo>
inline int
evaluate_map_helper_context(const map_helper_context & context,
                            const char * result_label,
                            double tol = 1.0e-13)
{
  constexpr unsigned int max_nodes = 27;

  Kokkos::View<double *> d_results(std::string(result_label), context.ref_values.size());
  const auto d_coords = context.d_coords;
  const auto d_xi = context.d_xi;
  const auto d_eta = context.d_eta;
  const auto d_zeta = context.d_zeta;
  const auto d_w = context.d_w;
  const unsigned int dim_ = context.dim;
  const unsigned int n_nodes_ = context.n_nodes;

  Kokkos::parallel_for(
    context.nqp,
    KOKKOS_LAMBDA(int q) {
      libMesh::Kokkos::RealVector nodes[max_nodes];
      for (unsigned int i = 0; i < n_nodes_; ++i)
        nodes[i] = libMesh::Kokkos::make_vector(
          d_coords(3 * i + 0), d_coords(3 * i + 1), d_coords(3 * i + 2));

      libMesh::Kokkos::RealVector xyz;
      libMesh::Kokkos::RealTensor J;
      libMesh::Kokkos::physical_point_and_jacobian<libMesh::LAGRANGE, Topo>(
        nodes, n_nodes_, d_xi(q), d_eta(q), d_zeta(q), xyz, J);

      const libMesh::Real jxw_q = libMesh::Kokkos::volume_jxw(J, dim_, d_w(q));
      const unsigned int base = 13 * static_cast<unsigned int>(q);

      d_results(base + 0) = vector_component(xyz, 0);
      d_results(base + 1) = vector_component(xyz, 1);
      d_results(base + 2) = vector_component(xyz, 2);
      d_results(base + 3) = tensor_component(J, 0, 0);
      d_results(base + 4) = tensor_component(J, 0, 1);
      d_results(base + 5) = tensor_component(J, 0, 2);
      d_results(base + 6) = tensor_component(J, 1, 0);
      d_results(base + 7) = tensor_component(J, 1, 1);
      d_results(base + 8) = tensor_component(J, 1, 2);
      d_results(base + 9) = tensor_component(J, 2, 0);
      d_results(base + 10) = tensor_component(J, 2, 1);
      d_results(base + 11) = tensor_component(J, 2, 2);
      d_results(base + 12) = jxw_q;
    });
  Kokkos::fence();

  return compare_device_values(d_results, context.ref_values, tol);
}

inline face_helper_context
build_face_helper_context(const element_fixture & fixture,
                          const libMesh::Elem & side,
                          unsigned int side_id,
                          const char * label_prefix)
{
  face_helper_context context;
  const unsigned int parent_dim = fixture.elem->dim();
  const libMesh::FEType fe_type(fixture.elem->default_order(),
                                libMesh::FEMap::map_fe_type(*fixture.elem));
  const unsigned int side_dim = side.dim();
  auto side_fe = libMesh::FEBase::build(parent_dim, fe_type);

  libMesh::QGauss qr(parent_dim - 1, libMesh::FOURTH);
  qr.allow_rules_with_negative_weights = true;
  qr.init(side.type());

  side_fe->attach_quadrature_rule(&qr);
  side_fe->get_JxW();
  side_fe->get_normals();
  side_fe->get_dxyzdxi();
  if (parent_dim >= 3)
    side_fe->get_dxyzdeta();
  side_fe->reinit(fixture.elem.get(), side_id);

  const unsigned int nqp = qr.n_points();
  const unsigned int n_parent_nodes = fixture.elem->n_nodes();
  const unsigned int n_face_nodes = side.n_nodes();

  std::vector<libMesh::Point> side_ref_points(nqp);
  for (unsigned int q = 0; q < nqp; ++q)
    side_ref_points[q] = qr.qp(q);

  std::vector<libMesh::Point> parent_ref_points;
  if (parent_dim == 2)
  {
    auto side_map_fe = libMesh::FEBase::build(parent_dim, fe_type);
    side_map_fe->get_xyz();
    side_map_fe->side_map(fixture.elem.get(), &side, side_id, side_ref_points, parent_ref_points);
  }

  context.ref_values.resize(13 * nqp);
  std::vector<libMesh::Real> xi_h(nqp), eta_h(nqp), zeta_h(nqp), w_h(nqp);
  std::vector<libMesh::Real> parent_xi_h(nqp, 0.0), parent_eta_h(nqp, 0.0), parent_zeta_h(nqp, 0.0);
  for (unsigned int q = 0; q < nqp; ++q)
  {
    const libMesh::Point row0 = libMesh::FEMap::map_deriv(side_dim, &side, 0, side_ref_points[q]);
    libMesh::Point row1(0.0);
    if (side_dim >= 2)
      row1 = libMesh::FEMap::map_deriv(side_dim, &side, 1, side_ref_points[q]);
    const auto & normal = side_fe->get_normals()[q];
    const unsigned int base = 13 * q;

    context.ref_values[base + 0] = row0(0);
    context.ref_values[base + 1] = row0(1);
    context.ref_values[base + 2] = row0(2);
    context.ref_values[base + 3] = row1(0);
    context.ref_values[base + 4] = row1(1);
    context.ref_values[base + 5] = row1(2);
    context.ref_values[base + 6] = 0.0;
    context.ref_values[base + 7] = 0.0;
    context.ref_values[base + 8] = 0.0;
    context.ref_values[base + 9] = side_fe->get_JxW()[q];
    context.ref_values[base + 10] = normal(0);
    context.ref_values[base + 11] = normal(1);
    context.ref_values[base + 12] = normal(2);

    xi_h[q] = qr.qp(q)(0);
    eta_h[q] = (parent_dim >= 3) ? qr.qp(q)(1) : libMesh::Real(0);
    zeta_h[q] = 0.0;
    w_h[q] = qr.w(q);

    if (parent_dim == 2)
    {
      parent_xi_h[q] = parent_ref_points[q](0);
      parent_eta_h[q] = parent_ref_points[q](1);
      parent_zeta_h[q] = parent_ref_points[q](2);
    }
  }

  context.d_face_coords = upload_point_coordinates(side, make_label(label_prefix, "_coords").c_str());
  context.d_parent_coords = upload_point_coordinates(*fixture.elem, make_label(label_prefix, "_parent_coords").c_str());
  context.d_xi = upload_real(xi_h, make_label(label_prefix, "_xi").c_str());
  context.d_eta = upload_real(eta_h, make_label(label_prefix, "_eta").c_str());
  context.d_zeta = upload_real(zeta_h, make_label(label_prefix, "_zeta").c_str());
  context.d_w = upload_real(w_h, make_label(label_prefix, "_w").c_str());
  context.d_parent_xi = upload_real(parent_xi_h, make_label(label_prefix, "_parent_xi").c_str());
  context.d_parent_eta = upload_real(parent_eta_h, make_label(label_prefix, "_parent_eta").c_str());
  context.d_parent_zeta = upload_real(parent_zeta_h, make_label(label_prefix, "_parent_zeta").c_str());
  context.nqp = nqp;
  context.parent_dim = parent_dim;
  context.n_parent_nodes = n_parent_nodes;
  context.n_face_nodes = n_face_nodes;

  return context;
}

template <libMesh::ElemType ParentTopo, libMesh::ElemType SideTopo>
inline int
evaluate_face_helper_context_2d(const face_helper_context & context,
                                const char * result_label,
                                double tol = 1.0e-13)
{
  constexpr unsigned int max_face_nodes = 9;
  constexpr unsigned int max_parent_nodes = 27;

  Kokkos::View<double *> d_results(std::string(result_label), context.ref_values.size());
  const auto d_face_coords = context.d_face_coords;
  const auto d_parent_coords = context.d_parent_coords;
  const auto d_xi = context.d_xi;
  const auto d_eta = context.d_eta;
  const auto d_zeta = context.d_zeta;
  const auto d_w = context.d_w;
  const auto d_parent_xi = context.d_parent_xi;
  const auto d_parent_eta = context.d_parent_eta;
  const auto d_parent_zeta = context.d_parent_zeta;
  const unsigned int n_parent_nodes_ = context.n_parent_nodes;
  const unsigned int n_face_nodes_ = context.n_face_nodes;

  Kokkos::parallel_for(
    context.nqp,
    KOKKOS_LAMBDA(int q) {
      libMesh::Kokkos::RealVector face_nodes[max_face_nodes];
      libMesh::Kokkos::RealVector parent_nodes[max_parent_nodes];
      for (unsigned int i = 0; i < n_face_nodes_; ++i)
        face_nodes[i] = libMesh::Kokkos::make_vector(
          d_face_coords(3 * i + 0), d_face_coords(3 * i + 1), d_face_coords(3 * i + 2));
      for (unsigned int i = 0; i < n_parent_nodes_; ++i)
        parent_nodes[i] = libMesh::Kokkos::make_vector(
          d_parent_coords(3 * i + 0), d_parent_coords(3 * i + 1), d_parent_coords(3 * i + 2));

      const libMesh::Kokkos::RealTensor J = libMesh::Kokkos::face_jacobian<libMesh::LAGRANGE, SideTopo>(
        face_nodes, n_face_nodes_, d_xi(q), d_eta(q), d_zeta(q));
      const libMesh::Kokkos::RealTensor parent_J = libMesh::Kokkos::jacobian<libMesh::LAGRANGE, ParentTopo>(
        parent_nodes, n_parent_nodes_, d_parent_xi(q), d_parent_eta(q), d_parent_zeta(q));
      const libMesh::Real jxw_q = libMesh::Kokkos::face_jxw(J, /*parent_dim=*/2u, d_w(q));
      const libMesh::Kokkos::RealVector normal_q = libMesh::Kokkos::edge_normal_on_parent_surface(J, parent_J);
      const unsigned int base = 13 * static_cast<unsigned int>(q);

      d_results(base + 0) = tensor_component(J, 0, 0);
      d_results(base + 1) = tensor_component(J, 0, 1);
      d_results(base + 2) = tensor_component(J, 0, 2);
      d_results(base + 3) = tensor_component(J, 1, 0);
      d_results(base + 4) = tensor_component(J, 1, 1);
      d_results(base + 5) = tensor_component(J, 1, 2);
      d_results(base + 6) = tensor_component(J, 2, 0);
      d_results(base + 7) = tensor_component(J, 2, 1);
      d_results(base + 8) = tensor_component(J, 2, 2);
      d_results(base + 9) = jxw_q;
      d_results(base + 10) = vector_component(normal_q, 0);
      d_results(base + 11) = vector_component(normal_q, 1);
      d_results(base + 12) = vector_component(normal_q, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_results, context.ref_values, tol);
}

template <libMesh::ElemType SideTopo>
inline int
evaluate_face_helper_context_3d(const face_helper_context & context,
                                const char * result_label,
                                double tol = 1.0e-13)
{
  constexpr unsigned int max_face_nodes = 9;

  Kokkos::View<double *> d_results(std::string(result_label), context.ref_values.size());
  const auto d_face_coords = context.d_face_coords;
  const auto d_xi = context.d_xi;
  const auto d_eta = context.d_eta;
  const auto d_zeta = context.d_zeta;
  const auto d_w = context.d_w;
  const unsigned int n_face_nodes_ = context.n_face_nodes;

  Kokkos::parallel_for(
    context.nqp,
    KOKKOS_LAMBDA(int q) {
      libMesh::Kokkos::RealVector face_nodes[max_face_nodes];
      for (unsigned int i = 0; i < n_face_nodes_; ++i)
        face_nodes[i] = libMesh::Kokkos::make_vector(
          d_face_coords(3 * i + 0), d_face_coords(3 * i + 1), d_face_coords(3 * i + 2));

      const libMesh::Kokkos::RealTensor J = libMesh::Kokkos::face_jacobian<libMesh::LAGRANGE, SideTopo>(
        face_nodes, n_face_nodes_, d_xi(q), d_eta(q), d_zeta(q));
      const libMesh::Real jxw_q = libMesh::Kokkos::face_jxw(J, /*parent_dim=*/3u, d_w(q));
      const libMesh::Kokkos::RealVector normal_q = libMesh::Kokkos::face_normal(J, /*parent_dim=*/3u);
      const unsigned int base = 13 * static_cast<unsigned int>(q);

      d_results(base + 0) = tensor_component(J, 0, 0);
      d_results(base + 1) = tensor_component(J, 0, 1);
      d_results(base + 2) = tensor_component(J, 0, 2);
      d_results(base + 3) = tensor_component(J, 1, 0);
      d_results(base + 4) = tensor_component(J, 1, 1);
      d_results(base + 5) = tensor_component(J, 1, 2);
      d_results(base + 6) = tensor_component(J, 2, 0);
      d_results(base + 7) = tensor_component(J, 2, 1);
      d_results(base + 8) = tensor_component(J, 2, 2);
      d_results(base + 9) = jxw_q;
      d_results(base + 10) = vector_component(normal_q, 0);
      d_results(base + 11) = vector_component(normal_q, 1);
      d_results(base + 12) = vector_component(normal_q, 2);
    });
  Kokkos::fence();

  return compare_device_values(d_results, context.ref_values, tol);
}

} // namespace kokkos_test_utils

#endif
