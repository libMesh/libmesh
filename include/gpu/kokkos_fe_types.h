// Kokkos FE type helpers.
//
// Defines the FEShapeKey aggregate and device-callable dispatch functions used
// by both host-side assembly setup and device-side evaluation.
//
// Uses libMesh's own ElemType, FEFamily, and FEElemClass enums directly —
// no wrapper enums are needed.

#ifndef LIBMESH_KOKKOS_FE_TYPES_H
#define LIBMESH_KOKKOS_FE_TYPES_H

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_fe_elem_class.h"
#include "libmesh/enum_order.h"
// ElemMappingType (LAGRANGE_MAP, RATIONAL_BERNSTEIN_MAP) is defined in enum_elem_type.h
#include "libmesh/libmesh_device.h"
#ifndef LIBMESH_KOKKOS_COMPILATION
#  include "libmesh/libmesh_common.h"
#endif

namespace libMesh::Kokkos
{

// Bring FEElemClass into this namespace so existing unqualified uses compile.
using libMesh::FEElemClass;

namespace detail
{

LIBMESH_DEVICE_INLINE void
abort_unsupported(const char * msg)
{
#ifdef LIBMESH_KOKKOS_COMPILATION
  ::Kokkos::abort(msg);
#else
  libmesh_error_msg(msg);
#endif
}

} // namespace detail

LIBMESH_DEVICE_INLINE bool
is_monomial_2d_elem_type(libMesh::ElemType elem_type)
{
  switch (elem_type)
  {
    case libMesh::C0POLYGON:
    case libMesh::TRI3:
    case libMesh::TRISHELL3:
    case libMesh::TRI6:
    case libMesh::TRI7:
    case libMesh::QUAD4:
    case libMesh::QUADSHELL4:
    case libMesh::QUAD8:
    case libMesh::QUADSHELL8:
    case libMesh::QUAD9:
    case libMesh::QUADSHELL9:
      return true;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
is_monomial_3d_elem_type(libMesh::ElemType elem_type,
                         bool include_pyramid18 = true)
{
  switch (elem_type)
  {
    case libMesh::TET4:
    case libMesh::TET10:
    case libMesh::TET14:
    case libMesh::HEX8:
    case libMesh::HEX20:
    case libMesh::HEX27:
    case libMesh::PRISM6:
    case libMesh::PRISM15:
    case libMesh::PRISM18:
    case libMesh::PRISM20:
    case libMesh::PRISM21:
    case libMesh::PYRAMID5:
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
    case libMesh::C0POLYHEDRON:
      return true;
    case libMesh::PYRAMID18:
      return include_pyramid18;
    default:
      return false;
  }
}

// ── Shape function space key ──────────────────────────────────────────────────
// Uniquely identifies a libMesh FE space, including the exact element topology.
// This must be exact for LAGRANGE spaces, since libMesh distinguishes e.g.
// QUAD8 from QUAD9 and HEX20 from HEX27 at the same polynomial order.
//
// Trivially copyable; fits in a register (enum + enum + enum, no heap).

struct FEShapeKey
{
  libMesh::FEFamily family;
  libMesh::ElemType elem_type;
  libMesh::Order    order;
};

// ── Device-callable conversion helpers ───────────────────────────────────────

/// Return the Kokkos side topology used for dispatch for any side of parent
/// element type \p parent.
/// This helper is valid only for elements whose side topology is uniform.
/// Mixed-face elements such as prisms and pyramids require side-specific logic.
/// In 1D, libMesh sides are NODEELEM objects; this helper returns EDGE2 as the
/// internal surrogate topology used by the Kokkos map/shape path.
LIBMESH_DEVICE_INLINE libMesh::ElemType
get_side_topology(libMesh::ElemType parent)
{
  switch (parent)
  {
    // 1D: libMesh sides are NodeElem, but Kokkos dispatches them through
    // a degenerate EDGE2 surrogate.
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::EDGE4:
      return libMesh::EDGE2;

    // 2D first-order: sides are linear edges
    case libMesh::TRI3:
    case libMesh::QUAD4:
      return libMesh::EDGE2;

    // 2D second-order: sides are quadratic edges
    case libMesh::TRI6:
    case libMesh::TRI7:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
      return libMesh::EDGE3;

    // 3D first-order: uniform-side-topology elements only
    case libMesh::TET4:
      return libMesh::TRI3;
    case libMesh::HEX8:
      return libMesh::QUAD4;

    // 3D second-order: uniform-side-topology elements only
    case libMesh::TET10:
      return libMesh::TRI6;
    case libMesh::TET14:
      return libMesh::TRI7;
    case libMesh::HEX20:
      return libMesh::QUAD8;
    case libMesh::HEX27:
      return libMesh::QUAD9;

    case libMesh::PRISM15:
    case libMesh::PRISM18:
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
    case libMesh::PRISM6:
    case libMesh::PRISM20:
    case libMesh::PRISM21:
    case libMesh::PYRAMID5:
    case libMesh::PYRAMID18:
      detail::abort_unsupported("get_side_topology(): mixed-face elements require side-specific topology");
      return libMesh::INVALID_ELEM;

    default:
      detail::abort_unsupported("get_side_topology(): unsupported element type");
      return libMesh::INVALID_ELEM; // unreachable after abort
  }
}

/// Map an ElemType to its base geometric class (order-independent).
/// e.g. QUAD4 / QUAD8 / QUAD9 all return FEElemClass::QUAD.
LIBMESH_DEVICE_INLINE libMesh::FEElemClass
class_from_topology(libMesh::ElemType topo)
{
  switch (topo)
  {
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::EDGE4:
      return libMesh::FEElemClass::EDGE;

    case libMesh::TRI3:
    case libMesh::TRI6:
    case libMesh::TRI7:
      return libMesh::FEElemClass::TRI;

    case libMesh::QUAD4:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
      return libMesh::FEElemClass::QUAD;

    case libMesh::TET4:
    case libMesh::TET10:
    case libMesh::TET14:
      return libMesh::FEElemClass::TET;

    case libMesh::HEX8:
    case libMesh::HEX20:
    case libMesh::HEX27:
      return libMesh::FEElemClass::HEX;

    case libMesh::PRISM6:
    case libMesh::PRISM15:
    case libMesh::PRISM18:
    case libMesh::PRISM20:
    case libMesh::PRISM21:
      return libMesh::FEElemClass::PRISM;

    case libMesh::PYRAMID5:
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
    case libMesh::PYRAMID18:
      return libMesh::FEElemClass::PYRAMID;

    default:
      detail::abort_unsupported("class_from_topology(): unsupported element type");
      return libMesh::FEElemClass::N_CLASSES; // unreachable after abort
  }
}

LIBMESH_DEVICE_INLINE libMesh::ElemType
lagrange_shape_topology_or_invalid(FEShapeKey key)
{
  switch (key.order)
  {
    case libMesh::FIRST:
      switch (key.elem_type)
      {
        case libMesh::EDGE2:
        case libMesh::EDGE3:
        case libMesh::EDGE4:
          return libMesh::EDGE2;

        case libMesh::TRI3:
        case libMesh::TRI6:
        case libMesh::TRI7:
          return libMesh::TRI3;

        case libMesh::QUAD4:
        case libMesh::QUAD8:
        case libMesh::QUAD9:
          return libMesh::QUAD4;

        case libMesh::TET4:
        case libMesh::TET10:
        case libMesh::TET14:
          return libMesh::TET4;

        case libMesh::HEX8:
        case libMesh::HEX20:
        case libMesh::HEX27:
          return libMesh::HEX8;

        default:
          return libMesh::INVALID_ELEM;
      }

    case libMesh::SECOND:
      switch (key.elem_type)
      {
        case libMesh::EDGE3:
          return libMesh::EDGE3;

        case libMesh::TRI6:
        case libMesh::TRI7:
          return libMesh::TRI6;

        case libMesh::QUAD8:
          return libMesh::QUAD8;

        case libMesh::QUAD9:
          return libMesh::QUAD9;

        case libMesh::TET10:
        case libMesh::TET14:
          return libMesh::TET10;

        case libMesh::HEX20:
          return libMesh::HEX20;

        case libMesh::HEX27:
          return libMesh::HEX27;

        default:
          return libMesh::INVALID_ELEM;
      }

    default:
      return libMesh::INVALID_ELEM;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
lagrange_exact_n_dofs_or_zero(libMesh::ElemType elem_type,
                              libMesh::Order order)
{
  switch (order)
  {
    case libMesh::CONSTANT:
      return (elem_type == libMesh::NODEELEM) ? 1u : 0u;

    case libMesh::FIRST:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
        case libMesh::EDGE4:
          return 2;

        case libMesh::TRI3:
        case libMesh::TRI6:
        case libMesh::TRI7:
          return 3;

        case libMesh::QUAD4:
        case libMesh::QUAD8:
        case libMesh::QUAD9:
          return 4;

        case libMesh::TET4:
        case libMesh::TET10:
        case libMesh::TET14:
          return 4;

        case libMesh::HEX8:
        case libMesh::HEX20:
        case libMesh::HEX27:
          return 8;

        case libMesh::PRISM6:
        case libMesh::PRISM15:
        case libMesh::PRISM18:
        case libMesh::PRISM20:
        case libMesh::PRISM21:
          return 6;

        case libMesh::PYRAMID5:
        case libMesh::PYRAMID13:
        case libMesh::PYRAMID14:
        case libMesh::PYRAMID18:
          return 5;

        default:
          return 0;
      }

    case libMesh::SECOND:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE3:
          return 3;

        case libMesh::TRI6:
        case libMesh::TRI7:
          return 6;

        case libMesh::QUAD8:
          return 8;

        case libMesh::QUAD9:
          return 9;

        case libMesh::TET10:
        case libMesh::TET14:
          return 10;

        case libMesh::HEX20:
          return 20;

        case libMesh::HEX27:
          return 27;

        case libMesh::PRISM15:
          return 15;

        case libMesh::PRISM18:
        case libMesh::PRISM20:
        case libMesh::PRISM21:
          return 18;

        case libMesh::PYRAMID13:
          return 13;

        case libMesh::PYRAMID14:
        case libMesh::PYRAMID18:
          return 14;

        default:
          return 0;
      }

    case libMesh::THIRD:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE4:
          return 4;

        case libMesh::TRI7:
          return 7;

        case libMesh::TET14:
          return 14;

        case libMesh::PRISM20:
          return 20;

        case libMesh::PRISM21:
          return 21;

        case libMesh::PYRAMID18:
          return 18;

        default:
          return 0;
      }

    default:
      return 0;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
monomial_exact_n_dofs_or_zero(libMesh::ElemType elem_type,
                              libMesh::Order order)
{
  if (elem_type == libMesh::INVALID_ELEM)
    return 0;
  if (order < libMesh::CONSTANT)
    return 0;

  switch (order)
  {
    case libMesh::CONSTANT:
      return 1;

    case libMesh::FIRST:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
        case libMesh::EDGE4:
          return 2;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 3;
      if (is_monomial_3d_elem_type(elem_type))
        return 4;
      return 0;

    case libMesh::SECOND:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
        case libMesh::EDGE4:
          return 3;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 6;
      if (is_monomial_3d_elem_type(elem_type))
        return 10;
      return 0;

    case libMesh::THIRD:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
        case libMesh::EDGE4:
          return 4;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 10;
      if (is_monomial_3d_elem_type(elem_type))
        return 20;
      return 0;

    case libMesh::FOURTH:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
          return 5;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 15;
      if (is_monomial_3d_elem_type(elem_type, false))
        return 35;
      return 0;

    case libMesh::FIFTH:
      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
          return 6;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 21;
      if (is_monomial_3d_elem_type(elem_type, false))
        return 56;
      return 0;

    default:
    {
      const unsigned int p = static_cast<unsigned int>(order);

      switch (elem_type)
      {
        case libMesh::NODEELEM:
          return 1;

        case libMesh::EDGE2:
        case libMesh::EDGE3:
          return p + 1;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return (p + 1) * (p + 2) / 2;
      if (is_monomial_3d_elem_type(elem_type, false))
        return (p + 1) * (p + 2) * (p + 3) / 6;
      return 0;
    }
  }
}

LIBMESH_DEVICE_INLINE unsigned int
monomial_evaluator_dim_or_zero(libMesh::ElemType elem_type)
{
  switch (elem_type)
  {
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::EDGE4:
      return 1;

    case libMesh::TRI3:
    case libMesh::TRI6:
    case libMesh::TRI7:
    case libMesh::QUAD4:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
      return 2;

    case libMesh::TET4:
    case libMesh::TET10:
    case libMesh::TET14:
    case libMesh::HEX8:
    case libMesh::HEX20:
    case libMesh::HEX27:
    case libMesh::PRISM6:
    case libMesh::PRISM15:
    case libMesh::PRISM18:
    case libMesh::PRISM20:
    case libMesh::PRISM21:
    case libMesh::PYRAMID5:
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
    case libMesh::PYRAMID18:
      return 3;

    default:
      return 0;
  }
}

/// Return true iff the current Kokkos physics evaluators can evaluate \p key.
/// This boundary is the intersection of:
/// 1. exact libMesh-valid (family, elem_type, order) keys, and
/// 2. currently implemented Kokkos evaluator topologies/orders.
LIBMESH_DEVICE_INLINE bool
supports_shape(FEShapeKey key)
{
  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return lagrange_exact_n_dofs_or_zero(key.elem_type, key.order) != 0 &&
             lagrange_shape_topology_or_invalid(key) != libMesh::INVALID_ELEM;

    case libMesh::MONOMIAL:
      return monomial_exact_n_dofs_or_zero(key.elem_type, key.order) != 0 &&
             monomial_evaluator_dim_or_zero(key.elem_type) != 0 &&
             key.order >= libMesh::CONSTANT &&
             key.order <= libMesh::FIFTH;

    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
supports_grad_shape(FEShapeKey key)
{
  return supports_shape(key);
}

LIBMESH_DEVICE_INLINE bool
supports_n_dofs(FEShapeKey key)
{
  return supports_shape(key);
}

/// Return the number of DOFs for a physics FE space described by \p key,
/// restricted to the current Kokkos evaluator support boundary.
LIBMESH_DEVICE_INLINE unsigned int
n_dofs(FEShapeKey key)
{
  if (!supports_n_dofs(key))
  {
    detail::abort_unsupported("n_dofs(FEShapeKey): unsupported FE key for current Kokkos evaluator support boundary");
    return 0;
  }

  switch (key.family)
  {
    case libMesh::LAGRANGE:
      return lagrange_exact_n_dofs_or_zero(key.elem_type, key.order);

    case libMesh::MONOMIAL:
      return monomial_exact_n_dofs_or_zero(key.elem_type, key.order);

    default:
      detail::abort_unsupported("n_dofs(FEShapeKey): unsupported FE family");
      return 0;
  }
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_TYPES_H
