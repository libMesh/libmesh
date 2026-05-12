#ifndef LIBMESH_FE_SHAPE_TRAITS_H
#define LIBMESH_FE_SHAPE_TRAITS_H

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_elem_class.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/libmesh_device.h"

namespace libMesh
{

struct FEShapeKey
{
  FEFamily family;
  ElemType elem_type;
  Order    order;
};

LIBMESH_DEVICE_INLINE bool
is_monomial_2d_elem_type(ElemType elem_type)
{
  switch (elem_type)
  {
    case C0POLYGON:
    case TRI3:
    case TRISHELL3:
    case TRI6:
    case TRI7:
    case QUAD4:
    case QUADSHELL4:
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
    case QUADSHELL9:
      return true;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE bool
is_monomial_3d_elem_type(ElemType elem_type,
                         bool include_pyramid18 = true)
{
  switch (elem_type)
  {
    case TET4:
    case TET10:
    case TET14:
    case HEX8:
    case HEX20:
    case HEX27:
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case C0POLYHEDRON:
      return true;
    case PYRAMID18:
      return include_pyramid18;
    default:
      return false;
  }
}

LIBMESH_DEVICE_INLINE ElemType
side_topology_or_invalid(ElemType parent)
{
  switch (parent)
  {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return EDGE2;

    case TRI3:
    case QUAD4:
      return EDGE2;

    case TRI6:
    case TRI7:
    case QUAD8:
    case QUAD9:
      return EDGE3;

    case TET4:
      return TRI3;
    case HEX8:
      return QUAD4;

    case TET10:
      return TRI6;
    case TET14:
      return TRI7;
    case HEX20:
      return QUAD8;
    case HEX27:
      return QUAD9;

    default:
      return INVALID_ELEM;
  }
}

LIBMESH_DEVICE_INLINE FEElemClass
class_from_topology_or_invalid(ElemType topo)
{
  switch (topo)
  {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return FEElemClass::EDGE;

    case TRI3:
    case TRI6:
    case TRI7:
      return FEElemClass::TRI;

    case QUAD4:
    case QUAD8:
    case QUAD9:
      return FEElemClass::QUAD;

    case TET4:
    case TET10:
    case TET14:
      return FEElemClass::TET;

    case HEX8:
    case HEX20:
    case HEX27:
      return FEElemClass::HEX;

    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
      return FEElemClass::PRISM;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      return FEElemClass::PYRAMID;

    default:
      return FEElemClass::N_CLASSES;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
elem_class_dim_or_zero(FEElemClass cls)
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

LIBMESH_DEVICE_INLINE unsigned int
topology_dim_or_zero(ElemType topo)
{
  return elem_class_dim_or_zero(class_from_topology_or_invalid(topo));
}

LIBMESH_DEVICE_INLINE ElemType
lagrange_shape_topology_or_invalid(FEShapeKey key)
{
  switch (key.order)
  {
    case CONSTANT:
    case FIRST:
      switch (key.elem_type)
      {
        case EDGE2:
        case EDGE3:
        case EDGE4:
          return EDGE2;

        case TRI3:
        case TRI6:
        case TRI7:
          return TRI3;

        case QUAD4:
        case QUAD8:
        case QUAD9:
          return QUAD4;

        case TET4:
        case TET10:
        case TET14:
          return TET4;

        case HEX8:
        case HEX20:
        case HEX27:
          return HEX8;

        default:
          return INVALID_ELEM;
      }

    case SECOND:
      switch (key.elem_type)
      {
        case EDGE3:
          return EDGE3;

        case TRI6:
        case TRI7:
          return TRI6;

        case QUAD8:
          return QUAD8;

        case QUAD9:
          return QUAD9;

        case TET10:
        case TET14:
          return TET10;

        case HEX20:
          return HEX20;

        case HEX27:
          return HEX27;

        default:
          return INVALID_ELEM;
      }

    default:
      return INVALID_ELEM;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
lagrange_exact_n_dofs_or_zero(ElemType elem_type,
                              Order order)
{
  switch (order)
  {
    case CONSTANT:
      return (elem_type == NODEELEM) ? 1u : 0u;

    case FIRST:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
        case EDGE4:
          return 2;

        case TRI3:
        case TRI6:
        case TRI7:
          return 3;

        case QUAD4:
        case QUAD8:
        case QUAD9:
          return 4;

        case TET4:
        case TET10:
        case TET14:
          return 4;

        case HEX8:
        case HEX20:
        case HEX27:
          return 8;

        case PRISM6:
        case PRISM15:
        case PRISM18:
        case PRISM20:
        case PRISM21:
          return 6;

        case PYRAMID5:
        case PYRAMID13:
        case PYRAMID14:
        case PYRAMID18:
          return 5;

        default:
          return 0;
      }

    case SECOND:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE3:
          return 3;

        case TRI6:
        case TRI7:
          return 6;

        case QUAD8:
          return 8;

        case QUAD9:
          return 9;

        case TET10:
        case TET14:
          return 10;

        case HEX20:
          return 20;

        case HEX27:
          return 27;

        case PRISM15:
          return 15;

        case PRISM18:
        case PRISM20:
        case PRISM21:
          return 18;

        case PYRAMID13:
          return 13;

        case PYRAMID14:
        case PYRAMID18:
          return 14;

        default:
          return 0;
      }

    case THIRD:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE4:
          return 4;

        case TRI7:
          return 7;

        case TET14:
          return 14;

        case PRISM20:
          return 20;

        case PRISM21:
          return 21;

        case PYRAMID18:
          return 18;

        default:
          return 0;
      }

    default:
      return 0;
  }
}

LIBMESH_DEVICE_INLINE unsigned int
monomial_exact_n_dofs_or_zero(ElemType elem_type,
                              Order order)
{
  if (elem_type == INVALID_ELEM)
    return 0;
  if (order < CONSTANT)
    return 0;

  switch (order)
  {
    case CONSTANT:
      return 1;

    case FIRST:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
        case EDGE4:
          return 2;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 3;
      if (is_monomial_3d_elem_type(elem_type))
        return 4;
      return 0;

    case SECOND:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
        case EDGE4:
          return 3;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 6;
      if (is_monomial_3d_elem_type(elem_type))
        return 10;
      return 0;

    case THIRD:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
        case EDGE4:
          return 4;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 10;
      if (is_monomial_3d_elem_type(elem_type))
        return 20;
      return 0;

    case FOURTH:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
          return 5;

        default:
          break;
      }

      if (is_monomial_2d_elem_type(elem_type))
        return 15;
      if (is_monomial_3d_elem_type(elem_type, false))
        return 35;
      return 0;

    case FIFTH:
      switch (elem_type)
      {
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
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
        case NODEELEM:
          return 1;

        case EDGE2:
        case EDGE3:
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
monomial_evaluator_dim_or_zero(ElemType elem_type)
{
  switch (elem_type)
  {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return 1;

    case TRI3:
    case TRI6:
    case TRI7:
    case QUAD4:
    case QUAD8:
    case QUAD9:
      return 2;

    case TET4:
    case TET10:
    case TET14:
    case HEX8:
    case HEX20:
    case HEX27:
    case PRISM6:
    case PRISM15:
    case PRISM18:
    case PRISM20:
    case PRISM21:
    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
    case PYRAMID18:
      return 3;

    default:
      return 0;
  }
}

LIBMESH_DEVICE_INLINE bool
supports_shape(FEShapeKey key)
{
  switch (key.family)
  {
    case LAGRANGE:
      return lagrange_exact_n_dofs_or_zero(key.elem_type, key.order) != 0 &&
             lagrange_shape_topology_or_invalid(key) != INVALID_ELEM;

    case MONOMIAL:
      return monomial_exact_n_dofs_or_zero(key.elem_type, key.order) != 0 &&
             monomial_evaluator_dim_or_zero(key.elem_type) != 0 &&
             key.order >= CONSTANT &&
             key.order <= FIFTH;

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

LIBMESH_DEVICE_INLINE unsigned int
n_dofs_or_zero(FEShapeKey key)
{
  switch (key.family)
  {
    case LAGRANGE:
      return lagrange_exact_n_dofs_or_zero(key.elem_type, key.order);

    case MONOMIAL:
      return monomial_exact_n_dofs_or_zero(key.elem_type, key.order);

    default:
      return 0;
  }
}

} // namespace libMesh

#endif // LIBMESH_FE_SHAPE_TRAITS_H
