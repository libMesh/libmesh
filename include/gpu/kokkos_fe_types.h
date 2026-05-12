// Kokkos FE type helpers.
//
// Shared FE topology/order support metadata lives in fe_shape_traits.h.
// This header keeps the Kokkos-facing hard-fail wrappers and namespace
// compatibility for existing FE device code.

#ifndef LIBMESH_KOKKOS_FE_TYPES_H
#define LIBMESH_KOKKOS_FE_TYPES_H

#include "libmesh/enum_fe_elem_class.h"
#include "libmesh/fe_reference_element_traits.h"
#include "libmesh/fe_shape_traits.h"
#include "libmesh/libmesh_device.h"
#ifndef LIBMESH_KOKKOS_COMPILATION
#  include "libmesh/libmesh_common.h"
#endif

namespace libMesh::Kokkos
{

using libMesh::FEElemClass;
using libMesh::FEShapeKey;
using libMesh::is_monomial_2d_elem_type;
using libMesh::is_monomial_3d_elem_type;
using libMesh::lagrange_shape_topology_or_invalid;
using libMesh::lagrange_exact_n_dofs_or_zero;
using libMesh::monomial_exact_n_dofs_or_zero;
using libMesh::monomial_evaluator_dim_or_zero;
using libMesh::supports_shape;
using libMesh::supports_grad_shape;
using libMesh::supports_n_dofs;

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

LIBMESH_DEVICE_INLINE libMesh::ElemType
get_side_topology(libMesh::ElemType parent)
{
  const libMesh::ElemType side_topology = libMesh::side_topology_or_invalid(parent);

  if (side_topology == libMesh::INVALID_ELEM)
  {
    if (requires_side_specific_topology(parent))
      detail::abort_unsupported("get_side_topology(): mixed-face elements require side-specific topology");
    else
      detail::abort_unsupported("get_side_topology(): unsupported element type");
    return libMesh::INVALID_ELEM;
  }

  return side_topology;
}

LIBMESH_DEVICE_INLINE libMesh::ElemType
get_side_topology(libMesh::ElemType parent,
                  unsigned int side)
{
  const libMesh::ElemType side_topology = libMesh::side_topology_or_invalid(parent, side);

  if (side_topology != libMesh::INVALID_ELEM)
    return side_topology;

  return get_side_topology(parent);
}

LIBMESH_DEVICE_INLINE libMesh::FEElemClass
class_from_topology(libMesh::ElemType topo)
{
  const libMesh::FEElemClass elem_class = libMesh::class_from_topology_or_invalid(topo);

  if (elem_class == libMesh::FEElemClass::N_CLASSES)
  {
    detail::abort_unsupported("class_from_topology(): unsupported element type");
    return libMesh::FEElemClass::N_CLASSES;
  }

  return elem_class;
}

LIBMESH_DEVICE_INLINE unsigned int
n_dofs(FEShapeKey key)
{
  if (!supports_n_dofs(key))
  {
    detail::abort_unsupported("n_dofs(FEShapeKey): unsupported FE key for current Kokkos evaluator support boundary");
    return 0;
  }

  return libMesh::n_dofs_or_zero(key);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_TYPES_H
