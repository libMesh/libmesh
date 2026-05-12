#ifndef LIBMESH_KOKKOS_FE_FACE_MAP_H
#define LIBMESH_KOKKOS_FE_FACE_MAP_H

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos_fe_evaluator.h"
#include "libmesh/elem.h"
#include "libmesh/fe_reference_element_traits.h"

namespace libMesh::Kokkos
{

LIBMESH_DEVICE_INLINE
RealVector point_to_real_vector(const libMesh::Point & pt)
{
#if LIBMESH_DIM == 1
  return make_vector(pt(0));
#elif LIBMESH_DIM == 2
  return make_vector(pt(0), pt(1));
#else
  return make_vector(pt(0), pt(1), pt(2));
#endif
}

inline unsigned int
parent_local_side_node(const libMesh::Elem & parent,
                       unsigned int side,
                       unsigned int side_node)
{
  unsigned int node = libMesh::invalid_uint;
  if (libMesh::try_local_side_node(parent.type(), side, side_node, node))
    return node;

  detail::abort_unsupported("map_face_qp_to_parent(): unsupported parent element type in local side-node lookup");
  return libMesh::invalid_uint;
}

inline unsigned int
recover_parent_side(const libMesh::Elem & parent,
                    const libMesh::Elem & side_in_parent)
{
  for (unsigned int side = 0; side < parent.n_sides(); ++side)
  {
    if (get_side_topology(parent.type(), side) != side_in_parent.type() ||
        ((libMesh::side_node_count_or_zero(parent.type(), side) &&
          libMesh::side_node_count_or_zero(parent.type(), side) != side_in_parent.n_nodes())))
      continue;

    bool same_side = true;
    for (unsigned int k = 0; k < side_in_parent.n_nodes(); ++k)
      if (parent.node_ptr(parent_local_side_node(parent, side, k)) != side_in_parent.node_ptr(k))
      {
        same_side = false;
        break;
      }

    if (same_side)
      return side;
  }

  return libMesh::invalid_uint;
}

inline libMesh::Point
parent_refspace_node(const libMesh::Elem & parent, unsigned int node)
{
  libMesh::Point pt;
  if (libMesh::try_reference_node(parent.type(), node, pt))
    return pt;

  detail::abort_unsupported("map_face_qp_to_parent(): unsupported parent element type in reference-node lookup");
  return libMesh::Point();
}

/**
 * Map a face quadrature point from the side element's reference coordinate system
 * to the parent element's reference coordinate system.
 *
 * side_in_parent must be obtained via build_side_ptr() (not side_ptr()), so that
 * second-order sides carry their midpoint nodes. Parent reference coordinates
 * are reconstructed from shared libMesh reference-element traits. They are not
 * reconstructed from side_in_parent.point(k), which lives in physical space.
 * Element types outside the Kokkos FE support boundary are rejected rather
 * than silently falling back to generic Elem runtime helpers.
 *
 * @param side_in_parent  The side element as embedded in the parent (from build_side_ptr())
 * @param mapping_type    Geometric mapping type (LAGRANGE_MAP, RATIONAL_BERNSTEIN_MAP)
 * @param side_topo       Topology of the side element (libMesh::ElemType)
 * @param face_qpt        Quadrature point in the side element's reference coordinates
 * @returns               Corresponding point in the parent element's reference coordinates
 */
inline RealVector
map_face_qp_to_parent(const libMesh::Elem & side_in_parent,
                  libMesh::ElemMappingType mapping_type,
                  libMesh::ElemType side_topo,
                  RealVector face_qpt)
{
  const libMesh::Elem * parent = side_in_parent.interior_parent();
  libmesh_error_msg_if(!parent,
                       "map_face_qp_to_parent(): side element must carry an interior_parent() from build_side_ptr()");

  const unsigned int side = recover_parent_side(*parent, side_in_parent);
  libmesh_error_msg_if(side == libMesh::invalid_uint,
                       "map_face_qp_to_parent(): could not recover parent side for the provided side element");

  const unsigned int n = side_in_parent.n_nodes();
  RealVector parent_pt = zero_vector();

  // 1-D elements: the "side" is a single vertex node. There is only one
  // point-side reference coordinate, (0,0,0), so we map directly to the
  // corresponding parent vertex in the parent reference element.
  if (n == 1)
  {
    const libMesh::Point pt = parent_refspace_node(*parent, parent_local_side_node(*parent, side, 0));
    return point_to_real_vector(pt);
  }

  for (unsigned int k = 0; k < n; ++k)
  {
    const Real s = face_qpt(0);
    const Real t = face_qpt(1);
    const Real psi = map_shape(mapping_type, side_topo, k, s, t, 0.0);

    const libMesh::Point pt = parent_refspace_node(*parent, parent_local_side_node(*parent, side, k));
    parent_pt.add_scaled(point_to_real_vector(pt), psi);
  }

  return parent_pt;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_HAVE_KOKKOS

#endif // LIBMESH_KOKKOS_FE_FACE_MAP_H
