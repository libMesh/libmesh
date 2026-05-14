#ifndef LIBMESH_KOKKOS_FE_FACE_MAP_H
#define LIBMESH_KOKKOS_FE_FACE_MAP_H

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos_fe_evaluator.h"
#include "libmesh/elem.h"

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
recover_parent_side(const libMesh::Elem & parent,
                    const libMesh::Elem & side_in_parent)
{
  for (unsigned int side = 0; side < parent.n_sides(); ++side)
  {
    auto candidate = parent.build_side_ptr(side);

    if (candidate->type() != side_in_parent.type() ||
        candidate->n_nodes() != side_in_parent.n_nodes())
      continue;

    bool same_side = true;
    for (unsigned int k = 0; k < candidate->n_nodes(); ++k)
      if (candidate->node_ptr(k) != side_in_parent.node_ptr(k))
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
  switch (parent.type())
  {
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
      switch (node)
      {
        case 9:
          return libMesh::Point(-0.5, -0.5, 0.5);
        case 10:
          return libMesh::Point(0.5, -0.5, 0.5);
        case 11:
          return libMesh::Point(0.5, 0.5, 0.5);
        case 12:
          return libMesh::Point(-0.5, 0.5, 0.5);
        default:
          return parent.master_point(node);
      }

    case libMesh::PYRAMID18:
      switch (node)
      {
        case 9:
          return libMesh::Point(-0.5, -0.5, 0.5);
        case 10:
          return libMesh::Point(0.5, -0.5, 0.5);
        case 11:
          return libMesh::Point(0.5, 0.5, 0.5);
        case 12:
          return libMesh::Point(-0.5, 0.5, 0.5);
        case 14:
          return libMesh::Point(-2. / 3., 0.0, 1. / 3.);
        case 15:
          return libMesh::Point(0.0, 2. / 3., 1. / 3.);
        case 16:
          return libMesh::Point(2. / 3., 0.0, 1. / 3.);
        case 17:
          return libMesh::Point(0.0, -2. / 3., 1. / 3.);
        default:
          return parent.master_point(node);
      }

    default:
      return parent.master_point(node);
  }
}

/**
 * Map a face quadrature point from the side element's reference coordinate system
 * to the parent element's reference coordinate system.
 *
 * side_in_parent must be obtained via build_side_ptr() (not side_ptr()), so that
 * second-order sides carry their midpoint nodes. Parent reference coordinates
 * are reconstructed from the FE reference-space node convention used by
 * FE::side_map(), not from side_in_parent.point(k), which lives in physical
 * space, and not from Elem::master_point() on pyramids, where those node
 * coordinates differ.
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
    const libMesh::Point pt = parent_refspace_node(*parent, parent->local_side_node(side, 0));
    return point_to_real_vector(pt);
  }

  for (unsigned int k = 0; k < n; ++k)
  {
    const Real s = face_qpt(0);
    const Real t = face_qpt(1);
    const Real psi = map_shape(mapping_type, side_topo, k, s, t, 0.0);

    const libMesh::Point pt = parent_refspace_node(*parent, parent->local_side_node(side, k));
    parent_pt.add_scaled(point_to_real_vector(pt), psi);
  }

  return parent_pt;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_HAVE_KOKKOS

#endif // LIBMESH_KOKKOS_FE_FACE_MAP_H
