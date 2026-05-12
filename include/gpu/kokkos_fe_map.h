// Kokkos device-compatible physical map evaluation.
//
// All functions are LIBMESH_DEVICE_INLINE — callable from both host and GPU.
//
// Two API levels:
//   1. Template on ElemType (preferred): eliminates the topology switch at
//      compile time, producing small inlined functions with no stack pressure.
//   2. Runtime ElemType dispatch: convenient but requires increased CUDA
//      stack size due to the large switch in map_shape.
//
// Given node coordinates and a reference-space point, these functions compute:
//   - Physical coordinates (xyz)
//   - Jacobian matrix (reference -> physical)
//   - Jacobian measures and JxW
//   - Outward normal helpers for face/edge integrals

#ifndef LIBMESH_KOKKOS_FE_MAP_H
#define LIBMESH_KOKKOS_FE_MAP_H

#include "kokkos_fe_evaluator.h"

namespace libMesh::Kokkos
{

template <libMesh::FEFamily family, libMesh::ElemType topo>
LIBMESH_DEVICE_INLINE RealVector
physical_point(const RealVector * nodes,
              unsigned int n_nodes,
              Real xi, Real eta, Real zeta)
{
  RealVector xyz = zero_vector();
  for (unsigned int i = 0; i < n_nodes; ++i)
    xyz += map_shape<family, topo>(i, xi, eta, zeta) * nodes[i];
  return xyz;
}

// =========================================================================
//  Compile-time dispatch (preferred for GPU — no switch overhead)
//
//  Template on FEFamily and ElemType so nvcc only instantiates the specific
//  FEEvaluator specialization. No topology switch means no stack pressure.
// =========================================================================

template <libMesh::FEFamily family, libMesh::ElemType topo>
LIBMESH_DEVICE_INLINE RealTensor
jacobian(const RealVector * nodes,
         unsigned int n_nodes,
         Real xi, Real eta, Real zeta)
{
  RealTensor J = zero_tensor();
  for (unsigned int k = 0; k < n_nodes; ++k)
    J += libMesh::outer_product(grad_map_shape<family, topo>(k, xi, eta, zeta), nodes[k]);
  return J;
}

template <libMesh::FEFamily family, libMesh::ElemType topo>
LIBMESH_DEVICE_INLINE void
physical_point_and_jacobian(const RealVector * nodes,
                         unsigned int n_nodes,
                         Real xi, Real eta, Real zeta,
                         RealVector & xyz,
                         RealTensor & J)
{
  xyz = zero_vector();
  J = zero_tensor();
  for (unsigned int k = 0; k < n_nodes; ++k)
  {
    const Real phi = map_shape<family, topo>(k, xi, eta, zeta);
    const RealVector grad = grad_map_shape<family, topo>(k, xi, eta, zeta);
    xyz += phi * nodes[k];
    J += libMesh::outer_product(grad, nodes[k]);
  }
}

template <libMesh::FEFamily family, libMesh::ElemType face_topo>
LIBMESH_DEVICE_INLINE RealTensor
face_jacobian(const RealVector * face_nodes,
             unsigned int n_face_nodes,
             Real xi, Real eta, Real zeta)
{
  RealTensor J = zero_tensor();
  for (unsigned int k = 0; k < n_face_nodes; ++k)
    J += libMesh::outer_product(grad_map_shape<family, face_topo>(k, xi, eta, zeta),
                                face_nodes[k]);
  return J;
}

// =========================================================================
//  Runtime topology dispatch (convenient, but larger GPU stack usage)
// =========================================================================

/// Compute physical coordinate (runtime topology).
LIBMESH_DEVICE_INLINE RealVector
physical_point(libMesh::ElemMappingType mapping_type,
              libMesh::ElemType topo,
              const RealVector * nodes,
              unsigned int n_nodes,
              Real xi, Real eta, Real zeta)
{
  RealVector xyz = zero_vector();
  for (unsigned int i = 0; i < n_nodes; ++i)
    xyz += map_shape(mapping_type, topo, i, xi, eta, zeta) * nodes[i];
  return xyz;
}

/// Compute Jacobian matrix (runtime topology), with rows d(x)/d(xi_r).
LIBMESH_DEVICE_INLINE RealTensor
jacobian(libMesh::ElemMappingType mapping_type,
         libMesh::ElemType topo,
         const RealVector * nodes,
         unsigned int n_nodes,
         Real xi, Real eta, Real zeta)
{
  RealTensor J = zero_tensor();
  for (unsigned int k = 0; k < n_nodes; ++k)
    J += libMesh::outer_product(grad_map_shape(mapping_type, topo, k, xi, eta, zeta),
                                nodes[k]);
  return J;
}

/// Compute physical point and Jacobian together (runtime topology).
LIBMESH_DEVICE_INLINE void
physical_point_and_jacobian(libMesh::ElemMappingType mapping_type,
                         libMesh::ElemType topo,
                         const RealVector * nodes,
                         unsigned int n_nodes,
                         Real xi, Real eta, Real zeta,
                         RealVector & xyz,
                         RealTensor & J)
{
  xyz = zero_vector();
  J = zero_tensor();
  for (unsigned int k = 0; k < n_nodes; ++k)
  {
    const Real phi = map_shape(mapping_type, topo, k, xi, eta, zeta);
    const RealVector grad = grad_map_shape(mapping_type, topo, k, xi, eta, zeta);
    xyz += phi * nodes[k];
    J += libMesh::outer_product(grad, nodes[k]);
  }
}

/// Face Jacobian (runtime topology).
LIBMESH_DEVICE_INLINE RealTensor
face_jacobian(libMesh::ElemMappingType mapping_type,
             libMesh::ElemType face_topo,
             const RealVector * face_nodes,
             unsigned int n_face_nodes,
             Real xi, Real eta, Real zeta)
{
  RealTensor J = zero_tensor();
  for (unsigned int k = 0; k < n_face_nodes; ++k)
    J += libMesh::outer_product(grad_map_shape(mapping_type, face_topo, k, xi, eta, zeta),
                                face_nodes[k]);
  return J;
}

// =========================================================================
//  Geometry helpers (topology-independent)
// =========================================================================

/// libMesh FEMap-compatible volume measure * quadrature_weight.
///   3D: det(J)                       * weight
///   2D: ||J_row0 x J_row1||          * weight
///   1D: ||J_row0||                   * weight
///   0D: weight
LIBMESH_DEVICE_INLINE Real
volume_jxw(const RealTensor & J, unsigned int dim, Real quad_weight)
{
  if (dim == 3)
    return detail::leading_determinant(J, 3) * quad_weight;
  else if (dim == 2)
    return J.row(0).cross(J.row(1)).norm() * quad_weight;
  else if (dim == 1)
    return J.row(0).norm() * quad_weight;
  else
    return quad_weight;
}

/// Face JxW: surface measure * quadrature_weight
///   3D: ||J_row0 x J_row1|| * weight
///   2D: ||J_row0||           * weight
///   1D: weight (face is a point)
LIBMESH_DEVICE_INLINE Real
face_jxw(const RealTensor & J, unsigned int parent_dim, Real quad_weight)
{
  if (parent_dim == 3)
    return J.row(0).cross(J.row(1)).norm() * quad_weight;
  else if (parent_dim == 2)
    return J.row(0).norm() * quad_weight;
  else
    return quad_weight;
}

/// Outward unit normal for a 3D face from the face Jacobian.
LIBMESH_DEVICE_INLINE RealVector
face_normal(const RealTensor & J, unsigned int parent_dim)
{
  if (parent_dim != 3)
  {
    detail::abort_unsupported("face_normal(): only 3D face normals are defined from face Jacobians alone; use edge_normal_on_parent_surface() for 2D parent elements");
    return zero_vector();
  }

  RealVector n = J.row(0).cross(J.row(1));

  const Real len = n.norm();
  if (len > 0.0)
    n *= 1.0 / len;
  return n;
}

/// Outward edge normal for a 2D parent element embedded in 3D.
/// Requires the edge Jacobian and the parent surface Jacobian at the mapped
/// parent-reference point.
LIBMESH_DEVICE_INLINE RealVector
edge_normal_on_parent_surface(const RealTensor & edge_J,
                              const RealTensor & parent_J)
{
  RealVector surface_normal = parent_J.row(0).cross(parent_J.row(1));
  const Real surface_len = surface_normal.norm();
  if (surface_len > 0.0)
    surface_normal *= 1.0 / surface_len;

  RealVector n = edge_J.row(0).cross(surface_normal);

  const Real len = n.norm();
  if (len > 0.0)
    n *= 1.0 / len;
  return n;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_MAP_H
