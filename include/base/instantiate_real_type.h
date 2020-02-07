// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_INSTANTIATE_REAL_TYPE_H
#define LIBMESH_INSTANTIATE_REAL_TYPE_H

#include "libmesh/default_coupling_impl.h"
#include "libmesh/ghost_point_neighbors_impl.h"
#include "libmesh/point_neighbor_coupling_impl.h"
#include "libmesh/sibling_coupling_impl.h"
#include "libmesh/fe_abstract_impl.h"
#include "libmesh/fe_base_impl.h"
#include "libmesh/fe_bernstein_impl.h"
#include "libmesh/fe_boundary_impl.h"
#include "libmesh/fe_clough_impl.h"
#include "libmesh/fe_hermite_impl.h"
#include "libmesh/fe_hierarchic_impl.h"
#include "libmesh/fe_impl.h"
#include "libmesh/fe_interface_impl.h"
#include "libmesh/fe_l2_hierarchic_impl.h"
#include "libmesh/fe_l2_lagrange_impl.h"
#include "libmesh/fe_lagrange_impl.h"
#include "libmesh/fe_lagrange_vec_impl.h"
#include "libmesh/fe_map_impl.h"
#include "libmesh/fe_monomial_impl.h"
#include "libmesh/fe_monomial_vec_impl.h"
#include "libmesh/fe_nedelec_one_impl.h"
#include "libmesh/fe_rational_impl.h"
#include "libmesh/fe_scalar_impl.h"
#include "libmesh/fe_subdivision_impl.h"
#include "libmesh/fe_szabab_impl.h"
#include "libmesh/fe_transformation_base_impl.h"
#include "libmesh/fe_xyz_boundary_impl.h"
#include "libmesh/fe_xyz_impl.h"
#include "libmesh/h1_fe_transformation_impl.h"
#include "libmesh/hcurl_fe_transformation_impl.h"
#include "libmesh/cell_hex20_impl.h"
#include "libmesh/cell_hex27_impl.h"
#include "libmesh/cell_hex8_impl.h"
#include "libmesh/cell_hex_impl.h"
#include "libmesh/cell_prism15_impl.h"
#include "libmesh/cell_prism18_impl.h"
#include "libmesh/cell_prism6_impl.h"
#include "libmesh/cell_prism_impl.h"
#include "libmesh/cell_pyramid5_impl.h"
#include "libmesh/cell_pyramid13_impl.h"
#include "libmesh/cell_pyramid14_impl.h"
#include "libmesh/cell_pyramid_impl.h"
#include "libmesh/cell_tet10_impl.h"
#include "libmesh/cell_tet4_impl.h"
#include "libmesh/cell_tet_impl.h"
#include "libmesh/edge_edge3_impl.h"
#include "libmesh/edge_edge4_impl.h"
#include "libmesh/elem_impl.h"
#include "libmesh/face_quad4_impl.h"
#include "libmesh/face_quad8_impl.h"
#include "libmesh/face_quad9_impl.h"
#include "libmesh/face_quad_impl.h"
#include "libmesh/face_tri3_impl.h"
#include "libmesh/face_tri6_impl.h"
#include "libmesh/face_tri_impl.h"
#include "libmesh/face_tri3_subdivision_impl.h"
#include "libmesh/remote_elem_impl.h"
#include "libmesh/boundary_info_impl.h"
#include "libmesh/distributed_mesh_impl.h"
#include "libmesh/mesh_base_impl.h"
#include "libmesh/point_locator_base_impl.h"
#include "libmesh/point_locator_tree_impl.h"
#include "libmesh/mesh_communication_global_indices_impl.h"
#include "libmesh/mesh_communication_impl.h"
#include "libmesh/mesh_generation_impl.h"
#include "libmesh/mesh_iterators_impl.h"
#include "libmesh/mesh_refinement_impl.h"
#include "libmesh/mesh_tools_impl.h"
#include "libmesh/replicated_mesh_impl.h"
#include "libmesh/unstructured_mesh_impl.h"
#include "libmesh/dense_matrix_base_impl.h"
#include "libmesh/dense_matrix_impl.h"
#include "libmesh/parallel_elem_impl.h"
#include "libmesh/topology_map_impl.h"
#include "libmesh/tree_node_impl.h"
#include "libmesh/tree_impl.h"
#include "libmesh/periodic_boundaries_impl.h"
#include "libmesh/periodic_boundary_base_impl.h"
#include "libmesh/partitioner_impl.h"
#include "libmesh/centroid_partitioner_impl.h"
#include "libmesh/linear_partitioner_impl.h"
#include "libmesh/mapped_subdomain_partitioner_impl.h"
#include "libmesh/metis_partitioner_impl.h"
#include "libmesh/parmetis_partitioner_impl.h"
#include "libmesh/sfc_partitioner_impl.h"
#include "libmesh/subdomain_partitioner_impl.h"

#define INSTANTIATE_ALL_REAL_TYPE0(RealType)                                                       \
  template class DefaultCouplingTempl<RealType>;                                                   \
  template class GhostPointNeighborsTempl<RealType>;                                               \
  template class PointNeighborCouplingTempl<RealType>;                                             \
  template class SiblingCouplingTempl<RealType>;                                                   \
  template class FEGenericBase<Real, RealType>;                                                    \
  template class FEGenericBase<RealGradient, RealType>;                                            \
  FACE_EDGE_SHAPE_ERROR_INSTANTIATE(init_face_shape_functions, RealType);                          \
  FACE_EDGE_SHAPE_ERROR_INSTANTIATE(init_edge_shape_functions, RealType);                          \
  INIT_FACE_EDGE_SHAPE_FUNCTIONS_INSTANTIATE(RealType);                                            \
  REINIT_AND_SIDE_MAPS_ALL(RealType);                                                              \
  ERRORS_IN_0D_INSTANTIATE_ALL(RealType);                                                          \
  ALL_FAMILY_1D_ERRORS_INSTANTIATE(RealType);                                                      \
  SUB_ERRORS_IN_0D_INSTANTIATE(XYZ, RealType);                                                     \
  template void FE<2, SUBDIVISION, RealType>::reinit(                                              \
      ElemTempl<RealType> const *,                                                                 \
      unsigned int,                                                                                \
      Real,                                                                                        \
      const std::vector<PointTempl<RealType>> * const,                                             \
      const std::vector<Real> * const);                                                            \
  template void FE<2, NEDELEC_ONE, RealType>::reinit(                                              \
      ElemTempl<RealType> const *,                                                                 \
      unsigned int,                                                                                \
      Real,                                                                                        \
      const std::vector<PointTempl<RealType>> * const,                                             \
      const std::vector<Real> * const);                                                            \
  template void FE<2, NEDELEC_ONE, RealType>::side_map(ElemTempl<RealType> const *,                \
                                                       ElemTempl<RealType> const *,                \
                                                       const unsigned int,                         \
                                                       const std::vector<PointTempl<RealType>> &,  \
                                                       std::vector<PointTempl<RealType>> &);       \
  template void FE<2, NEDELEC_ONE, RealType>::edge_reinit(                                         \
      ElemTempl<RealType> const *,                                                                 \
      unsigned int,                                                                                \
      Real,                                                                                        \
      const std::vector<PointTempl<RealType>> * const,                                             \
      const std::vector<Real> * const);                                                            \
  template void FE<3, NEDELEC_ONE, RealType>::reinit(                                              \
      ElemTempl<RealType> const *,                                                                 \
      unsigned int,                                                                                \
      Real,                                                                                        \
      const std::vector<PointTempl<RealType>> * const,                                             \
      const std::vector<Real> * const);                                                            \
  template void FE<3, NEDELEC_ONE, RealType>::side_map(ElemTempl<RealType> const *,                \
                                                       ElemTempl<RealType> const *,                \
                                                       const unsigned int,                         \
                                                       const std::vector<PointTempl<RealType>> &,  \
                                                       std::vector<PointTempl<RealType>> &);       \
  template void FE<3, NEDELEC_ONE, RealType>::edge_reinit(                                         \
      ElemTempl<RealType> const *,                                                                 \
      unsigned int,                                                                                \
      Real,                                                                                        \
      const std::vector<PointTempl<RealType>> * const,                                             \
      const std::vector<Real> * const);                                                            \
  template void FEMapTempl<RealType>::compute_face_map(                                            \
      int, const std::vector<Real> &, const ElemTempl<RealType> *);                                \
  template void FEMapTempl<RealType>::compute_edge_map(                                            \
      int, const std::vector<Real> &, const ElemTempl<RealType> *);                                \
  INSTANTIATE_FE(0, RealType);                                                                     \
  INSTANTIATE_FE(1, RealType);                                                                     \
  INSTANTIATE_FE(2, RealType);                                                                     \
  INSTANTIATE_FE(3, RealType);                                                                     \
  INSTANTIATE_SUBDIVISION_FE(RealType);                                                            \
  INSTANTIATE_FE_INTERFACE_METHODS(RealType);                                                      \
  template struct FEShim<0, CLOUGH, RealType>;                                                     \
  template struct FEShim<1, CLOUGH, RealType>;                                                     \
  template struct FEShim<2, CLOUGH, RealType>;                                                     \
  template struct FEShim<3, CLOUGH, RealType>;                                                     \
  template struct FEShim<0, HERMITE, RealType>;                                                    \
  template struct FEShim<1, HERMITE, RealType>;                                                    \
  template struct FEShim<2, HERMITE, RealType>;                                                    \
  template struct FEShim<3, HERMITE, RealType>;                                                    \
  template struct FEShim<0, HIERARCHIC, RealType>;                                                 \
  template struct FEShim<1, HIERARCHIC, RealType>;                                                 \
  template struct FEShim<2, HIERARCHIC, RealType>;                                                 \
  template struct FEShim<3, HIERARCHIC, RealType>;                                                 \
  template struct FEShim<0, LAGRANGE, RealType>;                                                   \
  template struct FEShim<1, LAGRANGE, RealType>;                                                   \
  template struct FEShim<2, LAGRANGE, RealType>;                                                   \
  template struct FEShim<3, LAGRANGE, RealType>;                                                   \
  template struct FEShim<0, MONOMIAL, RealType>;                                                   \
  template struct FEShim<1, MONOMIAL, RealType>;                                                   \
  template struct FEShim<2, MONOMIAL, RealType>;                                                   \
  template struct FEShim<3, MONOMIAL, RealType>;                                                   \
  template class MeshBaseTempl<RealType>;                                                          \
  template class UnstructuredMeshTempl<RealType>;                                                  \
  template class DistributedMeshTempl<RealType>;                                                   \
  template class ReplicatedMeshTempl<RealType>;                                                    \
  template class PointLocatorBaseTempl<RealType>; \
  template class PointLocatorTreeTempl<RealType>; \
  template class PeriodicBoundariesTempl<RealType>; \
  template class PeriodicBoundaryBaseTempl<RealType>; \
  template class ElemTempl<RealType>;                                                              \
  template class BoundaryInfoTempl<RealType>;                                                      \
  template class Edge3Templ<RealType>;                                                             \
  template class Edge4Templ<RealType>;                                                             \
  template class Hex8Templ<RealType>;                                                              \
  template class Hex20Templ<RealType>;                                                             \
  template class Hex27Templ<RealType>;                                                             \
  template class Prism6Templ<RealType>;                                                            \
  template class Prism15Templ<RealType>;                                                           \
  template class Prism18Templ<RealType>;                                                           \
  template class Pyramid5Templ<RealType>;                                                          \
  template class Tet4Templ<RealType>;                                                              \
  template class Tet10Templ<RealType>;                                                             \
  template class Tri3Templ<RealType>;                                                              \
  template class Tri6Templ<RealType>;                                                              \
  template class Quad4Templ<RealType>;                                                             \
  template class Quad8Templ<RealType>;                                                             \
  template class Quad9Templ<RealType>;                                                             \
  template class HexTempl<RealType>;                                                               \
  template class PrismTempl<RealType>;                                                             \
  template class PyramidTempl<RealType>;                                                           \
  template class RemoteElemTempl<RealType>;                                                        \
  template class TriTempl<RealType>; \
  template class Tri3SubdivisionTempl<RealType>; \
  template class Pyramid13Templ<RealType>; \
  template class Pyramid14Templ<RealType>; \
  template class TreeNode<2,RealType>; \
  template class TreeNode<4,RealType>; \
  template class TreeNode<8,RealType>; \
  template class Tree<2,RealType>; \
  template class Tree<4,RealType>; \
  template class Tree<8,RealType>; \
  template class PartitionerTempl<RealType>; \
  template class CentroidPartitionerTempl<RealType>; \
  template class LinearPartitionerTempl<RealType>; \
  template class MappedSubdomainPartitionerTempl<RealType>; \
  template class MetisPartitionerTempl<RealType>; \
  template class ParmetisPartitionerTempl<RealType>; \
  template class SFCPartitionerTempl<RealType>; \
  template class SubdomainPartitionerTempl<RealType>

#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
#define INSTANTIATE_ALL_REAL_TYPE(RealType)                                                        \
  INSTANTIATE_ALL_REAL_TYPE0(RealType);                                                            \
  template struct FEShim<0, BERNSTEIN, RealType>;                                                  \
  template struct FEShim<1, BERNSTEIN, RealType>;                                                  \
  template struct FEShim<2, BERNSTEIN, RealType>;                                                  \
  template struct FEShim<3, BERNSTEIN, RealType>;                                                  \
  template struct FEShim<0, RATIONAL_BERNSTEIN, RealType>;                                         \
  template struct FEShim<1, RATIONAL_BERNSTEIN, RealType>;                                         \
  template struct FEShim<2, RATIONAL_BERNSTEIN, RealType>;                                         \
  template struct FEShim<3, RATIONAL_BERNSTEIN, RealType>;                                         \
  template struct FEShim<0, SZABAB, RealType>;                                                     \
  template struct FEShim<1, SZABAB, RealType>;                                                     \
  template struct FEShim<2, SZABAB, RealType>;                                                     \
  template struct FEShim<3, SZABAB, RealType>

#else
#define INSTANTIATE_ALL_REAL_TYPE(RealType) INSTANTIATE_ALL_REAL_TYPE0(RealType)
#endif

#endif // LIBMESH_INSTANTIATE_REAL_TYPE_H
