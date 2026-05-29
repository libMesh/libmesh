// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_KOKKOS_HILBERT_SYSTEM_H
#define LIBMESH_KOKKOS_HILBERT_SYSTEM_H

#include "libmesh/libmesh_common.h"

#ifdef LIBMESH_HAVE_KOKKOS

#include "kokkos_hilbert_assembly.h"
#include "kokkos_parsed_function.h"

#include <type_traits>

namespace libMesh::Kokkos::detail
{

template <typename Storage>
LIBMESH_DEVICE_INLINE decltype(auto)
storage_at(const Storage & storage,
           const unsigned int i)
{
  return storage(i);
}

template <typename T>
LIBMESH_DEVICE_INLINE const T &
storage_at(const T * storage,
           const unsigned int i)
{
  return storage[i];
}

template <typename T, unsigned int N>
struct StaticArrayAccess
{
  using value_type = T;
  T values[N] = {};
  unsigned int size = 0;

  LIBMESH_DEVICE_INLINE
  const T & operator()(const unsigned int i) const
  {
    return values[i];
  }
};

template <typename NodeCoordinateView, typename ElemNodeIdView>
class ElementNodeAccess
{
public:
  LIBMESH_DEVICE_INLINE
  ElementNodeAccess(NodeCoordinateView node_coordinates,
                    ElemNodeIdView element_node_ids,
                    const unsigned int elem_index)
    : _node_coordinates(node_coordinates),
      _element_node_ids(element_node_ids),
      _elem_index(elem_index)
  {
  }

  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int node,
                            const unsigned int component) const
  {
    return _node_coordinates(_element_node_ids(_elem_index, node), component);
  }

private:
  NodeCoordinateView _node_coordinates;
  ElemNodeIdView _element_node_ids;
  unsigned int _elem_index;
};

template <typename NodeCoordinateView, typename ElemNodeIdView>
LIBMESH_DEVICE_INLINE auto
make_element_node_access(NodeCoordinateView node_coordinates,
                         ElemNodeIdView element_node_ids,
                         const unsigned int elem_index)
{
  return ElementNodeAccess<NodeCoordinateView, ElemNodeIdView>(node_coordinates,
                                                               element_node_ids,
                                                               elem_index);
}

template <typename GlobalCoeffView, typename LocalIndexView>
class GatheredCoeffAccess
{
public:
  LIBMESH_DEVICE_INLINE
  GatheredCoeffAccess(GlobalCoeffView global_coeffs,
                      LocalIndexView local_indices,
                      const unsigned int elem_index)
    : _global_coeffs(global_coeffs),
      _local_indices(local_indices),
      _elem_index(elem_index)
  {
  }

  LIBMESH_DEVICE_INLINE
  decltype(auto) operator()(const unsigned int i) const
  {
    return storage_at(_global_coeffs, _local_indices(_elem_index, i));
  }

private:
  GlobalCoeffView _global_coeffs;
  LocalIndexView _local_indices;
  unsigned int _elem_index;
};

template <typename GlobalCoeffView, typename LocalIndexView>
LIBMESH_DEVICE_INLINE auto
make_gathered_coeff_access(GlobalCoeffView global_coeffs,
                           LocalIndexView local_indices,
                           const unsigned int elem_index)
{
  return GatheredCoeffAccess<GlobalCoeffView, LocalIndexView>(global_coeffs,
                                                              local_indices,
                                                              elem_index);
}

template <typename ResidualView, typename JacobianView>
struct DenseElementOutputSink
{
  ResidualView residual;
  JacobianView jacobian;
  unsigned int n_dofs = 0;
  bool request_jacobian = false;

  template <typename Accumulator>
  LIBMESH_DEVICE_INLINE
  void write(const Accumulator & accum) const
  {
    for (unsigned int i = 0; i != n_dofs; ++i)
      {
        residual(i) = accum.residual(i);
        if (request_jacobian)
          for (unsigned int j = 0; j != n_dofs; ++j)
            jacobian(i, j) = accum.jacobian(i, j);
      }
  }
};

template <typename ResidualView, typename JacobianView>
struct FlatDeviceValueSink
{
  ResidualView residual;
  JacobianView jacobian;
  unsigned int n_dofs = 0;

  template <typename Accumulator>
  LIBMESH_DEVICE_INLINE
  void write(const Accumulator & accum) const
  {
    for (unsigned int i = 0; i != n_dofs; ++i)
      {
        residual(i) = -accum.residual(i);
        for (unsigned int j = 0; j != n_dofs; ++j)
          jacobian(i * n_dofs + j) = accum.jacobian(i, j);
      }
  }
};

template <typename View>
struct DirectScatterAccess
{
  using execution_space = typename std::decay_t<View>::execution_space;

  View values;

  LIBMESH_DEVICE_INLINE
  void add(const std::size_t slot,
           const Number value) const
  {
    ::Kokkos::atomic_add(&values(slot), value);
  }
};

template <typename LocalView, typename RemoteView>
struct SplitScatterAccess
{
  using execution_space = typename std::decay_t<LocalView>::execution_space;

  LocalView local_values;
  RemoteView remote_values;
  std::size_t local_size = 0;

  LIBMESH_DEVICE_INLINE
  void add(const std::size_t slot,
           const Number value) const
  {
    if (slot < local_size)
      ::Kokkos::atomic_add(&local_values(slot), value);
    else
      ::Kokkos::atomic_add(&remote_values(slot - local_size), value);
  }
};

template <typename DiagView, typename OffdiagView, typename RemoteView>
struct SplitMatrixScatterAccess
{
  using execution_space = typename std::decay_t<DiagView>::execution_space;

  DiagView diag_values;
  OffdiagView offdiag_values;
  RemoteView remote_values;
  std::size_t diag_size = 0;
  std::size_t offdiag_base = 0;

  LIBMESH_DEVICE_INLINE
  void add(const std::size_t slot,
           const Number value) const
  {
    if (slot < diag_size)
      ::Kokkos::atomic_add(&diag_values(slot), value);
    else if (slot < offdiag_base)
      ::Kokkos::atomic_add(&offdiag_values(slot - diag_size), value);
    else
      ::Kokkos::atomic_add(&remote_values(slot - offdiag_base), value);
  }
};

template <typename ResidualScatterAccess,
          typename JacobianScatterAccess,
          typename SlotStorage>
struct DirectHilbertScatterAccumulator
{
  ResidualScatterAccess rhs_scatter;
  JacobianScatterAccess mat_scatter;
  SlotStorage rhs_slots;
  SlotStorage mat_slots;
  std::size_t rhs_offset = 0;
  std::size_t mat_offset = 0;
  unsigned int dofs = 0;

  LIBMESH_DEVICE_INLINE
  void add_residual(const unsigned int i,
                    const Number value) const
  {
    rhs_scatter.add(rhs_slots(rhs_offset + i), -value);
  }

  LIBMESH_DEVICE_INLINE
  void add_jacobian(const unsigned int i,
                    const unsigned int j,
                    const Number value) const
  {
    mat_scatter.add(mat_slots(mat_offset + i * dofs + j), value);
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_dofs() const
  {
    return dofs;
  }
};

struct ZeroCoeffAccess
{
  static constexpr bool is_zero_coeff_access = true;

  LIBMESH_DEVICE_INLINE
  Number operator()(const unsigned int) const
  {
    return Number(0);
  }
};

template <typename FieldKeyStorage,
          typename FieldDofStorage,
          typename GlobalCoeffView,
          typename LocalIndexView,
          typename GoalFunction,
          unsigned int MaxFieldVariables = 16>
class GatheredParsedFEMGoalAccess
{
public:
  LIBMESH_DEVICE_INLINE
  GatheredParsedFEMGoalAccess(FieldKeyStorage field_keys,
                              FieldDofStorage field_dofs,
                              GlobalCoeffView global_coeffs,
                              const LocalIndexView * field_local_indices,
                              GoalFunction goal)
    : _field_keys(field_keys),
      _field_dofs(field_dofs),
      _global_coeffs(global_coeffs),
      _goal(goal)
  {
    for (unsigned int i = 0; i != MaxFieldVariables; ++i)
      _field_local_indices[i] = field_local_indices[i];
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Number value(const QpData & qp_data, const Point & xyz) const
  {
    Number vars[LIBMESH_DIM + 1 + MaxFieldVariables] = {};
    fill_variables(qp_data, xyz, vars);
    return _goal.value(vars);
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Gradient gradient(const QpData & qp_data, const Point & xyz) const
  {
    Number vars[LIBMESH_DIM + 1 + MaxFieldVariables] = {};
    Gradient field_gradients[MaxFieldVariables];
    fill_variables(qp_data, xyz, vars);

    for (unsigned int field = 0; field != _goal.n_field_variables(); ++field)
      field_gradients[field] = sample_field_gradient(qp_data, field);

    return _goal.gradient(vars, field_gradients);
  }

private:
  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  void fill_variables(const QpData & qp_data,
                      const Point & xyz,
                      Number * vars) const
  {
    vars[0] = xyz(0);
#if LIBMESH_DIM > 1
    vars[1] = xyz(1);
#endif
#if LIBMESH_DIM > 2
    vars[2] = xyz(2);
#endif
    vars[LIBMESH_DIM] = _goal.time();

    for (unsigned int field = 0; field != _goal.n_field_variables(); ++field)
      vars[LIBMESH_DIM + 1 + field] = sample_field_value(qp_data, field);
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Number sample_field_value(const QpData & qp_data,
                            const unsigned int field) const
  {
    const auto & qp_ref = qp_data.reference_point();
    const auto field_key = _field_keys(field);
    const unsigned int n_dofs = _field_dofs(field);
    const unsigned int elem_index = qp_data.elem_index();

    Number value = 0.;
    for (unsigned int i = 0; i != n_dofs; ++i)
      value += storage_at(_global_coeffs, _field_local_indices[field](elem_index, i)) *
               shape(field_key, i, qp_ref(0), qp_ref(1), qp_ref(2));

    return value;
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Gradient sample_field_gradient(const QpData & qp_data,
                                 const unsigned int field) const
  {
    Gradient grad;
    grad.zero();

    const auto & qp_ref = qp_data.reference_point();
    const auto & Jinv = qp_data.inverse_jacobian();
    const auto field_key = _field_keys(field);
    const unsigned int n_dofs = _field_dofs(field);
    const unsigned int elem_index = qp_data.elem_index();

    for (unsigned int i = 0; i != n_dofs; ++i)
      grad.add_scaled(Jinv * grad_shape(field_key, i, qp_ref(0), qp_ref(1), qp_ref(2)),
                      storage_at(_global_coeffs, _field_local_indices[field](elem_index, i)));

    return grad;
  }

  FieldKeyStorage _field_keys;
  FieldDofStorage _field_dofs;
  GlobalCoeffView _global_coeffs;
  LocalIndexView _field_local_indices[MaxFieldVariables] = {};
  GoalFunction _goal;
};

template <unsigned int MaxDofs,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename GoalAccess,
          typename CoeffStorage,
          typename Sink>
bool
run_hilbert_system_assembly(const libMesh::FEShapeKey key,
                            const libMesh::ElemMappingType mapping_type,
                            const NodeCoordinateStorage & node_coordinates,
                            const ElemNodeIdStorage & element_node_ids,
                            const unsigned int elem_index,
                            const unsigned int n_nodes,
                            const unsigned int quadrature_order,
                            const unsigned int hilbert_order,
                            const CoeffStorage & coeff,
                            const Number solution_derivative,
                            GoalAccess goal_access,
                            const bool request_jacobian,
                            const Sink & sink,
                            const char * const kernel_name)
{
  if (sink.n_dofs > MaxDofs)
    return false;

  if (!supports_hilbert_local_assembly(key, mapping_type, quadrature_order))
    return false;

  const auto elem_nodes =
    make_element_node_access(node_coordinates, element_node_ids, elem_index);
  const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(key,
                                                                          mapping_type,
                                                                          elem_nodes,
                                                                          n_nodes,
                                                                          quadrature_order,
                                                                          elem_index);

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<>(0, 1),
    KOKKOS_LAMBDA(const int) {
      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, coeff, solution_derivative);
      libMesh::Kokkos::detail::LocalHilbertAccumulator<MaxDofs> accum(sink.n_dofs);
      libMesh::detail::assemble_hilbert_element(fe,
                                                solution,
                                                goal_access,
                                                request_jacobian,
                                                hilbert_order,
                                                accum);
      sink.write(accum);
    });

  return true;
}

template <unsigned int MaxDofs,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage,
          typename ShapeKeyStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename QuadratureOrderStorage,
          typename OffsetStorage,
          typename GoalAccess,
          typename ResidualView,
          typename JacobianView>
void
run_hilbert_system_value_batch(const NodeCoordinateStorage & node_coordinates,
                               const ElemNodeIdStorage & element_node_ids,
                               const ElemMappingTypeStorage & element_mapping_types,
                               const ElemNodeCountStorage & element_n_nodes,
                               const ShapeKeyStorage & shape_keys,
                               const ElemIndexStorage & elem_indices,
                               const ElemDofCountStorage & elem_n_dofs,
                               const QuadratureOrderStorage & quadrature_orders,
                               const OffsetStorage & rhs_offsets,
                               const OffsetStorage & mat_offsets,
                               const unsigned int hilbert_order,
                               GoalAccess goal_access,
                               ResidualView rhs_values,
                               JacobianView mat_values,
                               const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);
  using ExecutionSpace = typename std::decay_t<ResidualView>::execution_space;

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<ExecutionSpace>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);
      const auto key = shape_keys(record_index);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key,
        element_mapping_types(elem_index),
        elem_nodes,
        element_n_nodes(elem_index),
        quadrature_orders(record_index),
        elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      libMesh::Kokkos::detail::LocalHilbertAccumulator<MaxDofs> accum(n_dofs);
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);

      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          rhs_values(rhs_offset + i) = -accum.residual(i);
          for (unsigned int j = 0; j != n_dofs; ++j)
            mat_values(mat_offset + i * n_dofs + j) = accum.jacobian(i, j);
        }
    });
}

template <unsigned int MaxDofs,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename OffsetStorage,
          typename SlotStorage,
          typename GoalAccess,
          typename ResidualScatterAccess,
          typename JacobianScatterAccess>
void
run_hilbert_system_bucket_scatter_batch(const libMesh::FEShapeKey key,
                                        const libMesh::ElemMappingType mapping_type,
                                        const unsigned int n_nodes,
                                        const unsigned int quadrature_order,
                                        const NodeCoordinateStorage & node_coordinates,
                                        const ElemNodeIdStorage & element_node_ids,
                                        const ElemIndexStorage & elem_indices,
                                        const ElemDofCountStorage & elem_n_dofs,
                                        const OffsetStorage & rhs_offsets,
                                        const OffsetStorage & mat_offsets,
                                        const SlotStorage & rhs_slots,
                                        const SlotStorage & mat_slots,
                                        const unsigned int hilbert_order,
                                        GoalAccess goal_access,
                                        ResidualScatterAccess rhs_scatter,
                                        JacobianScatterAccess mat_scatter,
                                        const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);
  using ExecutionSpace = typename std::decay_t<ResidualScatterAccess>::execution_space;

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<ExecutionSpace>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key, mapping_type, elem_nodes, n_nodes, quadrature_order, elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      DirectHilbertScatterAccumulator<ResidualScatterAccess,
                                      JacobianScatterAccess,
                                      SlotStorage> accum{
        rhs_scatter, mat_scatter, rhs_slots, mat_slots, rhs_offset, mat_offset, n_dofs};
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);
    });
}

template <unsigned int MaxDofs,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename OffsetStorage,
          typename GoalAccess,
          typename ResidualView,
          typename JacobianView>
void
run_hilbert_system_bucket_value_batch(const libMesh::FEShapeKey key,
                                      const libMesh::ElemMappingType mapping_type,
                                      const unsigned int n_nodes,
                                      const unsigned int quadrature_order,
                                      const NodeCoordinateStorage & node_coordinates,
                                      const ElemNodeIdStorage & element_node_ids,
                                      const ElemIndexStorage & elem_indices,
                                      const ElemDofCountStorage & elem_n_dofs,
                                      const OffsetStorage & rhs_offsets,
                                      const OffsetStorage & mat_offsets,
                                      const unsigned int hilbert_order,
                                      GoalAccess goal_access,
                                      ResidualView rhs_values,
                                      JacobianView mat_values,
                                      const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);
  using ExecutionSpace = typename std::decay_t<ResidualView>::execution_space;

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<ExecutionSpace>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key, mapping_type, elem_nodes, n_nodes, quadrature_order, elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      libMesh::Kokkos::detail::LocalHilbertAccumulator<MaxDofs> accum(n_dofs);
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);

      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          rhs_values(rhs_offset + i) = -accum.residual(i);
          for (unsigned int j = 0; j != n_dofs; ++j)
            mat_values(mat_offset + i * n_dofs + j) = accum.jacobian(i, j);
        }
    });
}

template <unsigned int MaxDofs,
          unsigned int MaxFieldVariables,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage,
          typename ShapeKeyStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename QuadratureOrderStorage,
          typename OffsetStorage,
          typename FieldKeyRecordStorage,
          typename FieldDofRecordStorage,
          typename FieldLocalIndexStorage,
          typename GlobalCoeffStorage,
          typename GoalFunction,
          typename ResidualView,
          typename JacobianView>
void
run_hilbert_system_fem_value_batch(const NodeCoordinateStorage & node_coordinates,
                                   const ElemNodeIdStorage & element_node_ids,
                                   const ElemMappingTypeStorage & element_mapping_types,
                                   const ElemNodeCountStorage & element_n_nodes,
                                   const ShapeKeyStorage & shape_keys,
                                   const ElemIndexStorage & elem_indices,
                                   const ElemDofCountStorage & elem_n_dofs,
                                   const QuadratureOrderStorage & quadrature_orders,
                                   const OffsetStorage & rhs_offsets,
                                   const OffsetStorage & mat_offsets,
                                   const FieldKeyRecordStorage & field_keys,
                                   const FieldDofRecordStorage & field_dofs,
                                   const FieldLocalIndexStorage & field_local_indices,
                                   const GlobalCoeffStorage & global_coeffs,
                                   GoalFunction goal_function,
                                   const unsigned int hilbert_order,
                                   ResidualView rhs_values,
                                   JacobianView mat_values,
                                   const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);
  using ExecutionSpace = typename std::decay_t<ResidualView>::execution_space;

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<ExecutionSpace>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);
      const auto key = shape_keys(record_index);

      StaticArrayAccess<libMesh::FEShapeKey, MaxFieldVariables> record_field_keys;
      StaticArrayAccess<unsigned int, MaxFieldVariables> record_field_dofs;
      record_field_keys.size = goal_function.n_field_variables();
      record_field_dofs.size = goal_function.n_field_variables();
      for (unsigned int field = 0; field != goal_function.n_field_variables(); ++field)
        {
          record_field_keys.values[field] = field_keys(field, record_index);
          record_field_dofs.values[field] = field_dofs(field, record_index);
        }

      const auto goal_access =
        GatheredParsedFEMGoalAccess<decltype(record_field_keys),
                                    decltype(record_field_dofs),
                                    GlobalCoeffStorage,
                                    typename FieldLocalIndexStorage::value_type,
                                    GoalFunction,
                                    MaxFieldVariables>(record_field_keys,
                                                       record_field_dofs,
                                                       global_coeffs,
                                                       field_local_indices.values,
                                                       goal_function);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key,
        element_mapping_types(elem_index),
        elem_nodes,
        element_n_nodes(elem_index),
        quadrature_orders(record_index),
        elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      libMesh::Kokkos::detail::LocalHilbertAccumulator<MaxDofs> accum(n_dofs);
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);

      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          rhs_values(rhs_offset + i) = -accum.residual(i);
          for (unsigned int j = 0; j != n_dofs; ++j)
            mat_values(mat_offset + i * n_dofs + j) = accum.jacobian(i, j);
        }
    });
}

template <unsigned int MaxDofs,
          unsigned int MaxFieldVariables,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename OffsetStorage,
          typename SlotStorage,
          typename FieldKeyStorage,
          typename FieldDofStorage,
          typename FieldLocalIndexStorage,
          typename GlobalCoeffStorage,
          typename GoalFunction,
          typename ResidualScatterAccess,
          typename JacobianScatterAccess>
void
run_hilbert_system_fem_bucket_scatter_batch(const libMesh::FEShapeKey key,
                                            const libMesh::ElemMappingType mapping_type,
                                            const unsigned int n_nodes,
                                            const unsigned int quadrature_order,
                                            const NodeCoordinateStorage & node_coordinates,
                                            const ElemNodeIdStorage & element_node_ids,
                                            const ElemIndexStorage & elem_indices,
                                            const ElemDofCountStorage & elem_n_dofs,
                                            const OffsetStorage & rhs_offsets,
                                            const OffsetStorage & mat_offsets,
                                            const SlotStorage & rhs_slots,
                                            const SlotStorage & mat_slots,
                                            FieldKeyStorage field_keys,
                                            FieldDofStorage field_dofs,
                                            const FieldLocalIndexStorage & field_local_indices,
                                            const GlobalCoeffStorage & global_coeffs,
                                            GoalFunction goal_function,
                                            const unsigned int hilbert_order,
                                            ResidualScatterAccess rhs_scatter,
                                            JacobianScatterAccess mat_scatter,
                                            const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);

      const auto goal_access =
        GatheredParsedFEMGoalAccess<FieldKeyStorage,
                                    FieldDofStorage,
                                    GlobalCoeffStorage,
                                    typename FieldLocalIndexStorage::value_type,
                                    GoalFunction,
                                    MaxFieldVariables>(field_keys,
                                                       field_dofs,
                                                       global_coeffs,
                                                       field_local_indices.values,
                                                       goal_function);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key, mapping_type, elem_nodes, n_nodes, quadrature_order, elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      DirectHilbertScatterAccumulator<ResidualScatterAccess,
                                      JacobianScatterAccess,
                                      SlotStorage> accum{
        rhs_scatter, mat_scatter, rhs_slots, mat_slots, rhs_offset, mat_offset, n_dofs};
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);
    });
}

template <unsigned int MaxDofs,
          unsigned int MaxFieldVariables,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemIndexStorage,
          typename ElemDofCountStorage,
          typename OffsetStorage,
          typename FieldKeyStorage,
          typename FieldDofStorage,
          typename FieldLocalIndexStorage,
          typename GlobalCoeffStorage,
          typename GoalFunction,
          typename ResidualView,
          typename JacobianView>
void
run_hilbert_system_fem_bucket_value_batch(const libMesh::FEShapeKey key,
                                          const libMesh::ElemMappingType mapping_type,
                                          const unsigned int n_nodes,
                                          const unsigned int quadrature_order,
                                          const NodeCoordinateStorage & node_coordinates,
                                          const ElemNodeIdStorage & element_node_ids,
                                          const ElemIndexStorage & elem_indices,
                                          const ElemDofCountStorage & elem_n_dofs,
                                          const OffsetStorage & rhs_offsets,
                                          const OffsetStorage & mat_offsets,
                                          FieldKeyStorage field_keys,
                                          FieldDofStorage field_dofs,
                                          const FieldLocalIndexStorage & field_local_indices,
                                          const GlobalCoeffStorage & global_coeffs,
                                          GoalFunction goal_function,
                                          const unsigned int hilbert_order,
                                          ResidualView rhs_values,
                                          JacobianView mat_values,
                                          const char * const kernel_name)
{
  const auto n_records = elem_indices.extent(0);

  ::Kokkos::parallel_for(
    kernel_name,
    ::Kokkos::RangePolicy<>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index) {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);

      const auto goal_access =
        GatheredParsedFEMGoalAccess<FieldKeyStorage,
                                    FieldDofStorage,
                                    GlobalCoeffStorage,
                                    typename FieldLocalIndexStorage::value_type,
                                    GoalFunction,
                                    MaxFieldVariables>(field_keys,
                                                       field_dofs,
                                                       global_coeffs,
                                                       field_local_indices.values,
                                                       goal_function);

      const auto elem_nodes =
        make_element_node_access(node_coordinates, element_node_ids, elem_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> fe(
        key, mapping_type, elem_nodes, n_nodes, quadrature_order, elem_index);

      const auto solution =
        libMesh::Kokkos::detail::make_hilbert_solution_access(fe, ZeroCoeffAccess{}, Number(1.));
      libMesh::Kokkos::detail::LocalHilbertAccumulator<MaxDofs> accum(n_dofs);
      libMesh::detail::assemble_hilbert_element(
        fe, solution, goal_access, true, hilbert_order, n_dofs, accum);

      const auto rhs_offset = rhs_offsets(record_index);
      const auto mat_offset = mat_offsets(record_index);
      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          rhs_values(rhs_offset + i) = -accum.residual(i);
          for (unsigned int j = 0; j != n_dofs; ++j)
            mat_values(mat_offset + i * n_dofs + j) = accum.jacobian(i, j);
        }
    });
}

} // namespace libMesh::Kokkos::detail

#endif // LIBMESH_HAVE_KOKKOS

#endif // LIBMESH_KOKKOS_HILBERT_SYSTEM_H
