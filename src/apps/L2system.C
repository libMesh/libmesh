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

#include "L2system.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/linear_solver.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

#include <array>
#include <chrono>
#include <type_traits>

#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_vector.h"
#endif

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
#include "../../include/gpu/kokkos_hilbert_system.h"
#include "libmesh/fe_shape_traits.h"

#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#undef __CUDACC_VER__
#endif

using namespace libMesh;

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
constexpr unsigned int kokkos_hilbert_max_dofs = 27;
constexpr unsigned int kokkos_parsed_fem_max_fields = 16;

template <typename T>
using KokkosScalarView = ::Kokkos::View<T *>;

using KokkosDenseJacobianView = ::Kokkos::View<Number **>;
using KokkosFlatJacobianView = ::Kokkos::View<Number *>;
using KokkosUnsignedIntView = ::Kokkos::View<unsigned int *>;
using KokkosSizeView = ::Kokkos::View<std::size_t *>;
using KokkosFieldKeyRecordView = ::Kokkos::View<FEShapeKey **>;
using KokkosFieldDofRecordView = ::Kokkos::View<unsigned int **>;
using KokkosFieldKeyStorage =
  libMesh::Kokkos::detail::StaticArrayAccess<FEShapeKey, kokkos_parsed_fem_max_fields>;
using KokkosFieldDofStorage =
  libMesh::Kokkos::detail::StaticArrayAccess<unsigned int, kokkos_parsed_fem_max_fields>;
using KokkosLocalIndexView = DofMap::KokkosLocalIndexCache::elem_local_index_view;

struct HilbertElementAssemblyRecord
{
  dof_id_type elem_id = DofObject::invalid_id;
  std::size_t rhs_offset = 0;
  std::size_t mat_offset = 0;
  unsigned int elem_index = libMesh::invalid_uint;
  unsigned int quadrature_order = 0;
  unsigned int n_dofs = 0;
};

struct KokkosHilbertAssemblyBucket
{
  std::size_t begin = 0;
  std::size_t end = 0;
  FEShapeKey key;
  ElemType elem_type = INVALID_ELEM;
  ElemMappingType mapping_type = LAGRANGE_MAP;
  unsigned int n_nodes = 0;
  unsigned int elem_p_level = 0;
  unsigned int quadrature_order = 0;
};

FEShapeKey
make_hilbert_shape_key(const Elem & elem,
                       const FEType & fe_type)
{
  return {fe_type.family,
          elem.type(),
          static_cast<Order>(fe_type.order.get_order() + cast_int<int>(elem.p_level()))};
}

FEShapeKey
make_hilbert_shape_key(const ElemType elem_type,
                       const unsigned int elem_p_level,
                       const FEType & fe_type)
{
  return {fe_type.family,
          elem_type,
          static_cast<Order>(fe_type.order.get_order() + cast_int<int>(elem_p_level))};
}

void
accumulate_hilbert_dense_outputs(const KokkosScalarView<Number> & d_F,
                                 const KokkosDenseJacobianView & d_K,
                                 const bool request_jacobian,
                                 DenseSubVector<Number> & F,
                                 DenseSubMatrix<Number> & K)
{
  auto h_F = ::Kokkos::create_mirror_view(d_F);
  ::Kokkos::deep_copy(h_F, d_F);

  for (unsigned int i = 0; i != F.size(); ++i)
    F(i) += h_F(i);

  if (!request_jacobian)
    return;

  auto h_K = ::Kokkos::create_mirror_view(d_K);
  ::Kokkos::deep_copy(h_K, d_K);

  for (unsigned int i = 0; i != F.size(); ++i)
    for (unsigned int j = 0; j != F.size(); ++j)
      K(i, j) += h_K(i, j);
}

void
build_hilbert_element_records(HilbertSystem & sys,
                              std::vector<HilbertElementAssemblyRecord> & records,
                              std::size_t & total_rhs_entries,
                              std::size_t & total_mat_entries)
{
  DofMap & dof_map = sys.get_dof_map();
  const auto * dof_index_cache = dof_map.get_kokkos_dof_index_cache(0);
  libmesh_assert(dof_index_cache);
  total_rhs_entries = 0;
  total_mat_entries = 0;

  for (auto elem_index : index_range(dof_index_cache->host_element_ids))
    {
      if (!sys.subdomains_list().empty() &&
          !sys.subdomains_list().count(dof_index_cache->host_element_subdomains[elem_index]))
        continue;

      const unsigned int n_dofs = dof_index_cache->host_element_n_dofs[elem_index];
      if (!n_dofs)
        continue;

      HilbertElementAssemblyRecord record;
      record.elem_id = dof_index_cache->host_element_ids[elem_index];
      record.rhs_offset = total_rhs_entries;
      record.mat_offset = total_mat_entries;
      record.elem_index = cast_int<unsigned int>(elem_index);
      record.n_dofs = n_dofs;
      total_rhs_entries += n_dofs;
      total_mat_entries += n_dofs * n_dofs;
      records.push_back(std::move(record));
    }
}

bool
build_hilbert_record_quadrature_orders(HilbertSystem & sys,
                                       std::vector<HilbertElementAssemblyRecord> & records)
{
  const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
  const FEType fe_type = sys.variable_type(0);
  const int quadrature_order =
    cast_int<int>(fe_type.default_quadrature_order()) + sys.extra_quadrature_order;

  for (auto & record : records)
    {
      libmesh_error_msg_if(quadrature_order < 0,
                           "Negative quadrature order is not supported for Kokkos Hilbert assembly");
      record.quadrature_order = cast_int<unsigned int>(quadrature_order);

      const FEShapeKey shape_key =
        make_hilbert_shape_key(geometry_cache.element_types(record.elem_index),
                               geometry_cache.element_p_levels(record.elem_index),
                               fe_type);
      if (!libMesh::Kokkos::detail::supports_hilbert_local_assembly(
             shape_key,
             geometry_cache.element_mapping_types(record.elem_index),
             record.quadrature_order) ||
          record.n_dofs > kokkos_hilbert_max_dofs)
        return false;
    }

  return true;
}

template <typename ParsedGoal, typename CachedGoal, typename HostGoalPtr, typename BuildCache>
const CachedGoal *
ensure_kokkos_goal_cache(std::unique_ptr<CachedGoal> & cache,
                         const HostGoalPtr & host_goal,
                         BuildCache && build_cache)
{
  if (cache)
    return cache.get();

  if (!host_goal)
    return nullptr;

  const auto * parsed_goal = dynamic_cast<const ParsedGoal *>(host_goal.get());
  if (!parsed_goal)
    return nullptr;

  cache = build_cache(*parsed_goal);
  return cache.get();
}

void
prewarm_kokkos_hilbert_entities(HilbertSystem & sys,
                                const libMesh::Kokkos::KokkosParsedFEMFunction<Number> * fem_goal)
{
  if (sys.current_local_solution)
    sys.get_dof_map().prepare_kokkos_local_index_cache(*sys.current_local_solution, 0);

  if (!fem_goal || !sys.input_system || !sys.input_system->current_local_solution)
    return;

  for (unsigned int field = 0; field != fem_goal->n_field_variables(); ++field)
    sys.input_system->get_dof_map().prepare_kokkos_local_index_cache(
      *sys.input_system->current_local_solution,
      fem_goal->field_variable_number(field));
}

#if defined(LIBMESH_HAVE_PETSC)
void
build_hilbert_coo_indices(const DofMap::KokkosDofIndexCache & dof_index_cache,
                          const std::vector<HilbertElementAssemblyRecord> & records,
                          std::vector<PetscInt> & rhs_rows,
                          std::vector<PetscInt> & mat_rows,
                          std::vector<PetscInt> & mat_cols)
{
  for (const auto & record : records)
    {
      const unsigned int n_dofs = record.n_dofs;

      for (unsigned int i = 0; i != n_dofs; ++i)
        rhs_rows[record.rhs_offset + i] =
          cast_int<PetscInt>(
            dof_index_cache.host_element_dof_indices[record.elem_index * dof_index_cache.max_dofs + i]);

      for (unsigned int i = 0; i != n_dofs; ++i)
        for (unsigned int j = 0; j != n_dofs; ++j)
          {
            const std::size_t offset = record.mat_offset + i * n_dofs + j;
            mat_rows[offset] = cast_int<PetscInt>(
              dof_index_cache.host_element_dof_indices[record.elem_index * dof_index_cache.max_dofs +
                                                       i]);
            mat_cols[offset] = cast_int<PetscInt>(
              dof_index_cache.host_element_dof_indices[record.elem_index * dof_index_cache.max_dofs +
                                                       j]);
          }
    }
}

class HostExactParsedFEMGoalAccess
{
public:
  HostExactParsedFEMGoalAccess(const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal,
                               FEMContext & input_context)
    : _goal(goal),
      _input_context(input_context)
  {
  }

  template <typename QpData>
  Number value(const QpData & qp_data, const Point & xyz) const
  {
    Number vars[LIBMESH_DIM + 1 + kokkos_parsed_fem_max_fields] = {};
    fill_variables(qp_data, xyz, vars);
    return _goal.value(vars);
  }

  template <typename QpData>
  Gradient gradient(const QpData & qp_data, const Point & xyz) const
  {
    Number vars[LIBMESH_DIM + 1 + kokkos_parsed_fem_max_fields] = {};
    Gradient field_gradients[kokkos_parsed_fem_max_fields];
    fill_variables(qp_data, xyz, vars);

    for (unsigned int field = 0; field != _goal.n_field_variables(); ++field)
      field_gradients[field] =
        _input_context.interior_gradient(_goal.field_variable_number(field), qp_data.qp_index());

    return _goal.gradient(vars, field_gradients);
  }

private:
  template <typename QpData>
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
      vars[LIBMESH_DIM + 1 + field] =
        _input_context.interior_value(_goal.field_variable_number(field), qp_data.qp_index());
  }

  const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & _goal;
  FEMContext & _input_context;
};

struct KokkosElementAssemblyState
{
  const MeshBase::KokkosGeometryCache * geometry_cache = nullptr;
  const DofMap::KokkosLocalIndexCache * solution_local_indices = nullptr;
  PetscVector<Number> * solution_vector = nullptr;
  unsigned int elem_index = libMesh::invalid_uint;
  ElemType elem_type = INVALID_ELEM;
  ElemMappingType mapping_type = LAGRANGE_MAP;
  unsigned int elem_n_nodes = 0;
  unsigned int elem_p_level = 0;
};

struct KokkosFEMGoalData
{
  KokkosFieldKeyStorage field_keys;
  KokkosFieldDofStorage field_dofs;
  std::array<KokkosLocalIndexView, kokkos_parsed_fem_max_fields> field_local_indices;
  PetscVector<Number> * input_vector = nullptr;
};

struct KokkosHilbertBatchData
{
  KokkosUnsignedIntView elem_indices;
  KokkosUnsignedIntView elem_n_dofs;
  KokkosUnsignedIntView quadrature_orders;
  KokkosSizeView rhs_offsets;
  KokkosSizeView mat_offsets;
};

struct KokkosFEMGoalBatchData
{
  std::vector<KokkosFieldKeyStorage> bucket_field_keys;
  std::vector<KokkosFieldDofStorage> bucket_field_dofs;
  libMesh::Kokkos::detail::StaticArrayAccess<KokkosLocalIndexView,
                                             kokkos_parsed_fem_max_fields>
    field_local_indices;
  PetscVector<Number> * input_vector = nullptr;
};

struct KokkosPetscAssemblyPlan
{
  std::vector<HilbertElementAssemblyRecord> records;
  std::vector<KokkosHilbertAssemblyBucket> buckets;
  std::vector<PetscInt> rhs_rows;
  std::vector<PetscInt> mat_rows;
  std::vector<PetscInt> mat_cols;
  KokkosHilbertBatchData batch_data;
  ::Kokkos::View<Number *> rhs_values;
  ::Kokkos::View<Number *> mat_values;
  KokkosFEMGoalBatchData fem_goal_batch_data;
  const void * geometry_cache_id = nullptr;
  const void * dof_index_cache_id = nullptr;
  FEType fe_type;
  unsigned int hilbert_order = 0;
  int extra_quadrature_order = 0;
  std::set<subdomain_id_type> subdomains;
  const void * matrix_target = nullptr;
  const void * rhs_target = nullptr;
  const void * fem_goal_target = nullptr;
  const void * input_vector_target = nullptr;
};

namespace
{

bool
build_kokkos_element_state(HilbertSystem & sys,
                           const Elem & elem,
                           KokkosElementAssemblyState & state)
{
  auto * solution_vector = dynamic_cast<PetscVector<Number> *>(sys.current_local_solution.get());
  if (!solution_vector || !solution_vector->supports_kokkos_access())
    return false;

  const DofMap & dof_map = sys.get_dof_map();
  state.geometry_cache = &sys.get_mesh().get_kokkos_geometry_cache();
  state.elem_index = sys.get_mesh().get_kokkos_elem_index(elem);
  if (state.elem_index == libMesh::invalid_uint)
    return false;

  state.elem_type = state.geometry_cache->element_types(state.elem_index);
  state.mapping_type = state.geometry_cache->element_mapping_types(state.elem_index);
  state.elem_n_nodes = state.geometry_cache->element_n_nodes(state.elem_index);
  state.elem_p_level = state.geometry_cache->element_p_levels(state.elem_index);

  state.solution_local_indices =
    dof_map.require_kokkos_local_index_cache(*sys.current_local_solution, 0);
  if (!state.solution_local_indices)
    return false;

  state.solution_vector = solution_vector;
  return true;
}

bool
build_kokkos_fem_goal_data(const ElemType elem_type,
                           const ElemMappingType mapping_type,
                           const unsigned int elem_p_level,
                           const unsigned int hilbert_order,
                           System & input_system,
                           const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                           KokkosFEMGoalData & goal_data)
{
  const unsigned int n_fields = goal_function.n_field_variables();
  if (n_fields > kokkos_parsed_fem_max_fields)
    return false;

  auto * input_vector = dynamic_cast<PetscVector<Number> *>(input_system.current_local_solution.get());
  if (!input_vector || !input_vector->supports_kokkos_access())
    return false;

  goal_data.field_keys.size = n_fields;
  goal_data.field_dofs.size = n_fields;
  goal_data.input_vector = input_vector;

  for (unsigned int field = 0; field != n_fields; ++field)
    {
      const unsigned int var_num = goal_function.field_variable_number(field);
      const FEType field_type = input_system.variable_type(var_num);
      const FEShapeKey field_key = make_hilbert_shape_key(elem_type, elem_p_level, field_type);

      if (mapping_type != LAGRANGE_MAP ||
          !supports_shape_with_lagrange_map(field_key) ||
          (hilbert_order > 0 && !supports_grad_shape(field_key)))
        return false;

      const unsigned int field_n_dofs =
        FEInterface::n_dofs(libMesh::Kokkos::dim_from_topology(elem_type), field_type, elem_type);
      if (field_n_dofs > kokkos_hilbert_max_dofs)
        return false;

      const auto * local_index_cache =
        input_system.get_dof_map().require_kokkos_local_index_cache(
          *input_system.current_local_solution, var_num);
      if (!local_index_cache)
        return false;

      goal_data.field_keys.values[field] = field_key;
      goal_data.field_dofs.values[field] = field_n_dofs;
      goal_data.field_local_indices[field] = local_index_cache->element_local_indices;
    }

  return true;
}

auto
make_hilbert_bucket_sort_key(const MeshBase::KokkosGeometryCache & geometry_cache,
                             const FEType & fe_type,
                             const HilbertElementAssemblyRecord & record)
{
  const unsigned int elem_index = record.elem_index;
  return std::make_tuple(cast_int<int>(geometry_cache.element_types(elem_index)),
                         cast_int<int>(geometry_cache.element_mapping_types(elem_index)),
                         cast_int<int>(fe_type.order.get_order() +
                                       cast_int<int>(geometry_cache.element_p_levels(elem_index))),
                         record.quadrature_order,
                         geometry_cache.element_n_nodes(elem_index));
}

void
sort_hilbert_element_records(const MeshBase::KokkosGeometryCache & geometry_cache,
                             const FEType & fe_type,
                             std::vector<HilbertElementAssemblyRecord> & records)
{
  std::sort(records.begin(),
            records.end(),
            [&geometry_cache, &fe_type](const auto & lhs, const auto & rhs)
            {
              return make_hilbert_bucket_sort_key(geometry_cache, fe_type, lhs) <
                     make_hilbert_bucket_sort_key(geometry_cache, fe_type, rhs);
            });
}

void
build_hilbert_assembly_buckets(const MeshBase::KokkosGeometryCache & geometry_cache,
                               const FEType & fe_type,
                               const std::vector<HilbertElementAssemblyRecord> & records,
                               std::vector<KokkosHilbertAssemblyBucket> & buckets)
{
  buckets.clear();
  if (records.empty())
    return;

  auto fill_bucket = [&geometry_cache, &fe_type, &records](KokkosHilbertAssemblyBucket & bucket,
                                                            const std::size_t begin,
                                                            const std::size_t end)
  {
    bucket.begin = begin;
    bucket.end = end;
    const auto & record = records[begin];
    const unsigned int elem_index = record.elem_index;
    bucket.key = make_hilbert_shape_key(geometry_cache.element_types(elem_index),
                                        geometry_cache.element_p_levels(elem_index),
                                        fe_type);
    bucket.elem_type = geometry_cache.element_types(elem_index);
    bucket.mapping_type = geometry_cache.element_mapping_types(elem_index);
    bucket.n_nodes = geometry_cache.element_n_nodes(elem_index);
    bucket.elem_p_level = geometry_cache.element_p_levels(elem_index);
    bucket.quadrature_order = record.quadrature_order;
  };

  std::size_t bucket_begin = 0;
  auto current_key = make_hilbert_bucket_sort_key(geometry_cache, fe_type, records.front());
  for (std::size_t i = 1; i != records.size(); ++i)
    {
      const auto next_key = make_hilbert_bucket_sort_key(geometry_cache, fe_type, records[i]);
      if (next_key == current_key)
        continue;

      buckets.emplace_back();
      fill_bucket(buckets.back(), bucket_begin, i);
      bucket_begin = i;
      current_key = next_key;
    }

  buckets.emplace_back();
  fill_bucket(buckets.back(), bucket_begin, records.size());
}

void
build_hilbert_batch_data(const std::vector<HilbertElementAssemblyRecord> & records,
                         KokkosHilbertBatchData & batch_data)
{
  batch_data.elem_indices = KokkosUnsignedIntView("hilbert_elem_indices", records.size());
  batch_data.elem_n_dofs = KokkosUnsignedIntView("hilbert_elem_n_dofs", records.size());
  batch_data.quadrature_orders =
    KokkosUnsignedIntView("hilbert_quadrature_orders", records.size());
  batch_data.rhs_offsets = KokkosSizeView("hilbert_rhs_offsets", records.size());
  batch_data.mat_offsets = KokkosSizeView("hilbert_mat_offsets", records.size());

  auto h_elem_indices = ::Kokkos::create_mirror_view(batch_data.elem_indices);
  auto h_elem_n_dofs = ::Kokkos::create_mirror_view(batch_data.elem_n_dofs);
  auto h_quadrature_orders = ::Kokkos::create_mirror_view(batch_data.quadrature_orders);
  auto h_rhs_offsets = ::Kokkos::create_mirror_view(batch_data.rhs_offsets);
  auto h_mat_offsets = ::Kokkos::create_mirror_view(batch_data.mat_offsets);

  for (auto record_index : index_range(records))
    {
      const auto & record = records[record_index];
      h_elem_indices(record_index) = record.elem_index;
      h_elem_n_dofs(record_index) = record.n_dofs;
      h_quadrature_orders(record_index) = record.quadrature_order;
      h_rhs_offsets(record_index) = record.rhs_offset;
      h_mat_offsets(record_index) = record.mat_offset;
    }

  ::Kokkos::deep_copy(batch_data.elem_indices, h_elem_indices);
  ::Kokkos::deep_copy(batch_data.elem_n_dofs, h_elem_n_dofs);
  ::Kokkos::deep_copy(batch_data.quadrature_orders, h_quadrature_orders);
  ::Kokkos::deep_copy(batch_data.rhs_offsets, h_rhs_offsets);
  ::Kokkos::deep_copy(batch_data.mat_offsets, h_mat_offsets);
}

bool
kokkos_petsc_plan_matches(const HilbertSystem & sys,
                          const MeshBase::KokkosGeometryCache & geometry_cache,
                          const DofMap::KokkosDofIndexCache & dof_index_cache,
                          const KokkosPetscAssemblyPlan & plan)
{
  return plan.geometry_cache_id == &geometry_cache &&
         plan.dof_index_cache_id == &dof_index_cache &&
         plan.fe_type == sys.variable_type(0) &&
         plan.hilbert_order == sys.hilbert_order() &&
         plan.extra_quadrature_order == sys.extra_quadrature_order &&
         plan.subdomains == sys.subdomains_list();
}

bool
build_kokkos_petsc_assembly_plan(HilbertSystem & sys,
                                 KokkosPetscAssemblyPlan & plan)
{
  const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
  const auto * dof_index_cache = sys.get_dof_map().get_kokkos_dof_index_cache(0);
  if (!dof_index_cache)
    return false;

  std::size_t total_rhs_entries = 0;
  std::size_t total_mat_entries = 0;
  std::vector<HilbertElementAssemblyRecord> records;
  build_hilbert_element_records(sys, records, total_rhs_entries, total_mat_entries);
  if (!build_hilbert_record_quadrature_orders(sys, records))
    return false;
  sort_hilbert_element_records(geometry_cache, sys.variable_type(0), records);

  std::vector<PetscInt> rhs_rows(total_rhs_entries);
  std::vector<PetscInt> mat_rows(total_mat_entries);
  std::vector<PetscInt> mat_cols(total_mat_entries);
  build_hilbert_coo_indices(*dof_index_cache, records, rhs_rows, mat_rows, mat_cols);

  KokkosHilbertBatchData batch_data;
  build_hilbert_batch_data(records, batch_data);
  std::vector<KokkosHilbertAssemblyBucket> buckets;
  build_hilbert_assembly_buckets(geometry_cache, sys.variable_type(0), records, buckets);

  plan.records = std::move(records);
  plan.buckets = std::move(buckets);
  plan.rhs_rows = std::move(rhs_rows);
  plan.mat_rows = std::move(mat_rows);
  plan.mat_cols = std::move(mat_cols);
  plan.batch_data = std::move(batch_data);
  plan.rhs_values = ::Kokkos::View<Number *>("hilbert_rhs_values", plan.rhs_rows.size());
  plan.mat_values = ::Kokkos::View<Number *>("hilbert_mat_values", plan.mat_rows.size());
  plan.geometry_cache_id = &geometry_cache;
  plan.dof_index_cache_id = dof_index_cache;
  plan.fe_type = sys.variable_type(0);
  plan.hilbert_order = sys.hilbert_order();
  plan.extra_quadrature_order = sys.extra_quadrature_order;
  plan.subdomains = sys.subdomains_list();
  plan.matrix_target = nullptr;
  plan.rhs_target = nullptr;
  plan.fem_goal_target = nullptr;
  plan.input_vector_target = nullptr;
  return true;
}

bool
build_kokkos_fem_goal_batch_data(HilbertSystem & sys,
                                 const KokkosPetscAssemblyPlan & plan,
                                 const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                 KokkosFEMGoalBatchData & batch_data)
{
  libmesh_assert(sys.input_system);
  System & input_system = *sys.input_system;
  const unsigned int n_fields = goal_function.n_field_variables();
  if (n_fields > kokkos_parsed_fem_max_fields)
    return false;

  auto * input_vector = dynamic_cast<PetscVector<Number> *>(input_system.current_local_solution.get());
  if (!input_vector || !input_vector->supports_kokkos_access())
    return false;

  batch_data.input_vector = input_vector;
  batch_data.field_local_indices.size = n_fields;
  batch_data.bucket_field_keys.assign(plan.buckets.size(), KokkosFieldKeyStorage{});
  batch_data.bucket_field_dofs.assign(plan.buckets.size(), KokkosFieldDofStorage{});

  for (unsigned int field = 0; field != n_fields; ++field)
    {
      const unsigned int var_num = goal_function.field_variable_number(field);
      const FEType field_type = input_system.variable_type(var_num);
      const auto * local_index_cache =
        input_system.get_dof_map().require_kokkos_local_index_cache(
          *input_system.current_local_solution, var_num);
      if (!local_index_cache)
        return false;

      batch_data.field_local_indices.values[field] = local_index_cache->element_local_indices;

      for (auto bucket_index : index_range(plan.buckets))
        {
          const auto & bucket = plan.buckets[bucket_index];
          const FEShapeKey field_key =
            make_hilbert_shape_key(bucket.elem_type, bucket.elem_p_level, field_type);

          if (bucket.mapping_type != LAGRANGE_MAP ||
              !supports_shape_with_lagrange_map(field_key) ||
              (sys.hilbert_order() > 0 && !supports_grad_shape(field_key)))
            return false;

          const unsigned int field_n_dofs =
            FEInterface::n_dofs(libMesh::Kokkos::dim_from_topology(bucket.elem_type),
                                field_type,
                                bucket.elem_type);
          if (field_n_dofs > kokkos_hilbert_max_dofs)
            return false;

          batch_data.bucket_field_keys[bucket_index].size = n_fields;
          batch_data.bucket_field_dofs[bucket_index].size = n_fields;
          batch_data.bucket_field_keys[bucket_index].values[field] = field_key;
          batch_data.bucket_field_dofs[bucket_index].values[field] = field_n_dofs;
        }
    }

  return true;
}

bool
ensure_kokkos_fem_goal_batch_data(HilbertSystem & sys,
                                  const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                  KokkosPetscAssemblyPlan & plan)
{
  libmesh_assert(sys.input_system);
  const void * input_vector_target = sys.input_system->current_local_solution.get();
  if (plan.fem_goal_target == &goal_function && plan.input_vector_target == input_vector_target)
    return true;

  KokkosFEMGoalBatchData batch_data;
  if (!build_kokkos_fem_goal_batch_data(sys, plan, goal_function, batch_data))
    return false;

  plan.fem_goal_batch_data = std::move(batch_data);
  plan.fem_goal_target = &goal_function;
  plan.input_vector_target = input_vector_target;
  return true;
}

bool
assemble_kokkos_hilbert_element(const FEShapeKey key,
                                const unsigned int quadrature_order,
                                const unsigned int hilbert_order,
                                const Number solution_derivative,
                                const libMesh::Kokkos::KokkosParsedFunction<Number> & goal_function,
                                const KokkosElementAssemblyState & state,
                                const bool request_jacobian,
                                DenseSubVector<Number> & F,
                                DenseSubMatrix<Number> & K)
{
  const unsigned int n_dofs = F.size();
  KokkosScalarView<Number> d_F("hilbert_residual", n_dofs);
  KokkosDenseJacobianView d_K("hilbert_jacobian", n_dofs, n_dofs);
  const auto goal_access =
    libMesh::Kokkos::detail::make_hilbert_analytic_goal_access(goal_function,
                                                               goal_function.gradient_function());
  const libMesh::Kokkos::detail::DenseElementOutputSink<decltype(d_F), decltype(d_K)> sink{
    d_F, d_K, n_dofs, request_jacobian};

  auto coeff_guard = state.solution_vector->make_kokkos_read_view_guard();
  const auto coeff =
    libMesh::Kokkos::detail::make_gathered_coeff_access(
      coeff_guard.view(), state.solution_local_indices->element_local_indices, state.elem_index);

  const bool success =
    libMesh::Kokkos::detail::run_hilbert_system_assembly<kokkos_hilbert_max_dofs>(
      key,
      state.mapping_type,
      state.geometry_cache->node_coordinates,
      state.geometry_cache->element_node_ids,
      state.elem_index,
      state.elem_n_nodes,
      quadrature_order,
      hilbert_order,
      coeff,
      solution_derivative,
      goal_access,
      request_jacobian,
      sink,
      "hilbert_local_assembly");

  if (!success)
    return false;

  accumulate_hilbert_dense_outputs(d_F, d_K, request_jacobian, F, K);
  return true;
}

bool
assemble_kokkos_hilbert_fem_goal_element(const FEShapeKey output_key,
                                         const unsigned int quadrature_order,
                                         const unsigned int hilbert_order,
                                         const Number solution_derivative,
                                         System & input_system,
                                         const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                         const KokkosElementAssemblyState & state,
                                         const bool request_jacobian,
                                         DenseSubVector<Number> & F,
                                         DenseSubMatrix<Number> & K)
{
  KokkosFEMGoalData goal_data;
  if (!build_kokkos_fem_goal_data(state.elem_type,
                                  state.mapping_type,
                                  state.elem_p_level,
                                  hilbert_order,
                                  input_system,
                                  goal_function,
                                  goal_data))
    return false;

  const unsigned int n_dofs = F.size();
  KokkosScalarView<Number> d_F("hilbert_residual", n_dofs);
  KokkosDenseJacobianView d_K("hilbert_jacobian", n_dofs, n_dofs);
  const libMesh::Kokkos::detail::DenseElementOutputSink<decltype(d_F), decltype(d_K)> sink{
    d_F, d_K, n_dofs, request_jacobian};

  const auto assemble_with_input_coeffs = [&](const auto & coeff_values,
                                              const auto & input_coeff_values)
  {
    const auto coeff =
      libMesh::Kokkos::detail::make_gathered_coeff_access(coeff_values,
                                                          state.solution_local_indices
                                                            ->element_local_indices,
                                                          state.elem_index);

    const auto goal_access =
      libMesh::Kokkos::detail::GatheredParsedFEMGoalAccess<KokkosFieldKeyStorage,
                                                           KokkosFieldDofStorage,
                                                           std::decay_t<decltype(input_coeff_values)>,
                                                           KokkosLocalIndexView,
                                                           libMesh::Kokkos::KokkosParsedFEMFunction<Number>,
                                                           kokkos_parsed_fem_max_fields>(
        goal_data.field_keys,
        goal_data.field_dofs,
        input_coeff_values,
        goal_data.field_local_indices.data(),
        goal_function);

    return libMesh::Kokkos::detail::run_hilbert_system_assembly<kokkos_hilbert_max_dofs>(
      output_key,
      state.mapping_type,
      state.geometry_cache->node_coordinates,
      state.geometry_cache->element_node_ids,
      state.elem_index,
      state.elem_n_nodes,
      quadrature_order,
      hilbert_order,
      coeff,
      solution_derivative,
      goal_access,
      request_jacobian,
      sink,
      "hilbert_local_fem_goal_assembly");
  };

  auto coeff_guard = state.solution_vector->make_kokkos_read_view_guard();
  const bool success =
    (goal_data.input_vector == state.solution_vector)
      ? assemble_with_input_coeffs(coeff_guard.view(), coeff_guard.view())
      : [&]()
        {
          auto input_guard = goal_data.input_vector->make_kokkos_read_view_guard();
          return assemble_with_input_coeffs(coeff_guard.view(), input_guard.view());
        }();

  if (!success)
    return false;

  accumulate_hilbert_dense_outputs(d_F, d_K, request_jacobian, F, K);
  return true;
}

#if defined(LIBMESH_HAVE_PETSC)
bool
assemble_host_exact_parsed_fem_goal_element(HilbertSystem & sys,
                                            FEMContext & c,
                                            const bool request_jacobian,
                                            DenseSubVector<Number> & F,
                                            DenseSubMatrix<Number> & K,
                                            const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function)
{
  if (!sys.input_system)
    return false;

  FEMContext * goal_context_ptr = sys.get_input_context(c);
  if (!goal_context_ptr)
    return false;

  FEMContext & goal_context = *goal_context_ptr;
  goal_context.pre_fe_reinit(*sys.input_system, &c.get_elem());
  goal_context.elem_fe_reinit();

  detail::HostHilbertFEAccess fe(c, 0, sys.hilbert_order());
  detail::HostHilbertAccumulator accum(F, K);
  auto solution =
    detail::make_hilbert_solution_access(fe,
                                         c.get_elem_solution(0),
                                         c.get_elem_solution_derivative());
  HostExactParsedFEMGoalAccess goal_access(goal_function, goal_context);
  detail::assemble_hilbert_element(fe,
                                   solution,
                                   goal_access,
                                   request_jacobian,
                                   sys.hilbert_order(),
                                   accum);
  return true;
}
#endif

template <typename ResidualView, typename JacobianView>
bool
assemble_kokkos_hilbert_element_device_values(const FEShapeKey key,
                                              const unsigned int quadrature_order,
                                              const unsigned int hilbert_order,
                                              const unsigned int elem_index,
                                              const MeshBase::KokkosGeometryCache & geometry_cache,
                                              const libMesh::Kokkos::KokkosParsedFunction<Number> & goal_function,
                                              ResidualView d_rhs_values,
                                              JacobianView d_mat_values)
{
  const unsigned int n_dofs = cast_int<unsigned int>(d_rhs_values.extent(0));
  const auto goal_access =
    libMesh::Kokkos::detail::make_hilbert_analytic_goal_access(goal_function,
                                                               goal_function.gradient_function());
  const libMesh::Kokkos::detail::FlatDeviceValueSink<ResidualView, JacobianView> sink{
    d_rhs_values, d_mat_values, n_dofs};

  return libMesh::Kokkos::detail::run_hilbert_system_assembly<kokkos_hilbert_max_dofs>(
    key,
    geometry_cache.element_mapping_types(elem_index),
    geometry_cache.node_coordinates,
    geometry_cache.element_node_ids,
    elem_index,
    geometry_cache.element_n_nodes(elem_index),
    quadrature_order,
    hilbert_order,
    libMesh::Kokkos::detail::ZeroCoeffAccess{},
    Number(1.),
    goal_access,
    true,
    sink,
    "hilbert_device_values");
}

template <typename ResidualView, typename JacobianView>
bool
assemble_kokkos_hilbert_fem_goal_device_values(const FEShapeKey output_key,
                                               const unsigned int quadrature_order,
                                               const unsigned int hilbert_order,
                                               const unsigned int elem_index,
                                               const MeshBase::KokkosGeometryCache & geometry_cache,
                                               System & input_system,
                                               const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                               ResidualView d_rhs_values,
                                               JacobianView d_mat_values)
{
  KokkosFEMGoalData goal_data;
  if (!build_kokkos_fem_goal_data(geometry_cache.element_types(elem_index),
                                  geometry_cache.element_mapping_types(elem_index),
                                  geometry_cache.element_p_levels(elem_index),
                                  hilbert_order,
                                  input_system,
                                  goal_function,
                                  goal_data))
    return false;

  const unsigned int n_dofs = cast_int<unsigned int>(d_rhs_values.extent(0));
  const libMesh::Kokkos::detail::FlatDeviceValueSink<ResidualView, JacobianView> sink{
    d_rhs_values, d_mat_values, n_dofs};

  auto input_guard = goal_data.input_vector->make_kokkos_read_view_guard();
  const auto goal_access =
    libMesh::Kokkos::detail::GatheredParsedFEMGoalAccess<KokkosFieldKeyStorage,
                                                         KokkosFieldDofStorage,
                                                         decltype(input_guard.view()),
                                                         KokkosLocalIndexView,
                                                         libMesh::Kokkos::KokkosParsedFEMFunction<Number>,
                                                         kokkos_parsed_fem_max_fields>(
      goal_data.field_keys,
      goal_data.field_dofs,
      input_guard.view(),
      goal_data.field_local_indices.data(),
      goal_function);

  return libMesh::Kokkos::detail::run_hilbert_system_assembly<kokkos_hilbert_max_dofs>(
    output_key,
    geometry_cache.element_mapping_types(elem_index),
    geometry_cache.node_coordinates,
    geometry_cache.element_node_ids,
    elem_index,
    geometry_cache.element_n_nodes(elem_index),
    quadrature_order,
    hilbert_order,
    libMesh::Kokkos::detail::ZeroCoeffAccess{},
    Number(1.),
    goal_access,
    true,
    sink,
    "hilbert_device_fem_goal_values");
}

template <typename ResidualView, typename JacobianView>
bool
assemble_kokkos_hilbert_record_values(HilbertSystem & sys,
                                      const HilbertElementAssemblyRecord & record,
                                      const unsigned int quadrature_order,
                                      const libMesh::Kokkos::KokkosParsedFunction<Number> * analytic_goal,
                                      const libMesh::Kokkos::KokkosParsedFEMFunction<Number> * fem_goal,
                                      ResidualView rhs_slice,
                                      JacobianView mat_slice)
{
  const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
  const FEShapeKey shape_key =
    make_hilbert_shape_key(geometry_cache.element_types(record.elem_index),
                           geometry_cache.element_p_levels(record.elem_index),
                           sys.variable_type(0));

  return analytic_goal ?
           assemble_kokkos_hilbert_element_device_values(shape_key,
                                                         quadrature_order,
                                                         sys.hilbert_order(),
                                                         record.elem_index,
                                                         geometry_cache,
                                                         analytic_goal->with_time(sys.time),
                                                         rhs_slice,
                                                         mat_slice) :
           assemble_kokkos_hilbert_fem_goal_device_values(shape_key,
                                                          quadrature_order,
                                                          sys.hilbert_order(),
                                                          record.elem_index,
                                                          geometry_cache,
                                                          *sys.input_system,
                                                          fem_goal->with_time(sys.time),
                                                          rhs_slice,
                                                          mat_slice);
}

bool
assemble_kokkos_petsc_global_system(HilbertSystem & sys,
                                    KokkosPetscAssemblyPlan & plan,
                                    const libMesh::Kokkos::KokkosParsedFunction<Number> * analytic_goal,
                                    const libMesh::Kokkos::KokkosParsedFEMFunction<Number> * fem_goal,
                                    PetscMatrixBase<Number> & system_matrix,
                                    PetscVector<Number> & system_rhs)
{
  if (sys.has_static_condensation() || sys.get_dof_map().n_constrained_dofs())
    return false;

  if (!analytic_goal && !fem_goal)
    return false;

  if (analytic_goal)
    {
      const auto timed_goal = analytic_goal->with_time(sys.time);
      const auto goal_access =
        libMesh::Kokkos::detail::make_hilbert_analytic_goal_access(timed_goal,
                                                                   timed_goal.gradient_function());

      const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
      for (const auto & bucket : plan.buckets)
        libMesh::Kokkos::detail::run_hilbert_system_bucket_value_batch<kokkos_hilbert_max_dofs>(
          bucket.key,
          bucket.mapping_type,
          bucket.n_nodes,
          bucket.quadrature_order,
          geometry_cache.node_coordinates,
          geometry_cache.element_node_ids,
          ::Kokkos::subview(plan.batch_data.elem_indices,
                            ::Kokkos::make_pair(bucket.begin, bucket.end)),
          ::Kokkos::subview(plan.batch_data.elem_n_dofs,
                            ::Kokkos::make_pair(bucket.begin, bucket.end)),
          ::Kokkos::subview(plan.batch_data.rhs_offsets,
                            ::Kokkos::make_pair(bucket.begin, bucket.end)),
          ::Kokkos::subview(plan.batch_data.mat_offsets,
                            ::Kokkos::make_pair(bucket.begin, bucket.end)),
          sys.hilbert_order(),
          goal_access,
          plan.rhs_values,
          plan.mat_values,
          "hilbert_value_bucket_batch");
    }
  else
    {
      const auto timed_goal = fem_goal->with_time(sys.time);
      if (!ensure_kokkos_fem_goal_batch_data(sys, *fem_goal, plan))
        return false;

      const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
      auto input_guard = plan.fem_goal_batch_data.input_vector->make_kokkos_read_view_guard();
      for (auto bucket_index : index_range(plan.buckets))
        {
          const auto & bucket = plan.buckets[bucket_index];
          libMesh::Kokkos::detail::run_hilbert_system_fem_bucket_value_batch<
            kokkos_hilbert_max_dofs,
            kokkos_parsed_fem_max_fields>(
            bucket.key,
            bucket.mapping_type,
            bucket.n_nodes,
            bucket.quadrature_order,
            geometry_cache.node_coordinates,
            geometry_cache.element_node_ids,
            ::Kokkos::subview(plan.batch_data.elem_indices,
                              ::Kokkos::make_pair(bucket.begin, bucket.end)),
            ::Kokkos::subview(plan.batch_data.elem_n_dofs,
                              ::Kokkos::make_pair(bucket.begin, bucket.end)),
            ::Kokkos::subview(plan.batch_data.rhs_offsets,
                              ::Kokkos::make_pair(bucket.begin, bucket.end)),
            ::Kokkos::subview(plan.batch_data.mat_offsets,
                              ::Kokkos::make_pair(bucket.begin, bucket.end)),
            plan.fem_goal_batch_data.bucket_field_keys[bucket_index],
            plan.fem_goal_batch_data.bucket_field_dofs[bucket_index],
            plan.fem_goal_batch_data.field_local_indices,
            input_guard.view(),
            timed_goal,
            sys.hilbert_order(),
            plan.rhs_values,
            plan.mat_values,
            "hilbert_fem_value_bucket_batch");
        }
    }

  if (plan.matrix_target != &system_matrix)
    {
      LibmeshPetscCall2(sys.comm(),
                        MatSetPreallocationCOO(system_matrix.mat(),
                                               static_cast<PetscCount>(plan.mat_rows.size()),
                                               plan.mat_rows.empty() ? nullptr : plan.mat_rows.data(),
                                               plan.mat_cols.empty() ? nullptr : plan.mat_cols.data()));
      plan.matrix_target = &system_matrix;
    }
  if (plan.rhs_target != &system_rhs)
    {
      LibmeshPetscCall2(sys.comm(),
                        VecSetPreallocationCOO(system_rhs.vec(),
                                               static_cast<PetscCount>(plan.rhs_rows.size()),
                                               plan.rhs_rows.empty() ? nullptr : plan.rhs_rows.data()));
      plan.rhs_target = &system_rhs;
    }
  LibmeshPetscCall2(sys.comm(),
                    MatSetValuesCOO(system_matrix.mat(),
                                    reinterpret_cast<const PetscScalar *>(plan.mat_values.data()),
                                    INSERT_VALUES));
  LibmeshPetscCall2(sys.comm(),
                    VecSetValuesCOO(system_rhs.vec(),
                                    reinterpret_cast<const PetscScalar *>(plan.rhs_values.data()),
                                    INSERT_VALUES));

  return true;
}
#endif

} // anonymous namespace
#endif

HilbertSystem::~HilbertSystem () = default;

HilbertSystem::HilbertSystem(libMesh::EquationSystems & es,
                             const std::string & name,
                             const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    input_system(nullptr),
    _fe_family("LAGRANGE"),
    _fe_order(1),
    _hilbert_order(0),
    _use_kokkos_backend(false),
    _fdm_eps(libMesh::TOLERANCE),
    _subdomains_list()
{
}

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
const libMesh::Kokkos::KokkosParsedFunction<libMesh::Number> *
HilbertSystem::ensure_kokkos_goal_func()
{
  return ensure_kokkos_goal_cache<libMesh::ParsedFunction<libMesh::Number>>(
    _kokkos_goal_func,
    _analytic_goal_func,
    [](const auto & parsed_goal) -> std::unique_ptr<libMesh::Kokkos::KokkosParsedFunction<libMesh::Number>>
    {
#ifdef LIBMESH_HAVE_FPARSER
      return std::make_unique<libMesh::Kokkos::KokkosParsedFunction<libMesh::Number>>(
        parsed_goal.build_program_bundle());
#else
      libmesh_ignore(parsed_goal);
      return nullptr;
#endif
    });
}

const libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number> *
HilbertSystem::ensure_kokkos_fem_goal_func()
{
  return ensure_kokkos_goal_cache<libMesh::ParsedFEMFunction<libMesh::Number>>(
    _kokkos_fem_goal_func,
    _goal_func,
    [](const auto & parsed_goal)
      -> std::unique_ptr<libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number>>
    {
#ifdef LIBMESH_HAVE_FPARSER
      const auto program_bundle = parsed_goal.build_program_bundle();
      if (!program_bundle.supports_kokkos_value_goal() ||
          program_bundle.value_variable_numbers.size() > kokkos_parsed_fem_max_fields)
        return nullptr;

      return std::make_unique<libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number>>(
        program_bundle);
#else
      libmesh_ignore(parsed_goal);
      return nullptr;
#endif
    });
}

void
HilbertSystem::reset_kokkos_goal_cache()
{
  _kokkos_goal_func.reset();
  _kokkos_fem_goal_func.reset();
#if defined(LIBMESH_HAVE_PETSC)
  _kokkos_petsc_plan.reset();
#endif
}

KokkosPetscAssemblyPlan *
HilbertSystem::ensure_kokkos_petsc_plan(bool * rebuilt)
{
  const auto & geometry_cache = this->get_mesh().get_kokkos_geometry_cache();
  const auto * dof_index_cache = this->get_dof_map().get_kokkos_dof_index_cache(0);
  if (!dof_index_cache)
    return nullptr;

  if (_kokkos_petsc_plan &&
      kokkos_petsc_plan_matches(*this, geometry_cache, *dof_index_cache, *_kokkos_petsc_plan))
    {
      if (rebuilt)
        *rebuilt = false;
      return _kokkos_petsc_plan.get();
    }

  auto plan = std::make_unique<KokkosPetscAssemblyPlan>();
  if (!build_kokkos_petsc_assembly_plan(*this, *plan))
    return nullptr;

  _kokkos_petsc_plan = std::move(plan);
  if (rebuilt)
    *rebuilt = true;
  return _kokkos_petsc_plan.get();
}

bool
HilbertSystem::try_kokkos_element_assembly(FEMContext & c,
                                           const bool request_jacobian,
                                           DenseSubVector<Number> & F,
                                           DenseSubMatrix<Number> & K)
{
#if defined(LIBMESH_HAVE_PETSC)
  const Elem & elem = c.get_elem();
  const unsigned int quadrature_order =
    cast_int<unsigned int>(c.get_element_qrule().get_order());
  KokkosElementAssemblyState state;
  if (!build_kokkos_element_state(*this, elem, state))
    return false;
  const FEShapeKey shape_key =
    make_hilbert_shape_key(state.elem_type, state.elem_p_level, this->variable_type(0));

  if (const auto * kokkos_goal = this->ensure_kokkos_goal_func();
      kokkos_goal &&
      assemble_kokkos_hilbert_element(shape_key,
                                      quadrature_order,
                                      _hilbert_order,
                                      c.get_elem_solution_derivative(),
                                      kokkos_goal->with_time(this->time),
                                      state,
                                      request_jacobian,
                                      F,
                                      K))
    return true;

  if (!input_system)
    return false;

  if (const auto * kokkos_goal = this->ensure_kokkos_fem_goal_func();
      kokkos_goal &&
      assemble_kokkos_hilbert_fem_goal_element(shape_key,
                                               quadrature_order,
                                               _hilbert_order,
                                               c.get_elem_solution_derivative(),
                                               *input_system,
                                               kokkos_goal->with_time(this->time),
                                               state,
                                               request_jacobian,
                                               F,
                                               K))
    return true;

  return false;
#else
  libmesh_ignore(c, request_jacobian, F, K);
  return false;
#endif
}

#if defined(LIBMESH_HAVE_PETSC)
bool
HilbertSystem::try_kokkos_petsc_solve()
{
  using clock = std::chrono::steady_clock;
  auto * petsc_matrix = dynamic_cast<PetscMatrixBase<Number> *>(this->matrix);
  auto * petsc_rhs = dynamic_cast<PetscVector<Number> *>(this->rhs);
  auto * petsc_solution = dynamic_cast<PetscVector<Number> *>(this->solution.get());

  const auto * analytic_goal = this->ensure_kokkos_goal_func();
  const auto * fem_goal = this->ensure_kokkos_fem_goal_func();

  if (!petsc_matrix || !petsc_rhs || !petsc_solution || !(analytic_goal || fem_goal))
    return false;

  prewarm_kokkos_hilbert_entities(*this, fem_goal);
  this->_last_kokkos_timing = {};
  const auto total_start = clock::now();
  petsc_matrix->zero();
  petsc_rhs->zero();
  petsc_solution->zero();

  bool rebuilt_plan = false;
  const auto plan_start = clock::now();
  auto * plan = this->ensure_kokkos_petsc_plan(&rebuilt_plan);
  if (!plan)
    return false;
  const auto plan_stop = clock::now();
  this->_last_kokkos_timing.plan_seconds =
    rebuilt_plan ?
      std::chrono::duration_cast<std::chrono::duration<Real>>(plan_stop - plan_start).count() :
      0.;

  const auto assembly_start = clock::now();
  if (!assemble_kokkos_petsc_global_system(*this,
                                           *plan,
                                           analytic_goal,
                                           fem_goal,
                                           *petsc_matrix,
                                           *petsc_rhs))
    return false;
  const auto assembly_stop = clock::now();
  this->_last_kokkos_timing.assembly_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(assembly_stop - assembly_start).count();

  petsc_matrix->close();
  petsc_rhs->close();
  petsc_solution->close();

  LinearSolver<Number> * solver = this->get_linear_solver();
  if (this->prefix_with_name())
    solver->init(this->prefix().c_str());
  else
    solver->init();

  const auto [maxlinearits, linear_tol] = this->get_linear_solve_parameters();
  const auto solve_start = clock::now();
  solver->solve(*this->matrix,
                *this->solution,
                *this->rhs,
                linear_tol,
                maxlinearits);
  const auto solve_stop = clock::now();
  this->_last_kokkos_timing.solve_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(solve_stop - solve_start).count();

  this->update();
  this->mesh_position_set();
  const auto total_stop = clock::now();
  this->_last_kokkos_timing.total_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(total_stop - total_start).count();
  return true;
}
#endif
#endif

void HilbertSystem::init_data ()
{
  this->add_variable ("u", static_cast<Order>(_fe_order),
                      Utility::string_to_enum<FEFamily>(_fe_family));

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();
}



void HilbertSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * my_fe = nullptr;

  // Now make sure we have requested all the data
  // we need to build the L2 system.

  // We might have a multi-dimensional mesh
  const std::set<unsigned char> & elem_dims =
    c.elem_dimensions();

  for (const auto & dim : elem_dims)
    {
      c.get_element_fe( 0, my_fe, dim );

      my_fe->get_JxW();
      my_fe->get_phi();
      my_fe->get_xyz();

      if (this->_hilbert_order > 0)
        my_fe->get_dphi();

      c.get_side_fe( 0, my_fe, dim );
      my_fe->get_nothing();
    }

  // Build a corresponding context for the input system if we haven't
  // already
  auto & input_context = input_contexts[&c];
  if (input_system && !input_context)
    {
      input_context = std::make_unique<FEMContext>(*input_system);
    }

  libmesh_assert(_goal_func || _analytic_goal_func);

  if (_goal_func)
    _goal_func->init_context(input_system ? *input_context : c);

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (input_system &&
      this->_hilbert_order > 0 &&
      dynamic_cast<ParsedFEMFunction<Number> *>(_goal_func.get()) &&
      this->ensure_kokkos_fem_goal_func())
    {
      for (const auto & dim : elem_dims)
        for (unsigned int var = 0; var != input_system->n_vars(); ++var)
          {
            input_context->get_element_fe(var, my_fe, dim);
            my_fe->get_dphi();
          }
    }
#endif

  FEMSystem::init_context(context);
}


bool HilbertSystem::element_time_derivative (bool request_jacobian,
                                             DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const Elem & elem = c.get_elem();

  if (!_subdomains_list.empty() &&
      !_subdomains_list.count(elem.subdomain_id()))
    return request_jacobian;

  // The subvectors and submatrices we need to fill:
  DenseSubMatrix<Number> & K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> & F = c.get_elem_residual(0);

#ifdef LIBMESH_HAVE_KOKKOS
  if (_use_kokkos_backend)
    {
#if !defined(LIBMESH_USE_COMPLEX_NUMBERS)
      if (this->try_kokkos_element_assembly(c, request_jacobian, F, K))
        return request_jacobian;
#else
      if (_analytic_goal_func &&
          dynamic_cast<ParsedFunction<Number> *>(_analytic_goal_func.get()))
        libmesh_error_msg("HilbertSystem Kokkos backend does not support ParsedFunction goals "
                          "when libMesh is built with complex Number.");
#endif
    }
#endif

  detail::HostHilbertFEAccess fe(c, 0, _hilbert_order);
  const auto assemble_with_goal = [&](auto & goal)
  {
    auto solution =
      detail::make_hilbert_solution_access(fe,
                                           c.get_elem_solution(0),
                                           c.get_elem_solution_derivative());
    detail::HostHilbertAccumulator accum(F, K);
    detail::assemble_hilbert_element(fe,
                                     solution,
                                     goal,
                                     request_jacobian,
                                     _hilbert_order,
                                     accum);
  };

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (const auto * kokkos_goal = this->ensure_kokkos_goal_func())
    {
      const auto parsed_goal = kokkos_goal->with_time(this->time);
      auto goal = detail::make_hilbert_analytic_goal_access(parsed_goal,
                                                            parsed_goal.gradient_function());
      assemble_with_goal(goal);
      return request_jacobian;
    }

#if defined(LIBMESH_HAVE_PETSC)
  if (input_system)
    if (const auto * kokkos_fem_goal = this->ensure_kokkos_fem_goal_func();
        kokkos_fem_goal &&
        assemble_host_exact_parsed_fem_goal_element(*this,
                                                    c,
                                                    request_jacobian,
                                                    F,
                                                    K,
                                                    kokkos_fem_goal->with_time(this->time)))
      return request_jacobian;
#endif
#endif

  if (_analytic_goal_func)
    {
      auto goal = detail::make_hilbert_analytic_goal_access(*_analytic_goal_func,
                                                            *_analytic_goal_grad);
      assemble_with_goal(goal);
    }
  else
    {
      FEMContext & goal_context =
        input_system ? *libmesh_map_find(input_contexts, &c) : c;

      if (input_system)
        {
          goal_context.pre_fe_reinit(*input_system, &elem);
          goal_context.elem_fe_reinit();
        }

      detail::HostHilbertGoalAccess goal(*_goal_func, _goal_grad.get(), goal_context);
      assemble_with_goal(goal);
    }

  return request_jacobian;
}

void HilbertSystem::solve()
{
  _last_kokkos_timing = {};
#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (_use_kokkos_backend && this->try_kokkos_petsc_solve())
    return;
#endif

  FEMSystem::solve();
}
