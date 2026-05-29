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

#include "../../include/systems/hilbert_assembly.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/linear_solver.h"
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/parsed_fem_function.h"
#include "libmesh/parsed_function.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef LIBMESH_HAVE_PETSC
#include "libmesh/petsc_matrix_base.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/petsc_vector.h"
#include <petscsf.h>
#endif

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
#include "../../include/gpu/kokkos_hilbert_system.h"
#include "libmesh/dof_map_kokkos.h"
#include "libmesh/fe_shape_traits.h"
#include "libmesh/mesh_base_kokkos.h"

#define PETSC_SKIP_CXX_COMPLEX_FIX 1
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#ifdef LIBMESH_HAVE_NVTX_API
#include <nvtx3/nvtx3.hpp>
#endif
#undef __CUDACC_VER__
#endif

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
namespace libMesh
{
constexpr unsigned int kokkos_hilbert_max_dofs = 27;
constexpr unsigned int kokkos_parsed_fem_max_fields = 16;

template <typename T>
using KokkosScalarView = ::Kokkos::View<T *>;

using KokkosDenseJacobianView = ::Kokkos::View<Number **>;
using KokkosGradientView = ::Kokkos::View<Real **>;
using KokkosUnsignedIntView = ::Kokkos::View<unsigned int *>;
using KokkosPetscIntView = ::Kokkos::View<PetscInt *>;
using KokkosSizeView = ::Kokkos::View<std::size_t *>;
using KokkosShapeKeyRecordView = ::Kokkos::View<FEShapeKey *>;
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

std::uint32_t
kokkos_region_color(const char * name)
{
  if (!name)
    return 0xff8f8f8f;
  if (std::strstr(name, "PETSc"))
    return 0xffe15759;
  if (std::strstr(name, "solve") || std::strstr(name, "KSPSolve"))
    return 0xffedc948;
  if (std::strstr(name, "direct storage"))
    return 0xfff28e2b;
  if (std::strstr(name, "remote"))
    return 0xff76b7b2;
  if (std::strstr(name, "analytic"))
    return 0xffb07aa1;
  if (std::strstr(name, "assemble") || std::strstr(name, "scatter") ||
      std::strstr(name, "entry values"))
    return 0xff4e79a7;
  if (std::strstr(name, "plan") || std::strstr(name, "batch") ||
      std::strstr(name, "graph") || std::strstr(name, "COO") ||
      std::strstr(name, "slot"))
    return 0xff59a14f;
  if (std::strstr(name, "zero") || std::strstr(name, "close") ||
      std::strstr(name, "cleanup"))
    return 0xffbab0ac;
  return 0xff9c755f;
}

struct KokkosProfilingRegion
{
  explicit KokkosProfilingRegion(const char * name)
  {
    ::Kokkos::Profiling::pushRegion(name);
#ifdef LIBMESH_HAVE_NVTX_API
    nvtxEventAttributes_t event_attributes = {};
    event_attributes.version = NVTX_VERSION;
    event_attributes.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    event_attributes.colorType = NVTX_COLOR_ARGB;
    event_attributes.color = kokkos_region_color(name);
    event_attributes.messageType = NVTX_MESSAGE_TYPE_ASCII;
    event_attributes.message.ascii = name;
    nvtxRangePushEx(&event_attributes);
#endif
  }

  ~KokkosProfilingRegion()
  {
#ifdef LIBMESH_HAVE_NVTX_API
    nvtxRangePop();
#endif
    ::Kokkos::Profiling::popRegion();
  }

  KokkosProfilingRegion(const KokkosProfilingRegion &) = delete;
  KokkosProfilingRegion & operator=(const KokkosProfilingRegion &) = delete;
};

void
accumulate_hilbert_dense_outputs(const KokkosScalarView<Number> & d_F,
                                 const KokkosDenseJacobianView & d_K,
                                 const bool request_jacobian,
                                 DenseSubVector<Number> & F,
                                 DenseSubMatrix<Number> & K)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::accumulate dense element outputs");

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
build_hilbert_element_records(::HilbertSystem & sys,
                              std::vector<HilbertElementAssemblyRecord> & records,
                              std::size_t & total_rhs_entries,
                              std::size_t & total_mat_entries)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build element records");

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
build_hilbert_record_quadrature_orders(::HilbertSystem & sys,
                                       std::vector<HilbertElementAssemblyRecord> & records)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build record quadrature orders");

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
        make_hilbert_shape_key(geometry_cache.host_element_types[record.elem_index],
                               geometry_cache.host_element_p_levels[record.elem_index],
                               fe_type);
      if (!libMesh::Kokkos::detail::supports_hilbert_local_assembly(
             shape_key,
             geometry_cache.host_element_mapping_types[record.elem_index],
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
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure parsed goal cache");

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

class SampledAnalyticGoalAccess
{
public:
  SampledAnalyticGoalAccess(KokkosSizeView qp_offsets,
                            KokkosScalarView<Number> values,
                            KokkosGradientView gradients)
    : _qp_offsets(qp_offsets),
      _values(values),
      _gradients(gradients)
  {
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Number value(const QpData & qp_data,
               const Point &) const
  {
    return _values(_qp_offsets(qp_data.elem_index()) + qp_data.qp_index());
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Gradient gradient(const QpData & qp_data,
                    const Point &) const
  {
    const auto offset = _qp_offsets(qp_data.elem_index()) + qp_data.qp_index();
    Gradient grad;
    for (unsigned int d = 0; d != LIBMESH_DIM; ++d)
      grad(d) = _gradients(offset, d);
    return grad;
  }

private:
  KokkosSizeView _qp_offsets;
  KokkosScalarView<Number> _values;
  KokkosGradientView _gradients;
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
  KokkosShapeKeyRecordView shape_keys;
};

struct KokkosFEMGoalBatchData
{
  KokkosFieldKeyRecordView field_keys;
  KokkosFieldDofRecordView field_dofs;
  libMesh::Kokkos::detail::StaticArrayAccess<KokkosLocalIndexView,
                                             kokkos_parsed_fem_max_fields>
    field_local_indices;
  PetscVector<Number> * input_vector = nullptr;
};

enum class KokkosDirectMatrixLayout
{
  none,
  seq_aij,
  mpi_aij
};

struct PetscIntPairHash
{
  std::size_t operator()(const std::pair<PetscInt, PetscInt> & values) const noexcept
  {
    const auto first_hash = std::hash<PetscInt>{}(values.first);
    const auto second_hash = std::hash<PetscInt>{}(values.second);
    return first_hash ^ (second_hash + 0x9e3779b97f4a7c15ULL + (first_hash << 6) + (first_hash >> 2));
  }
};

struct KokkosHilbertHomogeneousBucket
{
  FEShapeKey key{INVALID_FE, INVALID_ELEM, INVALID_ORDER};
  ElemMappingType mapping_type = INVALID_MAP;
  unsigned int n_nodes = 0;
  unsigned int quadrature_order = 0;
  unsigned int n_dofs = 0;
};

struct KokkosPetscAssemblyPlan
{
  std::vector<HilbertElementAssemblyRecord> records;
  std::size_t total_rhs_entries = 0;
  std::size_t total_mat_entries = 0;
  KokkosHilbertBatchData batch_data;
  std::vector<PetscInt> local_row_offsets;
  std::vector<PetscInt> local_column_indices;
  std::vector<PetscInt> diag_row_offsets;
  std::vector<PetscInt> diag_column_indices;
  std::vector<PetscInt> offdiag_row_offsets;
  std::vector<PetscInt> offdiag_column_indices;
  std::vector<PetscInt> offdiag_global_columns;
  std::vector<PetscInt> mat_coo_rows;
  std::vector<PetscInt> mat_coo_cols;
  KokkosPetscIntView local_row_offsets_view;
  KokkosPetscIntView local_column_indices_view;
  KokkosPetscIntView diag_row_offsets_view;
  KokkosPetscIntView diag_column_indices_view;
  KokkosPetscIntView offdiag_row_offsets_view;
  KokkosPetscIntView offdiag_column_indices_view;
  KokkosPetscIntView offdiag_global_columns_view;
  KokkosSizeView rhs_local_slots;
  KokkosSizeView mat_value_slots;
  KokkosScalarView<Number> rhs_entry_values;
  KokkosScalarView<Number> mat_entry_values;
  ::Kokkos::View<Number *> rhs_remote_values;
  ::Kokkos::View<Number *> mat_remote_values;
  ::Kokkos::View<Number *> rhs_remote_root_values;
  ::Kokkos::View<Number *> mat_remote_root_values;
  ::Kokkos::View<PetscScalar *> mat_coo_values;
  std::vector<PetscInt> rhs_remote_rows;
  std::vector<PetscInt> mat_remote_rows;
  std::vector<PetscInt> mat_remote_cols;
  std::vector<PetscInt> mat_remote_root_indices;
  KokkosPetscIntView rhs_remote_rows_view;
  KokkosPetscIntView mat_remote_rows_view;
  KokkosPetscIntView mat_remote_cols_view;
  KokkosPetscIntView mat_remote_root_indices_view;
  std::vector<processor_id_type> rhs_remote_owners;
  std::vector<processor_id_type> mat_remote_owners;
  PetscSF rhs_remote_sf = nullptr;
  PetscSF mat_remote_sf = nullptr;
  KokkosFEMGoalBatchData fem_goal_batch_data;
  KokkosSizeView analytic_goal_qp_offsets;
  KokkosScalarView<Number> analytic_goal_values;
  KokkosGradientView analytic_goal_gradients;
  const void * geometry_cache_id = nullptr;
  const void * dof_index_cache_id = nullptr;
  FEType fe_type;
  unsigned int hilbert_order = 0;
  int extra_quadrature_order = 0;
  std::set<subdomain_id_type> subdomains;
  const void * graph_matrix_target = nullptr;
  const void * direct_matrix_target = nullptr;
  const void * direct_rhs_target = nullptr;
  const void * fem_goal_target = nullptr;
  const void * input_vector_target = nullptr;
  const void * analytic_goal_target = nullptr;
  const void * analytic_goal_gradient_target = nullptr;
  Real analytic_goal_time = std::numeric_limits<Real>::quiet_NaN();
  KokkosHilbertHomogeneousBucket single_hilbert_bucket;
  bool has_single_hilbert_bucket = false;
  PetscInt row_start = 0;
  PetscInt row_stop = 0;
  PetscInt col_start = 0;
  PetscInt col_stop = 0;
  std::size_t rhs_local_size = 0;
  std::size_t mat_diag_size = 0;
  std::size_t mat_offdiag_size = 0;
  KokkosDirectMatrixLayout direct_matrix_layout = KokkosDirectMatrixLayout::none;
  KokkosDirectMatrixLayout direct_slot_layout = KokkosDirectMatrixLayout::none;
  numeric_index_type direct_slot_rhs_first_local = 0;
  numeric_index_type direct_slot_rhs_last_local = 0;
  PetscInt direct_slot_row_start = 0;
  PetscInt direct_slot_row_stop = 0;
  PetscInt direct_slot_col_start = 0;
  PetscInt direct_slot_col_stop = 0;
  bool direct_slot_cache_active = false;
  bool direct_storage_active = false;

  ~KokkosPetscAssemblyPlan()
  {
    if (rhs_remote_sf)
      {
        PetscErrorCode ierr = PetscSFDestroy(&rhs_remote_sf);
        libmesh_ignore(ierr);
      }
    if (mat_remote_sf)
      {
        PetscErrorCode ierr = PetscSFDestroy(&mat_remote_sf);
        libmesh_ignore(ierr);
      }
  }
};

bool
get_single_supported_hilbert_bucket(const KokkosPetscAssemblyPlan & plan,
                                    const MeshBase::KokkosGeometryCache & geometry_cache,
                                    KokkosHilbertHomogeneousBucket & bucket)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::detect homogeneous bucket");

  if (plan.records.empty())
    return false;

  const auto & first_record = plan.records.front();
  bucket.key = make_hilbert_shape_key(geometry_cache.host_element_types[first_record.elem_index],
                                      geometry_cache.host_element_p_levels[first_record.elem_index],
                                      plan.fe_type);
  bucket.mapping_type = geometry_cache.host_element_mapping_types[first_record.elem_index];
  bucket.n_nodes = geometry_cache.host_element_n_nodes[first_record.elem_index];
  bucket.quadrature_order = first_record.quadrature_order;
  bucket.n_dofs = first_record.n_dofs;

  if (bucket.n_dofs > kokkos_hilbert_max_dofs ||
      !libMesh::Kokkos::detail::supports_hilbert_local_assembly(bucket.key,
                                                               bucket.mapping_type,
                                                               bucket.quadrature_order))
    return false;

  for (const auto & record : plan.records)
    {
      const FEShapeKey key =
        make_hilbert_shape_key(geometry_cache.host_element_types[record.elem_index],
                               geometry_cache.host_element_p_levels[record.elem_index],
                               plan.fe_type);
      if (key.family != bucket.key.family ||
          key.elem_type != bucket.key.elem_type ||
          key.order != bucket.key.order ||
          geometry_cache.host_element_mapping_types[record.elem_index] != bucket.mapping_type ||
          geometry_cache.host_element_n_nodes[record.elem_index] != bucket.n_nodes ||
          record.quadrature_order != bucket.quadrature_order ||
          record.n_dofs != bucket.n_dofs)
        return false;
    }

  return true;
}

void
cache_single_supported_hilbert_bucket(KokkosPetscAssemblyPlan & plan,
                                      const MeshBase::KokkosGeometryCache & geometry_cache)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::cache homogeneous bucket");

  KokkosHilbertHomogeneousBucket bucket;
  plan.has_single_hilbert_bucket =
    get_single_supported_hilbert_bucket(plan, geometry_cache, bucket);
  plan.single_hilbert_bucket =
    plan.has_single_hilbert_bucket ? bucket : KokkosHilbertHomogeneousBucket{};
}

bool
get_cached_single_supported_hilbert_bucket(const KokkosPetscAssemblyPlan & plan,
                                           KokkosHilbertHomogeneousBucket & bucket)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::get homogeneous bucket cache");

  if (!plan.has_single_hilbert_bucket)
    return false;

  bucket = plan.single_hilbert_bucket;
  return true;
}

bool
debug_kokkos_bucket_value_batch()
{
  return libMesh::on_command_line("--kokkos-debug-bucket");
}

void
ensure_kokkos_entry_value_storage(KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure entry value storage");

  if (plan.rhs_entry_values.extent(0) != plan.total_rhs_entries)
    plan.rhs_entry_values = KokkosScalarView<Number>("hilbert_rhs_entry_values",
                                                     plan.total_rhs_entries);
  if (plan.mat_entry_values.extent(0) != plan.total_mat_entries)
    plan.mat_entry_values = KokkosScalarView<Number>("hilbert_mat_entry_values",
                                                     plan.total_mat_entries);
}

constexpr PetscMemType
kokkos_default_petsc_mem_type()
{
  return PETSC_MEMTYPE_KOKKOS;
}

bool
is_petsc_kokkos_vec(const Parallel::Communicator & comm, Vec vec)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::check PETSc Kokkos vector");

  const char * vec_type = nullptr;
  LibmeshPetscCall2(comm, VecGetType(vec, &vec_type));

  PetscBool is_seq_kokkos = PETSC_FALSE;
  PetscBool is_mpi_kokkos = PETSC_FALSE;
  LibmeshPetscCall2(comm, PetscStrcmp(vec_type, "seqkokkos", &is_seq_kokkos));
  LibmeshPetscCall2(comm, PetscStrcmp(vec_type, "mpikokkos", &is_mpi_kokkos));
  return is_seq_kokkos || is_mpi_kokkos;
}

void
clear_kokkos_petsc_remote_sf(KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::clear remote SF");

  if (plan.rhs_remote_sf)
    {
      const auto ierr = PetscSFDestroy(&plan.rhs_remote_sf);
      libmesh_ignore(ierr);
      plan.rhs_remote_sf = nullptr;
    }
  if (plan.mat_remote_sf)
    {
      const auto ierr = PetscSFDestroy(&plan.mat_remote_sf);
      libmesh_ignore(ierr);
      plan.mat_remote_sf = nullptr;
    }
  plan.rhs_remote_root_values = {};
  plan.mat_remote_root_values = {};
  plan.mat_remote_root_indices.clear();
  plan.rhs_remote_rows_view = {};
  plan.mat_remote_rows_view = {};
  plan.mat_remote_cols_view = {};
  plan.mat_remote_root_indices_view = {};
}

void
sync_kokkos_petsc_int_view(const std::vector<PetscInt> & host_values,
                           KokkosPetscIntView & device_view,
                           const std::string & name)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::sync PETSc int view");

  device_view = KokkosPetscIntView(name, host_values.size());
  auto host_view = ::Kokkos::create_mirror_view(device_view);
  for (auto i : index_range(host_values))
    host_view(i) = host_values[i];
  ::Kokkos::deep_copy(device_view, host_view);
}

void
sync_kokkos_petsc_owned_graph_views(KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::sync PETSc graph views");

  sync_kokkos_petsc_int_view(plan.local_row_offsets,
                             plan.local_row_offsets_view,
                             "hilbert_local_row_offsets");
  sync_kokkos_petsc_int_view(plan.local_column_indices,
                             plan.local_column_indices_view,
                             "hilbert_local_column_indices");
  sync_kokkos_petsc_int_view(plan.diag_row_offsets,
                             plan.diag_row_offsets_view,
                             "hilbert_diag_row_offsets");
  sync_kokkos_petsc_int_view(plan.diag_column_indices,
                             plan.diag_column_indices_view,
                             "hilbert_diag_column_indices");
  sync_kokkos_petsc_int_view(plan.offdiag_row_offsets,
                             plan.offdiag_row_offsets_view,
                             "hilbert_offdiag_row_offsets");
  sync_kokkos_petsc_int_view(plan.offdiag_column_indices,
                             plan.offdiag_column_indices_view,
                             "hilbert_offdiag_column_indices");
  sync_kokkos_petsc_int_view(plan.offdiag_global_columns,
                             plan.offdiag_global_columns_view,
                             "hilbert_offdiag_global_columns");
}

void
sync_kokkos_petsc_remote_slot_views(KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::sync remote slot views");

  sync_kokkos_petsc_int_view(plan.rhs_remote_rows,
                             plan.rhs_remote_rows_view,
                             "hilbert_rhs_remote_rows");
  sync_kokkos_petsc_int_view(plan.mat_remote_rows,
                             plan.mat_remote_rows_view,
                             "hilbert_mat_remote_rows");
  sync_kokkos_petsc_int_view(plan.mat_remote_cols,
                             plan.mat_remote_cols_view,
                             "hilbert_mat_remote_cols");
  sync_kokkos_petsc_int_view(plan.mat_remote_root_indices,
                             plan.mat_remote_root_indices_view,
                             "hilbert_mat_remote_root_indices");
}

bool
build_kokkos_petsc_owned_csr_graph(::HilbertSystem & sys,
                                   KokkosPetscAssemblyPlan & plan,
                                   const KokkosDirectMatrixLayout layout,
                                   const PetscInt row_start,
                                   const PetscInt row_stop,
                                   const PetscInt col_start,
                                   const PetscInt col_stop)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build owned CSR graph");

  const auto * sp = sys.get_dof_map().get_sparsity_pattern();
  if (!sp)
    return false;

  const auto & graph = sp->get_sparsity_pattern();
  plan.local_row_offsets.assign(graph.size() + 1, 0);
  plan.local_column_indices.clear();
  plan.diag_row_offsets.clear();
  plan.diag_column_indices.clear();
  plan.offdiag_row_offsets.clear();
  plan.offdiag_column_indices.clear();
  plan.offdiag_global_columns.clear();
  plan.mat_coo_rows.clear();
  plan.mat_coo_cols.clear();
  plan.row_start = row_start;
  plan.row_stop = row_stop;
  plan.col_start = col_start;
  plan.col_stop = col_stop;

  std::vector<std::vector<PetscInt>> row_columns(graph.size());
  std::vector<std::vector<PetscInt>> row_offdiag_globals;

  if (layout == KokkosDirectMatrixLayout::mpi_aij)
    row_offdiag_globals.resize(graph.size());

  std::vector<PetscInt> all_offdiag_globals;
  for (auto local_row : index_range(graph))
    {
      auto & cols = row_columns[local_row];
      cols.reserve(graph[local_row].size());
      for (const auto dof : graph[local_row])
        cols.push_back(cast_int<PetscInt>(dof));

      std::sort(cols.begin(), cols.end());
      cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
      plan.local_row_offsets[local_row + 1] =
        plan.local_row_offsets[local_row] + cast_int<PetscInt>(cols.size());
      plan.local_column_indices.insert(plan.local_column_indices.end(), cols.begin(), cols.end());

      if (layout != KokkosDirectMatrixLayout::mpi_aij)
        continue;

      auto & offdiag_globals = row_offdiag_globals[local_row];
      for (const auto col : cols)
        if (col < col_start || col >= col_stop)
          {
            offdiag_globals.push_back(col);
            all_offdiag_globals.push_back(col);
          }
    }

  if (layout != KokkosDirectMatrixLayout::mpi_aij)
    {
      sync_kokkos_petsc_owned_graph_views(plan);
      return true;
    }

  std::sort(all_offdiag_globals.begin(), all_offdiag_globals.end());
  all_offdiag_globals.erase(std::unique(all_offdiag_globals.begin(), all_offdiag_globals.end()),
                            all_offdiag_globals.end());
  plan.offdiag_global_columns = std::move(all_offdiag_globals);

  std::unordered_map<PetscInt, PetscInt> offdiag_column_map;
  offdiag_column_map.reserve(plan.offdiag_global_columns.size());
  for (auto i : index_range(plan.offdiag_global_columns))
    offdiag_column_map.emplace(plan.offdiag_global_columns[i], cast_int<PetscInt>(i));

  plan.diag_row_offsets.assign(graph.size() + 1, 0);
  plan.offdiag_row_offsets.assign(graph.size() + 1, 0);

  for (auto local_row : index_range(graph))
    {
      const auto & cols = row_columns[local_row];
      const auto & offdiag_globals = row_offdiag_globals[local_row];
      for (const auto col : cols)
        if (col >= col_start && col < col_stop)
          plan.diag_column_indices.push_back(col - col_start);

      for (const auto col : offdiag_globals)
        plan.offdiag_column_indices.push_back(offdiag_column_map.at(col));

      plan.diag_row_offsets[local_row + 1] =
        cast_int<PetscInt>(plan.diag_column_indices.size());
      plan.offdiag_row_offsets[local_row + 1] =
        cast_int<PetscInt>(plan.offdiag_column_indices.size());
    }

  sync_kokkos_petsc_owned_graph_views(plan);
  return true;
}

void
build_kokkos_petsc_owned_coo(KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build PETSc COO graph");

  plan.mat_coo_rows.clear();
  plan.mat_coo_cols.clear();
  plan.mat_coo_rows.reserve(plan.mat_diag_size + plan.mat_offdiag_size);
  plan.mat_coo_cols.reserve(plan.mat_diag_size + plan.mat_offdiag_size);

  if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij)
    {
      const PetscInt n_local_rows = plan.row_stop - plan.row_start;
      for (PetscInt local_row = 0; local_row != n_local_rows; ++local_row)
        {
          const PetscInt global_row = plan.row_start + local_row;
          for (PetscInt k = plan.local_row_offsets[local_row];
               k != plan.local_row_offsets[local_row + 1];
               ++k)
            {
              plan.mat_coo_rows.push_back(global_row);
              plan.mat_coo_cols.push_back(plan.local_column_indices[k]);
            }
        }
      return;
    }

  const PetscInt n_local_rows = plan.row_stop - plan.row_start;
  for (PetscInt local_row = 0; local_row != n_local_rows; ++local_row)
    {
      const PetscInt global_row = plan.row_start + local_row;
      for (PetscInt k = plan.diag_row_offsets[local_row];
           k != plan.diag_row_offsets[local_row + 1];
           ++k)
        {
          plan.mat_coo_rows.push_back(global_row);
          plan.mat_coo_cols.push_back(plan.col_start + plan.diag_column_indices[k]);
        }
    }

  for (PetscInt local_row = 0; local_row != n_local_rows; ++local_row)
    {
      const PetscInt global_row = plan.row_start + local_row;
      for (PetscInt k = plan.offdiag_row_offsets[local_row];
           k != plan.offdiag_row_offsets[local_row + 1];
           ++k)
        {
          plan.mat_coo_rows.push_back(global_row);
          plan.mat_coo_cols.push_back(plan.offdiag_global_columns[plan.offdiag_column_indices[k]]);
        }
    }
}

void
finalize_kokkos_direct_slots(::HilbertSystem & sys,
                             KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::finalize direct slots");

  plan.rhs_remote_owners.resize(plan.rhs_remote_rows.size());
  for (auto i : index_range(plan.rhs_remote_rows))
    plan.rhs_remote_owners[i] =
      sys.get_dof_map().dof_owner(cast_int<dof_id_type>(plan.rhs_remote_rows[i]));

  plan.mat_remote_owners.resize(plan.mat_remote_rows.size());
  for (auto i : index_range(plan.mat_remote_rows))
    plan.mat_remote_owners[i] =
      sys.get_dof_map().dof_owner(cast_int<dof_id_type>(plan.mat_remote_rows[i]));

  if (!plan.rhs_remote_rows.empty())
    plan.rhs_remote_values =
      ::Kokkos::View<Number *>("hilbert_rhs_remote_values", plan.rhs_remote_rows.size());

  if (!plan.mat_remote_rows.empty())
    plan.mat_remote_values =
      ::Kokkos::View<Number *>("hilbert_mat_remote_values", plan.mat_remote_rows.size());

  if (!plan.mat_coo_rows.empty())
    plan.mat_coo_values =
      ::Kokkos::View<PetscScalar *>("hilbert_mat_coo_values", plan.mat_coo_rows.size());
}

std::size_t
lookup_kokkos_owned_matrix_slot(const KokkosPetscAssemblyPlan & plan,
                                const PetscInt row,
                                const PetscInt col)
{
  const PetscInt local_row = row - plan.row_start;
  libmesh_error_msg_if(local_row < 0 || row >= plan.row_stop,
                       "HilbertSystem Kokkos remote matrix slot lookup received a nonlocal row.");

  if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij)
    {
      const PetscInt row_begin = plan.local_row_offsets[local_row];
      const PetscInt row_end = plan.local_row_offsets[local_row + 1];
      const auto slot_it = std::lower_bound(plan.local_column_indices.begin() + row_begin,
                                            plan.local_column_indices.begin() + row_end,
                                            col);
      libmesh_error_msg_if(slot_it == plan.local_column_indices.begin() + row_end || *slot_it != col,
                           "HilbertSystem Kokkos owned CSR graph is missing a sequential "
                           "remote matrix coupling.");
      return cast_int<std::size_t>(std::distance(plan.local_column_indices.begin(), slot_it));
    }

  if (col >= plan.col_start && col < plan.col_stop)
    {
      const PetscInt local_col = col - plan.col_start;
      const PetscInt row_begin = plan.diag_row_offsets[local_row];
      const PetscInt row_end = plan.diag_row_offsets[local_row + 1];
      const auto slot_it = std::lower_bound(plan.diag_column_indices.begin() + row_begin,
                                            plan.diag_column_indices.begin() + row_end,
                                            local_col);
      libmesh_error_msg_if(slot_it == plan.diag_column_indices.begin() + row_end || *slot_it != local_col,
                           "HilbertSystem Kokkos owned CSR graph is missing a diagonal MPI "
                           "remote matrix coupling.");
      return cast_int<std::size_t>(std::distance(plan.diag_column_indices.begin(), slot_it));
    }

  const auto offdiag_col_it =
    std::lower_bound(plan.offdiag_global_columns.begin(), plan.offdiag_global_columns.end(), col);
  libmesh_error_msg_if(offdiag_col_it == plan.offdiag_global_columns.end() || *offdiag_col_it != col,
                       "HilbertSystem Kokkos owned CSR graph is missing an off-diagonal MPI "
                       "remote matrix column.");
  const PetscInt offdiag_local_col =
    cast_int<PetscInt>(std::distance(plan.offdiag_global_columns.begin(), offdiag_col_it));
  const PetscInt row_begin = plan.offdiag_row_offsets[local_row];
  const PetscInt row_end = plan.offdiag_row_offsets[local_row + 1];
  const auto slot_it = std::lower_bound(plan.offdiag_column_indices.begin() + row_begin,
                                        plan.offdiag_column_indices.begin() + row_end,
                                        offdiag_local_col);
  libmesh_error_msg_if(slot_it == plan.offdiag_column_indices.begin() + row_end ||
                         *slot_it != offdiag_local_col,
                       "HilbertSystem Kokkos owned CSR graph is missing an off-diagonal MPI "
                       "remote matrix coupling.");
  return plan.mat_diag_size +
         cast_int<std::size_t>(std::distance(plan.offdiag_column_indices.begin(), slot_it));
}

bool
build_kokkos_petsc_remote_rhs_sf(::HilbertSystem & sys,
                                 KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build remote RHS SF");

  if (plan.direct_matrix_layout != KokkosDirectMatrixLayout::mpi_aij ||
      plan.rhs_remote_rows.empty())
    return true;

  std::vector<PetscSFNode> remote_nodes(plan.rhs_remote_rows.size());
  for (auto i : index_range(plan.rhs_remote_rows))
    {
      const processor_id_type owner = plan.rhs_remote_owners[i];
      remote_nodes[i].rank = cast_int<PetscMPIInt>(owner);
      remote_nodes[i].index =
        cast_int<PetscInt>(plan.rhs_remote_rows[i] - sys.get_dof_map().first_dof(owner));
    }

  LibmeshPetscCall2(sys.comm(), PetscSFCreate(sys.comm().get(), &plan.rhs_remote_sf));
  LibmeshPetscCall2(sys.comm(),
                    PetscSFSetGraph(plan.rhs_remote_sf,
                                    cast_int<PetscInt>(plan.rhs_local_size),
                                    cast_int<PetscInt>(plan.rhs_remote_rows.size()),
                                    nullptr,
                                    PETSC_COPY_VALUES,
                                    remote_nodes.data(),
                                    PETSC_COPY_VALUES));
  LibmeshPetscCall2(sys.comm(), PetscSFSetUp(plan.rhs_remote_sf));
  plan.rhs_remote_root_values =
    ::Kokkos::View<Number *>("hilbert_rhs_remote_root_values", plan.rhs_local_size);
  return true;
}

bool
build_kokkos_petsc_remote_mat_sf(::HilbertSystem & sys,
                                 KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build remote matrix SF");

  if (plan.direct_matrix_layout != KokkosDirectMatrixLayout::mpi_aij ||
      plan.mat_remote_rows.empty())
    return true;

  const PetscInt local_row_count = plan.row_stop - plan.row_start;
  std::vector<PetscSFNode> setup_nodes(plan.mat_remote_rows.size());
  for (auto i : index_range(plan.mat_remote_rows))
    {
      const processor_id_type owner = plan.mat_remote_owners[i];
      setup_nodes[i].rank = cast_int<PetscMPIInt>(owner);
      setup_nodes[i].index =
        cast_int<PetscInt>(plan.mat_remote_rows[i] - sys.get_dof_map().first_dof(owner));
    }

  PetscSF setup_sf = nullptr;
  LibmeshPetscCall2(sys.comm(), PetscSFCreate(sys.comm().get(), &setup_sf));
  LibmeshPetscCall2(sys.comm(),
                    PetscSFSetGraph(setup_sf,
                                    local_row_count,
                                    cast_int<PetscInt>(plan.mat_remote_rows.size()),
                                    nullptr,
                                    PETSC_COPY_VALUES,
                                    setup_nodes.data(),
                                    PETSC_COPY_VALUES));
  LibmeshPetscCall2(sys.comm(), PetscSFSetUp(setup_sf));

  const PetscInt * degrees = nullptr;
  LibmeshPetscCall2(sys.comm(), PetscSFComputeDegreeBegin(setup_sf, &degrees));
  LibmeshPetscCall2(sys.comm(), PetscSFComputeDegreeEnd(setup_sf, &degrees));

  std::vector<PetscInt> gathered_offsets(local_row_count + 1, 0);
  for (PetscInt local_row = 0; local_row != local_row_count; ++local_row)
    gathered_offsets[local_row + 1] = gathered_offsets[local_row] + degrees[local_row];

  std::vector<PetscInt> gathered_cols(gathered_offsets.back(), -1);
  if (!plan.mat_remote_cols.empty())
    {
      LibmeshPetscCall2(sys.comm(),
                        PetscSFGatherBegin(setup_sf,
                                           MPIU_INT,
                                           plan.mat_remote_cols.data(),
                                           gathered_cols.data()));
      LibmeshPetscCall2(sys.comm(),
                        PetscSFGatherEnd(setup_sf,
                                         MPIU_INT,
                                         plan.mat_remote_cols.data(),
                                         gathered_cols.data()));
    }

  std::vector<PetscInt> reply_slots(gathered_cols.size(), -1);
  for (PetscInt local_row = 0; local_row != local_row_count; ++local_row)
    {
      const PetscInt global_row = plan.row_start + local_row;
      for (PetscInt k = gathered_offsets[local_row]; k != gathered_offsets[local_row + 1]; ++k)
        reply_slots[k] =
          cast_int<PetscInt>(lookup_kokkos_owned_matrix_slot(plan, global_row, gathered_cols[k]));
    }

  std::vector<PetscInt> remote_root_indices(plan.mat_remote_rows.size(), -1);
  if (!remote_root_indices.empty())
    {
      LibmeshPetscCall2(sys.comm(),
                        PetscSFScatterBegin(setup_sf,
                                            MPIU_INT,
                                            reply_slots.data(),
                                            remote_root_indices.data()));
      LibmeshPetscCall2(sys.comm(),
                        PetscSFScatterEnd(setup_sf,
                                          MPIU_INT,
                                          reply_slots.data(),
                                          remote_root_indices.data()));
    }
  LibmeshPetscCall2(sys.comm(), PetscSFDestroy(&setup_sf));

  std::vector<PetscSFNode> remote_nodes(plan.mat_remote_rows.size());
  for (auto i : index_range(plan.mat_remote_rows))
    {
      libmesh_error_msg_if(remote_root_indices[i] < 0,
                           "HilbertSystem Kokkos remote matrix slot setup did not receive an "
                           "owner slot index.");
      remote_nodes[i].rank = cast_int<PetscMPIInt>(plan.mat_remote_owners[i]);
      remote_nodes[i].index = remote_root_indices[i];
    }

  LibmeshPetscCall2(sys.comm(), PetscSFCreate(sys.comm().get(), &plan.mat_remote_sf));
  LibmeshPetscCall2(sys.comm(),
                    PetscSFSetGraph(plan.mat_remote_sf,
                                    cast_int<PetscInt>(plan.mat_diag_size + plan.mat_offdiag_size),
                                    cast_int<PetscInt>(plan.mat_remote_rows.size()),
                                    nullptr,
                                    PETSC_COPY_VALUES,
                                    remote_nodes.data(),
                                    PETSC_COPY_VALUES));
  LibmeshPetscCall2(sys.comm(), PetscSFSetUp(plan.mat_remote_sf));
  plan.mat_remote_root_indices = std::move(remote_root_indices);
  plan.mat_remote_root_values =
    ::Kokkos::View<Number *>("hilbert_mat_remote_root_values",
                             plan.mat_diag_size + plan.mat_offdiag_size);
  sync_kokkos_petsc_remote_slot_views(plan);
  return true;
}

bool
build_kokkos_petsc_remote_sf(::HilbertSystem & sys,
                             KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build remote SF");

  clear_kokkos_petsc_remote_sf(plan);
  if (sys.n_processors() == 1)
    {
      sync_kokkos_petsc_remote_slot_views(plan);
      return true;
    }

  if (!build_kokkos_petsc_remote_rhs_sf(sys, plan))
    return false;
  if (!build_kokkos_petsc_remote_mat_sf(sys, plan))
    return false;
  sync_kokkos_petsc_remote_slot_views(plan);
  return true;
}

bool
kokkos_direct_slot_cache_matches(const KokkosPetscAssemblyPlan & plan,
                                 const KokkosDirectMatrixLayout layout,
                                 const numeric_index_type rhs_first_local,
                                 const numeric_index_type rhs_last_local,
                                 const PetscInt row_start,
                                 const PetscInt row_stop,
                                 const PetscInt col_start,
                                 const PetscInt col_stop)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::check direct slot cache");

  return plan.direct_slot_cache_active &&
         plan.direct_slot_layout == layout &&
         plan.direct_slot_rhs_first_local == rhs_first_local &&
         plan.direct_slot_rhs_last_local == rhs_last_local &&
         plan.direct_slot_row_start == row_start &&
         plan.direct_slot_row_stop == row_stop &&
         plan.direct_slot_col_start == col_start &&
         plan.direct_slot_col_stop == col_stop &&
         plan.rhs_local_slots.extent(0) == plan.total_rhs_entries &&
         plan.mat_value_slots.extent(0) == plan.total_mat_entries;
}

void
mark_kokkos_direct_slot_cache(KokkosPetscAssemblyPlan & plan,
                              const KokkosDirectMatrixLayout layout,
                              const numeric_index_type rhs_first_local,
                              const numeric_index_type rhs_last_local,
                              const PetscInt row_start,
                              const PetscInt row_stop,
                              const PetscInt col_start,
                              const PetscInt col_stop)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::mark direct slot cache");

  plan.direct_slot_layout = layout;
  plan.direct_slot_rhs_first_local = rhs_first_local;
  plan.direct_slot_rhs_last_local = rhs_last_local;
  plan.direct_slot_row_start = row_start;
  plan.direct_slot_row_stop = row_stop;
  plan.direct_slot_col_start = col_start;
  plan.direct_slot_col_stop = col_stop;
  plan.direct_slot_cache_active = true;
}

template <typename IntView>
KOKKOS_INLINE_FUNCTION
PetscInt
kokkos_lower_bound_petsc_index(const IntView values,
                               const PetscInt first,
                               const PetscInt last,
                               const PetscInt target)
{
  PetscInt count = last - first;
  PetscInt base = first;

  while (count > 0)
    {
      const PetscInt step = count / 2;
      const PetscInt middle = base + step;
      if (values(middle) < target)
        {
          base = middle + 1;
          count -= step + 1;
        }
      else
        count = step;
    }

  return base;
}

void
seed_kokkos_direct_remote_slots_from_host(::HilbertSystem & sys,
                                          const DofMap::KokkosDofIndexCache & dof_index_cache,
                                          KokkosPetscAssemblyPlan & plan,
                                          const numeric_index_type rhs_first_local,
                                          const numeric_index_type rhs_last_local)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::seed remote slots from host");

  if (sys.n_processors() == 1)
    return;

  auto h_rhs_slots = ::Kokkos::create_mirror_view(plan.rhs_local_slots);
  auto h_mat_slots = ::Kokkos::create_mirror_view(plan.mat_value_slots);
  std::unordered_map<PetscInt, std::size_t> rhs_remote_slot_map;
  std::unordered_map<std::pair<PetscInt, PetscInt>, std::size_t, PetscIntPairHash> mat_remote_slot_map;

  for (const auto & record : plan.records)
    {
      const auto elem_dofs =
        &dof_index_cache.host_element_dof_indices[record.elem_index * dof_index_cache.max_dofs];

      for (unsigned int i = 0; i != record.n_dofs; ++i)
        {
          const auto row = elem_dofs[i];
          if (row >= rhs_first_local && row < rhs_last_local)
            continue;

          const PetscInt petsc_row = cast_int<PetscInt>(row);
          const auto [it, inserted] =
            rhs_remote_slot_map.emplace(petsc_row, plan.rhs_remote_rows.size());
          if (inserted)
            plan.rhs_remote_rows.push_back(petsc_row);
          h_rhs_slots(record.rhs_offset + i) = (rhs_last_local - rhs_first_local) + it->second;
        }

      if (plan.direct_matrix_layout != KokkosDirectMatrixLayout::mpi_aij)
        continue;

      for (unsigned int i = 0; i != record.n_dofs; ++i)
        {
          const PetscInt row = cast_int<PetscInt>(elem_dofs[i]);
          if (row >= plan.row_start && row < plan.row_stop)
            continue;

          for (unsigned int j = 0; j != record.n_dofs; ++j)
            {
              const PetscInt col = cast_int<PetscInt>(elem_dofs[j]);
              const auto key = std::make_pair(row, col);
              const auto [it, inserted] =
                mat_remote_slot_map.emplace(key, plan.mat_remote_rows.size());
              if (inserted)
                {
                  plan.mat_remote_rows.push_back(row);
                  plan.mat_remote_cols.push_back(col);
                }
              h_mat_slots(record.mat_offset + i * record.n_dofs + j) =
                plan.mat_diag_size + plan.mat_offdiag_size + it->second;
            }
        }
    }

  ::Kokkos::deep_copy(plan.rhs_local_slots, h_rhs_slots);
  ::Kokkos::deep_copy(plan.mat_value_slots, h_mat_slots);
}

bool
bind_kokkos_direct_slots_from_plan_graph(::HilbertSystem & sys,
                                         const DofMap::KokkosDofIndexCache & dof_index_cache,
                                         KokkosPetscAssemblyPlan & plan,
                                         const numeric_index_type rhs_first_local,
                                         const numeric_index_type rhs_last_local)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::bind direct slots");

  plan.rhs_remote_rows.clear();
  plan.mat_remote_rows.clear();
  plan.mat_remote_cols.clear();
  plan.rhs_remote_owners.clear();
  plan.mat_remote_owners.clear();
  plan.rhs_remote_values = {};
  plan.mat_remote_values = {};
  clear_kokkos_petsc_remote_sf(plan);

  plan.rhs_local_slots = KokkosSizeView("hilbert_rhs_local_slots", plan.total_rhs_entries);
  plan.mat_value_slots = KokkosSizeView("hilbert_mat_value_slots", plan.total_mat_entries);
  seed_kokkos_direct_remote_slots_from_host(
    sys, dof_index_cache, plan, rhs_first_local, rhs_last_local);

  const bool is_seq_aij = plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij;
  const auto elem_indices = plan.batch_data.elem_indices;
  const auto elem_n_dofs = plan.batch_data.elem_n_dofs;
  const auto rhs_offsets = plan.batch_data.rhs_offsets;
  const auto mat_offsets = plan.batch_data.mat_offsets;
  const auto element_dof_indices = dof_index_cache.element_dof_indices;
  const auto rhs_slots = plan.rhs_local_slots;
  const auto mat_slots = plan.mat_value_slots;
  const auto local_row_offsets = plan.local_row_offsets_view;
  const auto local_column_indices = plan.local_column_indices_view;
  const auto diag_row_offsets = plan.diag_row_offsets_view;
  const auto diag_column_indices = plan.diag_column_indices_view;
  const auto offdiag_row_offsets = plan.offdiag_row_offsets_view;
  const auto offdiag_column_indices = plan.offdiag_column_indices_view;
  const auto offdiag_global_columns = plan.offdiag_global_columns_view;
  const auto first_local = static_cast<dof_id_type>(rhs_first_local);
  const auto last_local = static_cast<dof_id_type>(rhs_last_local);
  const PetscInt row_start = plan.row_start;
  const PetscInt row_stop = plan.row_stop;
  const PetscInt col_start = plan.col_start;
  const PetscInt col_stop = plan.col_stop;
  const PetscInt n_offdiag_global_columns =
    cast_int<PetscInt>(plan.offdiag_global_columns.size());
  const std::size_t mat_diag_size = plan.mat_diag_size;

  ::Kokkos::parallel_for(
    "hilbert_build_direct_slots",
    ::Kokkos::RangePolicy<>(0, plan.records.size()),
    KOKKOS_LAMBDA(const std::size_t record_index)
    {
      const unsigned int elem_index = elem_indices(record_index);
      const unsigned int n_dofs = elem_n_dofs(record_index);
      const std::size_t rhs_offset = rhs_offsets(record_index);
      const std::size_t mat_offset = mat_offsets(record_index);

      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          const dof_id_type row_id = element_dof_indices(elem_index, i);
          if (row_id >= first_local && row_id < last_local)
            rhs_slots(rhs_offset + i) = static_cast<std::size_t>(row_id - first_local);
        }

      for (unsigned int i = 0; i != n_dofs; ++i)
        {
          const PetscInt row = static_cast<PetscInt>(element_dof_indices(elem_index, i));
          if (!is_seq_aij && (row < row_start || row >= row_stop))
            continue;

          const PetscInt local_row = row - row_start;
          if (is_seq_aij)
            {
              const PetscInt row_begin = local_row_offsets(local_row);
              const PetscInt row_end = local_row_offsets(local_row + 1);
              for (unsigned int j = 0; j != n_dofs; ++j)
                {
                  const PetscInt col =
                    static_cast<PetscInt>(element_dof_indices(elem_index, j));
                  const PetscInt slot =
                    kokkos_lower_bound_petsc_index(local_column_indices, row_begin, row_end, col);
                  mat_slots(mat_offset + i * n_dofs + j) = static_cast<std::size_t>(slot);
                }
              continue;
            }

          const PetscInt diag_row_begin = diag_row_offsets(local_row);
          const PetscInt diag_row_end = diag_row_offsets(local_row + 1);
          const PetscInt offdiag_row_begin = offdiag_row_offsets(local_row);
          const PetscInt offdiag_row_end = offdiag_row_offsets(local_row + 1);
          for (unsigned int j = 0; j != n_dofs; ++j)
            {
              const PetscInt col = static_cast<PetscInt>(element_dof_indices(elem_index, j));
              if (col >= col_start && col < col_stop)
                {
                  const PetscInt local_col = col - col_start;
                  const PetscInt slot =
                    kokkos_lower_bound_petsc_index(diag_column_indices,
                                                   diag_row_begin,
                                                   diag_row_end,
                                                   local_col);
                  mat_slots(mat_offset + i * n_dofs + j) = static_cast<std::size_t>(slot);
                }
              else
                {
                  const PetscInt offdiag_local_col =
                    kokkos_lower_bound_petsc_index(offdiag_global_columns,
                                                   PetscInt(0),
                                                   n_offdiag_global_columns,
                                                   col);
                  const PetscInt offdiag_slot =
                    kokkos_lower_bound_petsc_index(offdiag_column_indices,
                                                   offdiag_row_begin,
                                                   offdiag_row_end,
                                                   offdiag_local_col);
                  mat_slots(mat_offset + i * n_dofs + j) =
                    mat_diag_size + static_cast<std::size_t>(offdiag_slot);
                }
            }
        }
    });
  ::Kokkos::fence();

  return true;
}

bool
ensure_kokkos_petsc_owned_matrix(::HilbertSystem & sys,
                                 KokkosPetscAssemblyPlan & plan,
                                 PetscMatrixBase<Number> & system_matrix)
{
  LOG_SCOPE("L2Kokkos::prepare PETSc COO matrix", "HilbertSystem");
  KokkosProfilingRegion kokkos_region("L2Kokkos::prepare PETSc COO matrix");

  if (plan.graph_matrix_target == &system_matrix)
    return true;

  const char * mat_type = nullptr;
  const char * options_prefix = nullptr;
  LibmeshPetscCall2(sys.comm(), MatGetType(system_matrix.mat(), &mat_type));
  LibmeshPetscCall2(sys.comm(), MatGetOptionsPrefix(system_matrix.mat(), &options_prefix));

  Mat new_mat = nullptr;
  LibmeshPetscCall2(sys.comm(), MatCreate(sys.comm().get(), &new_mat));
  LibmeshPetscCall2(sys.comm(),
                    MatSetSizes(new_mat,
                                cast_int<PetscInt>(system_matrix.local_m()),
                                cast_int<PetscInt>(system_matrix.local_n()),
                                cast_int<PetscInt>(system_matrix.m()),
                                cast_int<PetscInt>(system_matrix.n())));
  LibmeshPetscCall2(sys.comm(), MatSetBlockSize(new_mat, 1));
  if (options_prefix)
    LibmeshPetscCall2(sys.comm(), MatSetOptionsPrefix(new_mat, options_prefix));
  LibmeshPetscCall2(sys.comm(), MatSetType(new_mat, mat_type));
  LibmeshPetscCall2(sys.comm(), MatSetFromOptions(new_mat));

  if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::none)
    {
      LibmeshPetscCall2(sys.comm(), MatDestroy(&new_mat));
      return false;
    }

  {
    KokkosProfilingRegion preallocation_region("PETSc::MatSetPreallocationCOO");
    LibmeshPetscCall2(sys.comm(),
                      MatSetPreallocationCOO(new_mat,
                                             cast_int<PetscCount>(plan.mat_coo_rows.size()),
                                             plan.mat_coo_rows.data(),
                                             plan.mat_coo_cols.data()));
  }

  LibmeshPetscCall2(sys.comm(), MatSetOption(new_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
  {
    KokkosProfilingRegion assembly_region("PETSc::MatAssembly initial COO matrix");
    LibmeshPetscCall2(sys.comm(), MatAssemblyBegin(new_mat, MAT_FINAL_ASSEMBLY));
    LibmeshPetscCall2(sys.comm(), MatAssemblyEnd(new_mat, MAT_FINAL_ASSEMBLY));
  }
  system_matrix.reset_mat(new_mat, true);
  plan.graph_matrix_target = &system_matrix;
  return true;
}

bool
ensure_kokkos_petsc_direct_storage(::HilbertSystem & sys,
                                   const DofMap::KokkosDofIndexCache & dof_index_cache,
                                   KokkosPetscAssemblyPlan & plan,
                                   PetscMatrixBase<Number> & system_matrix,
                                   PetscVector<Number> & system_rhs)
{
  LOG_SCOPE("L2Kokkos::prepare direct PETSc storage", "HilbertSystem");
  KokkosProfilingRegion kokkos_region("L2Kokkos::prepare direct PETSc storage");

  if (plan.direct_matrix_target == &system_matrix &&
      plan.direct_rhs_target == &system_rhs &&
      plan.direct_storage_active)
    return true;

  plan.direct_matrix_target = &system_matrix;
  plan.direct_rhs_target = &system_rhs;
  plan.direct_storage_active = false;

  const char * mat_type = nullptr;
  LibmeshPetscCall2(sys.comm(), MatGetType(system_matrix.mat(), &mat_type));
  PetscBool is_seq_aijkokkos = PETSC_FALSE;
  PetscBool is_mpi_aijkokkos = PETSC_FALSE;
  LibmeshPetscCall2(sys.comm(), PetscStrcmp(mat_type, MATSEQAIJKOKKOS, &is_seq_aijkokkos));
  LibmeshPetscCall2(sys.comm(), PetscStrcmp(mat_type, MATMPIAIJKOKKOS, &is_mpi_aijkokkos));
  if (!is_seq_aijkokkos && !is_mpi_aijkokkos)
    libmesh_error_msg("HilbertSystem Kokkos direct PETSc storage requires PETSc AIJKokkos "
                      "matrices. Run with -mat_type aijkokkos.");
  const KokkosDirectMatrixLayout matrix_layout =
    is_seq_aijkokkos ? KokkosDirectMatrixLayout::seq_aij : KokkosDirectMatrixLayout::mpi_aij;
  plan.direct_matrix_layout = matrix_layout;

  PetscInt rhs_first_local_petsc = 0;
  PetscInt rhs_last_local_petsc = 0;
  LibmeshPetscCall2(sys.comm(),
                    VecGetOwnershipRange(system_rhs.vec(), &rhs_first_local_petsc, &rhs_last_local_petsc));
  const numeric_index_type rhs_first_local = cast_int<numeric_index_type>(rhs_first_local_petsc);
  const numeric_index_type rhs_last_local = cast_int<numeric_index_type>(rhs_last_local_petsc);
  PetscInt row_start = 0;
  PetscInt row_stop = 0;
  PetscInt col_start = 0;
  PetscInt col_stop = 0;
  LibmeshPetscCall2(sys.comm(), MatGetOwnershipRange(system_matrix.mat(), &row_start, &row_stop));
  LibmeshPetscCall2(sys.comm(),
                    MatGetOwnershipRangeColumn(system_matrix.mat(), &col_start, &col_stop));

  const bool slots_cached =
    kokkos_direct_slot_cache_matches(plan,
                                     matrix_layout,
                                     rhs_first_local,
                                     rhs_last_local,
                                     row_start,
                                     row_stop,
                                     col_start,
                                     col_stop) &&
    plan.graph_matrix_target == &system_matrix;
  if (!slots_cached)
    {
      plan.direct_matrix_layout = KokkosDirectMatrixLayout::none;
      plan.rhs_local_size = 0;
      plan.mat_diag_size = 0;
      plan.mat_offdiag_size = 0;
      plan.rhs_remote_rows.clear();
      plan.mat_remote_rows.clear();
      plan.mat_remote_cols.clear();
      plan.rhs_remote_owners.clear();
      plan.mat_remote_owners.clear();
      plan.rhs_remote_values = {};
      plan.mat_remote_values = {};
      plan.mat_coo_values = {};
      plan.rhs_local_slots = {};
      plan.mat_value_slots = {};
      plan.graph_matrix_target = nullptr;
      plan.direct_slot_cache_active = false;
      plan.direct_matrix_layout = matrix_layout;

      {
        LOG_SCOPE("L2Kokkos::direct storage graph setup", "HilbertSystem");
        KokkosProfilingRegion graph_region("L2Kokkos::direct storage graph setup");
        if (!build_kokkos_petsc_owned_csr_graph(
              sys, plan, matrix_layout, row_start, row_stop, col_start, col_stop))
          return false;
      }

      plan.rhs_local_size = cast_int<std::size_t>(rhs_last_local - rhs_first_local);
      plan.mat_diag_size =
        plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij ?
          cast_int<std::size_t>(plan.local_row_offsets.empty() ? 0 : plan.local_row_offsets.back()) :
          cast_int<std::size_t>(plan.diag_row_offsets.empty() ? 0 : plan.diag_row_offsets.back());
      plan.mat_offdiag_size =
        plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij ?
          0 :
          cast_int<std::size_t>(plan.offdiag_row_offsets.empty() ? 0 :
                                  plan.offdiag_row_offsets.back());
      {
        KokkosProfilingRegion coo_region("L2Kokkos::direct storage COO setup");
        build_kokkos_petsc_owned_coo(plan);
      }

      {
        LOG_SCOPE("L2Kokkos::direct storage slot setup", "HilbertSystem");
        KokkosProfilingRegion slot_region("L2Kokkos::direct storage slot setup");
        if (!bind_kokkos_direct_slots_from_plan_graph(
              sys, dof_index_cache, plan, rhs_first_local, rhs_last_local))
          return false;
        finalize_kokkos_direct_slots(sys, plan);
      }

      mark_kokkos_direct_slot_cache(plan,
                                    matrix_layout,
                                    rhs_first_local,
                                    rhs_last_local,
                                    row_start,
                                    row_stop,
                                    col_start,
                                    col_stop);
    }

  {
    LOG_SCOPE("L2Kokkos::direct storage PETSc matrix setup", "HilbertSystem");
    KokkosProfilingRegion matrix_region("L2Kokkos::direct storage PETSc matrix setup");
    if (!ensure_kokkos_petsc_owned_matrix(sys, plan, system_matrix))
      return false;
  }

  {
    LOG_SCOPE("L2Kokkos::direct storage remote setup", "HilbertSystem");
    KokkosProfilingRegion remote_region("L2Kokkos::direct storage remote setup");
    if (!build_kokkos_petsc_remote_sf(sys, plan))
      return false;
  }
  plan.direct_storage_active = true;
  return true;
}

namespace
{

bool
build_kokkos_element_state(::HilbertSystem & sys,
                           const Elem & elem,
                           KokkosElementAssemblyState & state)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build element state");

  auto * solution_vector = dynamic_cast<PetscVector<Number> *>(sys.current_local_solution.get());
  if (!solution_vector || !solution_vector->supports_kokkos_access())
    return false;

  const DofMap & dof_map = sys.get_dof_map();
  state.geometry_cache = &sys.get_mesh().get_kokkos_geometry_cache();
  state.elem_index = sys.get_mesh().get_kokkos_elem_index(elem);
  if (state.elem_index == libMesh::invalid_uint)
    return false;

  state.elem_type = state.geometry_cache->host_element_types[state.elem_index];
  state.mapping_type = state.geometry_cache->host_element_mapping_types[state.elem_index];
  state.elem_n_nodes = state.geometry_cache->host_element_n_nodes[state.elem_index];
  state.elem_p_level = state.geometry_cache->host_element_p_levels[state.elem_index];

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
  KokkosProfilingRegion kokkos_region("L2Kokkos::build FEM goal data");

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

void
build_hilbert_batch_data(const std::vector<HilbertElementAssemblyRecord> & records,
                         const FEType & fe_type,
                         const MeshBase::KokkosGeometryCache & geometry_cache,
                         KokkosHilbertBatchData & batch_data)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build batch data");

  batch_data.elem_indices = KokkosUnsignedIntView("hilbert_elem_indices", records.size());
  batch_data.elem_n_dofs = KokkosUnsignedIntView("hilbert_elem_n_dofs", records.size());
  batch_data.quadrature_orders =
    KokkosUnsignedIntView("hilbert_quadrature_orders", records.size());
  batch_data.rhs_offsets = KokkosSizeView("hilbert_rhs_offsets", records.size());
  batch_data.mat_offsets = KokkosSizeView("hilbert_mat_offsets", records.size());
  batch_data.shape_keys = KokkosShapeKeyRecordView("hilbert_shape_keys", records.size());

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

  const auto elem_indices = batch_data.elem_indices;
  const auto element_types = geometry_cache.element_types_dual.view_device();
  const auto element_p_levels = geometry_cache.element_p_levels_dual.view_device();
  const auto shape_keys = batch_data.shape_keys;
  const auto family = fe_type.family;
  const auto base_order = fe_type.order.get_order();
  const auto n_records = records.size();
  using GeometryExecutionSpace = typename std::decay_t<decltype(element_types)>::execution_space;

  ::Kokkos::parallel_for(
    "hilbert_build_shape_keys",
    ::Kokkos::RangePolicy<GeometryExecutionSpace>(0, cast_int<int>(n_records)),
    KOKKOS_LAMBDA(const int raw_record_index)
    {
      const unsigned int record_index = static_cast<unsigned int>(raw_record_index);
      const unsigned int elem_index = elem_indices(record_index);
      shape_keys(record_index) =
        FEShapeKey{family,
                   element_types(elem_index),
                   static_cast<Order>(base_order + static_cast<int>(element_p_levels(elem_index)))};
    });
  ::Kokkos::fence();
}

bool
kokkos_petsc_plan_matches(const ::HilbertSystem & sys,
                          const MeshBase::KokkosGeometryCache & geometry_cache,
                          const DofMap::KokkosDofIndexCache & dof_index_cache,
                          const KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::check assembly plan cache");

  return plan.geometry_cache_id == &geometry_cache &&
         plan.dof_index_cache_id == &dof_index_cache &&
         plan.fe_type == sys.variable_type(0) &&
         plan.hilbert_order == sys.hilbert_order() &&
         plan.extra_quadrature_order == sys.extra_quadrature_order &&
         plan.subdomains == sys.subdomains_list();
}

bool
build_kokkos_petsc_assembly_plan(::HilbertSystem & sys,
                                 KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build assembly plan");

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

  KokkosHilbertBatchData batch_data;
  build_hilbert_batch_data(records, sys.variable_type(0), geometry_cache, batch_data);

  plan.records = std::move(records);
  plan.total_rhs_entries = total_rhs_entries;
  plan.total_mat_entries = total_mat_entries;
  plan.batch_data = batch_data;
  plan.rhs_entry_values = {};
  plan.mat_entry_values = {};
  plan.mat_coo_rows.clear();
  plan.mat_coo_cols.clear();
  plan.rhs_local_slots = KokkosSizeView();
  plan.mat_value_slots = KokkosSizeView();
  plan.geometry_cache_id = &geometry_cache;
  plan.dof_index_cache_id = dof_index_cache;
  plan.fe_type = sys.variable_type(0);
  plan.hilbert_order = sys.hilbert_order();
  plan.extra_quadrature_order = sys.extra_quadrature_order;
  plan.subdomains = sys.subdomains_list();
  plan.graph_matrix_target = nullptr;
  plan.direct_matrix_target = nullptr;
  plan.direct_rhs_target = nullptr;
  plan.fem_goal_target = nullptr;
  plan.input_vector_target = nullptr;
  plan.analytic_goal_target = nullptr;
  plan.analytic_goal_gradient_target = nullptr;
  plan.analytic_goal_time = std::numeric_limits<Real>::quiet_NaN();
  cache_single_supported_hilbert_bucket(plan, geometry_cache);
  plan.direct_slot_layout = KokkosDirectMatrixLayout::none;
  plan.direct_slot_rhs_first_local = 0;
  plan.direct_slot_rhs_last_local = 0;
  plan.direct_slot_row_start = 0;
  plan.direct_slot_row_stop = 0;
  plan.direct_slot_col_start = 0;
  plan.direct_slot_col_stop = 0;
  plan.direct_slot_cache_active = false;
  plan.direct_storage_active = false;
  return true;
}

bool
build_kokkos_fem_goal_batch_data(::HilbertSystem & sys,
                                 const KokkosPetscAssemblyPlan & plan,
                                 const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                 KokkosFEMGoalBatchData & batch_data)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::build FEM goal batch data");

  libmesh_assert(sys.input_system);
  System & input_system = *sys.input_system;
  const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
  const unsigned int n_fields = goal_function.n_field_variables();
  if (n_fields > kokkos_parsed_fem_max_fields)
    return false;

  auto * input_vector = dynamic_cast<PetscVector<Number> *>(input_system.current_local_solution.get());
  if (!input_vector || !input_vector->supports_kokkos_access())
    return false;

  batch_data.input_vector = input_vector;
  batch_data.field_local_indices.size = n_fields;
  batch_data.field_keys =
    KokkosFieldKeyRecordView("hilbert_fem_goal_field_keys", n_fields, plan.records.size());
  batch_data.field_dofs =
    KokkosFieldDofRecordView("hilbert_fem_goal_field_dofs", n_fields, plan.records.size());
  auto h_field_keys = ::Kokkos::create_mirror_view(batch_data.field_keys);
  auto h_field_dofs = ::Kokkos::create_mirror_view(batch_data.field_dofs);

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

      for (auto record_index : index_range(plan.records))
        {
          const auto & record = plan.records[record_index];
          const ElemType elem_type = geometry_cache.host_element_types[record.elem_index];
          const unsigned int elem_p_level = geometry_cache.host_element_p_levels[record.elem_index];
          const ElemMappingType mapping_type =
            geometry_cache.host_element_mapping_types[record.elem_index];
          const FEShapeKey field_key = make_hilbert_shape_key(elem_type, elem_p_level, field_type);
          if (mapping_type != LAGRANGE_MAP ||
              !supports_shape_with_lagrange_map(field_key) ||
              (sys.hilbert_order() > 0 && !supports_grad_shape(field_key)))
            return false;

          const unsigned int field_n_dofs =
            FEInterface::n_dofs(libMesh::Kokkos::dim_from_topology(elem_type), field_type, elem_type);
          if (field_n_dofs > kokkos_hilbert_max_dofs)
            return false;

          h_field_keys(field, record_index) = field_key;
          h_field_dofs(field, record_index) = field_n_dofs;
        }
    }

  ::Kokkos::deep_copy(batch_data.field_keys, h_field_keys);
  ::Kokkos::deep_copy(batch_data.field_dofs, h_field_dofs);

  return true;
}

bool
ensure_kokkos_fem_goal_batch_data(::HilbertSystem & sys,
                                  const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function,
                                  KokkosPetscAssemblyPlan & plan)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure FEM goal batch data");

  libmesh_assert(sys.input_system);
  const void * input_vector_target = sys.input_system->current_local_solution.get();
  if (plan.fem_goal_target == &goal_function && plan.input_vector_target == input_vector_target)
    return true;

  KokkosFEMGoalBatchData batch_data;
  if (!build_kokkos_fem_goal_batch_data(sys, plan, goal_function, batch_data))
    return false;

  plan.fem_goal_batch_data = batch_data;
  plan.fem_goal_target = &goal_function;
  plan.input_vector_target = input_vector_target;
  return true;
}

bool
ensure_kokkos_sampled_analytic_goal_data(::HilbertSystem & sys,
                                         FunctionBase<Number> & goal_function,
                                         FunctionBase<Gradient> * goal_gradient,
                                         KokkosPetscAssemblyPlan & plan)
{
  LOG_SCOPE("L2Kokkos::sample analytic goal", "HilbertSystem");
  KokkosProfilingRegion kokkos_region("L2Kokkos::sample analytic goal");

  if (plan.analytic_goal_target == &goal_function &&
      plan.analytic_goal_gradient_target == goal_gradient &&
      plan.analytic_goal_time == sys.time)
    return true;

  if (sys.hilbert_order() > 0 && !goal_gradient)
    return false;

  const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
  const auto n_elements = geometry_cache.host_element_ids.size();
  std::vector<std::size_t> host_offsets(n_elements, std::numeric_limits<std::size_t>::max());
  std::vector<Number> host_values;
  std::vector<std::array<Real, LIBMESH_DIM>> host_gradients;

  {
    KokkosProfilingRegion host_sample_region("L2Kokkos::sample analytic goal on host");
    for (const auto & record : plan.records)
      {
        const Elem * elem = sys.get_mesh().query_elem_ptr(record.elem_id);
        if (!elem)
          return false;

        std::unique_ptr<QBase> qrule =
          QBase::build(QGAUSS, elem->dim(), static_cast<Order>(record.quadrature_order));
        std::unique_ptr<FEBase> fe = FEBase::build(elem->dim(), sys.variable_type(0));
        fe->attach_quadrature_rule(qrule.get());
        const auto & xyz = fe->get_xyz();
        fe->reinit(elem);

        host_offsets[record.elem_index] = host_values.size();

        for (const auto qp : index_range(xyz))
          {
            host_values.push_back(goal_function(xyz[qp], sys.time));

            std::array<Real, LIBMESH_DIM> grad_values = {};
            if (goal_gradient)
              {
                const Gradient grad = (*goal_gradient)(xyz[qp], sys.time);
                for (unsigned int d = 0; d != LIBMESH_DIM; ++d)
                  grad_values[d] = grad(d);
              }
            host_gradients.push_back(grad_values);
          }
      }
  }

  plan.analytic_goal_qp_offsets =
    KokkosSizeView("hilbert_analytic_goal_qp_offsets", host_offsets.size());
  plan.analytic_goal_values =
    KokkosScalarView<Number>("hilbert_analytic_goal_values", host_values.size());
  plan.analytic_goal_gradients =
    KokkosGradientView("hilbert_analytic_goal_gradients", host_gradients.size(), LIBMESH_DIM);

  KokkosProfilingRegion copy_region("L2Kokkos::copy analytic samples to device");
  auto h_offsets = ::Kokkos::create_mirror_view(plan.analytic_goal_qp_offsets);
  for (auto i : index_range(host_offsets))
    h_offsets(i) = host_offsets[i];
  ::Kokkos::deep_copy(plan.analytic_goal_qp_offsets, h_offsets);

  auto h_values = ::Kokkos::create_mirror_view(plan.analytic_goal_values);
  for (auto i : index_range(host_values))
    h_values(i) = host_values[i];
  ::Kokkos::deep_copy(plan.analytic_goal_values, h_values);

  auto h_gradients = ::Kokkos::create_mirror_view(plan.analytic_goal_gradients);
  for (auto i : index_range(host_gradients))
    for (unsigned int d = 0; d != LIBMESH_DIM; ++d)
      h_gradients(i, d) = host_gradients[i][d];
  ::Kokkos::deep_copy(plan.analytic_goal_gradients, h_gradients);

  plan.analytic_goal_target = &goal_function;
  plan.analytic_goal_gradient_target = goal_gradient;
  plan.analytic_goal_time = sys.time;
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
  KokkosProfilingRegion kokkos_region("L2Kokkos::assemble analytic element");

  const unsigned int n_dofs = F.size();
  KokkosScalarView<Number> d_F("hilbert_residual", n_dofs);
  KokkosDenseJacobianView d_K("hilbert_jacobian", n_dofs, n_dofs);
  auto goal_copy = goal_function;
  auto goal_grad = goal_copy.gradient_function();
  const auto goal_access =
    libMesh::Kokkos::detail::make_hilbert_analytic_goal_access(std::move(goal_copy),
                                                               std::move(goal_grad));
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
  KokkosProfilingRegion kokkos_region("L2Kokkos::assemble FEM goal element");

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
assemble_host_exact_parsed_fem_goal_element(::HilbertSystem & sys,
                                            FEMContext & c,
                                            const bool request_jacobian,
                                            DenseSubVector<Number> & F,
                                            DenseSubMatrix<Number> & K,
                                            const libMesh::Kokkos::KokkosParsedFEMFunction<Number> & goal_function)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::assemble host exact FEM goal element");

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

template <typename LocalRhsView,
          typename RemoteRhsView,
          typename DiagMatrixView,
          typename OffdiagMatrixView,
          typename RemoteMatrixView>
void
scatter_kokkos_hilbert_entry_values(const KokkosPetscAssemblyPlan & plan,
                                    const std::size_t mat_offdiag_base,
                                    LocalRhsView local_rhs_values,
                                    RemoteRhsView remote_rhs_values,
                                    DiagMatrixView diag_matrix_values,
                                    OffdiagMatrixView offdiag_matrix_values,
                                    RemoteMatrixView remote_matrix_values)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::launch scatter entry values");

  const std::size_t n_rhs_entries = plan.total_rhs_entries;
  const std::size_t n_mat_entries = plan.total_mat_entries;
  const std::size_t rhs_local_size = plan.rhs_local_size;
  const std::size_t mat_diag_size = plan.mat_diag_size;
  const auto rhs_slots = plan.rhs_local_slots;
  const auto mat_slots = plan.mat_value_slots;
  const auto rhs_entry_values = plan.rhs_entry_values;
  const auto mat_entry_values = plan.mat_entry_values;
  ::Kokkos::parallel_for(
    "hilbert_scatter_entry_values",
    ::Kokkos::RangePolicy<>(0, n_rhs_entries + n_mat_entries),
    KOKKOS_LAMBDA(const std::size_t entry_index)
    {
      if (entry_index < n_rhs_entries)
        {
          const std::size_t slot = rhs_slots(entry_index);
          const Number value = rhs_entry_values(entry_index);
          if (slot < rhs_local_size)
            ::Kokkos::atomic_add(&local_rhs_values(slot), value);
          else
            ::Kokkos::atomic_add(&remote_rhs_values(slot - rhs_local_size), value);
          return;
        }

      const std::size_t mat_entry = entry_index - n_rhs_entries;
      const std::size_t slot = mat_slots(mat_entry);
      const Number value = mat_entry_values(mat_entry);
      if (slot < mat_diag_size)
        ::Kokkos::atomic_add(&diag_matrix_values(slot), value);
      else if (slot < mat_offdiag_base)
        ::Kokkos::atomic_add(&offdiag_matrix_values(slot - mat_diag_size), value);
      else
        ::Kokkos::atomic_add(&remote_matrix_values(slot - mat_offdiag_base), value);
    });
}

template <unsigned int MaxDofs,
          typename GoalAccess,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage>
void
run_kokkos_hilbert_debug_bucket_value_batch(const KokkosHilbertHomogeneousBucket & bucket,
                                            const NodeCoordinateStorage & node_coordinates,
                                            const ElemNodeIdStorage & element_node_ids,
                                            KokkosPetscAssemblyPlan & plan,
                                            const unsigned int hilbert_order,
                                            GoalAccess goal_access,
                                            KokkosScalarView<Number> rhs_values,
                                            KokkosScalarView<Number> mat_values,
                                            const char * const kernel_name)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::debug launch bucket value batch");

  libMesh::Kokkos::detail::run_hilbert_system_bucket_value_batch<MaxDofs>(
    bucket.key,
    bucket.mapping_type,
    bucket.n_nodes,
    bucket.quadrature_order,
    node_coordinates,
    element_node_ids,
    plan.batch_data.elem_indices,
    plan.batch_data.elem_n_dofs,
    plan.batch_data.rhs_offsets,
    plan.batch_data.mat_offsets,
    hilbert_order,
    goal_access,
    rhs_values,
    mat_values,
    kernel_name);
}

template <typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage>
void
run_kokkos_hilbert_bucket_smoke_check(const KokkosHilbertHomogeneousBucket & bucket,
                                      const NodeCoordinateStorage & node_coordinates,
                                      const ElemNodeIdStorage & element_node_ids,
                                      const ElemMappingTypeStorage & element_mapping_types,
                                      const ElemNodeCountStorage & element_n_nodes,
                                      const KokkosPetscAssemblyPlan & plan,
                                      const unsigned int hilbert_order)
{
  if (plan.records.empty())
    return;

  KokkosProfilingRegion kokkos_region("L2Kokkos::debug bucket FE smoke check");

  KokkosUnsignedIntView smoke_uints("hilbert_debug_bucket_smoke_uints", 12);
  KokkosScalarView<Real> smoke_reals("hilbert_debug_bucket_smoke_reals", 4);
  ::Kokkos::deep_copy(smoke_uints, 0u);
  ::Kokkos::deep_copy(smoke_reals, Real(0));

  const auto elem_indices = plan.batch_data.elem_indices;
  const auto elem_n_dofs = plan.batch_data.elem_n_dofs;
  const auto shape_keys = plan.batch_data.shape_keys;
  const auto quadrature_orders = plan.batch_data.quadrature_orders;
  const FEShapeKey bucket_key = bucket.key;
  const ElemMappingType bucket_mapping_type = bucket.mapping_type;
  const unsigned int bucket_n_nodes = bucket.n_nodes;
  const unsigned int bucket_quadrature_order = bucket.quadrature_order;

  ::Kokkos::parallel_for(
    "hilbert_debug_bucket_fe_smoke",
    ::Kokkos::RangePolicy<>(0, 1),
    KOKKOS_LAMBDA(const int)
    {
      const unsigned int record_index = 0;
      const unsigned int elem_index = elem_indices(record_index);
      const auto elem_nodes =
        libMesh::Kokkos::detail::make_element_node_access(node_coordinates,
                                                          element_node_ids,
                                                          elem_index);

      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> bucket_fe(
        bucket_key,
        bucket_mapping_type,
        elem_nodes,
        bucket_n_nodes,
        bucket_quadrature_order,
        elem_index);

      const FEShapeKey record_key = shape_keys(record_index);
      const libMesh::Kokkos::detail::HilbertFEAccess<decltype(elem_nodes)> record_fe(
        record_key,
        element_mapping_types(elem_index),
        elem_nodes,
        element_n_nodes(elem_index),
        quadrature_orders(record_index),
        elem_index);

      smoke_uints(0) = elem_index;
      smoke_uints(1) = bucket_fe.n_qpoints();
      smoke_uints(2) = bucket_fe.n_dofs();
      smoke_uints(3) = record_fe.n_qpoints();
      smoke_uints(4) = record_fe.n_dofs();
      smoke_uints(5) = elem_n_dofs(record_index);
      smoke_uints(6) = bucket_quadrature_order;
      smoke_uints(7) = quadrature_orders(record_index);
      smoke_uints(8) = bucket_n_nodes;
      smoke_uints(9) = element_n_nodes(elem_index);
      smoke_uints(10) = static_cast<unsigned int>(bucket_key.elem_type);
      smoke_uints(11) = static_cast<unsigned int>(record_key.elem_type);

      if (bucket_fe.n_qpoints() && bucket_fe.n_dofs())
        {
          const auto qp_data = bucket_fe.qp_data(0, hilbert_order > 0);
          smoke_reals(0) = qp_data.JxW();
          smoke_reals(1) = qp_data.phi(0);
        }

      if (record_fe.n_qpoints() && record_fe.n_dofs())
        {
          const auto qp_data = record_fe.qp_data(0, hilbert_order > 0);
          smoke_reals(2) = qp_data.JxW();
          smoke_reals(3) = qp_data.phi(0);
        }
    });
  ::Kokkos::fence();

  auto h_smoke_uints = ::Kokkos::create_mirror_view(smoke_uints);
  auto h_smoke_reals = ::Kokkos::create_mirror_view(smoke_reals);
  ::Kokkos::deep_copy(h_smoke_uints, smoke_uints);
  ::Kokkos::deep_copy(h_smoke_reals, smoke_reals);

  libMesh::out << "L2Kokkos bucket smoke: "
               << "elem_index=" << h_smoke_uints(0)
               << " bucket_n_qpoints=" << h_smoke_uints(1)
               << " bucket_fe_n_dofs=" << h_smoke_uints(2)
               << " record_n_qpoints=" << h_smoke_uints(3)
               << " record_fe_n_dofs=" << h_smoke_uints(4)
               << " record_n_dofs=" << h_smoke_uints(5)
               << " bucket_quadrature_order=" << h_smoke_uints(6)
               << " record_quadrature_order=" << h_smoke_uints(7)
               << " bucket_n_nodes=" << h_smoke_uints(8)
               << " record_n_nodes=" << h_smoke_uints(9)
               << " bucket_elem_type=" << h_smoke_uints(10)
               << " record_elem_type=" << h_smoke_uints(11)
               << " bucket_JxW0=" << h_smoke_reals(0)
               << " bucket_phi00=" << h_smoke_reals(1)
               << " record_JxW0=" << h_smoke_reals(2)
               << " record_phi00=" << h_smoke_reals(3)
               << '\n';
}

template <typename ValueView>
Real
kokkos_abs_sum(const ValueView & values)
{
  Real sum = 0;
  ::Kokkos::parallel_reduce(
    "hilbert_debug_abs_sum",
    ::Kokkos::RangePolicy<>(0, values.extent(0)),
    KOKKOS_LAMBDA(const std::size_t i, Real & local_sum)
    {
      const Real value = static_cast<Real>(values(i));
      local_sum += value < Real(0) ? -value : value;
    },
    sum);
  return sum;
}

template <typename ValueView>
Real
kokkos_max_abs_diff(const ValueView & a,
                    const ValueView & b)
{
  Real max_diff = 0;
  ::Kokkos::parallel_reduce(
    "hilbert_debug_max_abs_diff",
    ::Kokkos::RangePolicy<>(0, a.extent(0)),
    KOKKOS_LAMBDA(const std::size_t i, Real & local_max)
    {
      Real diff = static_cast<Real>(a(i) - b(i));
      diff = diff < Real(0) ? -diff : diff;
      if (diff > local_max)
        local_max = diff;
    },
    ::Kokkos::Max<Real>(max_diff));
  return max_diff;
}

template <unsigned int MaxDofs,
          typename GoalAccess,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage>
void
run_kokkos_hilbert_record_value_batch(const NodeCoordinateStorage & node_coordinates,
                                      const ElemNodeIdStorage & element_node_ids,
                                      const ElemMappingTypeStorage & element_mapping_types,
                                      const ElemNodeCountStorage & element_n_nodes,
                                      KokkosPetscAssemblyPlan & plan,
                                      const unsigned int hilbert_order,
                                      GoalAccess goal_access,
                                      const char * const kernel_name)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::launch record value batch");

  libMesh::Kokkos::detail::run_hilbert_system_value_batch<MaxDofs>(
    node_coordinates,
    element_node_ids,
    element_mapping_types,
    element_n_nodes,
    plan.batch_data.shape_keys,
    plan.batch_data.elem_indices,
    plan.batch_data.elem_n_dofs,
    plan.batch_data.quadrature_orders,
    plan.batch_data.rhs_offsets,
    plan.batch_data.mat_offsets,
    hilbert_order,
    goal_access,
    plan.rhs_entry_values,
    plan.mat_entry_values,
    kernel_name);
}

template <typename GoalAccess,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage>
bool
try_kokkos_hilbert_specialized_record_value_batch(
  const KokkosPetscAssemblyPlan & const_plan,
  const NodeCoordinateStorage & node_coordinates,
  const ElemNodeIdStorage & element_node_ids,
  const ElemMappingTypeStorage & element_mapping_types,
  const ElemNodeCountStorage & element_n_nodes,
  KokkosPetscAssemblyPlan & plan,
  const unsigned int hilbert_order,
  GoalAccess goal_access)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::try specialized record value batch");

  KokkosHilbertHomogeneousBucket bucket;
  if (!get_cached_single_supported_hilbert_bucket(const_plan, bucket))
    return false;

  switch (bucket.n_dofs)
    {
    case 8:
      run_kokkos_hilbert_record_value_batch<8>(node_coordinates,
                                               element_node_ids,
                                               element_mapping_types,
                                               element_n_nodes,
                                               plan,
                                               hilbert_order,
                                               goal_access,
                                               "hilbert_record8_values");
      return true;
    case 27:
      run_kokkos_hilbert_record_value_batch<27>(node_coordinates,
                                                element_node_ids,
                                                element_mapping_types,
                                                element_n_nodes,
                                                plan,
                                                hilbert_order,
                                                goal_access,
                                                "hilbert_record27_values");
      return true;
    default:
      return false;
    }
}

template <typename GoalAccess,
          typename NodeCoordinateStorage,
          typename ElemNodeIdStorage,
          typename ElemMappingTypeStorage,
          typename ElemNodeCountStorage>
void
compare_kokkos_hilbert_bucket_value_batch(const KokkosPetscAssemblyPlan & const_plan,
                                          const NodeCoordinateStorage & node_coordinates,
                                          const ElemNodeIdStorage & element_node_ids,
                                          const ElemMappingTypeStorage & element_mapping_types,
                                          const ElemNodeCountStorage & element_n_nodes,
                                          KokkosPetscAssemblyPlan & plan,
                                          const unsigned int hilbert_order,
                                          GoalAccess goal_access)
{
  if (!debug_kokkos_bucket_value_batch())
    return;

  KokkosProfilingRegion kokkos_region("L2Kokkos::debug compare bucket value batch");

  KokkosHilbertHomogeneousBucket bucket;
  if (!get_cached_single_supported_hilbert_bucket(const_plan, bucket))
    {
      libMesh::out << "L2Kokkos bucket debug: no homogeneous bucket available\n";
      return;
    }

  run_kokkos_hilbert_bucket_smoke_check(bucket,
                                        node_coordinates,
                                        element_node_ids,
                                        element_mapping_types,
                                        element_n_nodes,
                                        plan,
                                        hilbert_order);

  KokkosScalarView<Number> bucket_rhs_values("hilbert_debug_bucket_rhs_values",
                                             plan.total_rhs_entries);
  KokkosScalarView<Number> bucket_mat_values("hilbert_debug_bucket_mat_values",
                                             plan.total_mat_entries);
  ::Kokkos::deep_copy(bucket_rhs_values, Number(0));
  ::Kokkos::deep_copy(bucket_mat_values, Number(0));

  switch (bucket.n_dofs)
    {
    case 8:
      run_kokkos_hilbert_debug_bucket_value_batch<8>(bucket,
                                                     node_coordinates,
                                                     element_node_ids,
                                                     plan,
                                                     hilbert_order,
                                                     goal_access,
                                                     bucket_rhs_values,
                                                     bucket_mat_values,
                                                     "hilbert_debug_bucket8_values");
      break;
    case 27:
      run_kokkos_hilbert_debug_bucket_value_batch<27>(bucket,
                                                      node_coordinates,
                                                      element_node_ids,
                                                      plan,
                                                      hilbert_order,
                                                      goal_access,
                                                      bucket_rhs_values,
                                                      bucket_mat_values,
                                                      "hilbert_debug_bucket27_values");
      break;
    default:
      libMesh::out << "L2Kokkos bucket debug: unsupported bucket n_dofs="
                   << bucket.n_dofs << '\n';
      return;
    }

  ::Kokkos::fence();

  const Real record_rhs_abs_sum = kokkos_abs_sum(plan.rhs_entry_values);
  const Real bucket_rhs_abs_sum = kokkos_abs_sum(bucket_rhs_values);
  const Real rhs_max_abs_diff = kokkos_max_abs_diff(plan.rhs_entry_values, bucket_rhs_values);
  const Real record_mat_abs_sum = kokkos_abs_sum(plan.mat_entry_values);
  const Real bucket_mat_abs_sum = kokkos_abs_sum(bucket_mat_values);
  const Real mat_max_abs_diff = kokkos_max_abs_diff(plan.mat_entry_values, bucket_mat_values);

  libMesh::out << "L2Kokkos bucket debug: "
               << "records=" << plan.records.size()
               << " family=" << static_cast<int>(bucket.key.family)
               << " elem_type=" << static_cast<int>(bucket.key.elem_type)
               << " order=" << static_cast<int>(bucket.key.order)
               << " mapping_type=" << static_cast<int>(bucket.mapping_type)
               << " n_nodes=" << bucket.n_nodes
               << " quadrature_order=" << bucket.quadrature_order
               << " n_dofs=" << bucket.n_dofs
               << " record_rhs_abs_sum=" << record_rhs_abs_sum
               << " bucket_rhs_abs_sum=" << bucket_rhs_abs_sum
               << " rhs_max_abs_diff=" << rhs_max_abs_diff
               << " record_mat_abs_sum=" << record_mat_abs_sum
               << " bucket_mat_abs_sum=" << bucket_mat_abs_sum
               << " mat_max_abs_diff=" << mat_max_abs_diff
               << '\n';
}

bool
assemble_kokkos_petsc_global_system(::HilbertSystem & sys,
                                    KokkosPetscAssemblyPlan & plan,
                                    const bool analytic_goal,
                                    const libMesh::Kokkos::KokkosParsedFunction<Number> *
                                      analytic_goal_func,
                                    const libMesh::Kokkos::KokkosParsedFEMFunction<Number> * fem_goal,
                                    PetscMatrixBase<Number> & system_matrix,
                                    PetscVector<Number> & system_rhs,
                                    ::HilbertSystem::KokkosAssemblyPath & assembly_path)
{
  LOG_SCOPE("L2Kokkos::assemble PETSc global system", "HilbertSystem");
  KokkosProfilingRegion kokkos_region("L2Kokkos::assemble PETSc global system");

  libmesh_error_msg_if(sys.has_static_condensation(),
                       "HilbertSystem Kokkos direct PETSc storage does not support static "
                       "condensation.");
  libmesh_error_msg_if(sys.get_dof_map().n_constrained_dofs(),
                       "HilbertSystem Kokkos direct PETSc storage does not yet support "
                       "constrained dofs.");
  libmesh_error_msg_if(!analytic_goal && !fem_goal,
                       "HilbertSystem Kokkos direct PETSc storage requires an analytic "
                       "goal or parsed FEM goal.");
  libmesh_ignore(analytic_goal_func);

  const auto * dof_index_cache = sys.get_dof_map().get_kokkos_dof_index_cache(0);
  libmesh_assert(dof_index_cache);

  {
    LOG_SCOPE("L2Kokkos::ensure direct storage", "HilbertSystem");
    KokkosProfilingRegion direct_region("L2Kokkos::ensure direct storage");
    if (!ensure_kokkos_petsc_direct_storage(
          sys, *dof_index_cache, plan, system_matrix, system_rhs))
      libmesh_error_msg("Failed to build the HilbertSystem Kokkos direct PETSc storage path. "
                        "The COO fallback path has been removed; this case must be supported "
                        "directly or fail.");
  }

  assembly_path = ::HilbertSystem::KokkosAssemblyPath::petsc_direct_storage;
  using KokkosVectorWriteView = typename PetscVector<Number>::kokkos_write_view;
  using KokkosMatrixWriteView =
    ::Kokkos::View<Number *,
                   typename ::Kokkos::DefaultExecutionSpace::memory_space,
                   ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;

  {
    LOG_SCOPE("L2Kokkos::write PETSc storage", "HilbertSystem");
    KokkosProfilingRegion write_region("L2Kokkos::write PETSc storage");

    auto rhs_guard = system_rhs.make_kokkos_write_view_guard();
    const KokkosVectorWriteView local_rhs_values = rhs_guard.view();
    KokkosVectorWriteView remote_rhs_values;
    KokkosMatrixWriteView diag_matrix_values;
    KokkosMatrixWriteView offdiag_matrix_values;
    KokkosMatrixWriteView remote_matrix_values;

    {
      KokkosProfilingRegion zero_region("L2Kokkos::zero direct storage values");
      ::Kokkos::deep_copy(local_rhs_values, Number(0));
      if (plan.rhs_remote_values.extent(0))
        {
          ::Kokkos::deep_copy(plan.rhs_remote_values, Number(0));
          remote_rhs_values =
            KokkosVectorWriteView(plan.rhs_remote_values.data(),
                                  plan.rhs_remote_values.extent(0));
        }

      if (plan.mat_coo_values.extent(0))
        ::Kokkos::deep_copy(plan.mat_coo_values, PetscScalar(0));
      if (plan.mat_remote_values.extent(0))
        {
          ::Kokkos::deep_copy(plan.mat_remote_values, Number(0));
          remote_matrix_values =
            KokkosMatrixWriteView(plan.mat_remote_values.data(),
                                  plan.mat_remote_values.extent(0));
        }
    }

    if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::seq_aij)
      diag_matrix_values =
        KokkosMatrixWriteView(plan.mat_coo_values.data(), plan.mat_diag_size);
    else if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::mpi_aij)
      {
        diag_matrix_values =
          KokkosMatrixWriteView(plan.mat_coo_values.data(), plan.mat_diag_size);
        offdiag_matrix_values =
          KokkosMatrixWriteView(plan.mat_coo_values.data() + plan.mat_diag_size,
                                plan.mat_offdiag_size);
      }
    else
      libmesh_error_msg("HilbertSystem Kokkos direct PETSc storage was built without a valid "
                        "PETSc AIJ storage layout.");

    {
      LOG_SCOPE("L2Kokkos::assemble all records", "HilbertSystem");
      KokkosProfilingRegion records_region("L2Kokkos::assemble all records");

      const auto & geometry_cache = sys.get_mesh().get_kokkos_geometry_cache();
      if (analytic_goal)
        {
          const SampledAnalyticGoalAccess goal_access(plan.analytic_goal_qp_offsets,
                                                      plan.analytic_goal_values,
                                                      plan.analytic_goal_gradients);
          ensure_kokkos_entry_value_storage(plan);
          const auto node_coordinates = geometry_cache.node_coordinates_dual.view_device();
          const auto element_node_ids = geometry_cache.element_node_ids_dual.view_device();
          const auto element_mapping_types =
            geometry_cache.element_mapping_types_dual.view_device();
          const auto element_n_nodes = geometry_cache.element_n_nodes_dual.view_device();

          if (!try_kokkos_hilbert_specialized_record_value_batch(plan,
                                                                 node_coordinates,
                                                                 element_node_ids,
                                                                 element_mapping_types,
                                                                 element_n_nodes,
                                                                 plan,
                                                                 sys.hilbert_order(),
                                                                 goal_access))
            libMesh::Kokkos::detail::run_hilbert_system_value_batch<kokkos_hilbert_max_dofs>(
              node_coordinates,
              element_node_ids,
              element_mapping_types,
              element_n_nodes,
              plan.batch_data.shape_keys,
              plan.batch_data.elem_indices,
              plan.batch_data.elem_n_dofs,
              plan.batch_data.quadrature_orders,
              plan.batch_data.rhs_offsets,
              plan.batch_data.mat_offsets,
              sys.hilbert_order(),
              goal_access,
              plan.rhs_entry_values,
              plan.mat_entry_values,
              "hilbert_all_record_values");

          compare_kokkos_hilbert_bucket_value_batch(plan,
                                                    node_coordinates,
                                                    element_node_ids,
                                                    element_mapping_types,
                                                    element_n_nodes,
                                                    plan,
                                                    sys.hilbert_order(),
                                                    goal_access);
        }
      else
        {
          if (!ensure_kokkos_fem_goal_batch_data(sys, *fem_goal, plan))
            return false;

          auto input_guard = plan.fem_goal_batch_data.input_vector->make_kokkos_read_view_guard();
          const auto timed_goal = fem_goal->with_time(sys.time);
          ensure_kokkos_entry_value_storage(plan);
          libMesh::Kokkos::detail::run_hilbert_system_fem_value_batch<
            kokkos_hilbert_max_dofs,
            kokkos_parsed_fem_max_fields>(geometry_cache.node_coordinates_dual.view_device(),
                                          geometry_cache.element_node_ids_dual.view_device(),
                                          geometry_cache.element_mapping_types_dual.view_device(),
                                          geometry_cache.element_n_nodes_dual.view_device(),
                                          plan.batch_data.shape_keys,
                                          plan.batch_data.elem_indices,
                                          plan.batch_data.elem_n_dofs,
                                          plan.batch_data.quadrature_orders,
                                          plan.batch_data.rhs_offsets,
                                          plan.batch_data.mat_offsets,
                                          plan.fem_goal_batch_data.field_keys,
                                          plan.fem_goal_batch_data.field_dofs,
                                          plan.fem_goal_batch_data.field_local_indices,
                                          input_guard.view(),
                                          timed_goal,
                                          sys.hilbert_order(),
                                          plan.rhs_entry_values,
                                          plan.mat_entry_values,
                                          "hilbert_fem_all_record_values");
        }
    }

    {
      LOG_SCOPE("L2Kokkos::scatter entry values", "HilbertSystem");
      KokkosProfilingRegion scatter_region("L2Kokkos::scatter entry values");
      scatter_kokkos_hilbert_entry_values(plan,
                                          plan.mat_diag_size + plan.mat_offdiag_size,
                                          local_rhs_values,
                                          remote_rhs_values,
                                          diag_matrix_values,
                                          offdiag_matrix_values,
                                          remote_matrix_values);
    }

    ::Kokkos::fence();

    if (plan.direct_matrix_layout == KokkosDirectMatrixLayout::mpi_aij)
      {
        LOG_SCOPE("L2Kokkos::apply remote reductions", "HilbertSystem");
        KokkosProfilingRegion remote_reduce_region("L2Kokkos::apply remote reductions");

        if (plan.rhs_remote_sf && plan.rhs_remote_values.extent(0))
          {
            ::Kokkos::deep_copy(plan.rhs_remote_root_values, Number(0));
            LibmeshPetscCall2(sys.comm(),
                              PetscSFReduceWithMemTypeBegin(plan.rhs_remote_sf,
                                                            MPIU_SCALAR,
                                                            kokkos_default_petsc_mem_type(),
                                                            plan.rhs_remote_values.data(),
                                                            kokkos_default_petsc_mem_type(),
                                                            plan.rhs_remote_root_values.data(),
                                                            MPIU_SUM));
            LibmeshPetscCall2(sys.comm(),
                              PetscSFReduceEnd(plan.rhs_remote_sf,
                                               MPIU_SCALAR,
                                               plan.rhs_remote_values.data(),
                                               plan.rhs_remote_root_values.data(),
                                               MPIU_SUM));
            ::Kokkos::parallel_for(
              "hilbert_apply_rhs_remote_reduce",
              ::Kokkos::RangePolicy<>(0, plan.rhs_remote_root_values.extent(0)),
              KOKKOS_LAMBDA(const std::size_t i)
              {
                local_rhs_values(i) += plan.rhs_remote_root_values(i);
              });
          }

        if (plan.mat_remote_sf && plan.mat_remote_values.extent(0))
          {
            ::Kokkos::deep_copy(plan.mat_remote_root_values, Number(0));
            LibmeshPetscCall2(sys.comm(),
                              PetscSFReduceWithMemTypeBegin(plan.mat_remote_sf,
                                                            MPIU_SCALAR,
                                                            kokkos_default_petsc_mem_type(),
                                                            plan.mat_remote_values.data(),
                                                            kokkos_default_petsc_mem_type(),
                                                            plan.mat_remote_root_values.data(),
                                                            MPIU_SUM));
            LibmeshPetscCall2(sys.comm(),
                              PetscSFReduceEnd(plan.mat_remote_sf,
                                               MPIU_SCALAR,
                                               plan.mat_remote_values.data(),
                                               plan.mat_remote_root_values.data(),
                                               MPIU_SUM));
            ::Kokkos::parallel_for(
              "hilbert_apply_mat_remote_reduce",
              ::Kokkos::RangePolicy<>(0, plan.mat_remote_root_values.extent(0)),
              KOKKOS_LAMBDA(const std::size_t i)
              {
                if (i < plan.mat_diag_size)
                  diag_matrix_values(i) += plan.mat_remote_root_values(i);
                else
                  offdiag_matrix_values(i - plan.mat_diag_size) += plan.mat_remote_root_values(i);
              });
          }

        ::Kokkos::fence();
      }

    {
      LOG_SCOPE("L2Kokkos::MatSetValuesCOO", "HilbertSystem");
      KokkosProfilingRegion mat_values_region("PETSc::MatSetValuesCOO");
      LibmeshPetscCall2(sys.comm(),
                        MatSetValuesCOO(system_matrix.mat(),
                                        plan.mat_coo_values.data(),
                                        INSERT_VALUES));
    }
  }

  return true;
}
} // anonymous namespace

} // namespace libMesh

using libMesh::DenseSubMatrix;
using libMesh::DenseSubVector;
using libMesh::Elem;
using libMesh::FEMContext;
using libMesh::FEShapeKey;
using libMesh::KokkosElementAssemblyState;
using libMesh::KokkosPetscAssemblyPlan;
using libMesh::KokkosProfilingRegion;
using libMesh::KokkosScalarView;
using libMesh::LinearSolver;
using libMesh::Number;
using libMesh::ParsedFEMFunction;
using libMesh::ParsedFunction;
using libMesh::PetscMatrixBase;
using libMesh::PetscVector;
using libMesh::Real;
using libMesh::assemble_host_exact_parsed_fem_goal_element;
using libMesh::assemble_kokkos_hilbert_element;
using libMesh::assemble_kokkos_hilbert_fem_goal_element;
using libMesh::assemble_kokkos_petsc_global_system;
using libMesh::build_kokkos_element_state;
using libMesh::build_kokkos_petsc_assembly_plan;
using libMesh::cast_int;
using libMesh::ensure_kokkos_goal_cache;
using libMesh::ensure_kokkos_sampled_analytic_goal_data;
using libMesh::is_petsc_kokkos_vec;
using libMesh::kokkos_petsc_plan_matches;
using libMesh::libmesh_ignore;
using libMesh::make_hilbert_shape_key;
using libMesh::kokkos_hilbert_max_dofs;
using libMesh::kokkos_parsed_fem_max_fields;
namespace detail = libMesh::detail;

struct HilbertSystemKokkosState
{
  std::unique_ptr<libMesh::Kokkos::KokkosParsedFunction<libMesh::Number, 64>> goal_func;
  std::unique_ptr<libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number, 64, 16>> fem_goal_func;
#if defined(LIBMESH_HAVE_PETSC)
  std::unique_ptr<libMesh::KokkosPetscAssemblyPlan> petsc_plan;
#endif
};

void
destroy_hilbert_system_kokkos_state(HilbertSystemKokkosState * state)
{
  delete state;
}

HilbertSystemKokkosState &
HilbertSystem::kokkos_state()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::get Kokkos state");

  if (!this->_kokkos_state)
    this->_kokkos_state.reset(new HilbertSystemKokkosState);

  return *this->_kokkos_state;
}

HilbertSystemKokkosState *
HilbertSystem::maybe_kokkos_state()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::maybe Kokkos state");

  return this->_kokkos_state.get();
}

const HilbertSystemKokkosState *
HilbertSystem::maybe_kokkos_state() const
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::maybe const Kokkos state");

  return this->_kokkos_state.get();
}

const libMesh::Kokkos::KokkosParsedFunction<libMesh::Number> *
HilbertSystem::ensure_kokkos_goal_func()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure analytic goal function");

  auto & state = this->kokkos_state();
  return ensure_kokkos_goal_cache<libMesh::ParsedFunction<libMesh::Number>>(
    state.goal_func,
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

bool
HilbertSystem::needs_exact_kokkos_fem_goal_context()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::check exact FEM goal context");

  return input_system &&
         dynamic_cast<ParsedFEMFunction<Number> *>(_goal_func.get()) &&
         this->ensure_kokkos_fem_goal_func();
}

bool
HilbertSystem::try_exact_kokkos_analytic_goal_host_assembly(FEMContext & c,
                                                            const bool request_jacobian,
                                                            DenseSubVector<Number> & F,
                                                            DenseSubMatrix<Number> & K)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::try exact analytic host assembly");

  if (const auto * kokkos_goal = this->ensure_kokkos_goal_func())
    {
      detail::HostHilbertFEAccess fe(c, 0, _hilbert_order);
      auto solution =
        detail::make_hilbert_solution_access(fe,
                                             c.get_elem_solution(0),
                                             c.get_elem_solution_derivative());
      detail::HostHilbertAccumulator accum(F, K);
      const auto parsed_goal = kokkos_goal->with_time(this->time);
      auto goal = detail::make_hilbert_analytic_goal_access(parsed_goal,
                                                            parsed_goal.gradient_function());
      detail::assemble_hilbert_element(fe,
                                       solution,
                                       goal,
                                       request_jacobian,
                                       _hilbert_order,
                                       accum);
      return true;
    }

  return false;
}

const libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number> *
HilbertSystem::ensure_kokkos_fem_goal_func()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure FEM goal function");

  auto & state = this->kokkos_state();
  return ensure_kokkos_goal_cache<libMesh::ParsedFEMFunction<libMesh::Number>>(
    state.fem_goal_func,
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

bool
HilbertSystem::try_exact_kokkos_fem_goal_host_assembly(FEMContext & c,
                                                       const bool request_jacobian,
                                                       DenseSubVector<Number> & F,
                                                       DenseSubMatrix<Number> & K)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::try exact FEM host assembly");

#if defined(LIBMESH_HAVE_PETSC)
  if (!input_system)
    return false;

  if (const auto * kokkos_fem_goal = this->ensure_kokkos_fem_goal_func())
    return assemble_host_exact_parsed_fem_goal_element(*this,
                                                       c,
                                                       request_jacobian,
                                                       F,
                                                       K,
                                                       kokkos_fem_goal->with_time(this->time));
#else
  libmesh_ignore(c, request_jacobian, F, K);
#endif

  return false;
}

void
HilbertSystem::reset_kokkos_goal_cache()
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::reset goal cache");

  if (auto * state = this->maybe_kokkos_state())
    {
      state->goal_func.reset();
      state->fem_goal_func.reset();
    }
}

KokkosPetscAssemblyPlan *
HilbertSystem::ensure_kokkos_petsc_plan(bool * rebuilt)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::ensure PETSc assembly plan");

  auto & state = this->kokkos_state();
  const auto & geometry_cache = this->get_mesh().get_kokkos_geometry_cache();
  const auto * dof_index_cache = this->get_dof_map().get_kokkos_dof_index_cache(0);
  if (!dof_index_cache)
    return nullptr;

  if (state.petsc_plan &&
      kokkos_petsc_plan_matches(*this, geometry_cache, *dof_index_cache, *state.petsc_plan))
    {
      if (rebuilt)
        *rebuilt = false;
      return state.petsc_plan.get();
    }

  auto plan = std::make_unique<KokkosPetscAssemblyPlan>();
  if (!build_kokkos_petsc_assembly_plan(*this, *plan))
    return nullptr;

  state.petsc_plan = std::move(plan);
  if (rebuilt)
    *rebuilt = true;
  return state.petsc_plan.get();
}

bool
HilbertSystem::try_kokkos_element_assembly(FEMContext & c,
                                           const bool request_jacobian,
                                           DenseSubVector<Number> & F,
                                           DenseSubMatrix<Number> & K)
{
  KokkosProfilingRegion kokkos_region("L2Kokkos::try element assembly");

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
  LOG_SCOPE("L2Kokkos::solve", "HilbertSystem");
  KokkosProfilingRegion kokkos_region("L2Kokkos::solve");

  using clock = std::chrono::steady_clock;
  auto * petsc_matrix = dynamic_cast<PetscMatrixBase<Number> *>(this->matrix);
  auto * petsc_rhs = dynamic_cast<PetscVector<Number> *>(this->rhs);
  auto * petsc_solution = dynamic_cast<PetscVector<Number> *>(this->solution.get());

  const bool analytic_goal = static_cast<bool>(_analytic_goal_func);
  const auto * fem_goal = this->ensure_kokkos_fem_goal_func();

  libmesh_error_msg_if(!petsc_matrix || !petsc_rhs || !petsc_solution,
                       "HilbertSystem Kokkos direct PETSc storage requires PETSc-backed matrix, "
                       "RHS, and solution objects.");
  libmesh_error_msg_if(!(analytic_goal || fem_goal),
                       "HilbertSystem Kokkos direct PETSc storage requires an analytic goal "
                       "or a parsed FEM goal with Kokkos support.");
  libmesh_error_msg_if(!is_petsc_kokkos_vec(this->comm(), petsc_rhs->vec()) ||
                         !is_petsc_kokkos_vec(this->comm(), petsc_solution->vec()),
                       "HilbertSystem Kokkos direct PETSc storage requires PETSc Kokkos vectors. "
                       "Run with -vec_type kokkos.");

  this->_last_kokkos_timing = {};
  const auto total_start = clock::now();
  {
    LOG_SCOPE("L2Kokkos::zero PETSc objects", "HilbertSystem");
    KokkosProfilingRegion zero_region("L2Kokkos::zero PETSc objects");
    petsc_matrix->zero();
    petsc_rhs->zero();
    petsc_solution->zero();
  }

  bool rebuilt_plan = false;
  const auto plan_start = clock::now();
  auto * plan = [&]()
  {
    LOG_SCOPE("L2Kokkos::ensure assembly plan", "HilbertSystem");
    KokkosProfilingRegion plan_region("L2Kokkos::ensure assembly plan");
    return this->ensure_kokkos_petsc_plan(&rebuilt_plan);
  }();
  libmesh_error_msg_if(!plan,
                       "HilbertSystem Kokkos direct PETSc storage could not build a supported "
                       "assembly plan for the current FE/mapping/quadrature configuration.");
  if (analytic_goal)
    {
      LOG_SCOPE("L2Kokkos::prepare analytic goal samples", "HilbertSystem");
      KokkosProfilingRegion analytic_region("L2Kokkos::prepare analytic goal samples");
      if (!ensure_kokkos_sampled_analytic_goal_data(*this,
                                                    *_analytic_goal_func,
                                                    _analytic_goal_grad.get(),
                                                    *plan))
        libmesh_error_msg("HilbertSystem Kokkos direct PETSc storage could not sample the "
                          "analytic goal on the host for device assembly.");
    }
  const auto plan_stop = clock::now();
  this->_last_kokkos_timing.plan_seconds =
    rebuilt_plan ?
      std::chrono::duration_cast<std::chrono::duration<Real>>(plan_stop - plan_start).count() :
      0.;

  const auto assembly_start = clock::now();
  auto assembly_path = HilbertSystem::KokkosAssemblyPath::none;
  {
    LOG_SCOPE("L2Kokkos::assemble global system", "HilbertSystem");
    KokkosProfilingRegion assembly_region("L2Kokkos::assemble global system");
    if (!assemble_kokkos_petsc_global_system(*this,
                                             *plan,
                                             analytic_goal,
                                             nullptr,
                                             fem_goal,
                                             *petsc_matrix,
                                             *petsc_rhs,
                                             assembly_path))
      libmesh_error_msg("HilbertSystem Kokkos direct PETSc storage assembly failed.");
  }
  const auto assembly_stop = clock::now();
  this->_last_kokkos_timing.assembly_seconds =
    std::chrono::duration_cast<std::chrono::duration<Real>>(assembly_stop - assembly_start).count();
  this->_last_kokkos_timing.assembly_path = assembly_path;

  {
    LOG_SCOPE("L2Kokkos::close PETSc objects", "HilbertSystem");
    KokkosProfilingRegion close_region("L2Kokkos::close PETSc objects");
    petsc_matrix->close();
    petsc_rhs->close();
    petsc_solution->close();
  }

  LinearSolver<Number> * solver = this->get_linear_solver();
  {
    LOG_SCOPE("L2Kokkos::init PETSc solver", "HilbertSystem");
    KokkosProfilingRegion init_region("L2Kokkos::init PETSc solver");
    if (this->prefix_with_name())
      solver->init(this->prefix().c_str());
    else
      solver->init();
  }

  const auto [maxlinearits, linear_tol] = this->get_linear_solve_parameters();
  const auto solve_start = clock::now();
  {
    LOG_SCOPE("L2Kokkos::PETSc solve", "HilbertSystem");
    KokkosProfilingRegion solve_region("PETSc::KSPSolve");
    solver->solve(*this->matrix,
                  *this->solution,
                  *this->rhs,
                  linear_tol,
                  maxlinearits);
  }
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

#endif // LIBMESH_HAVE_KOKKOS && LIBMESH_HAVE_PETSC && !LIBMESH_USE_COMPLEX_NUMBERS
