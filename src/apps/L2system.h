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

#ifndef LIBMESH_APPS_L2SYSTEM_H
#define LIBMESH_APPS_L2SYSTEM_H

// libMesh includes
#include "libmesh/enum_fe_family.h"
#include "libmesh/fdm_gradient.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"
#include "libmesh/function_base.h"
#include "libmesh/libmesh_common.h"

// C++ includes
#include <map>
#include <memory>
#include <set>
#include <string>

namespace libMesh
{
namespace Kokkos
{
  template <typename Output, unsigned int MaxStack>
  class KokkosParsedFunction;

  template <typename Output,
            unsigned int MaxStack,
            unsigned int MaxFieldVariables>
  class KokkosParsedFEMFunction;
}

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
struct KokkosPetscAssemblyPlan;
#endif
}

struct HilbertSystemKokkosState;

void libmesh_kokkos_initialize(int & argc,
                               char ** & argv,
                               bool enable);

void libmesh_kokkos_finalize(bool enable);

void destroy_hilbert_system_kokkos_state(HilbertSystemKokkosState *);

// FEMSystem, TimeSolver and NewtonSolver will handle most tasks,
// but we must specify element residuals.
class HilbertSystem : public libMesh::FEMSystem
{
public:
  enum class KokkosAssemblyPath
  {
    none,
    petsc_direct_storage
  };

  struct KokkosTimingInfo
  {
    libMesh::Real plan_seconds = 0.;
    libMesh::Real assembly_seconds = 0.;
    libMesh::Real solve_seconds = 0.;
    libMesh::Real total_seconds = 0.;
    KokkosAssemblyPath assembly_path = KokkosAssemblyPath::none;
  };

  using kokkos_state_ptr =
    std::unique_ptr<HilbertSystemKokkosState,
                    void (*)(HilbertSystemKokkosState *)>;

  HilbertSystem(libMesh::EquationSystems & es,
                const std::string & name,
                const unsigned int number);

  ~HilbertSystem() override;

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  std::set<libMesh::subdomain_id_type> & subdomains_list() { return _subdomains_list; }
  const std::set<libMesh::subdomain_id_type> & subdomains_list() const { return _subdomains_list; }

  unsigned int & hilbert_order() { return _hilbert_order; }
  unsigned int hilbert_order() const { return _hilbert_order; }

  void use_kokkos_backend(const bool use) { _use_kokkos_backend = use; }
  bool use_kokkos_backend() const { return _use_kokkos_backend; }

  void use_exact_parsed_fem_host_path(const bool use) { _use_exact_parsed_fem_host_path = use; }
  bool use_exact_parsed_fem_host_path() const { return _use_exact_parsed_fem_host_path; }

  const KokkosTimingInfo & last_kokkos_timing() const { return _last_kokkos_timing; }

  void solve() override;

  void set_fdm_eps(const libMesh::Real eps)
  {
    _fdm_eps = eps;
    rebuild_goal_gradient();
    rebuild_analytic_goal_gradient();
  }

  void set_goal_func(libMesh::FEMFunctionBase<libMesh::Number> & goal)
  {
    _goal_func = goal.clone();
    _analytic_goal_func.reset();
    _analytic_goal_grad.reset();
    reset_kokkos_goal_cache();
    rebuild_goal_gradient();
  }

  void set_goal_func(libMesh::FunctionBase<libMesh::Number> & goal)
  {
    _analytic_goal_func = goal.clone();
    _analytic_goal_func->init();
    _goal_func.reset();
    _goal_grad.reset();
    reset_kokkos_goal_cache();
    rebuild_analytic_goal_gradient();
  }

  // We want to be able to project functions based on *other* systems'
  // values. For that we need not only a FEMFunction but also a
  // reference to the system where it applies and a separate context
  // object (or multiple separate context objects, in the threaded
  // case) for that system.
  libMesh::System * input_system;

  libMesh::FEMContext * get_input_context(libMesh::FEMContext & c)
  {
    const auto it = input_contexts.find(&c);
    return (it == input_contexts.end()) ? nullptr : it->second.get();
  }

protected:
  std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number>> _goal_func;
  std::unique_ptr<libMesh::FunctionBase<libMesh::Number>> _analytic_goal_func;

  std::map<libMesh::FEMContext *, std::unique_ptr<libMesh::FEMContext>> input_contexts;

  void init_data() override;
  void init_context(libMesh::DiffContext & context) override;
  bool element_time_derivative(bool request_jacobian,
                               libMesh::DiffContext & context) override;

  std::string _fe_family;
  unsigned int _fe_order;
  unsigned int _hilbert_order;

  bool _use_kokkos_backend;
  bool _use_exact_parsed_fem_host_path;

  std::unique_ptr<libMesh::FDMGradient<libMesh::Gradient>> _goal_grad;
  std::unique_ptr<libMesh::FunctionBase<libMesh::Gradient>> _analytic_goal_grad;
  kokkos_state_ptr _kokkos_state;

  libMesh::Real _fdm_eps;
  std::set<libMesh::subdomain_id_type> _subdomains_list;
  KokkosTimingInfo _last_kokkos_timing;

private:
#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  HilbertSystemKokkosState & kokkos_state();
  HilbertSystemKokkosState * maybe_kokkos_state();
  const HilbertSystemKokkosState * maybe_kokkos_state() const;

  const libMesh::Kokkos::KokkosParsedFunction<libMesh::Number, 64> * ensure_kokkos_goal_func();
  const libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number, 64, 16> *
  ensure_kokkos_fem_goal_func();
  void reset_kokkos_goal_cache();
  bool needs_exact_kokkos_fem_goal_context();
  bool try_exact_kokkos_analytic_goal_host_assembly(libMesh::FEMContext & c,
                                                    bool request_jacobian,
                                                    libMesh::DenseSubVector<libMesh::Number> & F,
                                                    libMesh::DenseSubMatrix<libMesh::Number> & K);
  bool try_exact_kokkos_fem_goal_host_assembly(libMesh::FEMContext & c,
                                               bool request_jacobian,
                                               libMesh::DenseSubVector<libMesh::Number> & F,
                                               libMesh::DenseSubMatrix<libMesh::Number> & K);
  bool try_kokkos_element_assembly(libMesh::FEMContext & c,
                                   bool request_jacobian,
                                   libMesh::DenseSubVector<libMesh::Number> & F,
                                   libMesh::DenseSubMatrix<libMesh::Number> & K);
  libMesh::KokkosPetscAssemblyPlan * ensure_kokkos_petsc_plan(bool * rebuilt = nullptr);
  bool try_kokkos_petsc_solve();
#else
  void reset_kokkos_goal_cache() {}
  bool needs_exact_kokkos_fem_goal_context() { return false; }
#endif

  void rebuild_goal_gradient();
  void rebuild_analytic_goal_gradient();
};

#endif // LIBMESH_APPS_L2SYSTEM_H
