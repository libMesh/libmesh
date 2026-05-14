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

// libMesh includes
#include "../../include/systems/hilbert_assembly.h"
#include "../../include/numerics/parsed_fem_function.h"
#include "../../include/numerics/parsed_function.h"

#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
#include "../../include/gpu/kokkos_parsed_function.h"
#endif

#include "libmesh/enum_fe_family.h"
#include "libmesh/fdm_gradient.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fem_system.h"
#include "libmesh/function_base.h"
#include "libmesh/libmesh_common.h"

// C++ includes
#include <map>
#include <memory>


// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
struct KokkosPetscAssemblyPlan;
#endif

class HilbertSystem : public libMesh::FEMSystem
{
public:
  struct KokkosTimingInfo
  {
    libMesh::Real plan_seconds = 0.;
    libMesh::Real assembly_seconds = 0.;
    libMesh::Real solve_seconds = 0.;
    libMesh::Real total_seconds = 0.;
  };

  // Constructor
  HilbertSystem(libMesh::EquationSystems & es,
                const std::string & name,
                const unsigned int number);

  // Default destructor
  ~HilbertSystem();

  std::string & fe_family() { return _fe_family; }
  unsigned int & fe_order() { return _fe_order; }
  std::set<libMesh::subdomain_id_type> & subdomains_list() { return _subdomains_list; }
  const std::set<libMesh::subdomain_id_type> & subdomains_list() const { return _subdomains_list; }

  unsigned int & hilbert_order() { return _hilbert_order; }
  unsigned int hilbert_order() const { return _hilbert_order; }
  void use_kokkos_backend(bool use) { _use_kokkos_backend = use; }
  bool use_kokkos_backend() const { return _use_kokkos_backend; }
  const KokkosTimingInfo & last_kokkos_timing() const { return _last_kokkos_timing; }
  virtual void solve () override;
  void set_fdm_eps(libMesh::Real eps) {
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
  // values.  For that we need not only a FEMFunction but also a
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
  std::unique_ptr<libMesh::FEMFunctionBase<libMesh::Number> > _goal_func;
  std::unique_ptr<libMesh::FunctionBase<libMesh::Number>> _analytic_goal_func;

  std::map<libMesh::FEMContext *, std::unique_ptr<libMesh::FEMContext>>
    input_contexts;

  // System initialization
  virtual void init_data () override;

  // Context initialization
  virtual void init_context (libMesh::DiffContext & context) override;

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        libMesh::DiffContext & context) override;

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // The Hilbert order our subclass will project with
  unsigned int _hilbert_order;

  bool _use_kokkos_backend;

  // The function we will call to finite difference our goal
  // function
  std::unique_ptr<libMesh::FDMGradient<libMesh::Gradient>> _goal_grad;
  std::unique_ptr<libMesh::FunctionBase<libMesh::Gradient>> _analytic_goal_grad;
#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  std::unique_ptr<libMesh::Kokkos::KokkosParsedFunction<libMesh::Number>> _kokkos_goal_func;
  std::unique_ptr<libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number>> _kokkos_fem_goal_func;
#if defined(LIBMESH_HAVE_PETSC)
  std::unique_ptr<KokkosPetscAssemblyPlan> _kokkos_petsc_plan;
#endif
#endif

  // The perturbation we will use when finite differencing our goal
  // function
  libMesh::Real _fdm_eps;

  // Which subdomains to integrate on (all subdomains, if empty())
  std::set<libMesh::subdomain_id_type> _subdomains_list;
  KokkosTimingInfo _last_kokkos_timing;

private:
#if defined(LIBMESH_HAVE_KOKKOS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  const libMesh::Kokkos::KokkosParsedFunction<libMesh::Number> * ensure_kokkos_goal_func();
  const libMesh::Kokkos::KokkosParsedFEMFunction<libMesh::Number> * ensure_kokkos_fem_goal_func();
  bool try_kokkos_element_assembly(libMesh::FEMContext & c,
                                   bool request_jacobian,
                                   libMesh::DenseSubVector<libMesh::Number> & F,
                                   libMesh::DenseSubMatrix<libMesh::Number> & K);
#if defined(LIBMESH_HAVE_PETSC)
  bool try_kokkos_petsc_solve();
  KokkosPetscAssemblyPlan * ensure_kokkos_petsc_plan(bool * rebuilt = nullptr);
#endif

  void reset_kokkos_goal_cache();
#else
  void reset_kokkos_goal_cache() {}
#endif

  void rebuild_goal_gradient()
  {
    if (_goal_func)
      _goal_grad = std::make_unique<libMesh::FDMGradient<libMesh::Gradient>>(*_goal_func, _fdm_eps);
    else
      _goal_grad.reset();
  }

  void rebuild_analytic_goal_gradient()
  {
    if (_analytic_goal_func)
      {
        _analytic_goal_grad =
          std::make_unique<libMesh::detail::FunctionFDMGradient<libMesh::Gradient>>(*_analytic_goal_func,
                                                                                     _fdm_eps);
        _analytic_goal_grad->init();
      }
    else
      _analytic_goal_grad.reset();
  }
};
