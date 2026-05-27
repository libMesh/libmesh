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
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fem_context.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/parsed_function.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

using namespace libMesh;
namespace detail = libMesh::detail;

#if !defined(LIBMESH_HAVE_KOKKOS) || !defined(LIBMESH_HAVE_PETSC) || defined(LIBMESH_USE_COMPLEX_NUMBERS)
void
destroy_hilbert_system_kokkos_state(HilbertSystemKokkosState * state)
{
  libmesh_ignore(state);
}

#endif

HilbertSystem::HilbertSystem(libMesh::EquationSystems & es,
                             const std::string & name,
                             const unsigned int number)
  : libMesh::FEMSystem(es, name, number),
    input_system(nullptr),
    _fe_family("LAGRANGE"),
    _fe_order(1),
    _hilbert_order(0),
    _use_kokkos_backend(false),
    _use_exact_parsed_fem_host_path(false),
    _kokkos_state(nullptr, destroy_hilbert_system_kokkos_state),
    _fdm_eps(libMesh::TOLERANCE),
    _subdomains_list()
{
}

HilbertSystem::~HilbertSystem() = default;

void HilbertSystem::rebuild_goal_gradient()
{
  if (_goal_func)
    _goal_grad = std::make_unique<FDMGradient<Gradient>>(*_goal_func, _fdm_eps);
  else
    _goal_grad.reset();
}

void HilbertSystem::rebuild_analytic_goal_gradient()
{
  if (_analytic_goal_func)
    {
      _analytic_goal_grad =
        std::make_unique<detail::FunctionFDMGradient<Gradient>>(*_analytic_goal_func, _fdm_eps);
      _analytic_goal_grad->init();
    }
  else
    _analytic_goal_grad.reset();
}

void HilbertSystem::init_data()
{
  this->get_dof_map().full_sparsity_pattern_needed();
  this->add_variable("u",
                     static_cast<Order>(_fe_order),
                     Utility::string_to_enum<FEFamily>(_fe_family));

  // Do the parent's initialization after variables are defined.
  FEMSystem::init_data();
}

void HilbertSystem::init_context(DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  FEBase * my_fe = nullptr;

  // We might have a multi-dimensional mesh.
  const std::set<unsigned char> & elem_dims = c.elem_dimensions();

  for (const auto & dim : elem_dims)
    {
      c.get_element_fe(0, my_fe, dim);

      my_fe->get_JxW();
      my_fe->get_phi();
      my_fe->get_xyz();

      if (this->_hilbert_order > 0)
        my_fe->get_dphi();

      c.get_side_fe(0, my_fe, dim);
      my_fe->get_nothing();
    }

  // Build a corresponding context for the input system if we haven't already.
  auto & input_context = input_contexts[&c];
  if (input_system && !input_context)
    input_context = std::make_unique<FEMContext>(*input_system);

  libmesh_assert(_goal_func || _analytic_goal_func);

  if (_goal_func)
    _goal_func->init_context(input_system ? *input_context : c);

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (input_system &&
      this->_hilbert_order > 0 &&
      this->needs_exact_kokkos_fem_goal_context())
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

bool HilbertSystem::element_time_derivative(const bool request_jacobian,
                                            DiffContext & context)
{
  FEMContext & c = cast_ref<FEMContext &>(context);

  const Elem & elem = c.get_elem();

  if (!_subdomains_list.empty() &&
      !_subdomains_list.count(elem.subdomain_id()))
    return request_jacobian;

  DenseSubMatrix<Number> & K = c.get_elem_jacobian(0, 0);
  DenseSubVector<Number> & F = c.get_elem_residual(0);

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC)
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

#if defined(LIBMESH_HAVE_KOKKOS) && defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (this->try_exact_kokkos_analytic_goal_host_assembly(c, request_jacobian, F, K))
    return request_jacobian;

  if (_use_exact_parsed_fem_host_path && input_system)
    if (this->try_exact_kokkos_fem_goal_host_assembly(c, request_jacobian, F, K))
      return request_jacobian;
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
  if (_use_kokkos_backend)
    {
      if (this->try_kokkos_petsc_solve())
        return;

      libmesh_error_msg("HilbertSystem Kokkos backend did not complete the direct PETSc "
                        "storage solve path.");
    }
#else
  libmesh_error_msg_if(_use_kokkos_backend,
                       "HilbertSystem Kokkos backend requires a libMesh build with Kokkos, "
                       "PETSc, and real Number support.");
#endif

  FEMSystem::solve();
}
