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

#ifndef LIBMESH_HILBERT_ASSEMBLY_KERNEL_H
#define LIBMESH_HILBERT_ASSEMBLY_KERNEL_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"
#include "libmesh/vector_value.h"

#include <type_traits>
#include <utility>

namespace libMesh
{
namespace detail
{

template <typename T>
using hilbert_storage_t =
  std::conditional_t<std::is_lvalue_reference_v<T>, T, std::decay_t<T>>;

template <typename T, typename = void>
struct is_zero_hilbert_coeff_access : std::false_type {};

template <typename T>
struct is_zero_hilbert_coeff_access<
  T,
  std::void_t<decltype(std::decay_t<T>::is_zero_coeff_access)>>
  : std::bool_constant<std::decay_t<T>::is_zero_coeff_access> {};

template <typename CoeffStorage,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<CoeffStorage>> &&
                             !std::is_array_v<std::remove_reference_t<CoeffStorage>>,
                           int> = 0>
LIBMESH_DEVICE_INLINE decltype(auto)
coeff_at(const CoeffStorage & coeff, const unsigned int i)
{
  return coeff(i);
}

template <typename Scalar>
LIBMESH_DEVICE_INLINE const Scalar &
coeff_at(const Scalar * coeff, const unsigned int i)
{
  return coeff[i];
}

template <typename Scalar, std::size_t N>
LIBMESH_DEVICE_INLINE const Scalar &
coeff_at(const Scalar (&coeff)[N], const unsigned int i)
{
  libmesh_ignore(N);
  return coeff[i];
}

template <typename QpData, typename CoeffStorage>
LIBMESH_DEVICE_INLINE
Number interpolate_hilbert_value(const QpData & qp_data,
                                 const CoeffStorage & coeff,
                                 const unsigned int n_dofs)
{
  if constexpr (is_zero_hilbert_coeff_access<CoeffStorage>::value)
    return Number(0);

  Number u = 0.;

  for (unsigned int i = 0; i != n_dofs; ++i)
    u += coeff_at(coeff, i) * qp_data.phi(i);

  return u;
}

template <typename QpData, typename CoeffStorage>
LIBMESH_DEVICE_INLINE
Gradient interpolate_hilbert_gradient(const QpData & qp_data,
                                      const CoeffStorage & coeff,
                                      const unsigned int n_dofs)
{
  Gradient grad_u;
  grad_u.zero();

  if constexpr (is_zero_hilbert_coeff_access<CoeffStorage>::value)
    return grad_u;

  for (unsigned int i = 0; i != n_dofs; ++i)
    grad_u.add_scaled(qp_data.dphi(i), coeff_at(coeff, i));

  return grad_u;
}

template <typename FEAccess, typename CoeffStorage>
class HilbertSolutionAccess
{
public:
  LIBMESH_DEVICE_INLINE
  HilbertSolutionAccess(const FEAccess & fe,
                        CoeffStorage coeff,
                        const Number solution_derivative)
    : _fe(fe),
      _coeff(coeff),
      _solution_derivative(solution_derivative)
  {
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Number value(const QpData & qp_data) const
  {
    return interpolate_hilbert_value(qp_data, _coeff, _fe.n_dofs());
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Gradient gradient(const QpData & qp_data) const
  {
    return interpolate_hilbert_gradient(qp_data, _coeff, _fe.n_dofs());
  }

  LIBMESH_DEVICE_INLINE
  Number solution_derivative() const
  {
    return _solution_derivative;
  }

private:
  const FEAccess & _fe;
  CoeffStorage _coeff;
  Number _solution_derivative;
};

template <typename FEAccess, typename CoeffStorage>
LIBMESH_DEVICE_INLINE auto
make_hilbert_solution_access(const FEAccess & fe,
                             CoeffStorage && coeff,
                             const Number solution_derivative)
{
  return HilbertSolutionAccess<FEAccess, hilbert_storage_t<CoeffStorage &&>>(
    fe,
    std::forward<CoeffStorage>(coeff),
    solution_derivative);
}

template <typename GoalFunction, typename GoalGradient>
class HilbertAnalyticGoalAccess
{
public:
  LIBMESH_DEVICE_INLINE
  HilbertAnalyticGoalAccess(GoalFunction goal_func,
                            GoalGradient goal_grad)
    : _goal_func(goal_func),
      _goal_grad(goal_grad)
  {
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Number value(const QpData &, const Point & p) const
  {
    return _goal_func(p);
  }

  template <typename QpData>
  LIBMESH_DEVICE_INLINE
  Gradient gradient(const QpData &, const Point & p) const
  {
    return _goal_grad(p);
  }

private:
  GoalFunction _goal_func;
  GoalGradient _goal_grad;
};

template <typename GoalFunction, typename GoalGradient>
LIBMESH_DEVICE_INLINE auto
make_hilbert_analytic_goal_access(GoalFunction && goal_func,
                                  GoalGradient && goal_grad)
{
  return HilbertAnalyticGoalAccess<hilbert_storage_t<GoalFunction &&>,
                                   hilbert_storage_t<GoalGradient &&>>(
    std::forward<GoalFunction>(goal_func),
    std::forward<GoalGradient>(goal_grad));
}

template <typename FEAccess,
          typename SolutionAccess,
          typename GoalAccess,
          typename Accumulator>
LIBMESH_DEVICE_INLINE void
assemble_hilbert_element(const FEAccess & fe,
                         const SolutionAccess & solution,
                         GoalAccess & goal,
                         const bool request_jacobian,
                         const unsigned int hilbert_order,
                         const unsigned int n_u_dofs,
                         Accumulator & accum)
{
  const unsigned int n_qpoints = fe.n_qpoints();

  for (unsigned int qp = 0; qp != n_qpoints; qp++)
    {
      const auto qp_data = fe.qp_data(qp, hilbert_order > 0);
      const Point & xyz = qp_data.xyz();
      const Number err_u = solution.value(qp_data) - goal.value(qp_data, xyz);

      for (unsigned int i = 0; i != n_u_dofs; i++)
        accum.add_residual(i, qp_data.JxW() * (err_u * qp_data.phi(i)));

      if (hilbert_order > 0)
        {
          const Gradient err_grad_u =
            solution.gradient(qp_data) - goal.gradient(qp_data, xyz);

          for (unsigned int i = 0; i != n_u_dofs; i++)
            accum.add_residual(i, qp_data.JxW() * (err_grad_u * qp_data.dphi(i)));
        }

      if (request_jacobian)
        {
          const Number JxWxD = qp_data.JxW() * solution.solution_derivative();

          for (unsigned int i = 0; i != n_u_dofs; i++)
            for (unsigned int j = 0; j != n_u_dofs; ++j)
              accum.add_jacobian(i, j, JxWxD * (qp_data.phi(i) * qp_data.phi(j)));

          if (hilbert_order > 0)
            for (unsigned int i = 0; i != n_u_dofs; i++)
              for (unsigned int j = 0; j != n_u_dofs; ++j)
                accum.add_jacobian(i, j, JxWxD * (qp_data.dphi(i) * qp_data.dphi(j)));
        }
    }
}

template <typename FEAccess,
          typename SolutionAccess,
          typename GoalAccess,
          typename Accumulator>
LIBMESH_DEVICE_INLINE void
assemble_hilbert_element(const FEAccess & fe,
                         const SolutionAccess & solution,
                         GoalAccess & goal,
                         const bool request_jacobian,
                         const unsigned int hilbert_order,
                         Accumulator & accum)
{
  assemble_hilbert_element(fe,
                           solution,
                           goal,
                           request_jacobian,
                           hilbert_order,
                           accum.n_dofs(),
                           accum);
}

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_HILBERT_ASSEMBLY_KERNEL_H
