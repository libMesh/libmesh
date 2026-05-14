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

#ifndef LIBMESH_KOKKOS_HILBERT_ASSEMBLY_H
#define LIBMESH_KOKKOS_HILBERT_ASSEMBLY_H

#include "../systems/hilbert_assembly_kernel.h"

#include "kokkos_fe_evaluator.h"
#include "kokkos_fe_map.h"
#include "kokkos_quadrature.h"

#include <type_traits>
#include <utility>

namespace libMesh::Kokkos::detail
{

LIBMESH_DEVICE_INLINE bool
supports_hilbert_local_assembly(libMesh::FEShapeKey key,
                                libMesh::ElemMappingType mapping_type,
                                const unsigned int quadrature_order)
{
  return mapping_type == libMesh::LAGRANGE_MAP &&
         libMesh::supports_shape_with_lagrange_map(key) &&
         libMesh::supports_grad_shape(key) &&
         libMesh::Kokkos::GaussQuadrature::n_points(key.elem_type, quadrature_order) > 0;
}

template <typename NodeStorage>
class HilbertFEAccess
{
public:
  using node_storage_type = std::decay_t<NodeStorage>;
  class QpData
  {
  public:
    LIBMESH_DEVICE_INLINE
    QpData(const HilbertFEAccess & fe,
           const unsigned int qp,
           const bool need_gradients)
      : _fe(fe),
        _qp(qp),
        _qp_ref(GaussQuadrature::point(fe._key.elem_type, fe._quadrature_order, qp)),
        _JxW(0.),
        _need_gradients(need_gradients)
    {
      RealVector xyz = zero_vector();
      RealTensor J = zero_tensor();

      physical_point_and_jacobian(fe._mapping_type,
                                  fe._key.elem_type,
                                  fe._nodes,
                                  fe._n_nodes,
                                  _qp_ref(0),
                                  _qp_ref(1),
                                  _qp_ref(2),
                                  xyz,
                                  J);

      _xyz = Point(xyz(0), xyz(1), xyz(2));
      _JxW =
        volume_jxw(J,
                   fe._dim,
                   GaussQuadrature::weight(fe._key.elem_type, fe._quadrature_order, qp));

      if (_need_gradients)
        _Jinv = libMesh::Kokkos::inverse<libMesh::Kokkos::RealTensor>(J, fe._dim);
    }

    LIBMESH_DEVICE_INLINE
    Real JxW() const
    {
      return _JxW;
    }

    LIBMESH_DEVICE_INLINE
    Real phi(const unsigned int i) const
    {
      return shape(_fe._key, i, _qp_ref(0), _qp_ref(1), _qp_ref(2));
    }

    LIBMESH_DEVICE_INLINE
    Gradient dphi(const unsigned int i) const
    {
      libmesh_assert(_need_gradients);
      return _Jinv * grad_shape(_fe._key, i, _qp_ref(0), _qp_ref(1), _qp_ref(2));
    }

    LIBMESH_DEVICE_INLINE
    const Point & xyz() const
    {
      return _xyz;
    }

    LIBMESH_DEVICE_INLINE
    const RealVector & reference_point() const
    {
      return _qp_ref;
    }

    LIBMESH_DEVICE_INLINE
    const RealTensor & inverse_jacobian() const
    {
      libmesh_assert(_need_gradients);
      return _Jinv;
    }

    LIBMESH_DEVICE_INLINE
    unsigned int qp_index() const
    {
      return _qp;
    }

    LIBMESH_DEVICE_INLINE
    unsigned int elem_index() const
    {
      return _fe._elem_index;
    }

  private:
    const HilbertFEAccess & _fe;
    unsigned int _qp;
    RealVector _qp_ref;
    Point _xyz;
    Real _JxW;
    RealTensor _Jinv;
    bool _need_gradients;
  };

  LIBMESH_DEVICE_INLINE
  HilbertFEAccess(libMesh::FEShapeKey key,
                  libMesh::ElemMappingType mapping_type,
                  const NodeStorage & nodes,
                  const unsigned int n_nodes,
                  const unsigned int quadrature_order,
                  const unsigned int elem_index = 0)
    : _key(key),
      _mapping_type(mapping_type),
      _nodes(nodes),
      _n_nodes(n_nodes),
      _quadrature_order(quadrature_order),
      _dim(dim_from_topology(key.elem_type)),
      _elem_index(elem_index)
  {
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_qpoints() const
  {
    return GaussQuadrature::n_points(_key.elem_type, _quadrature_order);
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_dofs() const
  {
    return libMesh::Kokkos::n_dofs(_key);
  }

  LIBMESH_DEVICE_INLINE
  QpData qp_data(const unsigned int qp,
                 const bool need_gradients) const
  {
    return QpData(*this, qp, need_gradients);
  }

private:
  libMesh::FEShapeKey _key;
  libMesh::ElemMappingType _mapping_type;
  node_storage_type _nodes;
  unsigned int _n_nodes;
  unsigned int _quadrature_order;
  unsigned int _dim;
  unsigned int _elem_index;
};

template <typename FEAccess, typename CoeffStorage>
using HilbertSolutionAccess = libMesh::detail::HilbertSolutionAccess<FEAccess, CoeffStorage>;

template <typename FEAccess, typename CoeffStorage>
LIBMESH_DEVICE_INLINE auto
make_hilbert_solution_access(const FEAccess & fe,
                             CoeffStorage && coeff,
                             const Number solution_derivative)
{
  return libMesh::detail::make_hilbert_solution_access(
    fe,
    std::forward<CoeffStorage>(coeff),
    solution_derivative);
}

template <typename GoalFunction, typename GoalGradient>
using AnalyticHilbertGoalAccess =
  libMesh::detail::HilbertAnalyticGoalAccess<GoalFunction, GoalGradient>;

template <typename GoalFunction, typename GoalGradient>
LIBMESH_DEVICE_INLINE auto
make_hilbert_analytic_goal_access(GoalFunction && goal_func,
                                  GoalGradient && goal_grad)
{
  return libMesh::detail::make_hilbert_analytic_goal_access(
    std::forward<GoalFunction>(goal_func),
    std::forward<GoalGradient>(goal_grad));
}

template <unsigned int MaxDofs>
class LocalHilbertAccumulator
{
public:
  LIBMESH_DEVICE_INLINE
  explicit LocalHilbertAccumulator(const unsigned int n_dofs)
    : _n_dofs(n_dofs)
  {
    zero();
  }

  LIBMESH_DEVICE_INLINE
  void zero()
  {
    for (unsigned int i = 0; i != MaxDofs; ++i)
      {
        _F[i] = 0.;
        for (unsigned int j = 0; j != MaxDofs; ++j)
          _K[i][j] = 0.;
      }
  }

  LIBMESH_DEVICE_INLINE
  void add_residual(const unsigned int i,
                    const Number value)
  {
    _F[i] += value;
  }

  LIBMESH_DEVICE_INLINE
  void add_jacobian(const unsigned int i,
                    const unsigned int j,
                    const Number value)
  {
    _K[i][j] += value;
  }

  LIBMESH_DEVICE_INLINE
  Number residual(const unsigned int i) const
  {
    return _F[i];
  }

  LIBMESH_DEVICE_INLINE
  Number jacobian(const unsigned int i,
                  const unsigned int j) const
  {
    return _K[i][j];
  }

  LIBMESH_DEVICE_INLINE
  unsigned int n_dofs() const
  {
    return _n_dofs;
  }

private:
  Number _F[MaxDofs];
  Number _K[MaxDofs][MaxDofs];
  unsigned int _n_dofs;
};

} // namespace libMesh::Kokkos::detail

#endif // LIBMESH_KOKKOS_HILBERT_ASSEMBLY_H
