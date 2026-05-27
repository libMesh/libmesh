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

#ifndef LIBMESH_HILBERT_ASSEMBLY_H
#define LIBMESH_HILBERT_ASSEMBLY_H

#include "hilbert_assembly_kernel.h"

#include "libmesh/fdm_gradient.h"
#include "libmesh/fe_abstract.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/function_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/quadrature.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{
namespace detail
{

class HostHilbertFEAccess
{
public:
  class QpData
  {
  public:
    LIBMESH_DEVICE_INLINE
    QpData(const HostHilbertFEAccess & fe,
           const unsigned int qp)
      : _fe(fe),
        _qp(qp),
        _Jinv(fe.build_inverse_jacobian(qp))
    {
    }

    LIBMESH_DEVICE_INLINE
    Real JxW() const
    {
      return _fe.JxW(_qp);
    }

    LIBMESH_DEVICE_INLINE
    Real phi(const unsigned int i) const
    {
      return _fe.phi(i, _qp);
    }

    LIBMESH_DEVICE_INLINE
    const RealGradient & dphi(const unsigned int i) const
    {
      return _fe.dphi(i, _qp);
    }

    LIBMESH_DEVICE_INLINE
    const Point & xyz() const
    {
      return _fe.xyz(_qp);
    }

    LIBMESH_DEVICE_INLINE
    unsigned int qp_index() const
    {
      return _qp;
    }

    LIBMESH_DEVICE_INLINE
    const Point & reference_point() const
    {
      return _fe.reference_point(_qp);
    }

    LIBMESH_DEVICE_INLINE
    const RealTensor & inverse_jacobian() const
    {
      return _Jinv;
    }

    LIBMESH_DEVICE_INLINE
    unsigned int elem_index() const
    {
      return _fe.elem_index();
    }

  private:
    const HostHilbertFEAccess & _fe;
    unsigned int _qp;
    RealTensor _Jinv;
  };

  HostHilbertFEAccess(FEMContext & c,
                      const unsigned int var,
                      const unsigned int hilbert_order,
                      const unsigned int elem_index = libMesh::invalid_uint)
    : _n_dofs(c.n_dof_indices(var)),
      _elem_index(elem_index),
      _JxW(c.get_element_fe(var)->get_JxW()),
      _phi(c.get_element_fe(var)->get_phi()),
      _xyz(c.get_element_fe(var)->get_xyz()),
      _reference_points(c.get_element_qrule().get_points()),
      _fe_map(c.get_element_fe(var)->get_fe_map()),
      _dphi(hilbert_order > 0 ? &c.get_element_fe(var)->get_dphi() : nullptr)
  {
  }

  unsigned int n_qpoints() const
  {
    return cast_int<unsigned int>(_JxW.size());
  }

  unsigned int n_dofs() const
  {
    return _n_dofs;
  }

  Real JxW(const unsigned int qp) const
  {
    return _JxW[qp];
  }

  Real phi(const unsigned int i, const unsigned int qp) const
  {
    return _phi[i][qp];
  }

  const RealGradient & dphi(const unsigned int i, const unsigned int qp) const
  {
    libmesh_assert(_dphi);
    return (*_dphi)[i][qp];
  }

  const Point & xyz(const unsigned int qp) const
  {
    return _xyz[qp];
  }

  const Point & reference_point(const unsigned int qp) const
  {
    return _reference_points[qp];
  }

  unsigned int elem_index() const
  {
    return _elem_index;
  }

  RealTensor build_inverse_jacobian(const unsigned int qp) const
  {
    RealTensor Jinv;
    Jinv(0, 0) = _fe_map.get_dxidx()[qp];
#if LIBMESH_DIM > 1
    Jinv(0, 1) = _fe_map.get_dxidy()[qp];
    Jinv(1, 0) = _fe_map.get_detadx()[qp];
    Jinv(1, 1) = _fe_map.get_detady()[qp];
#endif
#if LIBMESH_DIM > 2
    Jinv(0, 2) = _fe_map.get_dxidz()[qp];
    Jinv(1, 2) = _fe_map.get_detadz()[qp];
    Jinv(2, 0) = _fe_map.get_dzetadx()[qp];
    Jinv(2, 1) = _fe_map.get_dzetady()[qp];
    Jinv(2, 2) = _fe_map.get_dzetadz()[qp];
#endif
    return Jinv;
  }

  LIBMESH_DEVICE_INLINE
  QpData qp_data(const unsigned int qp,
                 const bool) const
  {
    return QpData(*this, qp);
  }

private:
  const unsigned int _n_dofs;
  const unsigned int _elem_index;
  const std::vector<Real> & _JxW;
  const std::vector<std::vector<Real>> & _phi;
  const std::vector<Point> & _xyz;
  const std::vector<Point> & _reference_points;
  const FEMap & _fe_map;
  const std::vector<std::vector<RealGradient>> * _dphi;
};

class HostHilbertGoalAccess
{
public:
  HostHilbertGoalAccess(FEMFunctionBase<Number> & goal_func,
                        FDMGradient<Gradient> * goal_grad,
                        FEMContext & input_context)
    : _goal_func(goal_func),
      _goal_grad(goal_grad),
      _input_context(input_context)
  {
  }

  template <typename QpData>
  Number value(const QpData &, const Point & p)
  {
    return _goal_func(_input_context, p);
  }

  template <typename QpData>
  Gradient gradient(const QpData &, const Point & p)
  {
    libmesh_assert(_goal_grad);
    return (*_goal_grad)(_input_context, p);
  }

private:
  FEMFunctionBase<Number> & _goal_func;
  FDMGradient<Gradient> * const _goal_grad;
  FEMContext & _input_context;
};

template <typename GradType>
class FunctionFDMGradient : public FunctionBase<GradType>
{
public:
  typedef typename TensorTools::DecrementRank<GradType>::type ValType;

  FunctionFDMGradient(FunctionBase<ValType> & value_func,
                      const Real eps)
    : _val_func(value_func.clone()),
      _eps(eps)
  {
  }

  virtual std::unique_ptr<FunctionBase<GradType>> clone() const override
  {
    return std::make_unique<FunctionFDMGradient<GradType>>(*_val_func, _eps);
  }

  virtual GradType operator()(const Point & p,
                              const Real time = 0.) override
  {
    GradType g;

    auto & val = *_val_func;
    const Real one_over_dim = Real(0.5) / _eps;

    g(0) = (val(p + Point(_eps), time) -
            val(p + Point(-_eps), time)) * one_over_dim;
#if LIBMESH_DIM > 1
    g(1) = (val(p + Point(0, _eps), time) -
            val(p + Point(0, -_eps), time)) * one_over_dim;
#endif
#if LIBMESH_DIM > 2
    g(2) = (val(p + Point(0, 0, _eps), time) -
            val(p + Point(0, 0, -_eps), time)) * one_over_dim;
#endif

    return g;
  }

  virtual void operator()(const Point & p,
                          const Real time,
                          DenseVector<GradType> & output) override
  {
    const unsigned int sz = cast_int<unsigned int>(output.size());
    DenseVector<ValType> v(sz);

    auto & val = *_val_func;

    val(p + Point(_eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      output(i)(0) = v(i);

    val(p + Point(-_eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      {
        output(i)(0) -= v(i);
        output(i)(0) /= (2 * _eps);
      }

#if LIBMESH_DIM > 1
    val(p + Point(0, _eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      output(i)(1) = v(i);

    val(p + Point(0, -_eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      {
        output(i)(1) -= v(i);
        output(i)(1) /= (2 * _eps);
      }
#endif
#if LIBMESH_DIM > 2
    val(p + Point(0, 0, _eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      output(i)(2) = v(i);

    val(p + Point(0, 0, -_eps), time, v);
    for (unsigned int i = 0; i != sz; ++i)
      {
        output(i)(2) -= v(i);
        output(i)(2) /= (2 * _eps);
      }
#endif
  }

private:
  std::unique_ptr<FunctionBase<ValType>> _val_func;
  Real _eps;
};

class HostHilbertAccumulator
{
public:
  HostHilbertAccumulator(DenseSubVector<Number> & F,
                         DenseSubMatrix<Number> & K)
    : _F(F),
      _K(K)
  {
  }

  void add_residual(const unsigned int i,
                    const Number value)
  {
    _F(i) += value;
  }

  void add_jacobian(const unsigned int i,
                    const unsigned int j,
                    const Number value)
  {
    _K(i, j) += value;
  }

  unsigned int n_dofs() const
  {
    return _F.size();
  }

private:
  DenseSubVector<Number> & _F;
  DenseSubMatrix<Number> & _K;
};

} // namespace detail
} // namespace libMesh

#endif // LIBMESH_HILBERT_ASSEMBLY_H
