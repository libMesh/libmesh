// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FDM_GRADIENT_H
#define LIBMESH_FDM_GRADIENT_H

// libMesh includes
#include "libmesh/fem_function_base.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/tensor_tools.h"

namespace libMesh
{

template <typename GradType>
class FDMGradient : public FEMFunctionBase<GradType>
{
public:
  typedef typename TensorTools::DecrementRank<GradType>::type ValType;

  FDMGradient (FEMFunctionBase<ValType> & value_func,
               Real eps) :
    _val_func(value_func.clone()), _eps(eps)
  {}

  virtual void init_context (const FEMContext & c) override
  { _val_func->init_context(c); }

  virtual std::unique_ptr<FEMFunctionBase<GradType>> clone () const override
  { return std::make_unique<FDMGradient<GradType>>(*_val_func, _eps); }

  virtual GradType operator() (const FEMContext & c,
                               const Point & p,
                               const Real time = 0.) override
  {
    GradType g;

    auto & val = *_val_func;

    Real one_over_dim = Real(0.5) / _eps;

    g(0) = (val(c, p+Point(_eps), time) -
            val(c, p+Point(-_eps), time)) * one_over_dim;
#if LIBMESH_DIM > 1
    g(1) = (val(c, p+Point(0,_eps), time) -
            val(c, p+Point(0,-_eps), time)) * one_over_dim;
#endif
#if LIBMESH_DIM > 2
    g(2) = (val(c, p+Point(0,0,_eps), time) -
            val(c, p+Point(0,0,-_eps), time)) * one_over_dim;
#endif

    return g;
  }

  virtual void operator() (const FEMContext & c,
                           const Point & p,
                           const Real time,
                           DenseVector<GradType> & output) override
  {
    auto sz = output.size();
    DenseVector<ValType> v(sz);

    auto & val = *_val_func;

    val(c, p+Point(_eps), time, v);
    for (auto i : make_range(sz))
      output(i)(0) = v(i);

    val(c, p+Point(-_eps), time, v);
    for (auto i : make_range(sz))
      {
        output(i)(0) -= v(i);
        output(i)(0) /= 2;
      }

#if LIBMESH_DIM > 1
    val(c, p+Point(0,_eps), time, v);
    for (auto i : make_range(sz))
      output(i)(1) = v(i);

    val(c, p+Point(0,-_eps), time, v);
    for (auto i : make_range(sz))
      {
        output(i)(1) -= v(i);
        output(i)(1) /= 2;
      }
#endif
#if LIBMESH_DIM > 2
    val(c, p+Point(0,0,_eps), time, v);
    for (auto i : make_range(sz))
      output(i)(2) = v(i);

    val(c, p+Point(0,0,-_eps), time, v);
    for (auto i : make_range(sz))
      {
        output(i)(2) -= v(i);
        output(i)(2) /= 2;
      }
#endif
  }


  virtual GradType component (const FEMContext & c,
                              unsigned int i,
                              const Point & p,
                              Real time) override
  {
    GradType g;

    auto & val = *_val_func;

    Real one_over_dim = Real(0.5) / _eps;

    g(0) = (val.component(c, i, p+Point(_eps), time) -
            val.component(c, i, p+Point(-_eps), time)) * one_over_dim;
#if LIBMESH_DIM > 1
    g(1) = (val.component(c, i, p+Point(0,_eps), time) -
            val.component(c, i, p+Point(0,-_eps), time)) * one_over_dim;
#endif
#if LIBMESH_DIM > 2
    g(2) = (val.component(c, i, p+Point(0,0,_eps), time) -
            val.component(c, i, p+Point(0,0,-_eps), time)) * one_over_dim;
#endif

    return g;
  }

private:

  std::unique_ptr<FEMFunctionBase<ValType>> _val_func;

  Real _eps;
};


} // namespace libMesh

#endif // LIBMESH_FDM_GRADIENT_H
