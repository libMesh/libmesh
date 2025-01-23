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



#ifndef LIBMESH_WRAPPED_FUNCTOR_H
#define LIBMESH_WRAPPED_FUNCTOR_H

// Local Includes
#include "libmesh/fem_function_base.h"
#include "libmesh/function_base.h"
#include "libmesh/int_range.h"
#include "libmesh/point.h"
#include "libmesh/wrapped_function.h"

// C++ includes
#include <cstddef>
#include <memory>

namespace libMesh
{

/**
 * This class provides a wrapper with which to evaluate a
 * (libMesh-style) function pointer in a FunctionBase-compatible
 * interface. All overridden virtual functions are documented in
 * fem_function_base.h.
 *
 * \author Roy Stogner
 * \date 2015
 */
template <typename Output=Number>
class WrappedFunctor : public FEMFunctionBase<Output>
{
public:

  /**
   * Constructor to wrap FunctionBase functors in a FEMFunctionBase
   * compatible shim.
   */
  WrappedFunctor (const FunctionBase<Output> & func)
    : _func(func.clone())
  { }

  /**
   * Constructor to wrap scalar-valued function pointers.
   */
  WrappedFunctor (const System & sys,
                  Output fptr(const Point & p,
                              const Parameters & parameters,
                              const std::string & sys_name,
                              const std::string & unknown_name) = nullptr,
                  const Parameters * parameters = nullptr,
                  unsigned int varnum=0) :
    _func(std::make_unique<WrappedFunction<Output>>(sys, fptr, parameters, varnum)) {}

  /**
   * This class can't be copy constructed or assigned because it
   * contains a unique_ptr member.
   */
  WrappedFunctor (const WrappedFunctor &) = delete;
  WrappedFunctor & operator= (const WrappedFunctor &) = delete;

  /**
   * The remaining 5 special functions can be defaulted.
   */
  WrappedFunctor (WrappedFunctor &&) = default;
  WrappedFunctor & operator= (WrappedFunctor &&) = default;
  virtual ~WrappedFunctor () = default;

  /**
   * Any post-construction initialization
   */
  virtual void init () override { _func->init(); }

  /**
   * Tell the context we don't need anything from it
   */
  virtual void init_context (const FEMContext & c) override;

  virtual std::unique_ptr<FEMFunctionBase<Output>> clone () const override
  {
    return std::make_unique<WrappedFunctor<Output>>(*_func);
  }

  virtual Output operator() (const FEMContext &,
                             const Point & p,
                             const Real time = 0.) override
  { return _func->operator()(p, time); }

  virtual void operator() (const FEMContext &,
                           const Point & p,
                           const Real time,
                           DenseVector<Output> & output) override
  { _func->operator() (p, time, output); }

  virtual Output component (const FEMContext &,
                            unsigned int i,
                            const Point & p,
                            Real time=0.) override
  { return _func->component(i, p, time); }

protected:

  std::unique_ptr<FunctionBase<Output>> _func;
};


template <typename Output>
void WrappedFunctor<Output>::init_context (const FEMContext & c)
{
  for (auto dim : c.elem_dimensions())
    {
      for (auto v : make_range(c.n_vars()))
        {
          FEAbstract * fe;
          c.get_element_fe(v, fe, dim);
          fe->get_nothing();
        }
    }
}


} // namespace libMesh

#endif // LIBMESH_WRAPPED_FUNCTOR_H
