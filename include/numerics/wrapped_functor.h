// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/point.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

/**
 * This class provides a wrapper with which to evaluate a
 * (libMesh-style) function pointer in a FunctionBase-compatible
 * interface.
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
   * compatible shim
   */
  WrappedFunctor (const FunctionBase<Output> & func)
    : _func(func.clone())
  { }

  virtual UniquePtr<FEMFunctionBase<Output> > clone () const libmesh_override
  {
    return UniquePtr<FEMFunctionBase<Output> >
      (new WrappedFunctor<Output> (*_func));
  }

  /**
   * @returns the scalar value of variable varnum at coordinate \p p
   * and time \p time.
   */
  virtual Output operator() (const FEMContext &,
                             const Point & p,
                             const Real time = 0.) libmesh_override
  { return _func->operator()(p, time); }

  /**
   * Return function for vectors.
   * Returns in \p output the values of all system variables at the
   * coordinate \p p and for time \p time.
   */
  virtual void operator() (const FEMContext &,
                           const Point & p,
                           const Real time,
                           DenseVector<Output> & output) libmesh_override
  { _func->operator() (p, time, output); }

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   */
  virtual Output component (const FEMContext &,
                            unsigned int i,
                            const Point & p,
                            Real time=0.) libmesh_override
  { return _func->component(i, p, time); }

protected:

  UniquePtr<FunctionBase<Output> > _func;
};



} // namespace libMesh

#endif // LIBMESH_WRAPPED_FUNCTOR_H
