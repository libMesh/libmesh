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



#ifndef LIBMESH_FEM_FUNCTION_BASE_H
#define LIBMESH_FEM_FUNCTION_BASE_H

// C++ includes



// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_vector.h" // required to instantiate a DenseVector<> below
#include "libmesh/auto_ptr.h"
#include "libmesh/fem_context.h"

namespace libMesh
{

// Forward Declarations
class Point;

/**
 * FEMFunctionBase is a base class from which users can derive in
 * order to define "function-like" objects that can be used within
 * FEMSystem.
 *
 * \author Roy Stogner
 * \date 2012
 */
template <typename Output=Number>
class FEMFunctionBase
{
protected:

  /**
   * Constructor.
   */
  FEMFunctionBase () {}

public:

  /**
   * Destructor.
   */
  virtual ~FEMFunctionBase () {}


  /**
   * Prepares a context object for use.
   *
   * Most problems will want to reimplement this for efficiency, in
   * order to call FE::get_*() as their particular function requires.
   */
  virtual void init_context (const FEMContext &) {}

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  virtual UniquePtr<FEMFunctionBase<Output> > clone () const = 0;

  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual Output operator() (const FEMContext &,
                             const Point & p,
                             const Real time = 0.) = 0;


  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p.
   */
  void operator() (const FEMContext &,
                   const Point & p,
                   DenseVector<Output> & output);

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p and for time \p time.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual void operator() (const FEMContext &,
                           const Point & p,
                           const Real time,
                           DenseVector<Output> & output) = 0;

  /**
   * @returns the vector component \p i at coordinate
   * \p p and time \p time.
   * Subclasses aren't required to overload this, since the default
   * implementation is based on the full vector evaluation, which is
   * often correct.
   * Subclasses are recommended to overload this, since the default
   * implementation is based on a vector evaluation, which is usually
   * unnecessarily inefficient.
   */
  virtual Output component(const FEMContext &,
                           unsigned int i,
                           const Point & p,
                           Real time=0.);
};

template <typename Output>
inline
Output FEMFunctionBase<Output>::component (const FEMContext & context,
                                           unsigned int i,
                                           const Point & p,
                                           Real time)
{
  DenseVector<Output> outvec(i+1);
  (*this)(context, p, time, outvec);
  return outvec(i);
}

template <typename Output>
inline
void FEMFunctionBase<Output>::operator() (const FEMContext & context,
                                          const Point & p,
                                          DenseVector<Output> & output)
{
  // Call the time-dependent function with t=0.
  this->operator()(context, p, 0., output);
}

} // namespace libMesh

#endif // LIBMESH_FEM_FUNCTION_BASE_H
