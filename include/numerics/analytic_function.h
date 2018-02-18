// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ANALYTIC_FUNCTION_H
#define LIBMESH_ANALYTIC_FUNCTION_H

// Local Includes
#include "libmesh/function_base.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
template <typename T>
class DenseVector;

/**
 * This class provides function-like objects for which an analytical
 * expression can be provided.  The user may either provide
 * vector-return or number-return functions.  All overridden virtual
 * functions are documented in function_base.h.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
template <typename Output=Number>
class AnalyticFunction : public FunctionBase<Output>
{
public:

  /**
   * Constructor.  Takes a function pointer for scalar
   * return values.
   */
  typedef Output (*OutputFunction)(const Point & p, const Real time);
  AnalyticFunction (OutputFunction fptr);

  /**
   * Constructor.  Takes a function pointer for
   * vector valued functions.
   */
  typedef void (*OutputVectorFunction)(DenseVector<Output> & output,
                                       const Point & p,
                                       const Real time);
  AnalyticFunction (OutputVectorFunction fptr);

  /**
   * Destructor.
   */
  ~AnalyticFunction ();


  /**
   * Pointer to user-provided function that computes
   * the boundary values when an analytical expression
   * is available.
   */
  OutputFunction _number_fptr;

  /**
   * Pointer to user-provided vector valued function.
   */
  OutputVectorFunction _vector_fptr;

  virtual void init () libmesh_override;

  virtual void clear () libmesh_override;

  virtual std::unique_ptr<FunctionBase<Output>> clone () const libmesh_override;

  virtual Output operator() (const Point & p,
                             const Real time=0.) libmesh_override;

  virtual void operator() (const Point & p,
                           const Real time,
                           DenseVector<Output> & output) libmesh_override;
};



// ------------------------------------------------------------
// AnalyticFunction inline methods
template <typename Output>
inline
Output AnalyticFunction<Output>::operator() (const Point & p,
                                             const Real time)
{
  libmesh_assert (this->initialized());
  return (this->_number_fptr(p, time));
}



template <typename Output>
inline
void AnalyticFunction<Output>::operator() (const Point & p,
                                           const Real time,
                                           DenseVector<Output> & output)
{
  libmesh_assert (this->initialized());
  this->_vector_fptr(output, p, time);
  return;
}



template <typename Output>
AnalyticFunction<Output>::AnalyticFunction (OutputFunction fptr) :
  FunctionBase<Output> (),
  _number_fptr (fptr),
  _vector_fptr (libmesh_nullptr)
{
  libmesh_assert(fptr);
  this->_initialized = true;
}



template <typename Output>
inline
AnalyticFunction<Output>::AnalyticFunction (OutputVectorFunction fptr) :
  FunctionBase<Output> (),
  _number_fptr (libmesh_nullptr),
  _vector_fptr (fptr)
{
  libmesh_assert(fptr);
  this->_initialized = true;
}



template <typename Output>
inline
AnalyticFunction<Output>::~AnalyticFunction ()
{
}



template <typename Output>
void AnalyticFunction<Output>::init ()
{
  // dumb double-test
  libmesh_assert ((_number_fptr != libmesh_nullptr) || (_vector_fptr != libmesh_nullptr));

  // definitely ready
  this->_initialized = true;
}



template <typename Output>
inline
void AnalyticFunction<Output>::clear ()
{
  // We probably need a method to reset these later...
  _number_fptr = libmesh_nullptr;
  _vector_fptr = libmesh_nullptr;

  // definitely not ready
  this->_initialized = false;
}



template <typename Output>
inline
std::unique_ptr<FunctionBase<Output>>
AnalyticFunction<Output>::clone () const
{
  return std::unique_ptr<FunctionBase<Output>>
    ( _number_fptr ?
      new AnalyticFunction<Output>(_number_fptr) :
      new AnalyticFunction<Output>(_vector_fptr) );
}


} // namespace libMesh


#endif // LIBMESH_ANALYTIC_FUNCTION_H
