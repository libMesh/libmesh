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



#ifndef LIBMESH_FUNCTION_BASE_H
#define LIBMESH_FUNCTION_BASE_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_vector.h" // required to instantiate a DenseVector<> below
#include "libmesh/auto_ptr.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
class Point;

/**
 * This is the base class for functor-like classes.  These
 * entities are functions (in the mathematical sense) of time
 * and space, \f$ f(\mathbf{x},t) =  \mbox{\texttt{v}} \f$,
 * where \p v may be either a \p Number or a \p DenseVector<Number>.
 * Children of this base class implement different styles of
 * data retrieval for these functions.  Use the constructors
 * of the derived classes for creating new objects. The
 * required input of each derived class thwarts the effective
 * use of the commonly used \p build() member.  But afterwards
 * the virtual members allow the convenient and libMesh-common
 * usage through a \p FunctionBase *. Note that for functor objects
 * for vector-valued variables, it is assumed each component is indexed
 * contiguously; i.e. if u_var is index 3, then libMesh expects
 * the x-component of u_var is index 3, the y-component is index 4,
 * and the z-component is index 5. Note that for 2-D elements in 3
 * spatial dimensions, libMesh is expecting 2 components (i.e.
 * mesh_dimension() number of components).
 *
 * \author Daniel Dreyer
 * \date 2003
 */
template <typename Output=Number>
class FunctionBase
{
protected:

  /**
   * Constructor.  Optionally takes a master.
   */
  explicit
  FunctionBase (const FunctionBase * master = libmesh_nullptr);

public:

  /**
   * Destructor.
   */
  virtual ~FunctionBase ();

  /**
   * The actual initialization process.
   */
  virtual void init () {}

  /**
   * Clears the function.
   */
  virtual void clear () {}

  /**
   * Returns a new copy of the function.  The new copy should be as
   * ``deep'' as necessary to allow independent destruction and
   * simultaneous evaluations of the copies in different threads.
   */
  virtual UniquePtr<FunctionBase<Output> > clone () const = 0;


  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual Output operator() (const Point & p,
                             const Real time = 0.) = 0;

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p.
   */
  void operator() (const Point & p,
                   DenseVector<Output> & output);

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p and for time \p time.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   * Can optionally provide subdomain_ids which will restrict
   * the function to operate on elements with subdomain id contained
   * in the set. This is useful in cases where there are multiple
   * dimensioned elements, for example.
   */
  virtual void operator() (const Point & p,
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
  virtual Output component(unsigned int i,
                           const Point & p,
                           Real time=0.);


  /**
   * @returns \p true when this object is properly initialized
   * and ready for use, \p false otherwise.
   */
  bool initialized () const;

  /**
   * Function to set whether this is a time-dependent function or not.
   * This is intended to be only used by subclasses who cannot natively
   * determine time-dependence. In such a case, this function should
   * be used immediately following construction.
   */
  void set_is_time_dependent( bool is_time_dependent);

  /**
   * @returns \p true when the function this object represents
   * is actually time-dependent, \p false otherwise.
   */
  bool is_time_dependent() const;

protected:

  /**
   * Const pointer to our master, initialized to \p NULL.
   * There may be cases where multiple functions are required,
   * but to save memory, one master handles some centralized
   * data.
   */
  const FunctionBase * _master;

  /**
   * When \p init() was called so that everything is ready
   * for calls to \p operator() (...), then this \p bool is true.
   */
  bool _initialized;

  /**
   * Cache whether or not this function is actually time-dependent.
   */
  bool _is_time_dependent;

};


// ------------------------------------------------------------
// FunctionBase inline methods

template<typename Output>
inline
FunctionBase<Output>::FunctionBase (const FunctionBase * master) :
  _master             (master),
  _initialized        (false),
  _is_time_dependent  (true) // Assume we are time-dependent until the user says otherwise
{
}



template<typename Output>
inline
FunctionBase<Output>::~FunctionBase ()
{
}



template <typename Output>
inline
bool FunctionBase<Output>::initialized() const
{
  return (this->_initialized);
}

template <typename Output>
inline
void FunctionBase<Output>::set_is_time_dependent( bool is_time_dependent )
{
  this->_is_time_dependent = is_time_dependent;
}

template <typename Output>
inline
bool FunctionBase<Output>::is_time_dependent() const
{
  return (this->_is_time_dependent);
}


template <typename Output>
inline
Output FunctionBase<Output>::component (unsigned int i,
                                        const Point & p,
                                        Real time)
{
  DenseVector<Output> outvec(i+1);
  (*this)(p, time, outvec);
  return outvec(i);
}



template <typename Output>
inline
void FunctionBase<Output>::operator() (const Point & p,
                                       DenseVector<Output> & output)
{
  // Call the time-dependent function with t=0.
  this->operator()(p, 0., output);
}

} // namespace libMesh

#endif // LIBMESH_FUNCTION_BASE_H
