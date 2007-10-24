// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __function_base_h__
#define __function_base_h__

// C++ includes



// Local Includes
#include "libmesh_common.h"



// Forward Declarations
template <typename T> class DenseVector;
class FunctionBase;
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
 * usage through a \p FunctionBase*.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// FunctionBase class definition
class FunctionBase 
{
protected:

  /**
   * Constructor.  Optionally takes a master.
   */
  FunctionBase (const FunctionBase* master = NULL);

public:

  /**
   * Destructor.
   */
  virtual ~FunctionBase ();

  /**
   * The actual initialization process.
   */
  virtual void init () = 0;

  /**
   * Clears the function.
   */
  virtual void clear () = 0;


  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual Number operator() (const Point& p, 
			     const Real time = 0.) = 0;

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p.
   */
  void operator() (const Point& p,
		   DenseVector<Number>& output);

  /**
   * Return function for vectors.
   * Returns in \p output the values of the data at the
   * coordinate \p p and for time \p time.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
  virtual void operator() (const Point& p,
			   const Real time,
			   DenseVector<Number>& output) = 0;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use, \p false otherwise.
   */
  bool initialized () const;


protected:

  /**
   * Const pointer to our master, initialized to \p NULL.
   * There may be cases where multiple functions are required,
   * but to save memory, one master handles some centralized
   * data.
   */
  const FunctionBase* _master;

  /**
   * When \p init() was called so that everything is ready
   * for calls to \p operator() (...), then this \p bool is true.
   */
  bool _initialized;

};


// ------------------------------------------------------------
// FunctionBase inline methods
inline
bool FunctionBase::initialized() const
{
  return (this->_initialized);
}



inline
void FunctionBase::operator() (const Point& p,
			       DenseVector<Number>& output)
{
  // Call the time-dependent function with t=0.
  this->operator()(p, 0., output);
}

#endif
