// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __analytic_function_h__
#define __analytic_function_h__

// C++ includes


// Local Includes
#include "function_base.h"

namespace libMesh
{



// Forward Declarations
template <typename T>
class DenseVector;


/**
 * This class provides function-like objects for which an
 * analytical expression can be provided.  The user may
 * either provide vector-return or number-return functions.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// AnalyticFunction class definition
class AnalyticFunction : public FunctionBase
{
public:

  /**
   * Constructor.  Takes a function pointer for scalar
   * return values.
   */
  AnalyticFunction (Number fptr(const Point& p,
				const Real time));

  /**
   * Constructor.  Takes a function pointer for
   * vector valued functions.
   */
  AnalyticFunction (void fptr(DenseVector<Number>& output,
			      const Point& p,
			      const Real time));
  /**
   * Destructor.
   */
  ~AnalyticFunction ();


  /**
   * Pointer to user-provided function that computes
   * the boundary values when an analytical expression
   * is available.
   */
  Number (* _number_fptr) (const Point& p,
			   const Real time);

  /**
   * Pointer to user-provided vector valued function.
   */
  void (* _vector_fptr) (DenseVector<Number>& output,
			 const Point& p,
			 const Real time);

  /**
   * The actual initialization process.
   */
  void init ();

  /**
   * Clears the function.
   */
  void clear ();

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point& p,
		     const Real time=0.);

  /**
   * Like before, but returns the values in a
   * writable reference.
   */
  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output);

};



// ------------------------------------------------------------
// AnalyticFunction inline methods
inline
Number AnalyticFunction::operator() (const Point& p,
				     const Real time)
{
  libmesh_assert (this->initialized());
  return (this->_number_fptr(p, time));
}



inline
void AnalyticFunction::operator() (const Point& p,
				   const Real time,
				   DenseVector<Number>& output)
{
  libmesh_assert (this->initialized());
  this->_vector_fptr(output, p, time);
  return;
}


} // namespace libMesh


#endif

