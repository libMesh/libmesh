// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// ------------------------------------------------------------
// FEMFunctionBase class definition
template <typename Output=Number>
class FEMFunctionBase
{
protected:

  /**
   * Constructor.  Optionally takes a master.
   */
FEMFunctionBase () {};

public:

  /**
   * Destructor.
   */
virtual ~FEMFunctionBase () {};


  /**
   * Prepares a context object for use.
   * 
   * Most problems will want to reimplement this for efficiency, in
   * order to call FE::get_*() as their particular function requires.
   */
virtual void init_context (const FEMContext &) {}


  // ------------------------------------------------------
  // misc
  /**
   * @returns the scalar value at coordinate
   * \p p and time \p time, which defaults to zero.
   * Purely virtual, so you have to overload it.
   * Note that this cannot be a const method, check \p MeshFunction.
   */
virtual Output operator() (const FEMContext&, const Point& p,
			     const Real time = 0.) = 0;



protected:

/**
 * Variable index to decide which overloaded function should
 * be accessed
 */

unsigned int var_index ;

};

} // namespace libMesh

#endif // LIBMESH_FEM_FUNCTION_BASE_H
