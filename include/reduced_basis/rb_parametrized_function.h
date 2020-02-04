// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_RB_PARAMETRIZED_FUNCTION_H
#define LIBMESH_RB_PARAMETRIZED_FUNCTION_H



#include "libmesh/libmesh_common.h"


namespace libMesh
{

class RBParameters;
template <typename> class PointTempl;
typedef PointTempl<Real> Point;
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;

/**
 * A simple functor class that provides a RBParameter-dependent function.
 *
 * \author David Knezevic
 * \date 2012
 * \brief Provides a reduced basis parameterized function.
 */
class RBParametrizedFunction
{
public:

  /**
   * Virtual evaluate() gives us a vtable, so there's no cost in adding a
   * virtual destructor for safety's sake.
   */
  virtual ~RBParametrizedFunction() {}

  /**
   * Evaluate this parametrized function for the parameter value
   * \p mu at the point \p p.
   */
  virtual Number evaluate(const RBParameters &,
                          const Point &,
                          const Elem &) { return 0.; }
};

}

#endif // LIBMESH_RB_PARAMETRIZED_FUNCTION_H
