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

#ifndef LIBMESH_RB_THETA_H
#define LIBMESH_RB_THETA_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <vector>

namespace libMesh
{

class RBParameters;

/**
 * This class is part of the rbOOmit framework.
 *
 * RBTheta provides a default functor class with which
 * to define the parameter-dependent functions required
 * for the PDE decomposition employed by the Reduced
 * Basis method.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBTheta : public ReferenceCountedObject<RBTheta>
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  RBTheta () {}

  /**
   * Destructor.
   */
  virtual ~RBTheta () {}

  /**
   * Evaluate the functor object for the given parameter.
   * Default implementation is to return 1, overload
   * to provide problem dependent behavior.
   */
  virtual Number evaluate(const RBParameters &) { return 1.; }
};

}

#endif // LIBMESH_RB_THETA_H
