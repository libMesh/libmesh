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



#ifndef LIBMESH_REFINEMENT_SELECTOR_H
#define LIBMESH_REFINEMENT_SELECTOR_H

// Local Includes
#include "libmesh/libmesh_common.h"

// C++ includes
#include <vector>

namespace libMesh
{

// Forward Declarations
class EquationSystems;
class System;

/**
 * This abstract class provides an interface to methods for
 * selecting the type of refinement to be used on each element
 * in a given mesh.  Currently we assume that a set of elements
 * has already been flagged for h refinement, and the only concrete
 * subclass will change some of those elements to be flagged for p
 * refinement.  Future subclasses may handle anisotropic refinement
 * instead.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class RefinementSelector
{
public:

  /**
   * Constructor. Empty.
   */
  RefinementSelector() {}

  /**
   * Destructor.
   */
  virtual ~RefinementSelector() {}


  /**
   * This pure virtual function must be redefined
   * in derived classes to take a mesh flagged for h
   * refinement and potentially change the desired
   * refinement type.
   */
  virtual void select_refinement (const System & system);

  /**
   * This vector can be used to "scale" certain
   * variables in a system.
   * If the mask is not empty, the consideration given to each
   * component will be scaled by component_scale[c].
   */
  std::vector<float> component_scale;
};

} // namespace libMesh


#endif // LIBMESH_REFINEMENT_SELECTOR_H
