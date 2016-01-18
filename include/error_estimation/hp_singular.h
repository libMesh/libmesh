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



#ifndef LIBMESH_HP_SINGULAR_H
#define LIBMESH_HP_SINGULAR_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"

// C++ includes
#include <list>

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

// Forward Declarations
class System;


/**
 * This class uses a user-provided list of singularity locations
 * to choose between h refining and p elevation.
 * Currently we assume that a set of elements has already been flagged
 * for h refinement - any elements which do not contain a
 * user-provided singular point are instead flagged for p refinement.
 *
 * \author Roy H. Stogner
 * \date 2006
 */
class HPSingularity
{
public:

  /**
   * Constructor.
   */
  HPSingularity()
  {
    libmesh_experimental();
  }

  /**
   * Destructor.
   */
  virtual ~HPSingularity () {}


  /**
   * This pure virtual function must be redefined
   * in derived classes to take a mesh flagged for h
   * refinement and potentially change the desired
   * refinement type.
   */
  virtual void select_refinement (System & system);

  /**
   * This list, to be filled by the user, should include
   * all singular points in the solution.
   */
  std::list<Point> singular_points;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR

#endif // LIBMESH_HP_SINGULAR_H
