
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

#ifndef __periodic_boundary_h__
#define __periodic_boundary_h__

// Local Includes -----------------------------------
#include "libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes -----------------------------------
#include "point.h"
#include "vector_value.h" // RealVectorValue

// C++ Includes   -----------------------------------
#include <set>

namespace libMesh {

// Forward Declarations
class Elem;
class MeshBase;

/**
 * The definition of a periodic boundary.
 */
class PeriodicBoundary
{
public:
  /**
   * The boundary ID of this boundary and its counterpart
   */
  unsigned int myboundary,
    pairedboundary;

  /**
   * Constructor
   */
  PeriodicBoundary();

  /**
   * Destructor
   */
  virtual ~PeriodicBoundary() {}

  /**
   * Copy constructor
   */
  PeriodicBoundary(const PeriodicBoundary & o, bool inverse = false);

  /**
   * Constructor taking a reference to the translation vector.
   */
  PeriodicBoundary(const RealVectorValue & vector);

  /**
   * This function should be overloaded by derived classes to 
   * define how one finds corresponding nodes on the periodic 
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const;

  void set_variable(unsigned int var);

  void merge(const PeriodicBoundary & pb);

  bool is_my_variable(unsigned int var_num) const;

protected:
  // One of these days we'll support rotated boundaries
  // RealTensor rotation_matrix;

  // The vector which is added to points in myboundary
  // to produce corresponding points in pairedboundary
  RealVectorValue translation_vector;

  // Set of variables for this periodic boundary, empty means all varaibles possible
  std::set<unsigned int> variables;
};

} // namespace libmesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // __periodic_boundary_h__
