
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

#ifndef LIBMESH_PERIODIC_BOUNDARY_H
#define LIBMESH_PERIODIC_BOUNDARY_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes -----------------------------------
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/vector_value.h" // RealVectorValue

namespace libMesh {

// Forward Declarations
class Elem;
class MeshBase;

/**
 * The definition of a periodic boundary.
 */
class PeriodicBoundary : public PeriodicBoundaryBase
{
public:
  /**
   * Constructor
   */
  PeriodicBoundary();

  /**
   * Destructor
   */
  virtual ~PeriodicBoundary() {}

  /**
   * Copy constructor, with option for the copy to represent an inverse transformation.
   */
  PeriodicBoundary(const PeriodicBoundary & o, TransformationType t = FORWARD);

  /**
   * Constructor taking a reference to the translation vector.
   */
  PeriodicBoundary(const RealVectorValue & vector);

  /**
   * This function should be overloaded by derived classes to
   * define how one finds corresponding nodes on the periodic
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const libmesh_override;

  /**
   * If we want the DofMap to be able to make copies of references and
   * store them in the underlying map, this class must be clone'able,
   * i.e. have a kind of virtual construction mechanism.
   */
  virtual UniquePtr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const libmesh_override;

protected:
  // One of these days we'll support rotated boundaries
  // RealTensor rotation_matrix;

  // The vector which is added to points in myboundary
  // to produce corresponding points in pairedboundary
  RealVectorValue translation_vector;
};

} // namespace libmesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARY_H
