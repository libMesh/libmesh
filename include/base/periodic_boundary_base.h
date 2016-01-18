
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

#ifndef LIBMESH_PERIODIC_BOUNDARY_BASE_H
#define LIBMESH_PERIODIC_BOUNDARY_BASE_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes -----------------------------------
#include "libmesh/point.h"
#include "libmesh/auto_ptr.h"

// C++ Includes   -----------------------------------
#include <set>

namespace libMesh {

// Forward Declarations
class Elem;
class MeshBase;

/**
 * The base class for defining periodic boundaries.
 */
class PeriodicBoundaryBase
{
public:
  enum TransformationType
    { FORWARD=0,
      INVERSE=1 };

  /**
   * The boundary ID of this boundary and its counterpart
   */
  boundary_id_type myboundary, pairedboundary;

  /**
   * Constructor
   */
  PeriodicBoundaryBase();

  /**
   * Copy constructor
   */
  PeriodicBoundaryBase(const PeriodicBoundaryBase & other);

  /**
   * Destructor
   */
  virtual ~PeriodicBoundaryBase() {}

  /**
   * This function should be overloaded by derived classes to
   * define how one finds corresponding nodes on the periodic
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const = 0;

  /**
   * If we want the DofMap to be able to make copies of references and
   * store them in the underlying map, this class must be clone'able,
   * i.e. have a kind of virtual construction mechanism.  The user can
   * also pass a flag to enable an 'inverse transformation' to be cloned
   * from a forward transformation.  Note that not every transformation
   * needs to provide an automatic way to clone an inverse: you can simply
   * add a pair of PeriodicBoundaryBase objects using the appropriate
   * DofMap interface instead.  The simplest way to implement a clone
   * function like this is in terms of a copy constructor, see periodic_boundary.h.
   */
  virtual UniquePtr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const = 0;

  void set_variable(unsigned int var);

  void merge(const PeriodicBoundaryBase & pb);

  bool is_my_variable(unsigned int var_num) const;

protected:

  /**
   * Set of variables for this periodic boundary, empty means all variables possible
   */
  std::set<unsigned int> variables;
};

} // namespace libmesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARY_BASE_H
