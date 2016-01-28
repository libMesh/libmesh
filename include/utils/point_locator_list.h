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



#ifndef LIBMESH_POINT_LOCATOR_LIST_H
#define LIBMESH_POINT_LOCATOR_LIST_H

// Local Includes
#include "libmesh/point_locator_base.h"

// C++ includes
#include <cstddef>
#include <utility> // pair
#include <vector>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Point;
class Elem;

/**
 * This is a point locator.  It locates points in space
 * using a list of element centroids: given a mesh this
 * locator returns the element that is @e closest to the
 * given point in global coordinates.
 * Note that this may yield @e severe difficulties in
 * case of extremely distorted elements, e.g.\ infinite
 * elements.
 *
 * This list version is not efficient, but especially
 * @e reliable for the case of finding the @e closest
 * dim-1 element (nearest-surface-element, e.g. used
 * for projecting boundary conditions from a surface
 * mesh onto a volumetric mesh).
 * It should be noted that this class only works when
 * the element list in the associated mesh object is
 * @e not modified (like refinement etc).  Otherwise,
 * the point locator has to be cleared and re-initialized.
 * Use \p PointLocatorBase::build() to create objects
 * of this type at run time.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class PointLocatorList : public PointLocatorBase
{
public:

  /**
   * Constructor.  Needs the \p mesh which holds the
   * elements that should be identified as being close.
   * Optionally takes a master interpolator.  This master
   * helps in saving memory by simply only setting up @e
   * one list for all point locators.  Only the master
   * locator holds a list, the others simply  use the
   * master's list.
   */
  PointLocatorList (const MeshBase & mesh,
                    const PointLocatorBase * master = libmesh_nullptr);

  /**
   * Destructor.
   */
  ~PointLocatorList ();

  /**
   * Clears the locator.  Overloaded from base class.  This method
   * frees dynamic memory using "delete".
   */
  virtual void clear() libmesh_override;

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used.  Overloaded from base class.  This method allocates dynamic
   * memory using "new".
   */
  virtual void init() libmesh_override;

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located, optionally restricted to a set of allowed subdomains.
   * Overloaded from base class.
   */
  virtual const Elem * operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains = libmesh_nullptr) const libmesh_override;

  /**
   * Enables out-of-mesh mode.  In this mode, if asked to find a point
   * that is contained in no mesh at all, the point locator will
   * return a NULL pointer instead of crashing.  Per default, this
   * mode is off.
   */
  virtual void enable_out_of_mesh_mode () libmesh_override;

  /**
   * Disables out-of-mesh mode (default).  If asked to find a point
   * that is contained in no mesh at all, the point locator will now
   * crash.
   */
  virtual void disable_out_of_mesh_mode () libmesh_override;

  /**
   * Set a tolerance to use when determining
   * if a point is contained within the mesh.
   */
  virtual void set_close_to_point_tol(Real close_to_point_tol) libmesh_override;

  /**
   * Specify that we do not want to use a user-specified tolerance to
   * determine if a point is contained within the mesh.
   */
  virtual void unset_close_to_point_tol() libmesh_override;

protected:
  /**
   * Pointer to the list of element centroids.  Only the
   * master @e has such a list.  For servants, this
   * pointer points to the list of the master.   Note that
   * it's not a std::list as the name might suggest, but a std::vector.
   */
  std::vector<std::pair<Point, const Elem *> > * _list;
};

} // namespace libMesh

#endif // LIBMESH_POINT_LOCATOR_LIST_H
