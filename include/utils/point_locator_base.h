// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_POINT_LOCATOR_BASE_H
#define LIBMESH_POINT_LOCATOR_BASE_H

// Local Includes
#include "libmesh/reference_counted_object.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_point_locator_type.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{



// Forward Declarations
class PointLocatorBase;
class MeshBase;
class Point;
class TreeBase;
class Elem;
class Node;



/**
 * This is the base class for point locators.  They locate
 * points in space: given a mesh they return the element
 * and local coordinates for a given point in global coordinates.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class PointLocatorBase : public ReferenceCountedObject<PointLocatorBase>
{
protected:
  /**
   * Constructor.  Protected so that this base class
   * cannot be explicitly instantiated.  Takes a master
   * PointLocator that helps in saving memory.
   */
  PointLocatorBase (const MeshBase & mesh,
                    const PointLocatorBase * master);

public:
  /**
   * Destructor.
   */
  virtual ~PointLocatorBase ();

  /**
   * Builds an PointLocator for the mesh \p mesh.
   * Optionally takes a master PointLocator to save memory.
   * An \p UniquePtr<PointLocatorBase> is returned to prevent memory leak.
   * This way the user need not remember to delete the object.
   */
  static UniquePtr<PointLocatorBase> build (PointLocatorType t,
                                            const MeshBase & mesh,
                                            const PointLocatorBase * master = libmesh_nullptr);

  /**
   * Clears the \p PointLocator.
   */
  virtual void clear() = 0;

  /**
   * Initializes the point locator, so that the \p operator() methods can
   * be used.  Pure virtual.
   */
  virtual void init() = 0;

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located.  Pure virtual. Optionally allows the user to restrict
   * the subdomains searched.
   */
  virtual const Elem * operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains = libmesh_nullptr) const = 0;

  /**
   * Locates a set of elements in proximity to the point with global coordinates
   * \p p  Pure virtual. Optionally allows the user to restrict the subdomains searched.
   */
  virtual void operator() (const Point & p,
                           std::set<const Elem *> & candidate_elements,
                           const std::set<subdomain_id_type> * allowed_subdomains = libmesh_nullptr) const = 0;

  /**
   * Locates a Node with global coordinates \p p or returns NULL if no
   * such Node can be found.
   *
   * Virtual so subclasses can override for efficiency, but has a
   * default implementation that works based on element lookup.
   *
   * Optionally allows the user to restrict the subdomains searched;
   * with such a restriction, only a Node belonging to an element on
   * one or more of those subdomains will be returned.
   *
   * Will only return a Node whose distance from \p p is less than
   * \p tol multiplied by the size of a semilocal element which
   * contains \p p.
   */
  virtual const Node *
  locate_node (const Point & p,
               const std::set<subdomain_id_type> * allowed_subdomains = libmesh_nullptr,
               Real tol = TOLERANCE) const;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use, \p false otherwise.
   */
  bool initialized () const;

  /**
   * Enables out-of-mesh mode.  In this mode, if asked to find a point
   * that is contained in no mesh at all, the point locator will
   * return a NULL pointer instead of crashing.  Per default, this
   * mode is off.
   */
  virtual void enable_out_of_mesh_mode () = 0;

  /**
   * Disables out-of-mesh mode (default).  If asked to find a point
   * that is contained in no mesh at all, the point locator will now
   * crash.
   */
  virtual void disable_out_of_mesh_mode () = 0;

  /**
   * Set a tolerance to use when determining
   * if a point is contained within the mesh.
   */
  virtual void set_close_to_point_tol(Real close_to_point_tol);

  /**
   * Specify that we do not want to use a user-specified tolerance to
   * determine if a point is contained within the mesh.
   */
  virtual void unset_close_to_point_tol();

  /**
   * Boolean flag to indicate whether to print out extra info.
   */
  bool _verbose;

protected:
  /**
   * Const pointer to our master, initialized to \p NULL if none
   * given.  When using multiple PointLocators, one can be assigned
   * master and be in charge of something that all can have access to.
   */
  const PointLocatorBase * _master;

  /**
   * constant reference to the mesh in which the point is looked for.
   */
  const MeshBase & _mesh;

  /**
   * \p true when properly initialized, \p false otherwise.
   */
  bool _initialized;

  /**
   * \p true if we will use a user-specified tolerance for locating
   * the element.
   */
  bool _use_close_to_point_tol;

  /**
   * The tolerance to use.
   */
  Real _close_to_point_tol;
};

} // namespace libMesh

#endif // LIBMESH_POINT_LOCATOR_BASE_H
