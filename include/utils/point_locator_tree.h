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



#ifndef LIBMESH_POINT_LOCATOR_TREE_H
#define LIBMESH_POINT_LOCATOR_TREE_H

// Local Includes
#include "libmesh/point_locator_base.h"
#include "libmesh/tree_base.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Point;
class Elem;

/**
 * This is a point locator.  It locates points in space
 * using a tree: given a mesh they return the element
 * and local coordinates for a given point in global coordinates.
 * Use \p PointLocatorBase::build() to create objects of this
 * type at run time.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class PointLocatorTree : public PointLocatorBase
{
public:
  /**
   * Constructor.  Needs the \p mesh in which the points
   * should be located.  Optionally takes a master
   * interpolator.  This master helps in saving memory
   * by reducing the number of trees in use.  Only the
   * master locator holds a  tree, the others simply
   * use the master's tree.
   */
  PointLocatorTree (const MeshBase & mesh,
                    const PointLocatorBase * master = libmesh_nullptr);


  /**
   * Constructor.  Needs the \p mesh in which the points
   * should be located.  Allows the user to specify the
   * method to use when building the tree.
   * Optionally takes a master interpolator.
   * This master helps in saving memory
   * by reducing the number of trees in use.  Only the
   * master locator holds a  tree, the others simply
   * use the master's tree. Allows the user to specify
   * the build type.
   */
  PointLocatorTree (const MeshBase & mesh,
                    const Trees::BuildType build_type,
                    const PointLocatorBase * master = libmesh_nullptr);

  /**
   * Destructor.
   */
  ~PointLocatorTree ();

  /**
   * Clears the locator.  This function frees dynamic memory with "delete".
   */
  virtual void clear() libmesh_override;

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used.  This function allocates dynamic memory with "new".
   */
  void init(Trees::BuildType build_type);

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used.  This function allocates dynamic memory with "new".
   */
  virtual void init() libmesh_override;

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located, optionally restricted to a set of allowed subdomains.
   * The mutable _element member is used to cache
   * the result and allow it to be used during the next call to
   * operator().
   */
  virtual const Elem * operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains = libmesh_nullptr) const libmesh_override;

  /**
   * As a fallback option, it's helpful to be able to do a linear
   * search over the entire mesh. This can be used if operator()
   * fails to find an element that contains \p p, for example.
   * Optionally specify a "close to point" tolerance to use in
   * the linear search.
   * Return NULL if no element is found.
   */
  const Elem * perform_linear_search(const Point & p,
                                     const std::set<subdomain_id_type> * allowed_subdomains,
                                     bool use_close_to_point,
                                     Real close_to_point_tolerance=TOLERANCE) const;

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
   * Set the target bin size.
   */
  void set_target_bin_size(unsigned int target);

  /**
   * Get the target bin size.
   */
  unsigned int get_target_bin_size() const;

protected:
  /**
   * Pointer to our tree.  The tree is built at run-time
   * through \p init().  For servant PointLocators (not master),
   * this simply points to the tree of the master.
   */
  TreeBase * _tree;

  /**
   * Pointer to the last element that was found by the tree.
   * Chances are that this may be close to the next call to
   * \p operator()...
   */
  mutable const Elem * _element;

  /**
   * \p true if out-of-mesh mode is enabled.  See \p
   * enable_out_of_mesh_mode() for details.
   */
  bool _out_of_mesh_mode;

  /**
   * Target bin size, which gets passed to the constructor of _tree.
   */
  unsigned int _target_bin_size;

  /**
   * How the underlying tree is built.
   */
  Trees::BuildType _build_type;
};

} // namespace libMesh

#endif // LIBMESH_POINT_LOCATOR_TREE_H
