// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __point_locator_base_h__
#define __point_locator_base_h__

// C++ includes
#include <vector>



// Local Includes
#include "reference_counted_object.h"
#include "libmesh_common.h"
#include "auto_ptr.h"
#include "enum_point_locator_type.h"

namespace libMesh
{



// Forward Declarations
class PointLocatorBase;
class MeshBase;
class Point;
class TreeBase;
class Elem;



/**
 * This is the base class for point locators.  They locate
 * points in space: given a mesh they return the element
 * and local coordinates for a given point in global coordinates.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// PointLocatorBase class definition
class PointLocatorBase : public ReferenceCountedObject<PointLocatorBase>
{
protected:

  /**
   * Constructor.  Protected so that this base class 
   * cannot be explicitly instantiated.  Takes a master
   * PointLocator that helps in saving memory.
   */
  PointLocatorBase (const MeshBase& mesh,
		    const PointLocatorBase* master);


public:

  /**
   * Destructor.
   */
  virtual ~PointLocatorBase ();

  /**
   * Builds an PointLocator for the mesh \p mesh.  
   * Optionally takes a master PointLocator to save memory.
   * An \p AutoPtr<PointLocatorBase> is returned to prevent memory leak. 
   * This way the user need not remember to delete the object.
   */
  static AutoPtr<PointLocatorBase> build (const PointLocatorType t,
					  const MeshBase& mesh,
					  const PointLocatorBase* master = NULL);

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
   * \p p is located.  Pure virtual.  
   */
  virtual const Elem* operator() (const Point& p) const = 0;

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
  virtual void enable_out_of_mesh_mode (void) = 0;

  /**
   * Disables out-of-mesh mode (default).  If asked to find a point
   * that is contained in no mesh at all, the point locator will now
   * crash.
   */
  virtual void disable_out_of_mesh_mode (void) = 0;

protected:

  /**
   * Const pointer to our master, initialized to \p NULL if none
   * given.  When using multiple PointLocators, one can be assigned
   * master and be in charge of something that all can have access to.
   */
  const PointLocatorBase* _master;
  
  /**
   * constant reference to the mesh in which the point is looked for.
   */
  const MeshBase& _mesh;

  /**
   * \p true when properly initialized, \p false otherwise.
   */
  bool _initialized;

};


// ------------------------------------------------------------
// PointLocatorBase inline methods
inline
bool PointLocatorBase::initialized () const
{
  return (this->_initialized);
}

} // namespace libMesh

#endif


