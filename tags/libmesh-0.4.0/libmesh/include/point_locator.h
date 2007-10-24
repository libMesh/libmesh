// $Id: point_locator.h,v 1.2 2003-05-15 23:34:34 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __point_locator_h__
#define __point_locator_h__

// C++ includes



// Local Includes
#include "point_locator_base.h"
#include "enum_point_locator_type.h"



// Forward Declarations
class Mesh;
class Point;
class TreeBase;
class Elem;
template <PointLocatorType T> class PointLocator;


/**
 * This is a point locator class.  it locates
 * points in space: given a mesh they return the element
 * and local coordinates for a given point in global coordinates.
 * Currently, only the \p TREE version of a point locator is
 * available.  Use \p PointLocatorBase::build() to
 * create objects of this type at run time.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// PointLocator class definition
template <PointLocatorType T>
class PointLocator : public PointLocatorBase
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
  PointLocator (const Mesh& mesh,
		const PointLocatorBase* master = NULL);


public:

  /**
   * Destructor.
   */
  ~PointLocator ();

  /**
   * Clears the locator.
   */
  void clear();

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used.
   */
  void init();

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located.  This cannot be a const method, because the
   * resulting element is remembered for the next request.
   */
  virtual const Elem* operator() (const Point& p);

protected:

  /**
   * Pointer to our tree.  The tree is built at run-time
   * through \p init().  For servant PointLocators (not master),
   * this simply points to the tree of the master.
   */
  TreeBase* _tree;

  /**
   * Pointer to the last element that was found by the tree.
   * Chances are that this may be close to the next call to
   * \p operator()...
   */
  Elem* _element;

private:

  /**
   * Make all \p PointLocator<TREE> friend so that
   * the servants can access the master's tree.
   */
//  friend class PointLocator<TREE>;
};


// ------------------------------------------------------------
// PointLocator inline methods


#endif


