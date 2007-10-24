// $Id: point_locator_tree.h,v 1.6 2005-02-22 22:17:35 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __point_locator_tree_h__
#define __point_locator_tree_h__

// C++ includes



// Local Includes
#include "point_locator_base.h"



// Forward Declarations
class Mesh;
class Point;
class TreeBase;
class Elem;


/**
 * This is a point locator.  It locates points in space
 * using a tree: given a mesh they return the element
 * and local coordinates for a given point in global coordinates.
 * Use \p PointLocatorBase::build() to create objects of this 
 * type at run time.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// PointLocatorTree class definition
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
  PointLocatorTree (const Mesh& mesh,
		    const PointLocatorBase* master = NULL);


public:

  /**
   * Destructor.
   */
  ~PointLocatorTree ();

  /**
   * Clears the locator.  This function frees dynamic memory with "delete".
   */
  virtual void clear();

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used.  This function allocates dynamic memory with "new".
   */
  virtual void init();

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located.  The mutable _element member is used to cache
   * the result and allow it to be used during the next call to
   * operator().
   */
  virtual const Elem* operator() (const Point& p) const;

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
  mutable const Elem* _element;

};


// ------------------------------------------------------------
// PointLocatorTree inline methods


#endif


