// $Id: point_locator.C,v 1.2 2003-05-15 23:34:36 benkirk Exp $

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



// C++ includes


// Local Includes
#include "point_locator.h"
#include "mesh.h"
#include "point.h"
#include "tree_base.h"
#include "elem.h"




//------------------------------------------------------------------
// PointLocator methods
template <PointLocatorType T>
PointLocator<T>::PointLocator (const Mesh& mesh,
			       const PointLocatorBase* master) :
  PointLocatorBase (mesh, 
		    master),
  _tree            (NULL),
  _element         (NULL)
{
}




template <PointLocatorType T>
PointLocator<T>::~PointLocator ()
{
  this->clear ();
}




template <PointLocatorType T>
void PointLocator<T>::clear ()
{
  // only delete the tree when we are the master
  if (this->_tree != NULL)
    {
      if (this->_master == NULL)
	  // we own the tree
	  delete this->_tree;
      else
	  // someone else owns and therefore deletes the tree
	  this->_tree = NULL;
    }
}





template <PointLocatorType T>
void PointLocator<T>::init ()
{
  assert (this->_tree == NULL); 

  if (this->_initialized)
    {
      std::cerr << "ERROR: Already initialized!  Will ignore this call..."
		<< std::endl;
    }

  else

    {

      if (this->_master == NULL)
        {
	  /*
	   * We are the master, so we have to build the tree
	   */

	  switch (this->_mesh.mesh_dimension())
	    { 
	    case 2:
	      {
//TODO: What to do with the level of the tree?
		AutoPtr<TreeBase> ap(TreeBase::build(4, 
						     this->_mesh, 
						     100));
		_tree = ap.release();
		break;
	      }

	    case 3:
	      {
		AutoPtr<TreeBase> ap(TreeBase::build(8, 
						     this->_mesh, 
						     100));
		_tree = ap.release();
		break;
	      }

	    default:
	      {
		std::cerr << "ERROR: Bad dimension = " 
			  << this->_mesh.mesh_dimension() 
			  << std::endl;
		error();
	      }
	    }	

	}

      else
	  
        {
	  /*
	   * We are _not_ the master.  Let our Tree point to
	   * the master's tree.  But for this we first transform
	   * the master in a state for which we are friends
	   */
	  const PointLocator<T>* my_master =
	    dynamic_cast<const PointLocator<T>*>(this->_master);

	  this->_tree = my_master->_tree;
        }


      /*
       * Not all PointLocators may own a tree, but all of them
       * use their own element pointer.  Let the element pointer
       * be unique for every interpolator.
       * Suppose the interpolators are used concurrently
       * at different locations in the mesh, then it makes quite
       * sense to have unique start elements.
       */
      this->_element = this->_mesh.elem(0);

    }


  // ready for take-off
  this->_initialized = true;
}





template <PointLocatorType T>
const Elem* PointLocator<T>::operator() (const Point& p)
{
  assert (this->_initialized);

  /*
   * Check first the element from last time before asking the tree
   */
  if (!(this->_element->contains_point(p)))
    {
	// ask the tree
	this->_element = this->_tree->find_element (p);

	if (this->_element == NULL)
	  {
	    std::cerr << std::endl
		      << " ******** Serious Problem.  Could not find an Element "
		      << "in the Mesh" 
		      << std:: endl
		      << " ******** that contains the Point ";
	    p.print();
	    error();
	  }
    }

  /*
   * return the element
   */

  return this->_element;
}



// ------------------------------------------------------------
// Explicit Instantiations
template class PointLocator<TREE>;





