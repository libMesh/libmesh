// $Id: point_locator_tree.C,v 1.11 2005-02-22 22:17:43 jwpeterson Exp $

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



// C++ includes

// Local Includes
#include "mesh.h"
#include "elem.h"
#include "tree.h"
#include "point_locator_tree.h"




//------------------------------------------------------------------
// PointLocator methods
PointLocatorTree::PointLocatorTree (const Mesh& mesh,
				    const PointLocatorBase* master) :
  PointLocatorBase (mesh,master),
  _tree            (NULL),
  _element         (NULL)
{
  this->init();
}




PointLocatorTree::~PointLocatorTree ()
{
  this->clear ();
}




void PointLocatorTree::clear ()
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





void PointLocatorTree::init ()
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
	  // We are the master, so we have to build the tree
	  switch (this->_mesh.spatial_dimension())
	    { 
	    case 2:
	      {
//TODO: What to do with the level of the tree?
		
		_tree = new Trees::QuadTree (this->_mesh, 100,
					     Trees::QuadTree::ELEMENTS);
		break;
	      }

	    case 3:
	      {
		_tree = new Trees::OctTree (this->_mesh, 100,
					    Trees::OctTree::ELEMENTS);
		break;
	      }

	    default:
	      {
		std::cerr << "ERROR: Bad dimension = " 
			  << this->_mesh.spatial_dimension() 
			  << std::endl;
		error();
	      }
	    }	
	}

      else
	  
        {
	  // We are _not_ the master.  Let our Tree point to
	  // the master's tree.  But for this we first transform
	  // the master in a state for which we are friends.
	  // And make sure the master @e has a tree!
	  const PointLocatorTree* my_master =
	    dynamic_cast<const PointLocatorTree*>(this->_master);

	  if (my_master->initialized())
	    this->_tree = my_master->_tree;
	  else
	    {
	      std::cerr << "ERROR: Initialize master first, then servants!"
			<< std::endl;
	      error();
	    }
        }


      // Not all PointLocators may own a tree, but all of them
      // use their own element pointer.  Let the element pointer
      // be unique for every interpolator.
      // Suppose the interpolators are used concurrently
      // at different locations in the mesh, then it makes quite
      // sense to have unique start elements.
      this->_element = this->_mesh.elem(0);
    }


  // ready for take-off
  this->_initialized = true;
}





const Elem* PointLocatorTree::operator() (const Point& p) const
{
  assert (this->_initialized);
  
  // First check the element from last time before asking the tree
  if (!(this->_element->contains_point(p)))
    {
	// ask the tree
	this->_element = this->_tree->find_element (p);

	// Note that in some cases the tree may not find a point
	// e.g. when a point is located in an element that is not
	// entirely bounded by the tree node's bounding box.
	// In those cases we take a safe but slow way.
	if (this->_element == NULL)
	  {
// 	    const_elem_iterator pos           (this->_mesh.const_elements_begin());
// 	    const const_elem_iterator end_pos (this->_mesh.const_elements_end());

	    MeshBase::const_element_iterator       pos     = this->_mesh.elements_begin();
	    const MeshBase::const_element_iterator end_pos = this->_mesh.elements_end();

	    for ( ; pos != end_pos; ++pos)
	      if ((*pos)->contains_point(p))
		return this->_element = (*pos);
	  }

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

  // return the element
  return this->_element;
}

