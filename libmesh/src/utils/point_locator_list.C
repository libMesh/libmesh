// $Id: point_locator_list.C,v 1.2 2003-07-25 20:58:24 benkirk Exp $

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
#include "point_locator_list.h"
#include "mesh.h"
#include "elem_iterators.h"
//#include "point.h"
//#include "elem.h"


// typedefs
typedef std::vector<Point>::const_iterator   const_list_iterator;


//------------------------------------------------------------------
// PointLocator methods
PointLocatorList::PointLocatorList (const Mesh& mesh,
				    const PointLocatorBase* master) :
  PointLocatorBase (mesh, 
		    master),
  _list            (NULL)
{
}




PointLocatorList::~PointLocatorList ()
{
  this->clear ();
}




void PointLocatorList::clear ()
{
  // only delete the list when we are the master
  if (this->_list != NULL)
    {
      if (this->_master == NULL)
        {
	  // we own the list
	  this->_list->clear();
	  delete this->_list;
	}
      else
	  // someone else owns and therefore deletes the list
	  this->_list = NULL;
    }
}





void PointLocatorList::init ()
{
  assert (this->_list == NULL); 

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
	   * We are the master, so we have to build the list.
	   * First create it, then get a handy reference, and
	   * then try to speed up by reserving space...
	   */
	  this->_list = new std::vector<Point>;
	  std::vector<Point>& my_list = *(this->_list);

	  my_list.clear();
	  my_list.reserve(this->_mesh.n_elem());


	  /*
	   * fill our list with the centroids and element
	   * pointers of the mesh.  For this use the handy
	   * element iterators.
	   */
	  const_active_elem_iterator       el (this->_mesh.elements_begin());
	  const const_active_elem_iterator end(this->_mesh.elements_end()); 

	  for (; el!=end; ++el)
	      my_list.push_back((*el)->centroid());
	}

      else
	  
        {
	  /*
	   * We are _not_ the master.  Let our _list point to
	   * the master's list.  But for this we first transform
	   * the master in a state for which we are friends
	   * (this should also beware of a bad master pointer?).
	   * And make sure the master @e has a list!
	   */
	  const PointLocatorList* my_master =
	    dynamic_cast<const PointLocatorList*>(this->_master);

	  if (my_master->initialized())
	    this->_list = my_master->_list;
	  else
	    {
	      std::cerr << "ERROR: Initialize master first, then servants!"
			<< std::endl;
	      error();
	    }
        }

    }


  // ready for take-off
  this->_initialized = true;
}





const Elem* PointLocatorList::operator() (const Point& p)
{
  assert (this->_initialized);

  /*
   * Ask the list.  This is quite expensive, since
   * we loop through the whole list to try to find
   * the @e nearest element.
   * However, there is not much other to do: when
   * we would use bounding boxes like in a tree,
   * it may happen that a surface element is just
   * in plane with a bounding box face, and quite
   * close to it.  But when a point comes, this
   * point may belong to the bounding box (where the
   * coplanar element does @e not belong to).  Then
   * we would search through the elements in this 
   * bounding box, while the other bounding box'es
   * element is closer, but we simply don't consider
   * it!
   *
   * We _can_, however, use size_sq() instead of size()
   * here to avoid repeated calls to sqrt(), which is
   * pretty expensive.
   */
  {
    std::vector<Point>& my_list = *(this->_list);

    Real               last_distance_sq = Point(my_list[0] -p).size_sq();
    unsigned int       last_index       = 0;
    const unsigned int max_index        = my_list.size();


    for (unsigned int n=1; n<max_index; n++)
      {
	const Real current_distance_sq = Point(my_list[n] -p).size_sq();

	if (current_distance_sq < last_distance_sq)
	  {
	    last_distance_sq = current_distance_sq;
	    last_index       = n;
	  }
      }


//    std::cout << "Found element No. "<< last_index << std::endl;
//	      << std::endl << "with key: " << this->_mesh.elem(last_index).key() << std::endl;

    /*
     * return the element
     */
    return (this->_mesh.elem(last_index));

  }

}

