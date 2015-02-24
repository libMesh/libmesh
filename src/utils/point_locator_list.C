// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point_locator_list.h"

namespace libMesh
{


// typedefs
typedef std::vector<Point>::const_iterator   const_list_iterator;


//------------------------------------------------------------------
// PointLocator methods
PointLocatorList::PointLocatorList (const MeshBase& mesh,
                                    const PointLocatorBase* master) :
  PointLocatorBase (mesh,master),
  _list            (NULL)
{
  // This code will only work if your mesh is the Voroni mesh of it's
  // own elements' centroids.  If your mesh is that regular you might
  // as well hand-code an O(1) algorithm for locating points within
  // it. - RHS
  libmesh_experimental();

  this->init();
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
  libmesh_assert (!this->_list);

  if (this->_initialized)
    libMesh::err << "ERROR: Already initialized!  Will ignore this call..." << std::endl;

  else
    {
      if (this->_master == NULL)
        {
          START_LOG("init(no master)", "PointLocatorList");

          // We are the master, so we have to build the list.
          // First create it, then get a handy reference, and
          // then try to speed up by reserving space...
          this->_list = new std::vector<std::pair<Point, const Elem *> >;
          std::vector<std::pair<Point, const Elem *> >& my_list = *(this->_list);

          my_list.clear();
          my_list.reserve(this->_mesh.n_active_elem());

          // fill our list with the centroids and element
          // pointers of the mesh.  For this use the handy
          // element iterators.
          MeshBase::const_element_iterator       el  = _mesh.active_elements_begin();
          const MeshBase::const_element_iterator end = _mesh.active_elements_end();

          for (; el!=end; ++el)
            my_list.push_back(std::make_pair((*el)->centroid(), *el));

          STOP_LOG("init(no master)", "PointLocatorList");
        }

      else
        {
          // We are _not_ the master.  Let our _list point to
          // the master's list.  But for this we first transform
          // the master in a state for which we are friends
          // (this should also beware of a bad master pointer?).
          // And make sure the master @e has a list!
          const PointLocatorList* my_master =
            cast_ptr<const PointLocatorList*>(this->_master);

          if (my_master->initialized())
            this->_list = my_master->_list;
          else
            libmesh_error_msg("ERROR: Initialize master first, then servants!");
        }
    }

  // ready for take-off
  this->_initialized = true;
}



const Elem* PointLocatorList::operator() (const Point& p, const std::set<subdomain_id_type> *allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  START_LOG("operator()", "PointLocatorList");

  // Ask the list.  This is quite expensive, since
  // we loop through the whole list to try to find
  // the @e nearest element.
  // However, there is not much else to do: when
  // we would use bounding boxes like in a tree,
  // it may happen that a surface element is just
  // in plane with a bounding box face, and quite
  // close to it.  But when a point comes, this
  // point may belong to the bounding box (where the
  // coplanar element does @e not belong to).  Then
  // we would search through the elements in this
  // bounding box, while the other bounding box'es
  // element is closer, but we simply don't consider
  // it!
  //
  // We _can_, however, use size_sq() instead of size()
  // here to avoid repeated calls to std::sqrt(), which is
  // pretty expensive.
  {
    std::vector<std::pair<Point, const Elem *> >& my_list = *(this->_list);

    Real              last_distance_sq = std::numeric_limits<Real>::max();
    const Elem *      last_elem        = NULL;
    const std::size_t max_index        = my_list.size();

    for (std::size_t n=0; n<max_index; n++)
      {
        // Only consider elements in the allowed_subdomains list, if it exists
        if (!allowed_subdomains ||
            allowed_subdomains->count(my_list[n].second->subdomain_id()))
          {
            const Real current_distance_sq = Point(my_list[n].first -p).size_sq();

            if (current_distance_sq < last_distance_sq)
              {
                last_distance_sq = current_distance_sq;
                last_elem        = my_list[n].second;
              }
          }
      }

    // If we found an element, it should be active
    libmesh_assert (!last_elem || last_elem->active());

    // If we found an element and have a restriction list, they better match
    libmesh_assert (!last_elem || !allowed_subdomains || allowed_subdomains->count(last_elem->subdomain_id()));

    STOP_LOG("operator()", "PointLocatorList");

    // return the element
    return last_elem;
  }
}



void PointLocatorList::enable_out_of_mesh_mode ()
{
  // This functionality is not yet implemented for PointLocatorList.
  libmesh_not_implemented();
}



void PointLocatorList::disable_out_of_mesh_mode ()
{
  // This functionality is not yet implemented for PointLocatorList.
  libmesh_not_implemented();
}


void PointLocatorList::set_close_to_point_tol (Real /*close_to_point_tol*/)
{
  // This functionality is not yet implemented for PointLocatorList.
  libmesh_not_implemented();
}


void PointLocatorList::unset_close_to_point_tol ()
{
  // This functionality is not yet implemented for PointLocatorList.
  libmesh_not_implemented();
}

} // namespace libMesh
