// $Id: centroid_partitioner.C,v 1.2 2003-08-23 00:36:08 jwpeterson Exp $

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
#include <algorithm>

// Local includes
#include "centroid_partitioner.h"
#include "elem_iterators.h"
#include "mesh_base.h"








CentroidPartitioner::CentroidPartitioner(MeshBase& mesh,
					 CentroidSortMethod sm) : Partitioner(mesh),
								  _sort_method(sm)
{
  // construct list of centroids here?
}





void CentroidPartitioner::partition (const unsigned int n)
{

  // Possibly reconstruct centroids
  if (_mesh.n_elem() != _elem_centroids.size())
    this->compute_centroids();


  
  switch (_sort_method)
    {
    case X:
      {
	std::sort(_elem_centroids.begin(),
		  _elem_centroids.end(),
		  CentroidPartitioner::sort_x);
	
	break;
      }

      
    case Y:
      {
	std::sort(_elem_centroids.begin(),
		  _elem_centroids.end(),
		  CentroidPartitioner::sort_y);
	
	break;
	
      }

      
    case Z:
      {
	std::sort(_elem_centroids.begin(),
		  _elem_centroids.end(),
		  CentroidPartitioner::sort_z);
	
	break;
      }

      
     case RADIAL:
      {
	std::sort(_elem_centroids.begin(),
		  _elem_centroids.end(),
		  CentroidPartitioner::sort_radial);
	
	break;
      } 
    default:
      error();
    }

  
  // Make sure the user has not handed us an
  // invalid number of partitions.
  assert (n > 0);

  // the number of elements, e.g. 1000
  const unsigned int n_elem      = _mesh.n_elem();
  // the number of elements per processor, e.g 400
  const unsigned int target_size = n_elem / n;

  // Make sure the mesh hasn't changed since the
  // last time we computed the centroids.
  assert (_mesh.n_elem() == _elem_centroids.size());

  for (unsigned int i=0; i<n_elem; i++)
    {
      Elem* elem = _elem_centroids[i].second;

      elem->set_subdomain_id() =
	elem->set_processor_id() =

	std::min (i / target_size, n-1);
      
    }	
}




void CentroidPartitioner::repartition (const unsigned int n)
{
  this->partition (n);
}




void CentroidPartitioner::compute_centroids()
{
  _elem_centroids.clear();
  _elem_centroids.reserve(_mesh.n_elem());
  
  elem_iterator it(_mesh.elements_begin());
  const elem_iterator it_end(_mesh.elements_end());
  for (; it != it_end; ++it)
    {
      Elem* elem = *it;

      _elem_centroids.push_back(std::make_pair(elem->centroid(), elem));
    }
}




bool CentroidPartitioner::sort_x (const std::pair<Point, Elem*>& lhs,
				  const std::pair<Point, Elem*>& rhs)
{
  return (lhs.first(0) < rhs.first(0));
}




bool CentroidPartitioner::sort_y (const std::pair<Point, Elem*>& lhs,
				  const std::pair<Point, Elem*>& rhs)
{
  return (lhs.first(1) < rhs.first(1));
}





bool CentroidPartitioner::sort_z (const std::pair<Point, Elem*>& lhs,
				  const std::pair<Point, Elem*>& rhs)
{
  return (lhs.first(2) < rhs.first(2));
}



bool CentroidPartitioner::sort_radial (const std::pair<Point, Elem*>& lhs,
				       const std::pair<Point, Elem*>& rhs)
{
  return (lhs.first.size() < rhs.first.size());
}
