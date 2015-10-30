// $Id: centroid_partitioner.C,v 1.13 2005-02-22 22:17:42 jwpeterson Exp $

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
#include <algorithm> // for std::sort

// Local includes
#include "centroid_partitioner.h"
#include "mesh_base.h"
#include "elem.h"


//---------------------------------------------------------
// CentroidPartitioner methods
void CentroidPartitioner::_do_partition (MeshBase& mesh,
					 const unsigned int n)
{
  // Check for an easy return
  if (n == 1)
    {
      this->single_partition (mesh);
      return;
    }


  // Possibly reconstruct centroids
  if (mesh.n_elem() != _elem_centroids.size())
    this->compute_centroids (mesh);


  
  switch (this->sort_method())
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
  const unsigned int n_elem      = mesh.n_elem();
  // the number of elements per processor, e.g 400
  const unsigned int target_size = n_elem / n;

  // Make sure the mesh hasn't changed since the
  // last time we computed the centroids.
  assert (mesh.n_elem() == _elem_centroids.size());

  for (unsigned int i=0; i<n_elem; i++)
    {
      Elem* elem = _elem_centroids[i].second;

      elem->processor_id() = std::min (i / target_size, n-1);
    }	
}








void CentroidPartitioner::compute_centroids (MeshBase& mesh)
{
  _elem_centroids.clear();
  _elem_centroids.reserve(mesh.n_elem());
  
//   elem_iterator it(mesh.elements_begin());
//   const elem_iterator it_end(mesh.elements_end());

  MeshBase::element_iterator       it     = mesh.elements_begin();
  const MeshBase::element_iterator it_end = mesh.elements_end(); 

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
