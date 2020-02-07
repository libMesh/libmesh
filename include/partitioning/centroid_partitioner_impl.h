// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_CENTROID_PARTITIONER_IMPL_H
#define LIBMESH_CENTROID_PARTITIONER_IMPL_H

// C++ includes
#include <algorithm> // for std::sort

// Local includes
#include "libmesh/centroid_partitioner.h"
#include "libmesh/elem.h"
#include "libmesh/int_range.h"

namespace libMesh
{

template <typename RealType>
void CentroidPartitionerTempl<RealType>::partition_range(MeshBase & mesh,
                                          typename MeshBase::element_iterator it,
                                          typename MeshBase::element_iterator end,
                                          unsigned int n)
{
  // Check for easy returns
  if (it == end)
    return;

  if (n == 1)
    {
      this->single_partition_range (it, end);
      return;
    }

  // Make sure the user has not handed us an
  // invalid number of partitions.
  libmesh_assert_greater (n, 0);

  // We don't yet support distributed meshes with this Partitioner
  if (!mesh.is_serial())
    libmesh_not_implemented();

  // Compute the element centroids.  Note: we used to skip this step
  // if the number of elements was unchanged from the last call, but
  // that doesn't account for elements that have moved a lot since the
  // last time the Partitioner was called...
  this->compute_centroids (it, end);

  switch (this->sort_method())
    {
    case X:
      {
        std::sort(_elem_centroids.begin(),
                  _elem_centroids.end(),
                  CentroidPartitionerTempl<RealType>::sort_x);

        break;
      }


    case Y:
      {
        std::sort(_elem_centroids.begin(),
                  _elem_centroids.end(),
                  CentroidPartitionerTempl<RealType>::sort_y);

        break;

      }


    case Z:
      {
        std::sort(_elem_centroids.begin(),
                  _elem_centroids.end(),
                  CentroidPartitionerTempl<RealType>::sort_z);

        break;
      }


    case RADIAL:
      {
        std::sort(_elem_centroids.begin(),
                  _elem_centroids.end(),
                  CentroidPartitionerTempl<RealType>::sort_radial);

        break;
      }
    default:
      libmesh_error_msg("Unknown sort method: " << this->sort_method());
    }

  // Compute target_size, the approximate number of elements on each processor.
  const dof_id_type target_size = cast_int<dof_id_type>
    (_elem_centroids.size() / n);

  for (auto i : index_range(_elem_centroids))
    {
      Elem * elem = _elem_centroids[i].second;

      // FIXME: All "extra" elements go on the last processor... this
      // could probably be improved.
      elem->processor_id() =
        std::min (cast_int<processor_id_type>(i / target_size),
                  cast_int<processor_id_type>(n-1));
    }
}



template <typename RealType>
void CentroidPartitionerTempl<RealType>::_do_partition (MeshBase & mesh,
                                         const unsigned int n)
{
  this->partition_range(mesh,
                        mesh.elements_begin(),
                        mesh.elements_end(),
                        n);
}



template <typename RealType>
void CentroidPartitionerTempl<RealType>::compute_centroids (typename MeshBase::element_iterator it,
                                             typename MeshBase::element_iterator end)
{
  _elem_centroids.clear();

  for (auto & elem : as_range(it, end))
    _elem_centroids.push_back(std::make_pair(elem->centroid(), elem));
}




template <typename RealType>
bool CentroidPartitionerTempl<RealType>::sort_x (const std::pair<Point, Elem *> & lhs,
                                  const std::pair<Point, Elem *> & rhs)
{
  return (lhs.first(0) < rhs.first(0));
}




template <typename RealType>
bool CentroidPartitionerTempl<RealType>::sort_y (const std::pair<Point, Elem *> & lhs,
                                  const std::pair<Point, Elem *> & rhs)
{
  return (lhs.first(1) < rhs.first(1));
}





template <typename RealType>
bool CentroidPartitionerTempl<RealType>::sort_z (const std::pair<Point, Elem *> & lhs,
                                  const std::pair<Point, Elem *> & rhs)
{
  return (lhs.first(2) < rhs.first(2));
}



template <typename RealType>
bool CentroidPartitionerTempl<RealType>::sort_radial (const std::pair<Point, Elem *> & lhs,
                                       const std::pair<Point, Elem *> & rhs)
{
  return (lhs.first.norm() < rhs.first.norm());
}

} // namespace libMesh

#endif // LIBMESH_CENTROID_PARTITIONER_IMPL_H
