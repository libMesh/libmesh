// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_CENTROID_PARTITIONER_H
#define LIBMESH_CENTROID_PARTITIONER_H

// Local includes
#include "libmesh/partitioner.h"
#include "libmesh/point.h"

// C++ includes
#include <utility> // pair
#include <vector>

namespace libMesh
{


// Forward declarations
class Elem;


/**
 * The centroid partitioner partitions simply based on the
 * locations of element centroids.  You must define what
 * you mean by "less than" for the list of element centroids, e.g.
 * if you only care about distance in the z-direction, you would
 * define "less than" differently than if you cared about radial
 * distance.
 *
 * \author John W. Peterson
 * \author Benjamin S. Kirk
 * \date 2003
 */
class CentroidPartitioner : public Partitioner
{
public:


  /**
   * A typedef which is reserved only for use within
   * this class.  If \p X is chosen, then centroid locations
   * will be sorted according to their X-location, etc...
   */
  enum CentroidSortMethod {X=0,
                           Y,
                           Z,
                           RADIAL,
                           INVALID_METHOD};

  /**
   * Constructor.  Takes the \p CentroidSortMethod to use, which
   * defaults to \p X ordering.
   */
  explicit
  CentroidPartitioner (const CentroidSortMethod sm=X) : _sort_method(sm) {}

  /**
   * Creates a new partitioner of this type and returns it in
   * an \p UniquePtr.
   */
  virtual UniquePtr<Partitioner> clone () const libmesh_override
  {
    return UniquePtr<Partitioner>(new CentroidPartitioner(sort_method()));
  }

  /**
   * Specifies how the elements will be sorted.
   */
  CentroidSortMethod sort_method () const { return _sort_method; }

  /**
   * Change how the elements will be sorted.
   */
  void set_sort_method (const CentroidSortMethod sm) {_sort_method = sm; }


protected:
  /**
   * Partitions the mesh into n subdomains.  This is
   * a required interface for the class.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) libmesh_override;

private:

  /**
   * Computes a list of element centroids for the mesh.
   * This list will be kept around in case a repartition
   * is desired.
   */
  void compute_centroids (MeshBase & mesh);

  /**
   * Partition the list of centroids based on the
   * x-coordinate of the centroid.  This provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_x (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);

  /**
   * Partition the list of centroids based on the
   * y-coordinate of the centroid.  This  provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_y (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);

  /**
   * Partition the list of centroids based on the
   * z-coordinate of the centroid.  This provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_z (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);


  /**
   * Partition the list of centroids based on the
   * radial position of the centroid.  This provides
   * a function which may be passed to the std::sort
   * routine for sorting the elements by centroid.
   */
  static bool sort_radial (const std::pair<Point, Elem *> & lhs,
                           const std::pair<Point, Elem *> & rhs);

  /**
   * Store a flag which tells which type of
   * sort method we are using.
   */
  CentroidSortMethod _sort_method;

  /**
   * Vector which holds pairs of centroids and
   * their respective element pointers.
   */
  std::vector<std::pair<Point, Elem *> > _elem_centroids;
};

} // namespace libMesh


#endif // LIBMESH_CENTROID_PARTITIONER_H
