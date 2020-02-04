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


#ifndef LIBMESH_CENTROID_PARTITIONER_H
#define LIBMESH_CENTROID_PARTITIONER_H

// Local includes
#include "libmesh/partitioner.h"
#include "libmesh/point.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <utility> // pair
#include <vector>

namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;

/**
 * Partitions the Mesh based on the locations of element centroids.
 * You must define what you mean by "less than" for the list of
 * element centroids, e.g. if you only care about distance in the
 * z-direction, you would define "less than" differently than if you
 * cared about radial distance.
 *
 * \author John W. Peterson
 * \author Benjamin S. Kirk
 * \date 2003
 */
class CentroidPartitioner : public Partitioner
{
public:

  /**
   * A typedef which controls the sorting method used for ordering the
   * centroids. If e.g. \p X is chosen, then the centroids will be
   * sorted according to their x-coordinate.
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
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  CentroidPartitioner (const CentroidPartitioner &) = default;
  CentroidPartitioner (CentroidPartitioner &&) = default;
  CentroidPartitioner & operator= (const CentroidPartitioner &) = default;
  CentroidPartitioner & operator= (CentroidPartitioner &&) = default;
  virtual ~CentroidPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return libmesh_make_unique<CentroidPartitioner>(*this);
  }

  /**
   * Getter for the current sorting method.
   */
  CentroidSortMethod sort_method () const { return _sort_method; }

  /**
   * Setter for the current sorting method.
   */
  void set_sort_method (const CentroidSortMethod sm) { _sort_method = sm; }

  /**
   * Called by the SubdomainPartitioner to partition elements in the range (it, end).
   */
  virtual void partition_range(MeshBase & mesh,
                               MeshBase::element_iterator it,
                               MeshBase::element_iterator end,
                               const unsigned int n) override;

protected:

  /**
   * Partitions the mesh into n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override;

private:

  /**
   * Computes a list of element centroids for the mesh.
   */
  void compute_centroids (MeshBase::element_iterator it,
                          MeshBase::element_iterator end);

  /**
   * Helper function which sorts by the centroid's x-coordinate in the
   * internal std::sort call.
   */
  static bool sort_x (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);

  /**
   * Helper function which sorts by the centroid's y-coordinate in the
   * internal std::sort call.
   */
  static bool sort_y (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);

  /**
   * Helper function which sorts by the centroid's z-coordinate in the
   * internal std::sort call.
   */
  static bool sort_z (const std::pair<Point, Elem *> & lhs,
                      const std::pair<Point, Elem *> & rhs);

  /**
   * Helper function which sorts by the centroid's distance from the
   * origin in the internal std::sort call.
   */
  static bool sort_radial (const std::pair<Point, Elem *> & lhs,
                           const std::pair<Point, Elem *> & rhs);

  /**
   * Flag indicating the type of sort method we are using.
   */
  CentroidSortMethod _sort_method;

  /**
   * Vector which holds pairs of centroids and their respective
   * element pointers.
   */
  std::vector<std::pair<Point, Elem *>> _elem_centroids;
};

} // namespace libMesh

#endif // LIBMESH_CENTROID_PARTITIONER_H
