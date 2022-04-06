// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_POINT_LOCATOR_NANOFLANN_H
#define LIBMESH_POINT_LOCATOR_NANOFLANN_H

#include "libmesh/libmesh_config.h"

// This class is not defined unless libmesh is built with Nanoflann support
#ifdef LIBMESH_HAVE_NANOFLANN

// libmesh includes
#include "libmesh/point_locator_base.h"
#include "libmesh/point.h"
#include "libmesh/bounding_box.h"

// contrib includes
#include "libmesh/nanoflann.hpp"

// C++ includes
#include <vector>
#include <memory>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Point;
class Elem;

/**
 * This is a PointLocator that uses Nanoflann for its implementation.
 * Nanoflann is distributed with libmesh (see: contrib/nanoflann) and
 * libmesh must be built with nanoflann enabled for this class to
 * work.
 *
 * This PointLocator is still considered "experimental", so it is not
 * enabled in libMesh by default. In particular, there may be false
 * negatives (i.e. failures to find a containing Elem for a given
 * Point even though there actually is one in the Mesh) on adaptively
 * refined meshes, meshes with intersecting mixed-dimension manifolds,
 * and meshes with large aspect ratio changes across neighboring
 * elements. That said, the Nanoflann PointLocator did successfully
 * pass all the MOOSE CI testing that we threw at it (and even
 * uncovered some bugs in said tests) so it will likely work well for
 * most applications outside of the categories mentioned above.
 *
 * You must configure libmesh with --enable-nanoflann-pointlocator in
 * order to use this class.
 *
 * \author John W. Peterson
 * \date 2021
 */
class PointLocatorNanoflann : public PointLocatorBase
{
public:
  /**
   * Constructor. Needs the \p mesh in which the points should be
   * located. Optionally takes a pointer to a "master" PointLocator
   * object. If non-nullptr, this object simply forwards its calls
   * onto the master, so we can have multiple pointers that use the
   * same Nanoflann KD-Tree data structure.
   */
  PointLocatorNanoflann (const MeshBase & mesh,
                         const PointLocatorBase * master = nullptr);

  /**
   * Destructor.
   */
  virtual ~PointLocatorNanoflann ();

  /**
   * Restore to PointLocator to a just-constructed state.
   */
  virtual void clear() override final;

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used. This function allocates dynamic memory with "new".
   */
  virtual void init() override final;

  /**
   * Locates the element in which the point with global coordinates \p
   * p is located, optionally restricted to a set of allowed
   * subdomains.
   */
  virtual const Elem * operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override final;

  /**
   * Locates a set of elements in proximity to the point with global
   * coordinates \p p. Optionally allows the user to restrict the
   * subdomains searched.  The idea here is that if a Point lies on
   * the boundary between two elements, so that it is "in" both
   * elements (to within some tolerance) then the candidate_elements
   * set will contain both of these elements instead of just picking
   * one or the other.
   */
  virtual void operator() (const Point & p,
                           std::set<const Elem *> & candidate_elements,
                           const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override final;

  /**
   * Enables out-of-mesh mode. In this mode, if a searched-for Point
   * is not contained in any element of the Mesh, return nullptr
   * instead of throwing an error. By default, this mode is off.
   */
  virtual void enable_out_of_mesh_mode () override final;

  /**
   * Disables out-of-mesh mode (default). See above.
   */
  virtual void disable_out_of_mesh_mode () override final;

  /**
   * Set/get the number of results returned by each Nanoflann
   * findNeighbors() search.  If a containing Elem (to within
   * _contains_point_tol) is not found after linear searching the
   * first _num_results Elems returned by Nanoflann, then we
   * give up and return nullptr (or throw an error if not in
   * out-of-mesh-mode).
   *
   * Remarks:
   * 1.) Although we do a linear search through the Nanoflann results,
   * the Nanoflann results are always sorted by distance to the
   * searched-for Point, so it's likely that the containing Elem will
   * be found at the beginning of the linear search rather than the
   * end. For nicely shaped elements, the containing Elem is often the
   * first or second one found in the Nanoflann results.
   *
   * 2.) I'm not sure about the relative cost of requesting more
   * results from the Nanoflann search than one actually needs, but
   * presumably reducing _num_results will result in better
   * performance for your particular application up to a point. If, on
   * the other hand, _num_results is too small, then you risk
   * having a false negative result, so take this into account when
   * choosing the parameter for a particular application.
   */
  std::size_t get_num_results() const;
  void set_num_results(std::size_t val);

  //
  // Required Nanoflann typedefs and APIs
  //

  /**
   * Floating point type used for storing coordinates
   */
  typedef Real coord_t;

  /**
   * Must return the number of data points
   */
  std::size_t kdtree_get_point_count() const;

  /**
   * Returns the squared distance between the vector (p1[0], p1[1], p1[2])
   * and the vertex average of Elem "idx_p2" stored in _mesh
   */
  coord_t kdtree_distance(const coord_t * p1,
                          const std::size_t idx_p2,
                          std::size_t size) const;

  /**
   * Returns the dim'th component of the idx'th vertex average.
   */
  coord_t kdtree_get_pt(const std::size_t idx, int dim) const;

  /**
   * Optional bounding-box computation: return false to default to a
   * standard bbox computation loop.  Return true if the BBOX can
   * be computed more efficiently (and returned in "bb") than the
   * standard bounding box computation. The BBOX template parameter
   * must at least support bb.size() to find out the expected
   * dimensionality of the box (e.g. 2 or 3 for point clouds).
   */
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /*bb*/) const { return false; }

protected:

  /**
   * \p true if out-of-mesh mode is enabled.  See \p
   * enable_out_of_mesh_mode() for details.
   */
  bool _out_of_mesh_mode;

  /**
   * The number of results returned by Nanoflann when operator() is
   * called.
   */
  std::size_t _num_results;

  /**
   * Lists of Points and ids which make up the "point cloud" that is
   * to be searched via Nanoflann. The point cloud can be comprised of
   * Elem vertex averages or mesh Nodes, depending on the type of tree
   * created. We keep two separate vectors since the Points in the
   * cloud may not be numbered contiguously in general.
   *
   * These are shared_ptrs to vectors since, if we are not the "master"
   * PointLocator, they need to point at the master's vectors instead.
   */
  std::shared_ptr<std::vector<const Elem *>> _elems;
  std::shared_ptr<std::vector<Point>> _point_cloud;

  // kd_tree will be initialized during init() and then automatically
  // cleaned up by the destructor. We always create a LIBMESH_DIM
  // dimensional Nanoflann object.
  typedef nanoflann::L2_Simple_Adaptor<Real, PointLocatorNanoflann> adapter_t;
  typedef nanoflann::KDTreeSingleIndexAdaptor<adapter_t, PointLocatorNanoflann, LIBMESH_DIM> kd_tree_t;
  std::shared_ptr<kd_tree_t> _kd_tree;

  /**
   * The operator() functions on PointLocator-derived classes are
   * const, so to avoid re-allocating these result data structures every
   * time operator() is called, they have to be mutable.
   */
  mutable std::vector<std::size_t> _ret_index;
  mutable std::vector<Real> _out_dist_sqr;

  /**
   * Vector of indices used by indirect sort. Stored as a class member
   * so we don't have to allocate it from scratch for every search.
   */
  mutable std::vector<std::size_t> _b;

  /**
   * Helper function that wraps the call to the KDTree's
   * findNeighbors() routine.  Must be passed the Point to search for
   * and the number of results to return. Stores the results of the
   * search in the _ret_index and _out_dist_sqr class members and
   * returns a KNNResultSet, which has pointers to the index and
   * distance data.
   */
  nanoflann::KNNResultSet<Real>
  kd_tree_find_neighbors(const Point & p,
                         std::size_t num_results) const;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_NANOFLANN
#endif // LIBMESH_POINT_LOCATOR_NANOFLANN_H
