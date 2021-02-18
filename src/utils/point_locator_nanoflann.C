// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#include "libmesh/libmesh_config.h"

// This class is not defined unless libmesh is built with Nanoflann support
#ifdef LIBMESH_HAVE_NANOFLANN

// libmesh includes
#include "libmesh/point_locator_nanoflann.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/int_range.h" // make_range
#include "libmesh/mesh_tools.h" // create_local_bounding_box

// C++ includes
#include <array>

namespace libMesh
{

PointLocatorNanoflann::PointLocatorNanoflann (const MeshBase & mesh,
                                              const PointLocatorBase * master) :
  PointLocatorBase (mesh, master),
  _out_of_mesh_mode(false),
  _num_results(8)
{
  this->init();
}



PointLocatorNanoflann::~PointLocatorNanoflann () = default;



void
PointLocatorNanoflann::clear ()
{
  this->_initialized = false;
  this->_out_of_mesh_mode = false;

  // reset() actually frees the memory if we are master, otherwise it
  // just reduces the ref. count.
  _ids.reset();
  _point_cloud.reset();
  _local_bbox.reset();
  _kd_tree.reset();
}



void
PointLocatorNanoflann::init ()
{
  LOG_SCOPE("init()", "PointLocatorNanoflann");

  if (!_initialized)
    {
      // If _master == nullptr, then we _are_ the master, and thus
      // responsible for initializing.
      bool we_are_master = (_master == nullptr);

      // If we are the master PointLocator, fill in the _point_cloud
      // data structure with active, local element centroids.
      if (we_are_master)
        {
          _ids = std::make_shared<std::vector<dof_id_type>>();
          _point_cloud = std::make_shared<std::vector<Point>>();

          // Make the KD-Tree out of mesh element centroids.

          // We can either reserve exactly the right amount of space or
          // let push_back() take care of it, not sure what would be
          // faster actually, since it takes some time to count the number
          // of active+local elements.
          auto n_active_local_elem = _mesh.n_active_local_elem();
          _ids->reserve(n_active_local_elem);
          _point_cloud->reserve(n_active_local_elem);

          for (const auto & elem : _mesh.active_local_element_ptr_range())
            {
              _ids->push_back(elem->id());
              _point_cloud->push_back(elem->centroid());
            }

          // Construct the KD-Tree
          _kd_tree = std::make_shared<kd_tree_t>
            (LIBMESH_DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));

          _kd_tree->buildIndex();

          // First, create a BoundingBox for the local elements
          _local_bbox = std::make_shared<BoundingBox>(MeshTools::create_local_bounding_box(_mesh));
        }
      else // we are not master
        {
          // Cast the _master to the appropriate type, and point to
          // its data arrays.
          const auto my_master =
            cast_ptr<const PointLocatorNanoflann *>(this->_master);

          // Point our data structures at the master's
          _ids = my_master->_ids;
          _point_cloud = my_master->_point_cloud;
          _local_bbox = my_master->_local_bbox;
          _kd_tree = my_master->_kd_tree;
        }

      // We are initialized now
      this->_initialized = true;
    }
}

nanoflann::KNNResultSet<Real>
PointLocatorNanoflann::kd_tree_find_neighbors(const Point & p,
                                              std::size_t num_results) const
{
  // We are searching for the Point(s) closest to Point p.
  //
  // TODO: The kd_tree's findNeighbors() routine needs a pointer to
  // Real of length LIBMESH_DIM. It might be convenient if libMesh
  // Points had a data() member that provided this, for now we just
  // copy the coordinates into a std::array of the appropriate size.
  std::array<Real, LIBMESH_DIM> query_pt;
  for (int i=0; i<LIBMESH_DIM; ++i)
    query_pt[i] = p(i);

  // Allocate storage for the indices and distances
  _ret_index.resize(num_results);
  _out_dist_sqr.resize(num_results);

  // nanoflann::KNNResultSet cannot be resized/reused easily, I think
  // it is just meant to be re-created for each search.
  nanoflann::KNNResultSet<Real> result_set(num_results);

  // Initialize the result_set
  result_set.init(_ret_index.data(), _out_dist_sqr.data());

  // Do the search
  // We leave all the SearchParams() ctor args on their default values
  // (note that the first arg is ignored anyway).
  _kd_tree->findNeighbors(result_set, query_pt.data(), nanoflann::SearchParams());

  return result_set;
}



bool
PointLocatorNanoflann::search_local_bbox(const Point & p) const
{
  // If we are using a custom contains_point() tolerance, then we do a
  // BoundingBox intersection check with a tiny BoundingBox centered
  // on the search Point, otherwise we just do a simple
  // contains_point() check.

  if (_use_contains_point_tol)
    {
      Point min_copy = _local_bbox->min();
      Point max_copy = _local_bbox->max();

      // Compute absolute tolerance based on the bbox diagonal
      const Real abstol = (max_copy - min_copy).norm() * _contains_point_tol;

      // Inflate
      for (int i=0; i<LIBMESH_DIM; ++i)
        {
          min_copy(i) -= abstol;
          max_copy(i) += abstol;
        }

      BoundingBox p_bbox(min_copy, max_copy);
      return _local_bbox->intersects(p_bbox);
    }
  else
    return _local_bbox->contains_point(p);
}



const Elem *
PointLocatorNanoflann::operator() (const Point & p,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator()", "PointLocatorNanoflann");

  // Keep track of the number of elements checked in detail. This can
  // be useful for determining an optimal _num_results for
  // specific applications, by looking at the average and best/worse
  // case of each linear search.
  unsigned int n_elems_checked = 0;

  // We are not going to do any searching on this processor if the
  // Point doesn't fall in our (inflated) processor bounding box!
  bool point_in_local_bbox = search_local_bbox(p);

  // If a containing Elem is found locally, we will set this pointer.
  const Elem * found_elem = nullptr;

  // If the Point p is in our local bounding box, do a Nanoflann
  // findNeighbors() search for it. This returns a sorted list of the
  // closest _num_results elements which we then linearly
  // search.
  if (point_in_local_bbox)
    {
      auto result_set = this->kd_tree_find_neighbors(p, _num_results);

      // Note: we use result_set.size() here rather than
      // _num_results, in case the result returns fewer than
      // _num_results results!
      for (auto r : make_range(result_set.size()))
        {
          // Translate the Nanoflann index, which is from [0..n_points),
          // into the corresponding Elem id from the mesh.
          auto nanoflann_index = _ret_index[r];
          auto elem_id = (*_ids)[nanoflann_index];

          // Debugging: print the results
          // libMesh::err << "Centroid/Elem id = " << elem_id
          //              << ", dist^2 = " << _out_dist_sqr[r]
          //              << std::endl;

          const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

          // Before we even check whether the candidate Elem actually
          // contains the Point, we may need to check whether the
          // candidate Elem is from an allowed subdomain.  If the
          // candidate Elem is not from an allowed subdomain, we continue
          // to the next one.
          if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
            {
              // Debugging
              // libMesh::err << "Elem " << elem_id << " was not from an allowed subdomain, continuing search." << std::endl;
              continue;
            }

          // If we made it here, then the candidate Elem is from an
          // allowed subdomain, so let's next check whether it contains
          // the point. If the user set a custom tolerance, then we
          // actually check close_to_point() rather than contains_point(),
          // since this latter function warns about using non-default
          // tolerances, but otherwise does the same test.
          bool inside = _use_contains_point_tol ?
            candidate_elem->close_to_point(p, _contains_point_tol) :
            candidate_elem->contains_point(p);

          // Increment the number of elements checked
          n_elems_checked++;

          // If the point is inside an Elem from an allowed subdomain, we are done.
          if (inside)
            {
              // Debugging: report the number of Elems checked
              // libMesh::err << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;

              found_elem = candidate_elem;
              break;
            }
        } // end for(r)
    } // if (point_in_local_bbox)

  // If we made it here, then at least one of the following happened:
  // .) The search Point was not in the BoundingBox of local Elems.
  // .) All the candidate elements were from non-allowed subdomains.
  // .) The Point was not inside _any_ of the _num_results candidate Elems.
  // Thus, if we are not in _out_of_mesh_mode, throw an error,
  // otherwise return nullptr to indicate that no suitable element was
  // found.
  if (!_out_of_mesh_mode && !found_elem)
    libmesh_error_msg("Point was not contained within the closest " << n_elems_checked <<
                      " elems (by centroid distance), and _out_of_mesh_mode was not enabled.");

  return found_elem;
}


void
PointLocatorNanoflann::operator() (const Point & p,
                                   std::set<const Elem *> & candidate_elements,
                                   const std::set<subdomain_id_type> * allowed_subdomains) const
{
  libmesh_assert (this->_initialized);

  LOG_SCOPE("operator() returning set", "PointLocatorNanoflann");

  // Make sure the output set is initially empty
  candidate_elements.clear();

  // Keep track of the number of elements checked in detail
  unsigned int n_elems_checked = 0;

  // We are not going to do any searching on this processor if the
  // Point doesn't fall in our (inflated) processor bounding box!
  bool point_in_local_bbox = search_local_bbox(p);

  // If the Point p is in our local bounding box, do a Nanoflann
  // findNeighbors() search for it. We then store any of those
  // elements which also contain the Point in the candidate_elements
  // set.
  if (point_in_local_bbox)
    {
      // Do the KD-Tree search
      auto result_set = this->kd_tree_find_neighbors(p, _num_results);

      for (auto r : make_range(result_set.size()))
        {
          // Translate the Nanoflann index, which is from [0..n_points),
          // into the corresponding Elem id from the mesh.
          auto nanoflann_index = _ret_index[r];
          auto elem_id = (*_ids)[nanoflann_index];

          const Elem * candidate_elem = _mesh.elem_ptr(elem_id);

          // Before we even check whether the candidate Elem actually
          // contains the Point, we may need to check whether the
          // candidate Elem is from an allowed subdomain.  If the
          // candidate Elem is not from an allowed subdomain, we continue
          // to the next one.
          if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
            continue;

          // If we made it here, then the candidate Elem is from an
          // allowed subdomain, so let's next check whether it contains
          // the point. If the user set a custom tolerance, then we
          // actually check close_to_point() rather than contains_point(),
          // since this latter function warns about using non-default
          // tolerances, but otherwise does the same test.
          bool inside = _use_contains_point_tol ?
            candidate_elem->close_to_point(p, _contains_point_tol) :
            candidate_elem->contains_point(p);

          // Increment the number of elements checked
          n_elems_checked++;

          // If the point is contained in/close to an Elem from an
          // allowed subdomain, add it to the list.
          if (inside)
            candidate_elements.insert(candidate_elem);
        } // end for(r)
    } // if (point_in_local_bbox)

  // Debugging: for performance reasons, it may be useful to print the
  // number of Elems actually checked during the search.
  // libMesh::err << "Checked " << n_elems_checked << " nearby Elems before finding a containing Elem." << std::endl;
}


void
PointLocatorNanoflann::enable_out_of_mesh_mode ()
{
  // Out-of-mesh mode should now work properly even on meshes with
  // non-affine elements.
  _out_of_mesh_mode = true;
}


void
PointLocatorNanoflann::disable_out_of_mesh_mode ()
{
  _out_of_mesh_mode = false;
}



std::size_t
PointLocatorNanoflann::get_num_results() const
{
  return _num_results;
}



void
PointLocatorNanoflann::set_num_results(std::size_t val)
{
  // Must request at least 1 result
  _num_results = std::max(static_cast<std::size_t>(1), val);
}

//
// Required Nanoflann APIs
//

std::size_t PointLocatorNanoflann::kdtree_get_point_count() const
{
  return _point_cloud->size();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_distance(const coord_t * p1,
                                       const std::size_t idx_p2,
                                       std::size_t size) const
{
  // We only consider LIBMESH_DIM dimensional KD-Trees, so just make
  // sure we were called consistently.
  libmesh_assert(size == LIBMESH_DIM);

  // Construct a libmesh Point object from the LIBMESH_DIM-dimensional
  // input object, p1.
  Point point1;
  for (int i=0; i<LIBMESH_DIM; ++i)
    point1(i) = p1[i];

  // Compute Euclidean distance, squared
  return (point1 - (*_point_cloud)[idx_p2]).norm_sq();
}



PointLocatorNanoflann::coord_t
PointLocatorNanoflann::kdtree_get_pt(const std::size_t idx, int dim) const
{
  libmesh_assert_less (idx, _point_cloud->size());
  libmesh_assert_less (dim, LIBMESH_DIM);

  return (*_point_cloud)[idx](dim);
}

} // namespace libMesh

#endif
