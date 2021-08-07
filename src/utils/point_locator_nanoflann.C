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

// C++ includes
#include <array>
#include <numeric> // std::iota
#include <algorithm> // std::sort

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
  _elems.reset();
  _point_cloud.reset();
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
      // data structure with active, local element vertex averages.
      if (we_are_master)
        {
          _elems = std::make_shared<std::vector<Elem *>>();
          _point_cloud = std::make_shared<std::vector<Point>>();

          // Make the KD-Tree out of active element vertex averages.
          //
          // Note: we use active elements rather than active+local
          // elements since we sometimes need to be able to locate
          // points in ghosted elements (for example in Periodic BCs).
          // Active elements are also natural to use in the
          // DistributedMesh case. The size of the KD-Tree on each
          // processor therefore scales with the number of elements on
          // each processor in either case.
          //
          // Note 2: the approximate amount of space we should reserve
          // here is going to be different for ReplicatedMesh
          // (n_active_elem()) vs. DistributedMesh (n_active_local_elem())
          // so we just let the containers resize themselves
          // automatically.
          for (const auto & elem : _mesh.active_element_ptr_range())
            {
              _elems->push_back(elem);
              _point_cloud->push_back(elem->vertex_average());
            }

          // Construct the KD-Tree
          _kd_tree = std::make_shared<kd_tree_t>
            (LIBMESH_DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(/*max leaf=*/10));

          _kd_tree->buildIndex();
        }
      else // we are not master
        {
          // Cast the _master to the appropriate type, and point to
          // its data arrays.
          const auto my_master =
            cast_ptr<const PointLocatorNanoflann *>(this->_master);

          // Point our data structures at the master's
          _elems = my_master->_elems;
          _point_cloud = my_master->_point_cloud;
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

  // If a containing Elem is found locally, we will set this pointer.
  const Elem * found_elem = nullptr;

  auto result_set = this->kd_tree_find_neighbors(p, _num_results);

  // The results from Nanoflann are already sorted by (squared)
  // distance, but within that list of results, there may be some
  // vertex averages which are equidistant from the searched-for
  // Point. Therefore, we will now indirect_sort the results based on
  // elem id, so that the lowest-id Elem from the class of Elems which
  // are the same distance away from the search Point is always
  // selected.

  // operator< comparison lambda used in the indirect sort.  The
  // inputs to this vector are indices of the result_set array.
  auto comp = [this](std::size_t lhs, std::size_t rhs) -> bool
    {
      // First we sort by squared distance
      if (_out_dist_sqr[lhs] < _out_dist_sqr[rhs])
        return true;
      if (_out_dist_sqr[rhs] < _out_dist_sqr[lhs])
        return false;

      // If we made it here without returning, then the Points were
      // equidistant, so we sort based on Elem id instead.
      auto nanoflann_index_lhs = _ret_index[lhs];
      auto elem_id_lhs = (*_elems)[nanoflann_index_lhs]->id();

      auto nanoflann_index_rhs = _ret_index[rhs];
      auto elem_id_rhs = (*_elems)[nanoflann_index_rhs]->id();

      if (elem_id_lhs < elem_id_rhs)
        return true;

      return false;
    };

  // Set up the indices and do the indirect sort. The results are
  // stored in _b so that _b[0] is the first result to check, _b[1]
  // is second, etc.
  _b.resize(result_set.size());
  std::iota(_b.begin(), _b.end(), 0);
  std::sort(_b.begin(), _b.end(), comp);

  // Note: we use result_set.size() here rather than
  // _num_results, in case the result returns fewer than
  // _num_results results!
  for (auto r : make_range(result_set.size()))
    {
      // Translate the Nanoflann index, which is from [0..n_points),
      // into the corresponding Elem id from the mesh. Note that we
      // use index _b[r] instead of r, since we want to process the
      // Elems in order based on Elem id.
      auto nanoflann_index = _ret_index[_b[r]];
      const Elem * candidate_elem = (*_elems)[nanoflann_index];

      // Debugging: print the results
      // libMesh::err << "Vertex average/Elem id = " << candidate_elem->id()
      //              << ", dist^2 = " << _out_dist_sqr[r]
      //              << std::endl;

      // Before we even check whether the candidate Elem actually
      // contains the Point, we may need to check whether the
      // candidate Elem is from an allowed subdomain.  If the
      // candidate Elem is not from an allowed subdomain, we continue
      // to the next one.
      if (allowed_subdomains && !allowed_subdomains->count(candidate_elem->subdomain_id()))
        {
          // Debugging
          // libMesh::err << "Elem " << candidate_elem->id() << " was not from an allowed subdomain, continuing search." << std::endl;
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

  // If we made it here, then at least one of the following happened:
  // .) All the candidate elements were from non-allowed subdomains.
  // .) The Point was not inside _any_ of the _num_results candidate Elems.
  // Thus, if we are not in _out_of_mesh_mode, throw an error,
  // otherwise return nullptr to indicate that no suitable element was
  // found.
  if (!_out_of_mesh_mode && !found_elem)
    {
      // Debugging: we are about to throw an error, but before we do,
      // print information about the closest elements (by vertex average
      // distance) that the Point was not found in.
      // for (auto r : make_range(result_set.size()))
      //   {
      //     auto nanoflann_index = _ret_index[_b[r]];
      //     const Elem * candidate_elem = (*_elems)[nanoflann_index];
      //
      //     libMesh::err << "Vertex average/Elem id = " << candidate_elem->id()
      //                  << ", dist = " << std::sqrt(_out_dist_sqr[_b[r]])
      //                  << std::endl;
      //   } // end for(r)


      libmesh_error_msg("Point " << p << " was not contained within the closest " << n_elems_checked <<
                        " elems (by vertex average distance), and _out_of_mesh_mode was not enabled.");
    }

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

  // Do the KD-Tree search
  auto result_set = this->kd_tree_find_neighbors(p, _num_results);

  for (auto r : make_range(result_set.size()))
    {
      // Translate the Nanoflann index, which is from [0..n_points),
      // into the corresponding Elem id from the mesh.
      auto nanoflann_index = _ret_index[r];
      const Elem * candidate_elem = (*_elems)[nanoflann_index];

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
                                       std::size_t libmesh_dbg_var(size)) const
{
  // We only consider LIBMESH_DIM dimensional KD-Trees, so just make
  // sure we were called consistently.
  libmesh_assert_equal_to(size, LIBMESH_DIM);

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
