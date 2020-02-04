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

#ifndef LIBMESH_MESH_TOOLS_IMPL_H
#define LIBMESH_MESH_TOOLS_IMPL_H

// Local includes
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node_range.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/sphere.h"
#include "libmesh/threads.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h"

#ifdef DEBUG
#  include "libmesh/remote_elem.h"
#endif

// C++ includes
#include <limits>
#include <numeric> // for std::accumulate
#include <set>
#include <unordered_map>
#include <unordered_set>



// ------------------------------------------------------------
// anonymous namespace for helper classes
namespace {

using namespace libMesh;

/**
 * SumElemWeight(Range) sums the number of nodes per element
 * for each element in the provided range. The join() method
 * defines how to combine the reduction operation from two
 * distinct instances of this class which may be executed on
 * separate threads.
 */
class SumElemWeight
{
public:
  SumElemWeight () :
    _weight(0)
  {}

  SumElemWeight (SumElemWeight &, Threads::split) :
    _weight(0)
  {}

  template <typename RealType>
  void operator()(const ConstElemRangeTempl<RealType> & range)
  {
    for (const auto & elem : range)
      _weight += elem->n_nodes();
  }

  dof_id_type weight() const
  { return _weight; }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const SumElemWeight & other)
  { _weight += other.weight(); }
#endif

private:
  dof_id_type _weight;
};


/**
 * FindBBox(Range) computes the bounding box for the objects
 * in the specified range.  This class may be split and subranges
 * can be executed on separate threads.  The join() method
 * defines how the results from two separate threads are combined.
 */
template <typename RealType>
class FindBBox
{
public:
  FindBBox () : _bbox()
  {}

  FindBBox (FindBBox & other, Threads::split) :
    _bbox(other._bbox)
  {}

  void operator()(const ConstNodeRangeTempl<RealType> & range)
  {
    for (const auto & node : range)
      {
        libmesh_assert(node);
        _bbox.union_with(*node);
      }
  }

  void operator()(const ConstElemRangeTempl<RealType> & range)
  {
    for (const auto & elem : range)
      {
        libmesh_assert(elem);
        _bbox.union_with(elem->loose_bounding_box());
      }
  }

  PointTempl<RealType> & min() { return _bbox.min(); }

  PointTempl<RealType> & max() { return _bbox.max(); }

  // If we don't have threads we never need a join, and icpc yells a
  // warning if it sees an anonymous function that's never used
#if LIBMESH_USING_THREADS
  void join (const FindBBox & other)
  {
    _bbox.union_with(other._bbox);
  }
#endif

  libMesh::BoundingBoxTempl<RealType> & bbox ()
  {
    return _bbox;
  }

private:
  BoundingBoxTempl<RealType> _bbox;
};

#ifdef DEBUG
void assert_semiverify_dofobj(const Parallel::Communicator & communicator,
                              const DofObject * d,
                              unsigned int sysnum = libMesh::invalid_uint)
{
  if (d)
    {
      const unsigned int n_sys = d->n_systems();

      std::vector<unsigned int> n_vars (n_sys, 0);
      for (unsigned int s = 0; s != n_sys; ++s)
        if (sysnum == s ||
            sysnum == libMesh::invalid_uint)
          n_vars[s] = d->n_vars(s);

      const unsigned int tot_n_vars =
        std::accumulate(n_vars.begin(), n_vars.end(), 0);

      std::vector<unsigned int> n_comp (tot_n_vars, 0);
      std::vector<dof_id_type> first_dof (tot_n_vars, 0);

      for (unsigned int s = 0, i=0; s != n_sys; ++s)
        {
          if (sysnum != s &&
              sysnum != libMesh::invalid_uint)
            continue;

          for (unsigned int v = 0; v != n_vars[s]; ++v, ++i)
            {
              n_comp[i] = d->n_comp(s,v);
              first_dof[i] = n_comp[i] ? d->dof_number(s,v,0) : DofObject::invalid_id;
            }
        }

      libmesh_assert(communicator.semiverify(&n_sys));
      libmesh_assert(communicator.semiverify(&n_vars));
      libmesh_assert(communicator.semiverify(&n_comp));
      libmesh_assert(communicator.semiverify(&first_dof));
    }
  else
    {
      const unsigned int * p_ui = nullptr;
      const std::vector<unsigned int> * p_vui = nullptr;
      const std::vector<dof_id_type> * p_vdid = nullptr;

      libmesh_assert(communicator.semiverify(p_ui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vui));
      libmesh_assert(communicator.semiverify(p_vdid));
    }
}
#endif // DEBUG

}


namespace libMesh
{

// ------------------------------------------------------------
// MeshTools functions
template <typename RealType>
dof_id_type MeshTools::total_weight(const MeshBaseTempl<RealType> & mesh)
{
  if (!mesh.is_serial())
    {
      libmesh_parallel_only(mesh.comm());
      dof_id_type weight = MeshTools::weight (mesh, mesh.processor_id());
      mesh.comm().sum(weight);
      dof_id_type unpartitioned_weight =
        MeshTools::weight (mesh, DofObject::invalid_processor_id);
      return weight + unpartitioned_weight;
    }

  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.elements_begin(),
                                            mesh.elements_end()),
                            sew);
  return sew.weight();

}



template <typename RealType>
dof_id_type MeshTools::weight(const MeshBaseTempl<RealType> & mesh, const processor_id_type pid)
{
  SumElemWeight sew;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            sew);
  return sew.weight();
}



template <typename RealType>
void MeshTools::build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                                         std::vector<std::vector<dof_id_type>> & nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      {
        libmesh_assert_less (node.id(), nodes_to_elem_map.size());
        libmesh_assert_less (elem->id(), mesh.n_elem());

        nodes_to_elem_map[node.id()].push_back(elem->id());
      }
}



template <typename RealType>
void MeshTools::build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                                         std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      {
        libmesh_assert_less (node.id(), nodes_to_elem_map.size());

        nodes_to_elem_map[node.id()].push_back(elem);
      }
}



template <typename RealType>
void MeshTools::build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                                         std::unordered_map<dof_id_type, std::vector<dof_id_type>> & nodes_to_elem_map)
{
  nodes_to_elem_map.clear();

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      nodes_to_elem_map[node.id()].push_back(elem->id());
}



template <typename RealType>
void MeshTools::build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                                         std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map)
{
  nodes_to_elem_map.clear();

  for (const auto & elem : mesh.element_ptr_range())
    for (auto & node : elem->node_ref_range())
      nodes_to_elem_map[node.id()].push_back(elem);
}



#ifdef LIBMESH_ENABLE_DEPRECATED
template <typename RealType>
void MeshTools::find_boundary_nodes (const MeshBaseTempl<RealType> & mesh,
                                     std::vector<bool> & on_boundary)
{
  libmesh_deprecated();

  // Resize the vector which holds boundary nodes and fill with false.
  on_boundary.resize(mesh.max_node_id());
  std::fill(on_boundary.begin(),
            on_boundary.end(),
            false);

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr) // on the boundary
        {
          std::unique_ptr<const Elem> side = elem->build_side_ptr(s);

          auto nodes_on_side = elem->nodes_on_side(s);

          for (auto & node_id : nodes_on_side)
            on_boundary[node_id] = true;
        }
}
#endif

template <typename RealType>
std::unordered_set<dof_id_type>
MeshTools::find_boundary_nodes(const MeshBaseTempl<RealType> & mesh)
{
  std::unordered_set<dof_id_type> boundary_nodes;

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) == nullptr) // on the boundary
        {
          auto nodes_on_side = elem->nodes_on_side(s);

          for (auto & local_id : nodes_on_side)
            boundary_nodes.insert(elem->node_ptr(local_id)->id());
        }

  return boundary_nodes;
}

template <typename RealType>
std::unordered_set<dof_id_type>
MeshTools::find_block_boundary_nodes(const MeshBaseTempl<RealType> & mesh)
{
  std::unordered_set<dof_id_type> block_boundary_nodes;

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto s : elem->side_index_range())
      if (elem->neighbor_ptr(s) && (elem->neighbor_ptr(s)->subdomain_id() != elem->subdomain_id()))
        {
          auto nodes_on_side = elem->nodes_on_side(s);

          for (auto & local_id : nodes_on_side)
            block_boundary_nodes.insert(elem->node_ptr(local_id)->id());
        }

  return block_boundary_nodes;
}


#ifdef LIBMESH_ENABLE_DEPRECATED
inline
MeshTools::BoundingBox
MeshTools::bounding_box(const MeshBaseTempl<Real> & mesh)
{
  // This function is deprecated.  It simply calls
  // create_bounding_box() and converts the result to a
  // MeshTools::BoundingBox.
  libmesh_deprecated();
  return MeshTools::create_bounding_box(mesh);
}
#endif



template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
MeshTools::create_bounding_box (const MeshBaseTempl<RealType> & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  FindBBox<RealType> find_bbox;

  // Start with any unpartitioned elements we know about locally
  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(DofObject::invalid_processor_id),
                                            mesh.pid_elements_end(DofObject::invalid_processor_id)),
                            find_bbox);

  // And combine with our local elements
  find_bbox.bbox().union_with(MeshTools::create_local_bounding_box(mesh));

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
MeshTools::create_nodal_bounding_box (const MeshBaseTempl<RealType> & mesh)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  FindBBox<RealType> find_bbox;

  // Start with any unpartitioned nodes we know about locally
  Threads::parallel_reduce (ConstNodeRange (mesh.pid_nodes_begin(DofObject::invalid_processor_id),
                                            mesh.pid_nodes_end(DofObject::invalid_processor_id)),
                            find_bbox);

  // Add our local nodes
  Threads::parallel_reduce (ConstNodeRange (mesh.local_nodes_begin(),
                                            mesh.local_nodes_end()),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



template <typename RealType>
Sphere
MeshTools::bounding_sphere(const MeshBaseTempl<RealType> & mesh)
{
  libMesh::BoundingBoxTempl<RealType> bbox = MeshTools::create_bounding_box(mesh);

  const auto diag = (bbox.second - bbox.first).norm();
  const auto cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
MeshTools::create_local_bounding_box (const MeshBaseTempl<RealType> & mesh)
{
  FindBBox<RealType> find_bbox;

  Threads::parallel_reduce (ConstElemRange (mesh.local_elements_begin(),
                                            mesh.local_elements_end()),
                            find_bbox);

  return find_bbox.bbox();
}



#ifdef LIBMESH_ENABLE_DEPRECATED
inline
MeshTools::BoundingBox
MeshTools::processor_bounding_box (const MeshBaseTempl<Real> & mesh,
                                   const processor_id_type pid)
{
  libmesh_deprecated();
  return MeshTools::create_processor_bounding_box(mesh, pid);
}
#endif



template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
MeshTools::create_processor_bounding_box (const MeshBaseTempl<RealType> & mesh,
                                          const processor_id_type pid)
{
  // This can only be run in parallel, with consistent arguments.
  libmesh_parallel_only(mesh.comm());
  libmesh_assert(mesh.comm().verify(pid));

  libmesh_assert_less (pid, mesh.n_processors());

  FindBBox<RealType> find_bbox;

  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
                                            mesh.pid_elements_end(pid)),
                            find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



template <typename RealType>
Sphere
MeshTools::processor_bounding_sphere (const MeshBaseTempl<RealType> & mesh,
                                      const processor_id_type pid)
{
  auto bbox =
    MeshTools::create_processor_bounding_box(mesh, pid);

  const auto diag = (bbox.second - bbox.first).norm();
  const auto cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



#ifdef LIBMESH_ENABLE_DEPRECATED
inline
MeshTools::BoundingBox
MeshTools::subdomain_bounding_box (const MeshBaseTempl<Real> & mesh,
                                   const subdomain_id_type sid)
{
  libmesh_deprecated();
  return MeshTools::create_subdomain_bounding_box(mesh, sid);
}
#endif



template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
MeshTools::create_subdomain_bounding_box (const MeshBaseTempl<RealType> & mesh,
                                          const subdomain_id_type sid)
{
  // This can only be run in parallel, with consistent arguments.
  libmesh_parallel_only(mesh.comm());
  libmesh_assert(mesh.comm().verify(sid));

  FindBBox<RealType> find_bbox;

  Threads::parallel_reduce
    (ConstElemRange (mesh.active_local_subdomain_elements_begin(sid),
                     mesh.active_local_subdomain_elements_end(sid)),
     find_bbox);

  // Compare the bounding boxes across processors
  mesh.comm().min(find_bbox.min());
  mesh.comm().max(find_bbox.max());

  return find_bbox.bbox();
}



template <typename RealType>
Sphere
MeshTools::subdomain_bounding_sphere (const MeshBaseTempl<RealType> & mesh,
                                      const subdomain_id_type sid)
{
  auto bbox =
    MeshTools::create_subdomain_bounding_box(mesh, sid);

  const auto diag = (bbox.second - bbox.first).norm();
  const auto cent = (bbox.second + bbox.first)/2;

  return Sphere (cent, .5*diag);
}



template <typename RealType>
void MeshTools::elem_types (const MeshBaseTempl<RealType> & mesh,
                            std::vector<ElemType> & et)
{
  // Loop over the the elements.  If the current element type isn't in
  // the vector, insert it.
  for (const auto & elem : mesh.element_ptr_range())
    if (!std::count(et.begin(), et.end(), elem->type()))
      et.push_back(elem->type());
}



template <typename RealType>
dof_id_type MeshTools::n_elem_of_type (const MeshBaseTempl<RealType> & mesh,
                                       const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.type_elements_begin(type),
                                                mesh.type_elements_end  (type)));
}



template <typename RealType>
dof_id_type MeshTools::n_active_elem_of_type (const MeshBaseTempl<RealType> & mesh,
                                              const ElemType type)
{
  return static_cast<dof_id_type>(std::distance(mesh.active_type_elements_begin(type),
                                                mesh.active_type_elements_end  (type)));
}

template <typename RealType>
dof_id_type MeshTools::n_non_subactive_elem_of_type_at_level(const MeshBaseTempl<RealType> & mesh,
                                                             const ElemType type,
                                                             const unsigned int level)
{
  dof_id_type cnt = 0;

  // iterate over the elements of the specified type
  for (const auto & elem : as_range(mesh.type_elements_begin(type),
                                    mesh.type_elements_end(type)))
    if ((elem->level() == level) && !elem->subactive())
      cnt++;

  return cnt;
}


template <typename RealType>
unsigned int MeshTools::n_active_local_levels(const MeshBaseTempl<RealType> & mesh)
{
  unsigned int nl = 0;

  for (auto & elem : mesh.active_local_element_ptr_range())
    nl = std::max(elem->level() + 1, nl);

  return nl;
}



template <typename RealType>
unsigned int MeshTools::n_active_levels(const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = MeshTools::n_active_local_levels(mesh);

  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    if (elem->active())
      nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



template <typename RealType>
unsigned int MeshTools::n_local_levels(const MeshBaseTempl<RealType> & mesh)
{
  unsigned int nl = 0;

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    nl = std::max(elem->level() + 1, nl);

  return nl;
}



template <typename RealType>
unsigned int MeshTools::n_levels(const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = MeshTools::n_local_levels(mesh);

  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);

  // n_levels() is only valid and should only be called in cases where
  // the mesh is validly distributed (or serialized).  Let's run an
  // expensive test in debug mode to make sure this is such a case.
#ifdef DEBUG
  const unsigned int paranoid_nl = MeshTools::paranoid_n_levels(mesh);
  libmesh_assert_equal_to(nl, paranoid_nl);
#endif
  return nl;
}



template <typename RealType>
unsigned int MeshTools::paranoid_n_levels(const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int nl = 0;
  for (const auto & elem : mesh.element_ptr_range())
    nl = std::max(elem->level() + 1, nl);

  mesh.comm().max(nl);
  return nl;
}



template <typename RealType>
void MeshTools::get_not_subactive_node_ids(const MeshBaseTempl<RealType> & mesh,
                                           std::set<dof_id_type> & not_subactive_node_ids)
{
  for (const auto & elem : mesh.element_ptr_range())
    if (!elem->subactive())
      for (auto & n : elem->node_ref_range())
        not_subactive_node_ids.insert(n.id());
}



template <typename IteratorType>
dof_id_type MeshTools::n_elem (const IteratorType & begin,
                               const IteratorType & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



template <typename IteratorType>
dof_id_type MeshTools::n_nodes (const IteratorType & begin,
                                const IteratorType & end)
{
  return cast_int<dof_id_type>(std::distance(begin, end));
}



template <typename RealType>
unsigned int MeshTools::n_p_levels (const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());

  unsigned int max_p_level = 0;

  // first my local elements
  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    max_p_level = std::max(elem->p_level(), max_p_level);

  // then any unpartitioned objects
  for (const auto & elem : as_range(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()))
    max_p_level = std::max(elem->p_level(), max_p_level);

  mesh.comm().max(max_p_level);
  return max_p_level + 1;
}



template <typename RealType>
void MeshTools::find_nodal_neighbors(const MeshBaseTempl<RealType> &,
                                     const NodeTempl<RealType> & node,
                                     const std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map,
                                     std::vector<const NodeTempl<RealType> *> & neighbors)
{
  // We'll refer back to the Node ID several times
  dof_id_type global_id = node.id();

  // We'll construct a std::set<const NodeTempl<RealType> *> for more efficient
  // searching while finding the nodal neighbors, and return it to the
  // user in a std::vector.
  std::set<const NodeTempl<RealType> *> neighbor_set;

  // Look through the elements that contain this node
  // find the local node id... then find the side that
  // node lives on in the element
  // next, look for the _other_ node on that side
  // That other node is a "nodal_neighbor"... save it
  for (const auto & elem : nodes_to_elem_map[global_id])
    {
      // We only care about active elements...
      if (elem->active())
        {
          // Which local node number is global_id?
          unsigned local_node_number = elem->local_node(global_id);

          // Make sure it was found
          libmesh_assert_not_equal_to(local_node_number, libMesh::invalid_uint);

          const unsigned short n_edges = elem->n_edges();

          // If this element has no edges, the edge-based algorithm below doesn't make sense.
          if (!n_edges)
            {
              switch (elem->type())
                {
                case EDGE2:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      case 1:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(0));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE3:
                  {
                    switch (local_node_number)
                      {
                        // The outside nodes have node 2 as a neighbor
                      case 0:
                      case 1:
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                        // The middle node has the outer nodes as neighbors
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE4:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The left-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      case 1:
                        // The right-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The left-middle node
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The right-middle node
                      case 3:
                        neighbor_set.insert(elem->node_ptr(1));
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                default:
                  libmesh_error_msg("Unrecognized ElemType: " << Utility::enum_to_string(elem->type()) << std::endl);
                }
            }

          // Index of the current edge
          unsigned current_edge = 0;

          const unsigned short n_nodes = elem->n_nodes();

          while (current_edge < n_edges)
            {
              // Find the edge the node is on
              bool found_edge = false;
              for (; current_edge<n_edges; ++current_edge)
                if (elem->is_node_on_edge(local_node_number, current_edge))
                  {
                    found_edge = true;
                    break;
                  }

              // Did we find one?
              if (found_edge)
                {
                  const NodeTempl<RealType> * node_to_save = nullptr;

                  // Find another node in this element on this edge
                  for (unsigned other_node_this_edge = 0; other_node_this_edge != n_nodes; other_node_this_edge++)
                    if ( (elem->is_node_on_edge(other_node_this_edge, current_edge)) && // On the current edge
                         (elem->node_id(other_node_this_edge) != global_id))               // But not the original node
                      {
                        // We've found a nodal neighbor!  Save a pointer to it..
                        node_to_save = elem->node_ptr(other_node_this_edge);
                        break;
                      }

                  // Make sure we found something
                  libmesh_assert(node_to_save != nullptr);

                  neighbor_set.insert(node_to_save);
                }

              // Keep looking for edges, node may be on more than one edge
              current_edge++;
            }
        } // if (elem->active())
    } // for

  // Assign the entries from the set to the vector.  Note: this
  // replaces any existing contents in neighbors and modifies its size
  // accordingly.
  neighbors.assign(neighbor_set.begin(), neighbor_set.end());
}



template <typename RealType>
void MeshTools::find_nodal_neighbors(const MeshBaseTempl<RealType> &,
                                     const NodeTempl<RealType> & node,
                                     const std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map,
                                     std::vector<const NodeTempl<RealType> *> & neighbors)
{
  // We'll refer back to the NodeTempl<RealType> ID several times
  dof_id_type global_id = node.id();

  // We'll construct a std::set<const NodeTempl<RealType> *> for more efficient
  // searching while finding the nodal neighbors, and return it to the
  // user in a std::vector.
  std::set<const NodeTempl<RealType> *> neighbor_set;

  // List of Elems attached to this node.
  const auto & elem_vec = libmesh_map_find(nodes_to_elem_map, global_id);

  // Look through the elements that contain this node
  // find the local node id... then find the side that
  // node lives on in the element
  // next, look for the _other_ node on that side
  // That other node is a "nodal_neighbor"... save it
  for (const auto & elem : elem_vec)
    {
      // We only care about active elements...
      if (elem->active())
        {
          // Which local node number is global_id?
          unsigned local_node_number = elem->local_node(global_id);

          // Make sure it was found
          libmesh_assert_not_equal_to(local_node_number, libMesh::invalid_uint);

          const unsigned short n_edges = elem->n_edges();

          // If this element has no edges, the edge-based algorithm below doesn't make sense.
          if (!n_edges)
            {
              switch (elem->type())
                {
                case EDGE2:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      case 1:
                        // The other node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(0));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE3:
                  {
                    switch (local_node_number)
                      {
                        // The outside nodes have node 2 as a neighbor
                      case 0:
                      case 1:
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                        // The middle node has the outer nodes as neighbors
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(1));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                case EDGE4:
                  {
                    switch (local_node_number)
                      {
                      case 0:
                        // The left-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      case 1:
                        // The right-middle node is a nodal neighbor
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The left-middle node
                      case 2:
                        neighbor_set.insert(elem->node_ptr(0));
                        neighbor_set.insert(elem->node_ptr(3));
                        break;

                        // The right-middle node
                      case 3:
                        neighbor_set.insert(elem->node_ptr(1));
                        neighbor_set.insert(elem->node_ptr(2));
                        break;

                      default:
                        libmesh_error_msg("Invalid local node number: " << local_node_number << " found." << std::endl);
                      }
                    break;
                  }

                default:
                  libmesh_error_msg("Unrecognized ElemType: " << Utility::enum_to_string(elem->type()) << std::endl);
                }
            }

          // Index of the current edge
          unsigned current_edge = 0;

          const unsigned short n_nodes = elem->n_nodes();

          while (current_edge < n_edges)
            {
              // Find the edge the node is on
              bool found_edge = false;
              for (; current_edge<n_edges; ++current_edge)
                if (elem->is_node_on_edge(local_node_number, current_edge))
                  {
                    found_edge = true;
                    break;
                  }

              // Did we find one?
              if (found_edge)
                {
                  const NodeTempl<RealType> * node_to_save = nullptr;

                  // Find another node in this element on this edge
                  for (unsigned other_node_this_edge = 0; other_node_this_edge != n_nodes; other_node_this_edge++)
                    if ( (elem->is_node_on_edge(other_node_this_edge, current_edge)) && // On the current edge
                         (elem->node_id(other_node_this_edge) != global_id))               // But not the original node
                      {
                        // We've found a nodal neighbor!  Save a pointer to it..
                        node_to_save = elem->node_ptr(other_node_this_edge);
                        break;
                      }

                  // Make sure we found something
                  libmesh_assert(node_to_save != nullptr);

                  neighbor_set.insert(node_to_save);
                }

              // Keep looking for edges, node may be on more than one edge
              current_edge++;
            }
        } // if (elem->active())
    } // for

  // Assign the entries from the set to the vector.  Note: this
  // replaces any existing contents in neighbors and modifies its size
  // accordingly.
  neighbors.assign(neighbor_set.begin(), neighbor_set.end());
}



template <typename RealType>
void MeshTools::find_hanging_nodes_and_parents(const MeshBaseTempl<RealType> & mesh,
                                               std::map<dof_id_type, std::vector<dof_id_type>> & hanging_nodes)
{
  // Loop through all the elements
  for (auto & elem : mesh.active_local_element_ptr_range())
    if (elem->type() == QUAD4)
      for (auto s : elem->side_index_range())
        {
          // Loop over the sides looking for sides that have hanging nodes
          // This code is inspired by compute_proj_constraints()
          const ElemTempl<RealType> * neigh = elem->neighbor_ptr(s);

          // If not a boundary side
          if (neigh != nullptr)
            {
              // Is there a coarser element next to this one?
              if (neigh->level() < elem->level())
                {
                  const ElemTempl<RealType> * ancestor = elem;
                  while (neigh->level() < ancestor->level())
                    ancestor = ancestor->parent();
                  unsigned int s_neigh = neigh->which_neighbor_am_i(ancestor);
                  libmesh_assert_less (s_neigh, neigh->n_neighbors());

                  // Couple of helper uints...
                  unsigned int local_node1=0;
                  unsigned int local_node2=0;

                  bool found_in_neighbor = false;

                  // Find the two vertices that make up this side
                  while (!elem->is_node_on_side(local_node1++,s)) { }
                  local_node1--;

                  // Start looking for the second one with the next node
                  local_node2=local_node1+1;

                  // Find the other one
                  while (!elem->is_node_on_side(local_node2++,s)) { }
                  local_node2--;

                  //Pull out their global ids:
                  dof_id_type node1 = elem->node_id(local_node1);
                  dof_id_type node2 = elem->node_id(local_node2);

                  // Now find which node is present in the neighbor
                  // FIXME This assumes a level one rule!
                  // The _other_ one is the hanging node

                  // First look for the first one
                  // FIXME could be streamlined a bit
                  for (unsigned int n=0;n<neigh->n_sides();n++)
                    if (neigh->node_id(n) == node1)
                      found_in_neighbor=true;

                  dof_id_type hanging_node=0;

                  if (!found_in_neighbor)
                    hanging_node=node1;
                  else // If it wasn't node1 then it must be node2!
                    hanging_node=node2;

                  // Reset these for reuse
                  local_node1=0;
                  local_node2=0;

                  // Find the first node that makes up the side in the neighbor (these should be the parent nodes)
                  while (!neigh->is_node_on_side(local_node1++,s_neigh)) { }
                  local_node1--;

                  local_node2=local_node1+1;

                  // Find the second node...
                  while (!neigh->is_node_on_side(local_node2++,s_neigh)) { }
                  local_node2--;

                  // Save them if we haven't already found the parents for this one
                  if (hanging_nodes[hanging_node].size()<2)
                    {
                      hanging_nodes[hanging_node].push_back(neigh->node_id(local_node1));
                      hanging_nodes[hanging_node].push_back(neigh->node_id(local_node2));
                    }
                }
            }
        }
}



#ifdef DEBUG
template <typename RealType>
void MeshTools::libmesh_assert_equal_n_systems (const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_equal_n_systems()", "MeshTools");

  unsigned int n_sys = libMesh::invalid_uint;

  for (const auto & elem : mesh.element_ptr_range())
    {
      if (n_sys == libMesh::invalid_uint)
        n_sys = elem->n_systems();
      else
        libmesh_assert_equal_to (elem->n_systems(), n_sys);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      if (n_sys == libMesh::invalid_uint)
        n_sys = node->n_systems();
      else
        libmesh_assert_equal_to (node->n_systems(), n_sys);
    }
}



#ifdef LIBMESH_ENABLE_AMR
template <typename RealType>
void MeshTools::libmesh_assert_old_dof_objects (const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_old_dof_objects()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      if (elem->refinement_flag() == Elem::JUST_REFINED ||
          elem->refinement_flag() == Elem::INACTIVE)
        continue;

      if (elem->has_dofs())
        libmesh_assert(elem->old_dof_object);

      for (auto & node : elem->node_ref_range())
        if (node.has_dofs())
          libmesh_assert(node.old_dof_object);
    }
}
#else
template <typename RealType>
void MeshTools::libmesh_assert_old_dof_objects (const MeshBaseTempl<RealType> &) {}
#endif // LIBMESH_ENABLE_AMR



template <typename RealType>
void MeshTools::libmesh_assert_valid_node_pointers(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_node_pointers()", "MeshTools");

  // Here we specifically do not want "auto &" because we need to
  // reseat the (temporary) pointer variable in the loop below,
  // without modifying the original.
  for (const ElemTempl<RealType> * elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      while (elem)
        {
          elem->libmesh_assert_valid_node_pointers();
          for (auto n : elem->neighbor_ptr_range())
            if (n && n != RemoteElem::get_instance())
              n->libmesh_assert_valid_node_pointers();

          libmesh_assert_not_equal_to (elem->parent(), RemoteElem::get_instance());
          elem = elem->parent();
        }
    }
}


template <typename RealType>
void MeshTools::libmesh_assert_valid_remote_elems(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_remote_elems()", "MeshTools");

  for (const auto & elem : as_range(mesh.local_elements_begin(),
                                    mesh.local_elements_end()))
    {
      libmesh_assert (elem);

      // We currently don't allow active_local_elements to have
      // remote_elem neighbors
      if (elem->active())
        for (auto n : elem->neighbor_ptr_range())
          libmesh_assert_not_equal_to (n, RemoteElem::get_instance());

#ifdef LIBMESH_ENABLE_AMR
      const ElemTempl<RealType> * parent = elem->parent();
      if (parent)
        libmesh_assert_not_equal_to (parent, RemoteElem::get_instance());

      // We can only be strict about active elements' subactive
      // children
      if (elem->active() && elem->has_children())
        for (auto & child : elem->child_ref_range())
          libmesh_assert_not_equal_to (&child, RemoteElem::get_instance());
#endif
    }
}


template <typename RealType>
void MeshTools::libmesh_assert_no_links_to_elem(const MeshBaseTempl<RealType> & mesh,
                                                const ElemTempl<RealType> * bad_elem)
{
  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      libmesh_assert_not_equal_to (elem->parent(), bad_elem);
      for (auto n : elem->neighbor_ptr_range())
        libmesh_assert_not_equal_to (n, bad_elem);

#ifdef LIBMESH_ENABLE_AMR
      if (elem->has_children())
        for (auto & child : elem->child_ref_range())
          libmesh_assert_not_equal_to (&child, bad_elem);
#endif
    }
}



template <typename RealType>
void MeshTools::libmesh_assert_valid_elem_ids(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_elem_ids()", "MeshTools");

  processor_id_type lastprocid = 0;
  dof_id_type lastelemid = 0;

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      libmesh_assert (elem);
      processor_id_type elemprocid = elem->processor_id();
      dof_id_type elemid = elem->id();

      libmesh_assert_greater_equal (elemid, lastelemid);
      libmesh_assert_greater_equal (elemprocid, lastprocid);

      lastelemid = elemid;
      lastprocid = elemprocid;
    }
}



template <typename RealType>
void MeshTools::libmesh_assert_valid_amr_elem_ids(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_amr_elem_ids()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      const ElemTempl<RealType> * parent = elem->parent();

      if (parent)
        {
          libmesh_assert_greater_equal (elem->id(), parent->id());
          libmesh_assert_greater_equal (elem->processor_id(), parent->processor_id());
        }
    }
}



template <typename RealType>
void MeshTools::libmesh_assert_valid_amr_interior_parents(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_amr_interior_parents()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      // We can skip to the next element if we're full-dimension
      // and therefore don't have any interior parents
      if (elem->dim() >= LIBMESH_DIM)
        continue;

      const ElemTempl<RealType> * ip = elem->interior_parent();

      const ElemTempl<RealType> * parent = elem->parent();

      if (ip && (ip != RemoteElem::get_instance()) && parent)
        {
          libmesh_assert_equal_to (ip->top_parent(),
                                   elem->top_parent()->interior_parent());

          if (ip->level() == elem->level())
            libmesh_assert_equal_to (ip->parent(),
                                     parent->interior_parent());
          else
            {
              libmesh_assert_less (ip->level(), elem->level());
              libmesh_assert_equal_to (ip, parent->interior_parent());
            }
        }
    }
}



template <typename RealType>
void MeshTools::libmesh_assert_connected_nodes (const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_connected_nodes()", "MeshTools");

  std::set<const NodeTempl<RealType> *> used_nodes;

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);

      for (auto & n : elem->node_ref_range())
        used_nodes.insert(&n);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      libmesh_assert(node);
      libmesh_assert(used_nodes.count(node));
    }
}



namespace MeshTools {

template <typename RealType>
void libmesh_assert_valid_boundary_ids(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_boundary_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  const BoundaryInfo & boundary_info = mesh.get_boundary_info();

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);
      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      const unsigned int my_n_edges = elem ? elem->n_edges() : 0;
      const unsigned int my_n_sides = elem ? elem->n_sides() : 0;
      unsigned int
        n_nodes = my_n_nodes,
        n_edges = my_n_edges,
        n_sides = my_n_sides;

      mesh.comm().max(n_nodes);
      mesh.comm().max(n_edges);
      mesh.comm().max(n_sides);

      if (elem)
        {
          libmesh_assert_equal_to(my_n_nodes, n_nodes);
          libmesh_assert_equal_to(my_n_edges, n_edges);
          libmesh_assert_equal_to(my_n_sides, n_sides);
        }

      // Let's test all IDs on the element with one communication
      // rather than n_nodes + n_edges + n_sides communications, to
      // cut down on latency in dbg modes.
      std::vector<boundary_id_type> all_bcids;

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          std::vector<boundary_id_type> bcids;
          if (elem)
            {
              boundary_info.boundary_ids(elem->node_ptr(n), bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }
          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          all_bcids.insert(all_bcids.end(), bcids.begin(),
                           bcids.end());
          // Separator
          all_bcids.push_back(BoundaryInfo::invalid_id);
        }

      for (unsigned short e=0; e != n_edges; ++e)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          all_bcids.insert(all_bcids.end(), bcids.begin(),
                           bcids.end());
          // Separator
          all_bcids.push_back(BoundaryInfo::invalid_id);

          if (elem)
            {
              boundary_info.raw_edge_boundary_ids(elem, e, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      for (unsigned short s=0; s != n_sides; ++s)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          if (elem)
            {
              boundary_info.raw_boundary_ids(elem, s, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      for (unsigned short sf=0; sf != 2; ++sf)
        {
          std::vector<boundary_id_type> bcids;

          if (elem)
            {
              boundary_info.shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));

          if (elem)
            {
              boundary_info.raw_shellface_boundary_ids(elem, sf, bcids);

              // Ordering of boundary ids shouldn't matter
              std::sort(bcids.begin(), bcids.end());

              all_bcids.insert(all_bcids.end(), bcids.begin(),
                               bcids.end());
              // Separator
              all_bcids.push_back(BoundaryInfo::invalid_id);
            }

          // libmesh_assert(mesh.comm().semiverify (elem ? &bcids : nullptr));
        }

      libmesh_assert(mesh.comm().semiverify
                     (elem ? &all_bcids : nullptr));
    }
}


template <typename RealType>
void libmesh_assert_valid_dof_ids(const MeshBaseTempl<RealType> & mesh, unsigned int sysnum)
{
  LOG_SCOPE("libmesh_assert_valid_dof_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_elem_ptr(i),
                             sysnum);

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    assert_semiverify_dofobj(mesh.comm(),
                             mesh.query_node_ptr(i),
                             sysnum);
}


template <typename RealType>
void libmesh_assert_contiguous_dof_ids(const MeshBaseTempl<RealType> & mesh, unsigned int sysnum)
{
  LOG_SCOPE("libmesh_assert_contiguous_dof_ids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type min_dof_id = std::numeric_limits<dof_id_type>::max(),
              max_dof_id = std::numeric_limits<dof_id_type>::min();

  // Figure out what our local dof id range is
  for (const auto * node : mesh.local_node_ptr_range())
    {
      for (auto v : IntRange<unsigned int>(0, node->n_vars(sysnum)))
        for (auto c : IntRange<unsigned int>(0, node->n_comp(sysnum, v)))
          {
            dof_id_type id = node->dof_number(sysnum, v, c);
            min_dof_id = std::min (min_dof_id, id);
            max_dof_id = std::max (max_dof_id, id);
          }
    }

  // Make sure no other processors' ids are inside it
  for (const auto * node : mesh.node_ptr_range())
    {
      if (node->processor_id() == mesh.processor_id())
        continue;
      for (auto v : IntRange<unsigned int>(0, node->n_vars(sysnum)))
        for (auto c : IntRange<unsigned int>(0, node->n_comp(sysnum, v)))
          {
            dof_id_type id = node->dof_number(sysnum, v, c);
            libmesh_assert (id < min_dof_id ||
                            id > max_dof_id);
          }
    }
}


#ifdef LIBMESH_ENABLE_UNIQUE_ID
template <typename RealType>
void libmesh_assert_valid_unique_ids(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_unique_ids()", "MeshTools");

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);
      const unique_id_type unique_id = elem ? elem->unique_id() : 0;
      const unique_id_type * uid_ptr = elem ? &unique_id : nullptr;
      libmesh_assert(mesh.comm().semiverify(uid_ptr));
    }

  dof_id_type pmax_node_id = mesh.max_node_id();
  mesh.comm().max(pmax_node_id);

  for (dof_id_type i=0; i != pmax_node_id; ++i)
    {
      const NodeTempl<RealType> * node = mesh.query_node_ptr(i);
      const unique_id_type unique_id = node ? node->unique_id() : 0;
      const unique_id_type * uid_ptr = node ? &unique_id : nullptr;
      libmesh_assert(mesh.comm().semiverify(uid_ptr));
    }
}
#endif

template <typename RealType>
void libmesh_assert_consistent_distributed(const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());

  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);
      processor_id_type pid =
        elem ? elem->processor_id() : DofObject::invalid_processor_id;
      mesh.comm().min(pid);
      libmesh_assert(elem || pid != mesh.processor_id());
    }

  dof_id_type parallel_max_node_id = mesh.max_node_id();
  mesh.comm().max(parallel_max_node_id);

  for (dof_id_type i=0; i != parallel_max_node_id; ++i)
    {
      const NodeTempl<RealType> * node = mesh.query_node_ptr(i);
      processor_id_type pid =
        node ? node->processor_id() : DofObject::invalid_processor_id;
      mesh.comm().min(pid);
      libmesh_assert(node || pid != mesh.processor_id());
    }
}


template <typename RealType>
void libmesh_assert_consistent_distributed_nodes(const MeshBaseTempl<RealType> & mesh)
{
  libmesh_parallel_only(mesh.comm());
  auto locator = mesh.sub_point_locator();

  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      unsigned int n_nodes = my_n_nodes;
      mesh.comm().max(n_nodes);

      if (n_nodes)
        libmesh_assert(mesh.comm().semiverify(elem ? &my_n_nodes : nullptr));

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          const NodeTempl<RealType> * node = elem ? elem->node_ptr(n) : nullptr;
          processor_id_type pid =
            node ? node->processor_id() : DofObject::invalid_processor_id;
          mesh.comm().min(pid);
          libmesh_assert(node || pid != mesh.processor_id());
        }
    }
}

namespace
{
template <typename>
struct LibMeshAssertShim;

template <typename RealType>
struct LibMeshAssertShim<ElemTempl<RealType>>
{
  static void topology_consistent_procids(const MeshBaseTempl<RealType> & mesh)
    {
      LOG_SCOPE("libmesh_assert_topology_consistent_procids()", "MeshTools");

      // This parameter is not used when !LIBMESH_ENABLE_AMR
      libmesh_ignore(mesh);

      // If we're adaptively refining, check processor ids for consistency
      // between parents and children.
#ifdef LIBMESH_ENABLE_AMR

      // Ancestor elements we won't worry about, but subactive and active
      // elements ought to have parents with consistent processor ids
      for (const auto & elem : mesh.element_ptr_range())
      {
        libmesh_assert(elem);

        if (!elem->active() && !elem->subactive())
          continue;

        const ElemTempl<RealType> * parent = elem->parent();

        if (parent)
        {
          libmesh_assert(parent->has_children());
          processor_id_type parent_procid = parent->processor_id();
          bool matching_child_id = false;
          // If we've got a remote_elem then we don't know whether
          // it's responsible for the parent's processor id; all
          // we can do is assume it is and let its processor fail
          // an assert if there's something wrong.
          for (auto & child : parent->child_ref_range())
            if (&child == RemoteElem::get_instance() ||
                child.processor_id() == parent_procid)
              matching_child_id = true;
          libmesh_assert(matching_child_id);
        }
      }
#endif
    }

  static void parallel_consistent_procids(const MeshBaseTempl<RealType> & mesh)
    {
      LOG_SCOPE("libmesh_assert_parallel_consistent_procids()", "MeshTools");

      if (mesh.n_processors() == 1)
        return;

      libmesh_parallel_only(mesh.comm());

      // Some code (looking at you, stitch_meshes) modifies DofObject ids
      // without keeping max_elem_id()/max_node_id() consistent, but
      // that's done in a safe way for performance reasons, so we'll play
      // along and just figure out new max ids ourselves.
      dof_id_type parallel_max_elem_id = 0;
      for (const auto & elem : mesh.element_ptr_range())
        parallel_max_elem_id = std::max<dof_id_type>(parallel_max_elem_id,
                                                     elem->id()+1);
      mesh.comm().max(parallel_max_elem_id);

      // Check processor ids for consistency between processors

      for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
      {
        const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);

        processor_id_type min_id =
          elem ? elem->processor_id() :
          std::numeric_limits<processor_id_type>::max();
        mesh.comm().min(min_id);

        processor_id_type max_id =
          elem ? elem->processor_id() :
          std::numeric_limits<processor_id_type>::min();
        mesh.comm().max(max_id);

        if (elem)
        {
          libmesh_assert_equal_to (min_id, elem->processor_id());
          libmesh_assert_equal_to (max_id, elem->processor_id());
        }

        if (min_id == mesh.processor_id())
          libmesh_assert(elem);
      }
    }

};

template <typename RealType>
struct LibMeshAssertShim<NodeTempl<RealType>>
{
  static void topology_consistent_procids(const MeshBaseTempl<RealType> & mesh)
    {
      LOG_SCOPE("libmesh_assert_topology_consistent_procids()", "MeshTools");

      if (mesh.n_processors() == 1)
        return;

      libmesh_parallel_only(mesh.comm());

      // We want this test to be valid even when called even after nodes
      // have been added asynchronously but before they're renumbered.
      //
      // Plus, some code (looking at you, stitch_meshes) modifies
      // DofObject ids without keeping max_elem_id()/max_node_id()
      // consistent, but that's done in a safe way for performance
      // reasons, so we'll play along and just figure out new max ids
      // ourselves.
      dof_id_type parallel_max_node_id = 0;
      for (const auto & node : mesh.node_ptr_range())
        parallel_max_node_id = std::max<dof_id_type>(parallel_max_node_id,
                                                     node->id()+1);
      mesh.comm().max(parallel_max_node_id);


      std::vector<bool> node_touched_by_me(parallel_max_node_id, false);

      for (const auto & elem : as_range(mesh.local_elements_begin(),
                                        mesh.local_elements_end()))
      {
        libmesh_assert (elem);

        for (auto & node : elem->node_ref_range())
        {
          dof_id_type nodeid = node.id();
          node_touched_by_me[nodeid] = true;
        }
      }
      std::vector<bool> node_touched_by_anyone(node_touched_by_me);
      mesh.comm().max(node_touched_by_anyone);

      for (const auto & node : mesh.local_node_ptr_range())
      {
        libmesh_assert(node);
        dof_id_type nodeid = node->id();
        libmesh_assert(!node_touched_by_anyone[nodeid] ||
                       node_touched_by_me[nodeid]);
      }
    }

  static void parallel_consistent_procids(const MeshBaseTempl<RealType> & mesh)
    {
      LOG_SCOPE("libmesh_assert_parallel_consistent_procids()", "MeshTools");

      if (mesh.n_processors() == 1)
        return;

      libmesh_parallel_only(mesh.comm());

      // We want this test to be valid even when called even after nodes
      // have been added asynchronously but before they're renumbered
      //
      // Plus, some code (looking at you, stitch_meshes) modifies
      // DofObject ids without keeping max_elem_id()/max_node_id()
      // consistent, but that's done in a safe way for performance
      // reasons, so we'll play along and just figure out new max ids
      // ourselves.
      dof_id_type parallel_max_node_id = 0;
      for (const auto & node : mesh.node_ptr_range())
        parallel_max_node_id = std::max<dof_id_type>(parallel_max_node_id,
                                                     node->id()+1);
      mesh.comm().max(parallel_max_node_id);

      std::vector<bool> node_touched_by_anyone(parallel_max_node_id, false);

      for (const auto & elem : as_range(mesh.local_elements_begin(),
                                        mesh.local_elements_end()))
      {
        libmesh_assert (elem);

        for (auto & node : elem->node_ref_range())
        {
          dof_id_type nodeid = node.id();
          node_touched_by_anyone[nodeid] = true;
        }
      }
      mesh.comm().max(node_touched_by_anyone);

      // Check processor ids for consistency between processors
      // on any node an element touches
      for (dof_id_type i=0; i != parallel_max_node_id; ++i)
      {
        if (!node_touched_by_anyone[i])
          continue;

        const NodeTempl<RealType> * node = mesh.query_node_ptr(i);
        const processor_id_type pid = node ? node->processor_id() : 0;

        libmesh_assert(mesh.comm().semiverify (node ? &pid : nullptr));
      }
    }
};
}


template <typename DofObjectSubclass, typename RealType>
void libmesh_assert_topology_consistent_procids(const MeshBaseTempl<RealType> & mesh)
{
  LibMeshAssertShim<DofObjectSubclass>::topology_consistent_procids(mesh);
}



template <typename DofObjectSubclass, typename RealType>
void libmesh_assert_parallel_consistent_procids(const MeshBaseTempl<RealType> & mesh)
{
  LibMeshAssertShim<DofObjectSubclass>::parallel_consistent_procids(mesh);
}

template <typename RealType>
void libmesh_assert_parallel_consistent_new_node_procids(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_parallel_consistent_new_node_procids()", "MeshTools");

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  // We want this test to hit every node when called even after nodes
  // have been added asynchronously but before everything has been
  // renumbered.
  dof_id_type parallel_max_elem_id = mesh.max_elem_id();
  mesh.comm().max(parallel_max_elem_id);

  std::vector<bool> elem_touched_by_anyone(parallel_max_elem_id, false);

  for (dof_id_type i=0; i != parallel_max_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_nodes = elem ? elem->n_nodes() : 0;
      unsigned int n_nodes = my_n_nodes;
      mesh.comm().max(n_nodes);

      if (n_nodes)
        libmesh_assert(mesh.comm().semiverify(elem ? &my_n_nodes : nullptr));

      for (unsigned int n=0; n != n_nodes; ++n)
        {
          const NodeTempl<RealType> * node = elem ? elem->node_ptr(n) : nullptr;
          const processor_id_type pid = node ? node->processor_id() : 0;
          libmesh_assert(mesh.comm().semiverify (node ? &pid : nullptr));
        }
    }
}

template <typename RealType>
void libmesh_assert_canonical_node_procids (const MeshBaseTempl<RealType> & mesh)
{
  for (const auto & elem : mesh.active_element_ptr_range())
    for (auto & node : elem->node_ref_range())
      libmesh_assert_equal_to
        (node.processor_id(),
         node.choose_processor_id(node.processor_id(),
                                  elem->processor_id()));
}



} // namespace MeshTools



#ifdef LIBMESH_ENABLE_AMR
template <typename RealType>
void MeshTools::libmesh_assert_valid_refinement_flags(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_refinement_flags()", "MeshTools");

  libmesh_parallel_only(mesh.comm());
  if (mesh.n_processors() == 1)
    return;

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  std::vector<unsigned char> my_elem_h_state(pmax_elem_id, 255);
  std::vector<unsigned char> my_elem_p_state(pmax_elem_id, 255);

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      dof_id_type elemid = elem->id();

      my_elem_h_state[elemid] =
        static_cast<unsigned char>(elem->refinement_flag());

      my_elem_p_state[elemid] =
        static_cast<unsigned char>(elem->p_refinement_flag());
    }
  std::vector<unsigned char> min_elem_h_state(my_elem_h_state);
  mesh.comm().min(min_elem_h_state);

  std::vector<unsigned char> min_elem_p_state(my_elem_p_state);
  mesh.comm().min(min_elem_p_state);

  for (dof_id_type i=0; i!= pmax_elem_id; ++i)
    {
      libmesh_assert(my_elem_h_state[i] == 255 ||
                     my_elem_h_state[i] == min_elem_h_state[i]);
      libmesh_assert(my_elem_p_state[i] == 255 ||
                     my_elem_p_state[i] == min_elem_p_state[i]);
    }
}
#else
template <typename RealType>
void MeshTools::libmesh_assert_valid_refinement_flags(const MeshBaseTempl<RealType> &)
{
}
#endif // LIBMESH_ENABLE_AMR



#ifdef LIBMESH_ENABLE_AMR
template <typename RealType>
void MeshTools::libmesh_assert_valid_refinement_tree(const MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("libmesh_assert_valid_refinement_tree()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert(elem);
      if (elem->has_children())
        for (auto & child : elem->child_ref_range())
          if (&child != RemoteElem::get_instance())
            libmesh_assert_equal_to (child.parent(), elem);
      if (elem->active())
        {
          libmesh_assert(!elem->ancestor());
          libmesh_assert(!elem->subactive());
        }
      else if (elem->ancestor())
        {
          libmesh_assert(!elem->subactive());
        }
      else
        libmesh_assert(elem->subactive());

      if (elem->p_refinement_flag() == Elem::JUST_REFINED)
        libmesh_assert_greater(elem->p_level(), 0);
    }
}
#else
template <typename RealType>
void MeshTools::libmesh_assert_valid_refinement_tree(const MeshBaseTempl<RealType> &)
{
}
#endif // LIBMESH_ENABLE_AMR



template <typename RealType>
void MeshTools::libmesh_assert_valid_neighbors(const MeshBaseTempl<RealType> & mesh,
                                               bool assert_valid_remote_elems)
{
  LOG_SCOPE("libmesh_assert_valid_neighbors()", "MeshTools");

  for (const auto & elem : mesh.element_ptr_range())
    {
      libmesh_assert (elem);
      elem->libmesh_assert_valid_neighbors();
    }

  if (mesh.n_processors() == 1)
    return;

  libmesh_parallel_only(mesh.comm());

  dof_id_type pmax_elem_id = mesh.max_elem_id();
  mesh.comm().max(pmax_elem_id);

  for (dof_id_type i=0; i != pmax_elem_id; ++i)
    {
      const ElemTempl<RealType> * elem = mesh.query_elem_ptr(i);

      const unsigned int my_n_neigh = elem ? elem->n_neighbors() : 0;
      unsigned int n_neigh = my_n_neigh;
      mesh.comm().max(n_neigh);
      if (elem)
        libmesh_assert_equal_to (my_n_neigh, n_neigh);

      for (unsigned int n = 0; n != n_neigh; ++n)
        {
          dof_id_type my_neighbor = DofObject::invalid_id;
          dof_id_type * p_my_neighbor = nullptr;

          // If we have a non-remote_elem neighbor link, then we can
          // verify it.
          if (elem && elem->neighbor_ptr(n) != RemoteElem::get_instance())
            {
              p_my_neighbor = &my_neighbor;
              if (elem->neighbor_ptr(n))
                my_neighbor = elem->neighbor_ptr(n)->id();

              // But wait - if we haven't set remote_elem links yet then
              // some nullptr links on ghost elements might be
              // future-remote_elem links, so we can't verify those.
              if (!assert_valid_remote_elems &&
                  !elem->neighbor_ptr(n) &&
                  elem->processor_id() != mesh.processor_id())
                p_my_neighbor = nullptr;
            }
          libmesh_assert(mesh.comm().semiverify(p_my_neighbor));
        }
    }
}



#endif // DEBUG



// Functors for correct_node_proc_ids
namespace {

typedef std::unordered_map<dof_id_type, processor_id_type> proc_id_map_type;

template <typename RealType = Real>
struct SyncNodeSet
{
  typedef unsigned char datum; // bool but without bit twiddling issues
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef NodeTempl<RealType> Node;

  SyncNodeSet(std::unordered_set<const Node *> & _set,
              MeshBase & _mesh) :
    node_set(_set), mesh(_mesh) {}

  std::unordered_set<const Node *> & node_set;

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find whether each requested node belongs in the set
    data.resize(ids.size());

    for (auto i : index_range(ids))
      {
        const dof_id_type id = ids[i];

        // We'd better have every node we're asked for
        Node * node = mesh.node_ptr(id);

        // Return if the node is in the set.
        data[i] = (node_set.find(node) != node_set.end());
      }
  }

  // ------------------------------------------------------------
  bool act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> in_set)
  {
    bool data_changed = false;

    // Add nodes we've been informed of to our own set
    for (auto i : index_range(ids))
      {
        if (in_set[i])
          {
            Node * node = mesh.node_ptr(ids[i]);
            if (!node_set.count(node))
              {
                node_set.insert(node);
                data_changed = true;
              }
          }
      }

    return data_changed;
  }
};

template <typename RealType = Real>
struct NodesNotInSet
{
  typedef NodeTempl<RealType> Node;

  NodesNotInSet(const std::unordered_set<const Node *> _set)
    : node_set(_set) {}

  bool operator() (const Node * node) const
  {
    if (node_set.count(node))
      return false;
    return true;
  }

  const std::unordered_set<const Node *> node_set;
};

template <typename RealType = Real>
struct SyncProcIdsFromMap
{
  typedef processor_id_type datum;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef NodeTempl<RealType> Node;

  SyncProcIdsFromMap(const proc_id_map_type & _map,
                     MeshBase & _mesh) :
    new_proc_ids(_map), mesh(_mesh) {}

  const proc_id_map_type & new_proc_ids;

  MeshBase & mesh;

  // ------------------------------------------------------------
  void gather_data (const std::vector<dof_id_type> & ids,
                    std::vector<datum> & data)
  {
    // Find the new processor id of each requested node
    data.resize(ids.size());

    for (auto i : index_range(ids))
      {
        const dof_id_type id = ids[i];
        const proc_id_map_type::const_iterator it = new_proc_ids.find(id);

        // Return the node's new processor id if it has one, or its
        // old processor id if not.
        if (it != new_proc_ids.end())
          data[i] = it->second;
        else
          {
            // We'd better find every node we're asked for
            const Node & node = mesh.node_ref(id);
            data[i] = node.processor_id();
          }
      }
  }

  // ------------------------------------------------------------
  void act_on_data (const std::vector<dof_id_type> & ids,
                    const std::vector<datum> proc_ids)
  {
    // Set the node processor ids we've now been informed of
    for (auto i : index_range(ids))
      {
        Node & node = mesh.node_ref(ids[i]);
        node.processor_id() = proc_ids[i];
      }
  }
};
}



template <typename RealType>
void MeshTools::correct_node_proc_ids (MeshBaseTempl<RealType> & mesh)
{
  LOG_SCOPE("correct_node_proc_ids()","MeshTools");

  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  // We require all processors to agree on nodal processor ids before
  // going into this algorithm.
#ifdef DEBUG
  MeshTools::libmesh_assert_parallel_consistent_procids<Node>(mesh);
#endif

  // If we have any unpartitioned elements at this
  // stage there is a problem
  libmesh_assert (MeshTools::n_elem(mesh.unpartitioned_elements_begin(),
                                    mesh.unpartitioned_elements_end()) == 0);

  // Fix nodes' processor ids.  Coarsening may have left us with nodes
  // which are no longer touched by any elements of the same processor
  // id, and for DofMap to work we need to fix that.

  // This is harder now that libMesh no longer requires a distributed
  // mesh to ghost all nodal neighbors: it is possible for two active
  // elements on two different processors to share the same node in
  // such a way that neither processor knows the others' element
  // exists!

  // While we're at it, if this mesh is configured to allow
  // repartitioning, we'll repartition *all* the nodes' processor ids
  // using the canonical Node heuristic, to try and improve DoF load
  // balancing.  But if the mesh is disallowing repartitioning, we
  // won't touch processor_id on any node where it's valid, regardless
  // of whether or not it's canonical.
  bool repartition_all_nodes = !mesh.skip_noncritical_partitioning();
  std::unordered_set<const NodeTempl<RealType> *> valid_nodes;

  // If we aren't allowed to repartition, then we're going to leave
  // every node we can at its current processor_id, and *only*
  // repartition the nodes whose current processor id is incompatible
  // with DoFMap (because it doesn't touch an active element, e.g. due
  // to coarsening)
  if (!repartition_all_nodes)
    {
      for (const auto & elem : mesh.active_element_ptr_range())
        for (const auto & node : elem->node_ref_range())
          if (elem->processor_id() == node.processor_id())
            valid_nodes.insert(&node);

      SyncNodeSet<RealType> syncv(valid_nodes, mesh);

      Parallel::sync_dofobject_data_by_id
        (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), syncv);
    }

  // We build up a set of compatible processor ids for each node
  proc_id_map_type new_proc_ids;

  for (auto & elem : mesh.active_element_ptr_range())
    {
      processor_id_type pid = elem->processor_id();

      for (auto & node : elem->node_ref_range())
        {
          const dof_id_type id = node.id();
          const proc_id_map_type::iterator it = new_proc_ids.find(id);
          if (it == new_proc_ids.end())
            new_proc_ids.insert(std::make_pair(id,pid));
          else
            it->second = node.choose_processor_id(it->second, pid);
        }
    }

  // Sort the new pids to push to each processor
  std::map<processor_id_type, std::vector<std::pair<dof_id_type, processor_id_type>>>
    ids_to_push;

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type id = node->id();
      const proc_id_map_type::iterator it = new_proc_ids.find(id);
      if (it == new_proc_ids.end())
        continue;
      const processor_id_type pid = it->second;
      if (node->processor_id() != DofObject::invalid_processor_id)
        ids_to_push[node->processor_id()].push_back(std::make_pair(id, pid));
    }

  auto action_functor =
    [& mesh, & new_proc_ids]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, processor_id_type>> & data)
    {
      for (auto & p : data)
        {
          const dof_id_type id = p.first;
          const processor_id_type pid = p.second;
          const proc_id_map_type::iterator it = new_proc_ids.find(id);
          if (it == new_proc_ids.end())
            new_proc_ids.insert(std::make_pair(id,pid));
          else
            {
              const NodeTempl<RealType> & node = mesh.node_ref(id);
              it->second = node.choose_processor_id(it->second, pid);
            }
        }
    };

  Parallel::push_parallel_vector_data
    (mesh.comm(), ids_to_push, action_functor);

  // Now new_proc_ids is correct for every node we used to own.  Let's
  // ask every other processor about the nodes they used to own.  But
  // first we'll need to keep track of which nodes we used to own,
  // lest we get them confused with nodes we newly own.
  std::unordered_set<NodeTempl<RealType> *> ex_local_nodes;
  for (auto & node : mesh.local_node_ptr_range())
    {
      const proc_id_map_type::iterator it = new_proc_ids.find(node->id());
      if (it != new_proc_ids.end() && it->second != mesh.processor_id())
        ex_local_nodes.insert(node);
    }

  SyncProcIdsFromMap<RealType> sync(new_proc_ids, mesh);
  if (repartition_all_nodes)
    Parallel::sync_dofobject_data_by_id
      (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), sync);
  else
    {
      NodesNotInSet<RealType> nnis(valid_nodes);

      Parallel::sync_dofobject_data_by_id
        (mesh.comm(), mesh.nodes_begin(), mesh.nodes_end(), nnis, sync);
    }

  // And finally let's update the nodes we used to own.
  for (const auto & node : ex_local_nodes)
    {
      if (valid_nodes.count(node))
        continue;

      const dof_id_type id = node->id();
      const proc_id_map_type::iterator it = new_proc_ids.find(id);
      libmesh_assert(it != new_proc_ids.end());
      node->processor_id() = it->second;
    }

  // We should still have consistent nodal processor ids coming out of
  // this algorithm, but if we're allowed to repartition the mesh then
  // they should be canonically correct too.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_procids<Node>(mesh);
  //if (repartition_all_nodes)
  //  MeshTools::libmesh_assert_canonical_node_procids(mesh);
#endif
}



template <typename RealType>
void MeshTools::Private::globally_renumber_nodes_and_elements (MeshBaseTempl<RealType> & mesh)
{
  MeshCommunication().assign_global_indices(mesh);
}

} // namespace libMesh

#define INSTANTIATE_MESH_TOOLS0(RealType)                                                          \
  template dof_id_type total_weight(const MeshBaseTempl<RealType> & mesh);                         \
  template dof_id_type weight(const MeshBaseTempl<RealType> & mesh, const processor_id_type pid);  \
  template void build_nodes_to_elem_map(                                                           \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      std::vector<std::vector<dof_id_type>> & nodes_to_elem_map);                                  \
  template void build_nodes_to_elem_map(                                                           \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map);                  \
  template void build_nodes_to_elem_map(                                                           \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      std::unordered_map<dof_id_type, std::vector<dof_id_type>> & nodes_to_elem_map);              \
  template void build_nodes_to_elem_map(                                                           \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> &                  \
          nodes_to_elem_map);                                                                      \
  template std::unordered_set<dof_id_type> find_boundary_nodes(                                    \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template std::unordered_set<dof_id_type> find_block_boundary_nodes(                              \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template libMesh::BoundingBoxTempl<RealType> create_bounding_box(                                \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template Sphere bounding_sphere(const MeshBaseTempl<RealType> & mesh);                           \
  template libMesh::BoundingBoxTempl<RealType> create_nodal_bounding_box(                          \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template libMesh::BoundingBoxTempl<RealType> create_local_bounding_box(                          \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template libMesh::BoundingBoxTempl<RealType> create_processor_bounding_box(                      \
      const MeshBaseTempl<RealType> & mesh, const processor_id_type pid);                          \
  template Sphere processor_bounding_sphere(const MeshBaseTempl<RealType> & mesh,                  \
                                            const processor_id_type pid);                          \
  template libMesh::BoundingBoxTempl<RealType> create_subdomain_bounding_box(                      \
      const MeshBaseTempl<RealType> & mesh, const subdomain_id_type sid);                          \
  template Sphere subdomain_bounding_sphere(const MeshBaseTempl<RealType> & mesh,                  \
                                            const subdomain_id_type sid);                          \
  template void elem_types(const MeshBaseTempl<RealType> & mesh, std::vector<ElemType> & et);      \
  template dof_id_type n_elem_of_type(const MeshBaseTempl<RealType> & mesh, const ElemType type);  \
  template dof_id_type n_active_elem_of_type(const MeshBaseTempl<RealType> & mesh,                 \
                                             const ElemType type);                                 \
  template dof_id_type n_non_subactive_elem_of_type_at_level(                                      \
      const MeshBaseTempl<RealType> & mesh, const ElemType type, const unsigned int level);        \
  template unsigned int n_levels(const MeshBaseTempl<RealType> & mesh);                            \
  template unsigned int n_local_levels(const MeshBaseTempl<RealType> & mesh);                      \
  template unsigned int n_active_levels(const MeshBaseTempl<RealType> & mesh);                     \
  template unsigned int n_active_local_levels(const MeshBaseTempl<RealType> & mesh);               \
  template unsigned int n_p_levels(const MeshBaseTempl<RealType> & mesh);                          \
  template unsigned int paranoid_n_levels(const MeshBaseTempl<RealType> & mesh);                   \
  template void get_not_subactive_node_ids(const MeshBaseTempl<RealType> & mesh,                   \
                                           std::set<dof_id_type> & not_subactive_node_ids);        \
  template dof_id_type n_elem(                                                                     \
      const typename MeshBaseTempl<RealType>::const_element_iterator & begin,                      \
      const typename MeshBaseTempl<RealType>::const_element_iterator & end);                       \
  template dof_id_type n_nodes(                                                                    \
      const typename MeshBaseTempl<RealType>::const_node_iterator & begin,                         \
      const typename MeshBaseTempl<RealType>::const_node_iterator & end);                          \
  template void find_nodal_neighbors(                                                              \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      const NodeTempl<RealType> & n,                                                               \
      const std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map,             \
      std::vector<const NodeTempl<RealType> *> & neighbors);                                       \
  template void find_nodal_neighbors(                                                              \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      const NodeTempl<RealType> & n,                                                               \
      const std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> &            \
          nodes_to_elem_map,                                                                       \
      std::vector<const NodeTempl<RealType> *> & neighbors);                                       \
  template void find_hanging_nodes_and_parents(                                                    \
      const MeshBaseTempl<RealType> & mesh,                                                        \
      std::map<dof_id_type, std::vector<dof_id_type>> & hanging_nodes);                            \
  template void correct_node_proc_ids(MeshBaseTempl<RealType> &);                                  \
  template void Private::globally_renumber_nodes_and_elements(MeshBaseTempl<RealType> &)

#ifdef LIBMESH_ENABLE_DEPRECATED
#define INSTANTIATE_MESH_TOOLS1(RealType)                                                          \
  INSTANTIATE_MESH_TOOLS0(RealType);                                                               \
  template void find_boundary_nodes(const MeshBaseTempl<RealType> & mesh,                          \
                                    std::vector<bool> & on_boundary)
#else
#define INSTANTIATE_MESH_TOOLS1(RealType) INSTANTIATE_MESH_TOOLS0(RealType)
#endif

#ifdef DEBUG
#define INSTANTIATE_MESH_TOOLS2(RealType)                                                          \
  INSTANTIATE_MESH_TOOLS1(RealType);                                                               \
  template void libmesh_assert_no_links_to_elem(const MeshBaseTempl<RealType> & mesh,              \
                                                const ElemTempl<RealType> * bad_elem);             \
  template void libmesh_assert_equal_n_systems(const MeshBaseTempl<RealType> & mesh);              \
  template void libmesh_assert_old_dof_objects(const MeshBaseTempl<RealType> & mesh);              \
  template void libmesh_assert_valid_node_pointers(const MeshBaseTempl<RealType> & mesh);          \
  template void libmesh_assert_valid_remote_elems(const MeshBaseTempl<RealType> & mesh);           \
  template void libmesh_assert_valid_elem_ids(const MeshBaseTempl<RealType> & mesh);               \
  template void libmesh_assert_valid_amr_elem_ids(const MeshBaseTempl<RealType> & mesh);           \
  template void libmesh_assert_valid_amr_interior_parents(const MeshBaseTempl<RealType> & mesh);   \
  template void libmesh_assert_connected_nodes(const MeshBaseTempl<RealType> & mesh);              \
  template void libmesh_assert_valid_boundary_ids(const MeshBaseTempl<RealType> & mesh);           \
  template void libmesh_assert_valid_dof_ids(const MeshBaseTempl<RealType> & mesh,                 \
                                             unsigned int sysnum = libMesh::invalid_uint);         \
  template void libmesh_assert_contiguous_dof_ids(const MeshBaseTempl<RealType> & mesh,            \
                                                  unsigned int sysnum);                            \
  template void libmesh_assert_consistent_distributed(const MeshBaseTempl<RealType> & mesh);       \
  template void libmesh_assert_consistent_distributed_nodes(const MeshBaseTempl<RealType> & mesh); \
  template void libmesh_assert_parallel_consistent_new_node_procids(                               \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_parallel_consistent_procids<ElemTempl<RealType>>(                   \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_parallel_consistent_procids<NodeTempl<RealType>>(                   \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_topology_consistent_procids<ElemTempl<RealType>>(                   \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_topology_consistent_procids<NodeTempl<RealType>>(                   \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_valid_procids<ElemTempl<RealType>>(                                 \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_valid_procids<NodeTempl<RealType>>(                                 \
      const MeshBaseTempl<RealType> & mesh);                                                       \
  template void libmesh_assert_canonical_node_procids(const MeshBaseTempl<RealType> & mesh);       \
  template void libmesh_assert_valid_refinement_flags(const MeshBaseTempl<RealType> & mesh);       \
  template void libmesh_assert_valid_refinement_tree(const MeshBaseTempl<RealType> & mesh);        \
  template void libmesh_assert_valid_neighbors(const MeshBaseTempl<RealType> & mesh,               \
                                               bool assert_valid_remote_elems = true)
#ifdef LIBMESH_ENABLE_UNIQUE_ID
#define INSTANTIATE_MESH_TOOLS(RealType)                                                           \
  INSTANTIATE_MESH_TOOLS2(RealType);                                                               \
  template void libmesh_assert_valid_unique_ids(const MeshBaseTempl<RealType> & mesh)
#else
#define INSTANTIATE_MESH_TOOLS(RealType) INSTANTIATE_MESH_TOOLS2(RealType)
#endif // LIBMESH_ENABLE_UNIQUE_ID

#else

#define INSTANTIATE_MESH_TOOLS(RealType) INSTANTIATE_MESH_TOOLS1(RealType)

#endif // DEBUG

#endif // LIBMESH_MESH_TOOLS_IMPL_H
