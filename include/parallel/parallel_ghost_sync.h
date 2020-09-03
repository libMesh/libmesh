// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARALLEL_GHOST_SYNC_H
#define LIBMESH_PARALLEL_GHOST_SYNC_H

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/int_range.h"
#include "libmesh/location_maps.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_algebra.h"

// TIMPI includes
#include "timpi/communicator.h"
#include "timpi/parallel_sync.h"

// C++ includes
#include <map> // FIXME - pid > comm.size() breaks with unordered_map
#include <vector>


namespace libMesh
{



//--------------------------------------------------------------------------
namespace Parallel {

//------------------------------------------------------------------------
/**
 * Request data about a range of ghost nodes uniquely identified by
 * their xyz location or a range of active ghost elements uniquely
 * identified by their centroids' xyz location.  Fulfill requests
 * with
 * sync.gather_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 * The user-provided location_map will be used and left unchanged
 * if it is provided, or filled and cleared if it is empty.
 */
template <typename Iterator,
          typename DofObjType,
          typename SyncFunctor>
void sync_dofobject_data_by_xyz(const Communicator &      comm,
                                const Iterator &          range_begin,
                                const Iterator &          range_end,
                                LocationMap<DofObjType> * location_map,
                                SyncFunctor &             sync);

//------------------------------------------------------------------------
/**
 * Request data about a range of ghost dofobjects uniquely
 * identified by their id.  Fulfill requests with
 * sync.gather_data(const std::vector<dof_id_type> & ids,
 *                  std::vector<sync::datum> & data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<dof_id_type> & ids,
 *                  std::vector<sync::datum> & data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 */
template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator &     range_begin,
                               const Iterator &     range_end,
                               SyncFunctor &        sync);

/**
 * Request data about a range of ghost dofobjects uniquely
 * identified by their id.
 *
 * Elements within the range can be excluded from the request by
 * returning false from dofobj_check(dof_object)
 */
template <typename Iterator,
          typename DofObjectCheckFunctor,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator & range_begin,
                               const Iterator & range_end,
                               const DofObjectCheckFunctor & dofobj_check,
                               SyncFunctor &    sync);

//------------------------------------------------------------------------
/**
 * Request data about a range of ghost elements uniquely
 * identified by their parent id and which child they are.
 * Fulfill requests with
 * sync.gather_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 */
template <typename Iterator,
          typename SyncFunctor>
void sync_element_data_by_parent_id(MeshBase &       mesh,
                                    const Iterator & range_begin,
                                    const Iterator & range_end,
                                    SyncFunctor &    sync);

//------------------------------------------------------------------------
/**
 * Synchronize data about a range of ghost nodes uniquely identified
 * by an element id and local node id, assuming a single
 * synchronization pass is necessary.
 *
 * Data for all nodes connected to elements in the given range of
 * *element* iterators will be requested.
 *
 * Elements can be further excluded from the request by returning
 * false from element_check(elem)
 *
 * Nodes can be further excluded from the request by returning false
 * from node_check(elem, local_node_num)
 *
 * Fulfill requests with
 * sync.gather_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * bool sync.act_on_data(const std::vector<unsigned int> & ids,
 *                       std::vector<sync::datum> & data)
 * and return true iff the response changed any data.
 *
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 *
 * This method returns true iff the sync pass changed any data on any
 * processor.
 */
template <typename ElemCheckFunctor,
          typename NodeCheckFunctor,
          typename SyncFunctor>
bool sync_node_data_by_element_id_once(MeshBase & mesh,
                                       const MeshBase::const_element_iterator & range_begin,
                                       const MeshBase::const_element_iterator & range_end,
                                       const ElemCheckFunctor & elem_check,
                                       const NodeCheckFunctor & node_check,
                                       SyncFunctor & sync);



//------------------------------------------------------------------------
/**
 * Synchronize data about a range of ghost nodes uniquely identified
 * by an element id and local node id, iterating until data is
 * completely in sync and futher synchronization passes cause no
 * changes.
 *
 * Imagine a vertex surrounded by triangles, each on a different
 * processor, with a ghosting policy that include only face neighbors
 * and not point neighbors.  Then the only way for authoritative
 * information to trickle out from that vertex is by being passed
 * along, one neighbor at a time, to processors who mostly don't even
 * see the node's true owner!
 *
 * Data for all nodes connected to elements in the given range of
 * *element* iterators will be requested.
 *
 * Elements can be further excluded from the request by returning
 * false from element_check(elem)
 *
 * Nodes can be further excluded from the request by returning false
 * from node_check(elem, local_node_num)
 *
 * Fulfill requests with
 * sync.gather_data(const std::vector<unsigned int> & ids,
 *                  std::vector<sync::datum> & data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * bool sync.act_on_data(const std::vector<unsigned int> & ids,
 *                       std::vector<sync::datum> & data)
 * and return true iff the response changed any data.
 *
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 */
template <typename ElemCheckFunctor,
          typename NodeCheckFunctor,
          typename SyncFunctor>
void sync_node_data_by_element_id(MeshBase & mesh,
                                  const MeshBase::const_element_iterator & range_begin,
                                  const MeshBase::const_element_iterator & range_end,
                                  const ElemCheckFunctor & elem_check,
                                  const NodeCheckFunctor & node_check,
                                  SyncFunctor & sync);


//------------------------------------------------------------------------
// Parallel members


// "Check" Functor to perform sync operations with no exclusions
struct SyncEverything
{
  SyncEverything() {}

  bool operator() (const DofObject *) const { return true; }

  bool operator() (const Elem *, unsigned int) const
  { return true; }
};



template <typename Iterator,
          typename DofObjType,
          typename SyncFunctor>
void sync_dofobject_data_by_xyz(const Communicator & comm,
                                const Iterator & range_begin,
                                const Iterator & range_end,
                                LocationMap<DofObjType> & location_map,
                                SyncFunctor & sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  // We need a valid location_map
#ifdef DEBUG
  bool need_map_update = (range_begin != range_end && location_map.empty());
  comm.max(need_map_update);
  libmesh_assert(!need_map_update);
#endif

  // Count the objects to ask each processor about
  std::map<processor_id_type, dof_id_type>
    ghost_objects_from_proc;

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType * obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::map<processor_id_type, std::vector<Point>>
    requested_objs_pt;
  // Corresponding ids to keep track of
  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_objs_id;

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (auto pair : ghost_objects_from_proc)
    {
      const processor_id_type p = pair.first;
      if (p != comm.rank())
        {
          requested_objs_pt[p].reserve(pair.second);
          requested_objs_id[p].reserve(pair.second);
        }
    }

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType * obj = *it;
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      Point p = location_map.point_of(*obj);
      requested_objs_pt[obj_procid].push_back(p);
      requested_objs_id[obj_procid].push_back(obj->id());
    }

  std::map<const std::vector<Point> *, processor_id_type>
    requested_objs_pt_inv;
  for (auto & pair : requested_objs_pt)
    requested_objs_pt_inv[&pair.second] = pair.first;

  auto gather_functor =
    [&location_map, &sync]
    (processor_id_type /*pid*/, const std::vector<Point> & pts,
     std::vector<typename SyncFunctor::datum> & data)
    {
      // Find the local id of each requested object
      std::size_t query_size = pts.size();
      std::vector<dof_id_type> query_id(query_size);
      for (std::size_t i=0; i != query_size; ++i)
        {
          Point pt = pts[i];

          // Look for this object in the multimap
          DofObjType * obj = location_map.find(pt);

          // We'd better find every object we're asked for
          libmesh_assert (obj);

          // Return the object's correct processor id,
          // and our (correct if it's local) id for it.
          query_id[i] = obj->id();
        }

      // Gather whatever data the user wants
      sync.gather_data(query_id, data);
    };

  auto action_functor =
    [&sync, &requested_objs_id,
     &requested_objs_pt_inv]
    (processor_id_type /* pid */, const std::vector<Point> & point_request,
     const std::vector<typename SyncFunctor::datum> & data)
    {
      // With splits working on more pids than ranks, query_pid may not equal pid
      const processor_id_type query_pid =
        requested_objs_pt_inv[&point_request];

      // Let the user process the results
      sync.act_on_data(requested_objs_id[query_pid], data);
    };

  // Trade requests with other processors
  typename SyncFunctor::datum * ex = nullptr;
  pull_parallel_vector_data
    (comm, requested_objs_pt, gather_functor, action_functor, ex);
}



template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator & range_begin,
                               const Iterator & range_end,
                               SyncFunctor &    sync)
{
  sync_dofobject_data_by_id(comm, range_begin, range_end, SyncEverything(), sync);
}

template <typename Iterator,
          typename DofObjectCheckFunctor,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator & range_begin,
                               const Iterator & range_end,
                               const DofObjectCheckFunctor & dofobj_check,
                               SyncFunctor &    sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  // Count the objects to ask each processor about
  std::map<processor_id_type, dof_id_type>
    ghost_objects_from_proc;

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject * obj = *it;
      libmesh_assert (obj);

      // We may want to pass Elem* or Node* to the check function, not
      // just DofObject*
      if (!dofobj_check(*it))
        continue;

      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_objs_id;

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (auto pair : ghost_objects_from_proc)
    {
      const processor_id_type p = pair.first;
      if (p != comm.rank())
        requested_objs_id[p].reserve(pair.second);
    }

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject * obj = *it;

      if (!dofobj_check(*it))
        continue;

      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      requested_objs_id[obj_procid].push_back(obj->id());
    }

  auto gather_functor =
    [&sync]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     std::vector<typename SyncFunctor::datum> & data)
    {
      sync.gather_data(ids, data);
    };

  auto action_functor =
    [&sync]
    (processor_id_type, const std::vector<dof_id_type> & ids,
     const std::vector<typename SyncFunctor::datum> & data)
    {
      // Let the user process the results
      sync.act_on_data(ids, data);
    };

  // Trade requests with other processors
  typename SyncFunctor::datum * ex = nullptr;
  pull_parallel_vector_data
    (comm, requested_objs_id, gather_functor, action_functor, ex);
}



// If there's no refined elements, there's nothing to sync
#ifdef LIBMESH_ENABLE_AMR
template <typename Iterator,
          typename SyncFunctor>
void sync_element_data_by_parent_id(MeshBase &       mesh,
                                    const Iterator & range_begin,
                                    const Iterator & range_end,
                                    SyncFunctor &    sync)
{
  const Communicator & comm (mesh.comm());

  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  // Count the objects to ask each processor about
  std::map<processor_id_type, dof_id_type>
    ghost_objects_from_proc;

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      Elem * elem = *it;
      processor_id_type obj_procid = elem->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;
      const Elem * parent = elem->parent();
      if (!parent || !elem->active())
        continue;

      ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::map<processor_id_type, std::vector<dof_id_type>>
    requested_objs_id;
  std::map<processor_id_type, std::vector<std::pair<dof_id_type,unsigned char>>>
    requested_objs_parent_id_child_num;

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (auto pair : ghost_objects_from_proc)
    {
      const processor_id_type p = pair.first;
      if (p != comm.rank())
        {
          requested_objs_id[p].reserve(pair.second);
          requested_objs_parent_id_child_num[p].reserve(pair.second);
        }
    }

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      Elem * elem = *it;
      processor_id_type obj_procid = elem->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;
      const Elem * parent = elem->parent();
      if (!parent || !elem->active())
        continue;

      requested_objs_id[obj_procid].push_back(elem->id());
      requested_objs_parent_id_child_num[obj_procid].emplace_back
        (parent->id(), cast_int<unsigned char>(parent->which_child_am_i(elem)));
    }

  std::map<const std::vector<std::pair<dof_id_type,unsigned char>> *, processor_id_type>
    requested_objs_parent_id_child_num_inv;
  for (auto & pair : requested_objs_parent_id_child_num)
    requested_objs_parent_id_child_num_inv[&pair.second] = pair.first;

  auto gather_functor =
    [&mesh, &sync]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, unsigned char>> & parent_id_child_num,
     std::vector<typename SyncFunctor::datum> & data)
    {
      // Find the id of each requested element
      std::size_t query_size = parent_id_child_num.size();
      std::vector<dof_id_type> query_id(query_size);
      for (std::size_t i=0; i != query_size; ++i)
        {
          Elem & parent = mesh.elem_ref(parent_id_child_num[i].first);
          libmesh_assert(parent.has_children());
          Elem * child = parent.child_ptr(parent_id_child_num[i].second);
          libmesh_assert(child);
          libmesh_assert(child->active());
          query_id[i] = child->id();
        }

      // Gather whatever data the user wants
      sync.gather_data(query_id, data);
    };

  auto action_functor =
    [&sync, &requested_objs_id,
     &requested_objs_parent_id_child_num_inv]
    (processor_id_type /* pid */,
     const std::vector<std::pair<dof_id_type, unsigned char>> & parent_id_child_num_request,
     const std::vector<typename SyncFunctor::datum> & data)
    {
      // With splits working on more pids than ranks, query_pid may not equal pid
      const processor_id_type query_pid =
        requested_objs_parent_id_child_num_inv[&parent_id_child_num_request];

      // Let the user process the results
      sync.act_on_data(requested_objs_id[query_pid], data);
    };

  // Trade requests with other processors
  typename SyncFunctor::datum * ex = nullptr;
  pull_parallel_vector_data
    (comm, requested_objs_parent_id_child_num, gather_functor,
     action_functor, ex);
}
#else
template <typename Iterator,
          typename SyncFunctor>
void sync_element_data_by_parent_id(MeshBase &,
                                    const Iterator &,
                                    const Iterator &,
                                    SyncFunctor &)
{
}
#endif // LIBMESH_ENABLE_AMR



template <typename ElemCheckFunctor,
          typename NodeCheckFunctor,
          typename SyncFunctor>
bool sync_node_data_by_element_id_once(MeshBase & mesh,
                                       const MeshBase::const_element_iterator & range_begin,
                                       const MeshBase::const_element_iterator & range_end,
                                       const ElemCheckFunctor & elem_check,
                                       const NodeCheckFunctor & node_check,
                                       SyncFunctor & sync)
{
  const Communicator & comm (mesh.comm());

  // Count the objects to ask each processor about
  std::map<processor_id_type, dof_id_type>
    ghost_objects_from_proc;

  for (const auto & elem : as_range(range_begin, range_end))
    {
      libmesh_assert (elem);

      if (!elem_check(elem))
        continue;

      const processor_id_type proc_id = elem->processor_id();

      bool i_have_elem =
        (proc_id == comm.rank() ||
         proc_id == DofObject::invalid_processor_id);

      if (elem->active() && i_have_elem)
        continue;

      for (auto n : elem->node_index_range())
        {
          if (!node_check(elem, n))
            continue;

          const processor_id_type node_pid =
            elem->node_ref(n).processor_id();

          if (i_have_elem && (node_pid == comm.rank()))
            continue;

          if (i_have_elem)
            {
              libmesh_assert_not_equal_to
                (node_pid, DofObject::invalid_processor_id);
              ghost_objects_from_proc[node_pid]++;
            }
          else
            {
              const processor_id_type request_pid =
                (node_pid == DofObject::invalid_processor_id) ?
                 proc_id : node_pid;
              ghost_objects_from_proc[request_pid]++;
            }
        }
    }

  // Now repeat that iteration, filling request sets this time.

  // Request sets to send to each processor
  std::map<processor_id_type, std::vector<std::pair<dof_id_type, unsigned char>>>
    requested_objs_elem_id_node_num;

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (auto pair : ghost_objects_from_proc)
    {
      const processor_id_type p = pair.first;
      if (p != comm.rank())
        requested_objs_elem_id_node_num[p].reserve(ghost_objects_from_proc[p]);
    }

  for (const auto & elem : as_range(range_begin, range_end))
    {
      libmesh_assert (elem);

      if (!elem_check(elem))
        continue;

      const processor_id_type proc_id = elem->processor_id();

      bool i_have_elem =
        (proc_id == comm.rank() ||
         proc_id == DofObject::invalid_processor_id);

      if (elem->active() && i_have_elem)
        continue;

      const dof_id_type elem_id = elem->id();

      for (auto n : elem->node_index_range())
        {
          if (!node_check(elem, n))
            continue;

          const Node & node = elem->node_ref(n);
          const processor_id_type node_pid = node.processor_id();

          if (i_have_elem && (node_pid == comm.rank()))
            continue;

          if (i_have_elem)
            {
              libmesh_assert_not_equal_to
                (node_pid, DofObject::invalid_processor_id);
              requested_objs_elem_id_node_num[node_pid].emplace_back
                (elem_id, cast_int<unsigned char>(n));
            }
          else
            {
              const processor_id_type request_pid =
                (node_pid == DofObject::invalid_processor_id) ?
                 proc_id : node_pid;
              requested_objs_elem_id_node_num[request_pid].emplace_back
                (elem_id,cast_int<unsigned char>(n));
            }
        }
    }

  auto gather_functor =
    [&mesh, &sync]
    (processor_id_type,
     const std::vector<std::pair<dof_id_type, unsigned char>> & elem_id_node_num,
     std::vector<typename SyncFunctor::datum> & data)
    {
      // Find the id of each requested element
      std::size_t request_size = elem_id_node_num.size();
      std::vector<dof_id_type> query_id(request_size);
      for (std::size_t i=0; i != request_size; ++i)
        {
          // We might now get queries about remote elements, in which
          // case we'll have to ignore them and wait for the query
          // answer to filter to the querier via another source.
          const Elem * elem = mesh.query_elem_ptr(elem_id_node_num[i].first);

          if (elem)
            {
              const unsigned int n = elem_id_node_num[i].second;
              libmesh_assert_less (n, elem->n_nodes());

              const Node & node = elem->node_ref(n);

              // This isn't a safe assertion in the case where we're
              // syncing processor ids
              // libmesh_assert_equal_to (node->processor_id(), comm.rank());

              query_id[i] = node.id();
            }
          else
            query_id[i] = DofObject::invalid_id;
        }

      // Gather whatever data the user wants
      sync.gather_data(query_id, data);
    };

  bool data_changed = false;

  auto action_functor =
    [&sync, &mesh, &data_changed]
    (processor_id_type /* pid */,
     const std::vector<std::pair<dof_id_type, unsigned char>> & elem_id_node_num,
     const std::vector<typename SyncFunctor::datum> & data)
    {
      const std::size_t data_size = data.size();

      libmesh_assert_equal_to(elem_id_node_num.size(), data_size);

      std::vector<dof_id_type> requested_objs_id(data.size());

      for (auto i : IntRange<std::size_t>(0,data_size))
        {
          const Elem & elem = mesh.elem_ref(elem_id_node_num[i].first);
          const Node & node = elem.node_ref(elem_id_node_num[i].second);
          requested_objs_id[i] = node.id();
        }

      // Let the user process the results.  If any of the results
      // were different than what the user expected, then we may
      // need to sync again just in case this processor has to
      // pass on the changes to yet another processor.
      if (sync.act_on_data(requested_objs_id, data))
        data_changed = true;
    };

  // Trade requests with other processors
  typename SyncFunctor::datum * ex = nullptr;
  pull_parallel_vector_data
    (comm, requested_objs_elem_id_node_num, gather_functor,
     action_functor, ex);

  comm.max(data_changed);

  return data_changed;
}



template <typename ElemCheckFunctor,
          typename NodeCheckFunctor,
          typename SyncFunctor>
void sync_node_data_by_element_id(MeshBase & mesh,
                                  const MeshBase::const_element_iterator & range_begin,
                                  const MeshBase::const_element_iterator & range_end,
                                  const ElemCheckFunctor & elem_check,
                                  const NodeCheckFunctor & node_check,
                                  SyncFunctor & sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(mesh.comm());

  bool need_sync = false;

  do
    {
      need_sync =
        sync_node_data_by_element_id_once
          (mesh, range_begin, range_end, elem_check, node_check,
           sync);
    } while (need_sync);
}


}



// This struct can be created and passed to the
// Parallel::sync_dofobject_data_by_id() function.
struct SyncNodalPositions
{
  // The constructor.  You need a reference to the mesh where you will
  // be setting/getting nodal positions.
  explicit
  SyncNodalPositions(MeshBase & m);

  // The datum typedef is required of this functor, so that the
  // Parallel::sync_dofobject_data_by_id() function can create e.g.
  // std::vector<datum>.
  typedef Point datum;

  // First required interface.  This function must fill up the data vector for the
  // ids specified in the ids vector.
  void gather_data (const std::vector<dof_id_type> & ids, std::vector<datum> & data) const;

  // Second required interface.  This function must do something with the data in
  // the data vector for the ids in the ids vector.
  void act_on_data (const std::vector<dof_id_type> & ids, const std::vector<datum> & data) const;

  MeshBase & mesh;
};

// This struct can be created and passed to the
// Parallel::sync_dofobject_data_by_id() function.
struct SyncSubdomainIds
{
  // The constructor.  You need a reference to the mesh where you will
  // be setting/getting element subdomain IDs.
  explicit
  SyncSubdomainIds(MeshBase & m);

  // The datum typedef is required of this functor, so that the
  // Parallel::sync_dofobject_data_by_id() function can create e.g.
  // std::vector<datum>.
  typedef subdomain_id_type datum;

  // First required interface.  This function must fill up the data vector for the
  // ids specified in the ids vector.
  void gather_data (const std::vector<dof_id_type> & ids, std::vector<datum> & data) const;

  // Second required interface.  This function must do something with the data in
  // the data vector for the ids in the ids vector.
  void act_on_data (const std::vector<dof_id_type> & ids, const std::vector<datum> & data) const;

  MeshBase & mesh;
};


} // namespace libMesh

#endif // LIBMESH_PARALLEL_GHOST_SYNC_H
