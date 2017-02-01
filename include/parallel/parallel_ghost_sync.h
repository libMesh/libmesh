// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/location_maps.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"

// C++ Includes
#include LIBMESH_INCLUDE_UNORDERED_SET


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
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator &     range_begin,
                               const Iterator &     range_end,
                               SyncFunctor &        sync);

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
 * Request data about a range of ghost nodes uniquely identified by
 * an element id and local node id.
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

  // Count the objectss to ask each processor about
  std::vector<dof_id_type>
    ghost_objects_from_proc(comm.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType * obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<Real> >
    requested_objs_x(comm.size()),
    requested_objs_y(comm.size()),
    requested_objs_z(comm.size());
  // Corresponding ids to keep track of
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(comm.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != comm.size(); ++p)
    if (p != comm.rank())
      {
        requested_objs_x[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_y[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_z[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType * obj = *it;
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      Point p = location_map.point_of(*obj);
      requested_objs_x[obj_procid].push_back(p(0));
      requested_objs_y[obj_procid].push_back(p(1));
      requested_objs_z[obj_procid].push_back(p(2));
      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != comm.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        cast_int<processor_id_type>
        ((comm.rank() + p) % comm.size());
      const processor_id_type procdown =
        cast_int<processor_id_type>
        ((comm.size() + comm.rank() - p) %
         comm.size());
      std::vector<Real> request_to_fill_x,
        request_to_fill_y,
        request_to_fill_z;
      comm.send_receive(procup, requested_objs_x[procup],
                        procdown, request_to_fill_x);
      comm.send_receive(procup, requested_objs_y[procup],
                        procdown, request_to_fill_y);
      comm.send_receive(procup, requested_objs_z[procup],
                        procdown, request_to_fill_z);

      // Find the local id of each requested object
      std::vector<dof_id_type> request_to_fill_id(request_to_fill_x.size());
      for (std::size_t i=0; i != request_to_fill_x.size(); ++i)
        {
          Point pt(request_to_fill_x[i],
                   request_to_fill_y[i],
                   request_to_fill_z[i]);

          // Look for this object in the multimap
          DofObjType * obj = location_map.find(pt);

          // We'd better find every object we're asked for
          libmesh_assert (obj);

          // Return the object's correct processor id,
          // and our (correct if it's local) id for it.
          request_to_fill_id[i] = obj->id();
        }

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      comm.send_receive(procdown, data,
                        procup, received_data);
      libmesh_assert_equal_to (requested_objs_x[procup].size(),
                               received_data.size());

      // Let the user process the results
      sync.act_on_data(requested_objs_id[procup], received_data);
    }
}



template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator & comm,
                               const Iterator & range_begin,
                               const Iterator & range_end,
                               SyncFunctor &    sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  // Count the objects to ask each processor about
  std::vector<dof_id_type>
    ghost_objects_from_proc(comm.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject * obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(comm.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != comm.size(); ++p)
    if (p != comm.rank())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject * obj = *it;
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == comm.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != comm.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        cast_int<processor_id_type>
        ((comm.rank() + p) % comm.size());
      const processor_id_type procdown =
        cast_int<processor_id_type>
        ((comm.size() + comm.rank() - p) %
         comm.size());
      std::vector<dof_id_type> request_to_fill_id;
      comm.send_receive(procup, requested_objs_id[procup],
                        procdown, request_to_fill_id);

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      comm.send_receive(procdown, data,
                        procup, received_data);
      libmesh_assert_equal_to (requested_objs_id[procup].size(),
                               received_data.size());

      // Let the user process the results
      sync.act_on_data(requested_objs_id[procup], received_data);
    }
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
  std::vector<dof_id_type>
    ghost_objects_from_proc(comm.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject * obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(comm.size()),
    requested_objs_parent_id(comm.size());
  std::vector<std::vector<unsigned char> >
    requested_objs_child_num(comm.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != comm.size(); ++p)
    if (p != comm.rank())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_parent_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_child_num[p].reserve(ghost_objects_from_proc[p]);
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
      requested_objs_parent_id[obj_procid].push_back(parent->id());
      requested_objs_child_num[obj_procid].push_back
        (cast_int<unsigned char>
         (parent->which_child_am_i(elem)));
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != comm.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        cast_int<processor_id_type>
        ((comm.rank() + p) % comm.size());
      const processor_id_type procdown =
        cast_int<processor_id_type>
        ((comm.size() + comm.rank() - p) %
         comm.size());
      std::vector<dof_id_type>   request_to_fill_parent_id;
      std::vector<unsigned char> request_to_fill_child_num;
      comm.send_receive(procup, requested_objs_parent_id[procup],
                        procdown, request_to_fill_parent_id);
      comm.send_receive(procup, requested_objs_child_num[procup],
                        procdown, request_to_fill_child_num);

      // Find the id of each requested element
      std::size_t request_size = request_to_fill_parent_id.size();
      std::vector<dof_id_type> request_to_fill_id(request_size);
      for (std::size_t i=0; i != request_size; ++i)
        {
          Elem & parent = mesh.elem_ref(request_to_fill_parent_id[i]);
          libmesh_assert(parent.has_children());
          Elem * child = parent.child_ptr(request_to_fill_child_num[i]);
          libmesh_assert(child);
          libmesh_assert(child->active());
          request_to_fill_id[i] = child->id();
        }

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      comm.send_receive(procdown, data,
                        procup, received_data);
      libmesh_assert_equal_to (requested_objs_id[procup].size(),
                               received_data.size());

      // Let the user process the results
      sync.act_on_data(requested_objs_id[procup], received_data);
    }
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


struct SyncEverything
{
  SyncEverything() {}

  bool operator() (const Elem *) const { return true; }

  bool operator() (const Elem *, unsigned int) const
  { return true; }
};



template <typename ElemCheckFunctor,
          typename NodeCheckFunctor,
          typename SyncFunctor>
void sync_node_data_by_element_id(MeshBase &       mesh,
                                  const MeshBase::const_element_iterator & range_begin,
                                  const MeshBase::const_element_iterator & range_end,
                                  const ElemCheckFunctor & elem_check,
                                  const NodeCheckFunctor & node_check,
                                  SyncFunctor & sync)
{
  const Communicator & comm (mesh.comm());

  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  // Keep track of which nodes we've asked about, so we only hit each
  // once?
  // LIBMESH_BEST_UNORDERED_SET<dof_id_type> queried_nodes;

  // No.  We need to ask every neighboring processor about every node,
  // probably repeatedly.  Imagine a vertex surrounded by triangles,
  // each on a different processor, with a ghosting policy that
  // include only face neighbors and not point neighbors.  Then the
  // only way for authoritative information to trickle out from that
  // vertex is by being passed along, one neighbor at a time, to
  // processors who mostly don't even see the node's true owner!

  bool need_sync = false;

  do
    {
      // This is the last sync we need, unless we later discover
      // otherwise
      need_sync = false;

      // Count the objects to ask each processor about
      std::vector<dof_id_type>
        ghost_objects_from_proc(comm.size(), 0);

      for (MeshBase::const_element_iterator it = range_begin;
           it != range_end; ++it)
        {
          const Elem * elem = *it;
          libmesh_assert (elem);

          if (!elem_check(elem))
            continue;

          const processor_id_type proc_id = elem->processor_id();
          if (proc_id == comm.rank() ||
              proc_id == DofObject::invalid_processor_id)
            continue;

          for (unsigned int n=0; n != elem->n_nodes(); ++n)
            {
              if (!node_check(elem, n))
                continue;

              ghost_objects_from_proc[proc_id]++;
            }
        }

      // Now repeat that iteration, filling request sets this time.

      // Request sets to send to each processor
      std::vector<std::vector<dof_id_type> >
        requested_objs_elem_id(comm.size());
      std::vector<std::vector<unsigned char> >
        requested_objs_node_num(comm.size());

      // Keep track of current local ids for each too
      std::vector<std::vector<dof_id_type> >
        requested_objs_id(comm.size());

      // We know how many objects live on each processor, so reserve()
      // space for each.
      for (processor_id_type p=0; p != comm.size(); ++p)
        if (p != comm.rank())
          {
            requested_objs_elem_id[p].reserve(ghost_objects_from_proc[p]);
            requested_objs_node_num[p].reserve(ghost_objects_from_proc[p]);
            requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
          }

      for (MeshBase::const_element_iterator it = range_begin;
           it != range_end; ++it)
        {
          const Elem * elem = *it;
          libmesh_assert (elem);

          if (!elem_check(elem))
            continue;

          const processor_id_type proc_id = elem->processor_id();
          if (proc_id == comm.rank() ||
              proc_id == DofObject::invalid_processor_id)
            continue;

          const dof_id_type elem_id = elem->id();

          for (unsigned int n=0; n != elem->n_nodes(); ++n)
            {
              if (!node_check(elem, n))
                continue;

              const Node & node = elem->node_ref(n);
              const dof_id_type node_id = node.id();

              requested_objs_elem_id[proc_id].push_back(elem_id);
              requested_objs_node_num[proc_id].push_back
                (cast_int<unsigned char>(n));
              requested_objs_id[proc_id].push_back(node_id);
            }
        }

      // Trade requests with other processors
      for (processor_id_type p=1; p != comm.size(); ++p)
        {
          // Trade my requests with processor procup and procdown
          const processor_id_type procup =
            cast_int<processor_id_type>
            ((comm.rank() + p) % comm.size());
          const processor_id_type procdown =
            cast_int<processor_id_type>
            ((comm.size() + comm.rank() - p) %
             comm.size());

          libmesh_assert_equal_to (requested_objs_id[procup].size(),
                                   ghost_objects_from_proc[procup]);
          libmesh_assert_equal_to (requested_objs_elem_id[procup].size(),
                                   ghost_objects_from_proc[procup]);
          libmesh_assert_equal_to (requested_objs_node_num[procup].size(),
                                   ghost_objects_from_proc[procup]);

          std::vector<dof_id_type>   request_to_fill_elem_id;
          std::vector<unsigned char> request_to_fill_node_num;
          comm.send_receive(procup, requested_objs_elem_id[procup],
                            procdown, request_to_fill_elem_id);
          comm.send_receive(procup, requested_objs_node_num[procup],
                            procdown, request_to_fill_node_num);

          // Find the id of each requested element
          std::size_t request_size = request_to_fill_elem_id.size();
          std::vector<dof_id_type> request_to_fill_id(request_size);
          for (std::size_t i=0; i != request_size; ++i)
            {
              const Elem & elem = mesh.elem_ref(request_to_fill_elem_id[i]);

              const unsigned int n = request_to_fill_node_num[i];
              libmesh_assert_less (n, elem.n_nodes());

              const Node & node = elem.node_ref(n);

              // This isn't a safe assertion in the case where we're
              // synching processor ids
              // libmesh_assert_equal_to (node->processor_id(), comm.rank());

              request_to_fill_id[i] = node.id();
            }

          // Gather whatever data the user wants
          std::vector<typename SyncFunctor::datum> data;
          sync.gather_data(request_to_fill_id, data);

          // Trade back the results
          std::vector<typename SyncFunctor::datum> received_data;
          comm.send_receive(procdown, data,
                            procup, received_data);
          libmesh_assert_equal_to (requested_objs_elem_id[procup].size(),
                                   received_data.size());

          // Let the user process the results.  If any of the results
          // were different than what the user expected, then we'll
          // need to sync again just in case this processor has to
          // pass on the changes to yet another processor.
          bool data_changed =
            sync.act_on_data(requested_objs_id[procup], received_data);

          if (data_changed)
            need_sync = true;
        }
      comm.max(need_sync);
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
  void act_on_data (const std::vector<dof_id_type> & ids, std::vector<datum> & data) const;

  MeshBase & mesh;
};


} // namespace libMesh

#endif // LIBMESH_PARALLEL_GHOST_SYNC_H
