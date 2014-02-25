// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes -----------------------------------
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/location_maps.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"

// C++ Includes   -----------------------------------

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
 * sync.gather_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 * The user-provided location_map will be used and left unchanged
 * if it is provided, or filled and cleared if it is empty.
 */
template <typename Iterator,
          typename DofObjType,
          typename SyncFunctor>
void sync_dofobject_data_by_xyz(const Communicator&      communicator,
                                const Iterator&          range_begin,
                                const Iterator&          range_end,
                                LocationMap<DofObjType>* location_map,
                                SyncFunctor&             sync);

//------------------------------------------------------------------------
/**
 * Request data about a range of ghost dofobjects uniquely
 * identified by their id.  Fulfill requests with
 * sync.gather_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 */
template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator& communicator,
                               const Iterator&     range_begin,
                               const Iterator&     range_end,
                               SyncFunctor&        sync);

//------------------------------------------------------------------------
/**
 * Request data about a range of ghost elements uniquely
 * identified by their parent id and which child they are.
 * Fulfill requests with
 * sync.gather_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data),
 * by resizing and setting the values of the data vector.
 * Respond to fulfillment with
 * sync.act_on_data(const std::vector<unsigned int>& ids,
 *                  std::vector<sync::datum>& data)
 * The user must define Parallel::StandardType<sync::datum> if
 * sync::datum isn't a built-in type.
 */
template <typename Iterator,
          typename SyncFunctor>
void sync_element_data_by_parent_id(MeshBase&       mesh,
                                    const Iterator& range_begin,
                                    const Iterator& range_end,
                                    SyncFunctor&    sync);

//------------------------------------------------------------------------
// Parallel members

template <typename Iterator,
          typename DofObjType,
          typename SyncFunctor>
void sync_dofobject_data_by_xyz(const Communicator&      communicator,
                                const Iterator&          range_begin,
                                const Iterator&          range_end,
                                LocationMap<DofObjType>& location_map,
                                SyncFunctor&             sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(communicator);

  // We need a valid location_map
#ifdef DEBUG
  bool need_map_update = (range_begin != range_end && location_map.empty());
  communicator.max(need_map_update);
  libmesh_assert(!need_map_update);
#endif

  // Count the objectss to ask each processor about
  std::vector<dof_id_type>
    ghost_objects_from_proc(communicator.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType *obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<Real> >
    requested_objs_x(communicator.size()),
    requested_objs_y(communicator.size()),
    requested_objs_z(communicator.size());
  // Corresponding ids to keep track of
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(communicator.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != communicator.size(); ++p)
    if (p != communicator.rank())
      {
        requested_objs_x[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_y[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_z[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType *obj = *it;
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == communicator.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      Point p = location_map.point_of(*obj);
      requested_objs_x[obj_procid].push_back(p(0));
      requested_objs_y[obj_procid].push_back(p(1));
      requested_objs_z[obj_procid].push_back(p(2));
      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != communicator.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        libmesh_cast_int<processor_id_type>
        ((communicator.rank() + p) % communicator.size());
      const processor_id_type procdown =
        libmesh_cast_int<processor_id_type>
        ((communicator.size() + communicator.rank() - p) %
         communicator.size());
      std::vector<Real> request_to_fill_x,
        request_to_fill_y,
        request_to_fill_z;
      communicator.send_receive(procup, requested_objs_x[procup],
                                procdown, request_to_fill_x);
      communicator.send_receive(procup, requested_objs_y[procup],
                                procdown, request_to_fill_y);
      communicator.send_receive(procup, requested_objs_z[procup],
                                procdown, request_to_fill_z);

      // Find the local id of each requested object
      std::vector<dof_id_type> request_to_fill_id(request_to_fill_x.size());
      for (std::size_t i=0; i != request_to_fill_x.size(); ++i)
        {
          Point pt(request_to_fill_x[i],
                   request_to_fill_y[i],
                   request_to_fill_z[i]);

          // Look for this object in the multimap
          DofObjType *obj = location_map.find(pt);

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
      communicator.send_receive(procdown, data,
                                procup, received_data);
      libmesh_assert_equal_to (requested_objs_x[procup].size(),
                               received_data.size());

      // Let the user process the results
      sync.act_on_data(requested_objs_id[procup], received_data);
    }
}



template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Communicator& communicator,
                               const Iterator& range_begin,
                               const Iterator& range_end,
                               SyncFunctor&    sync)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(communicator);

  // Count the objects to ask each processor about
  std::vector<dof_id_type>
    ghost_objects_from_proc(communicator.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(communicator.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != communicator.size(); ++p)
    if (p != communicator.rank())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid == communicator.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;

      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != communicator.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        libmesh_cast_int<processor_id_type>
        (communicator.rank() + p) % communicator.size();
      const processor_id_type procdown =
        libmesh_cast_int<processor_id_type>
        ((communicator.size() + communicator.rank() - p) %
         communicator.size());
      std::vector<dof_id_type> request_to_fill_id;
      communicator.send_receive(procup, requested_objs_id[procup],
                                procdown, request_to_fill_id);

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      communicator.send_receive(procdown, data,
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
void sync_element_data_by_parent_id(MeshBase&       mesh,
                                    const Iterator& range_begin,
                                    const Iterator& range_end,
                                    SyncFunctor&    sync)
{
  const Communicator &communicator (mesh.comm());

  // This function must be run on all processors at once
  libmesh_parallel_only(communicator);

  // Count the objects to ask each processor about
  std::vector<dof_id_type>
    ghost_objects_from_proc(communicator.size(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      libmesh_assert (obj);
      processor_id_type obj_procid = obj->processor_id();
      if (obj_procid != DofObject::invalid_processor_id)
        ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<dof_id_type> >
    requested_objs_id(communicator.size()),
    requested_objs_parent_id(communicator.size());
  std::vector<std::vector<unsigned char> >
    requested_objs_child_num(communicator.size());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (processor_id_type p=0; p != communicator.size(); ++p)
    if (p != communicator.rank())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_parent_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_child_num[p].reserve(ghost_objects_from_proc[p]);
      }

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      Elem *elem = *it;
      processor_id_type obj_procid = elem->processor_id();
      if (obj_procid == communicator.rank() ||
          obj_procid == DofObject::invalid_processor_id)
        continue;
      const Elem *parent = elem->parent();
      if (!parent || !elem->active())
        continue;

      requested_objs_id[obj_procid].push_back(elem->id());
      requested_objs_parent_id[obj_procid].push_back(parent->id());
      requested_objs_child_num[obj_procid].push_back
        (libmesh_cast_int<unsigned char>
         (parent->which_child_am_i(elem)));
    }

  // Trade requests with other processors
  for (processor_id_type p=1; p != communicator.size(); ++p)
    {
      // Trade my requests with processor procup and procdown
      const processor_id_type procup =
        libmesh_cast_int<processor_id_type>
        (communicator.rank() + p) % communicator.size();
      const processor_id_type procdown =
        libmesh_cast_int<processor_id_type>
        ((communicator.size() + communicator.rank() - p) %
         communicator.size());
      std::vector<dof_id_type>   request_to_fill_parent_id;
      std::vector<unsigned char> request_to_fill_child_num;
      communicator.send_receive(procup, requested_objs_parent_id[procup],
                                procdown, request_to_fill_parent_id);
      communicator.send_receive(procup, requested_objs_child_num[procup],
                                procdown, request_to_fill_child_num);

      // Find the id of each requested element
      std::size_t request_size = request_to_fill_parent_id.size();
      std::vector<dof_id_type> request_to_fill_id(request_size);
      for (std::size_t i=0; i != request_size; ++i)
        {
          Elem *parent = mesh.elem(request_to_fill_parent_id[i]);
          libmesh_assert(parent);
          libmesh_assert(parent->has_children());
          Elem *child = parent->child(request_to_fill_child_num[i]);
          libmesh_assert(child);
          libmesh_assert(child->active());
          request_to_fill_id[i] = child->id();
        }

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      communicator.send_receive(procdown, data,
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
void sync_element_data_by_parent_id(MeshBase&,
                                    const Iterator&,
                                    const Iterator&,
                                    SyncFunctor&)
{
}
#endif // LIBMESH_ENABLE_AMR




}


// This struct can be created and passed to the
// Parallel::sync_dofobject_data_by_id() function.
struct SyncNodalPositions
{
  // The constructor.  You need a reference to the mesh where you will
  // be setting/getting nodal positions.
  explicit
  SyncNodalPositions(MeshBase& m);

  // The datum typedef is required of this functor, so that the
  // Parallel::sync_dofobject_data_by_id() function can create e.g.
  // std::vector<datum>.
  typedef Point datum;

  // First required interface.  This function must fill up the data vector for the
  // ids specified in the ids vector.
  void gather_data (const std::vector<dof_id_type>& ids, std::vector<datum>& data);

  // Second required interface.  This function must do something with the data in
  // the data vector for the ids in the ids vector.
  void act_on_data (const std::vector<dof_id_type>& ids, std::vector<datum>& data);

  MeshBase &mesh;
};


} // namespace libMesh

#endif // LIBMESH_PARALLEL_GHOST_SYNC_H
