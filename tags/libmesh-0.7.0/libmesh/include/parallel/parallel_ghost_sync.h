
// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __parallel_ghost_sync_h__
#define __parallel_ghost_sync_h__
// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "auto_ptr.h"
#include "elem.h"
#include "location_maps.h"
#include "mesh_base.h"
#include "parallel.h"

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
   * The user must define Parallel::datatype<sync::datum> if
   * sync::datum isn't a built-in type.
   * The user-provided location_map will be used and left unchanged
   * if it is provided, or filled and cleared if it is empty.
   */
  template <typename Iterator,
            typename DofObjType,
            typename SyncFunctor>
  void sync_dofobject_data_by_xyz(const Iterator&          range_begin,
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
   * The user must define Parallel::datatype<sync::datum> if
   * sync::datum isn't a built-in type.
   */
  template <typename Iterator,
            typename SyncFunctor>
  void sync_dofobject_data_by_id(const Iterator& range_begin,
                                 const Iterator& range_end,
                                 SyncFunctor&    sync);

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
   * The user must define Parallel::datatype<sync::datum> if
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
void sync_dofobject_data_by_xyz(const Iterator&          range_begin,
                                const Iterator&          range_end,
                                LocationMap<DofObjType>& location_map,
                                SyncFunctor&             sync)
{
  // This function must be run on all processors at once
  parallel_only();

  // We need a valid location_map
#ifdef DEBUG
  bool need_map_update = (range_begin != range_end && location_map.empty());
  Parallel::max(need_map_update);
  libmesh_assert(!need_map_update);
#endif

  // Count the objectss to ask each processor about
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType *obj = *it;
      libmesh_assert (obj);
      unsigned int obj_procid = obj->processor_id();
      libmesh_assert (obj_procid != DofObject::invalid_processor_id);

      ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<Real> >
    requested_objs_x(libMesh::n_processors()),
    requested_objs_y(libMesh::n_processors()),
    requested_objs_z(libMesh::n_processors());
  // Corresponding ids to keep track of
  std::vector<std::vector<unsigned int> >
    requested_objs_id(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        requested_objs_x[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_y[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_z[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObjType *obj = *it;
      unsigned int obj_procid = obj->processor_id();
      if (obj_procid == libMesh::processor_id())
        continue;

      Point p = location_map.point_of(*obj);
      requested_objs_x[obj_procid].push_back(p(0));
      requested_objs_y[obj_procid].push_back(p(1));
      requested_objs_z[obj_procid].push_back(p(2));
      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<Real> request_to_fill_x,
                        request_to_fill_y,
                        request_to_fill_z;
      Parallel::send_receive(procup, requested_objs_x[procup],
                             procdown, request_to_fill_x);
      Parallel::send_receive(procup, requested_objs_y[procup],
                             procdown, request_to_fill_y);
      Parallel::send_receive(procup, requested_objs_z[procup],
                             procdown, request_to_fill_z);

      // Find the local id of each requested object
      std::vector<unsigned int> request_to_fill_id(request_to_fill_x.size());
      for (unsigned int i=0; i != request_to_fill_x.size(); ++i)
        {
          Point p(request_to_fill_x[i],
                  request_to_fill_y[i],
                  request_to_fill_z[i]);

          // Look for this object in the multimap
          DofObjType *obj = location_map.find(p);

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
      Parallel::send_receive(procdown, data,
                             procup, received_data);
      libmesh_assert(requested_objs_x[procup].size() ==
                     received_data.size());

      // Let the user process the results
      sync.act_on_data(requested_objs_id[procup], received_data);
    }
}



template <typename Iterator,
          typename SyncFunctor>
void sync_dofobject_data_by_id(const Iterator& range_begin,
                               const Iterator& range_end,
                               SyncFunctor&    sync)
{
  // This function must be run on all processors at once
  parallel_only();

  // Count the objects to ask each processor about
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      libmesh_assert (obj);
      unsigned int obj_procid = obj->processor_id();
      libmesh_assert (obj_procid != DofObject::invalid_processor_id);

      ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_objs_id(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
      }
  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      unsigned int obj_procid = obj->processor_id();
      if (obj_procid == libMesh::processor_id())
        continue;

      requested_objs_id[obj_procid].push_back(obj->id());
    }

  // Trade requests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill_id;
      Parallel::send_receive(procup, requested_objs_id[procup],
                             procdown, request_to_fill_id);

      // Gather whatever data the user wants
      std::vector<typename SyncFunctor::datum> data;
      sync.gather_data(request_to_fill_id, data);

      // Trade back the results
      std::vector<typename SyncFunctor::datum> received_data;
      Parallel::send_receive(procdown, data,
                             procup, received_data);
      libmesh_assert(requested_objs_id[procup].size() ==
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
  // This function must be run on all processors at once
  parallel_only();

  // Count the objects to ask each processor about
  std::vector<unsigned int>
    ghost_objects_from_proc(libMesh::n_processors(), 0);

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      DofObject *obj = *it;
      libmesh_assert (obj);
      unsigned int obj_procid = obj->processor_id();
      libmesh_assert (obj_procid != DofObject::invalid_processor_id);

      ghost_objects_from_proc[obj_procid]++;
    }

  // Request sets to send to each processor
  std::vector<std::vector<unsigned int> >
    requested_objs_id(libMesh::n_processors()),
    requested_objs_parent_id(libMesh::n_processors()),
    requested_objs_child_num(libMesh::n_processors());

  // We know how many objects live on each processor, so reserve()
  // space for each.
  for (unsigned int p=0; p != libMesh::n_processors(); ++p)
    if (p != libMesh::processor_id())
      {
        requested_objs_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_parent_id[p].reserve(ghost_objects_from_proc[p]);
        requested_objs_child_num[p].reserve(ghost_objects_from_proc[p]);
      }

  for (Iterator it = range_begin; it != range_end; ++it)
    {
      Elem *elem = *it;
      unsigned int obj_procid = elem->processor_id();
      if (obj_procid == libMesh::processor_id())
        continue;
      const Elem *parent = elem->parent();
      if (!parent || !elem->active())
        continue;

      requested_objs_id[obj_procid].push_back(elem->id());
      requested_objs_parent_id[obj_procid].push_back(parent->id());
      requested_objs_child_num[obj_procid].push_back
	(parent->which_child_am_i(elem));
    }

  // Trade requests with other processors
  for (unsigned int p=1; p != libMesh::n_processors(); ++p)
    {
      // Trade my requests with processor procup and procdown
      unsigned int procup = (libMesh::processor_id() + p) %
                             libMesh::n_processors();
      unsigned int procdown = (libMesh::n_processors() +
                               libMesh::processor_id() - p) %
                               libMesh::n_processors();
      std::vector<unsigned int> request_to_fill_parent_id,
                                request_to_fill_child_num;
      Parallel::send_receive(procup, requested_objs_parent_id[procup],
                             procdown, request_to_fill_parent_id);
      Parallel::send_receive(procup, requested_objs_child_num[procup],
                             procdown, request_to_fill_child_num);

      // Find the id of each requested element
      unsigned int request_size = request_to_fill_parent_id.size();
      std::vector<unsigned int> request_to_fill_id(request_size);
      for (unsigned int i=0; i != request_size; ++i)
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
      Parallel::send_receive(procdown, data,
                             procup, received_data);
      libmesh_assert(requested_objs_id[procup].size() ==
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

} // namespace libMesh

#endif // #define __parallel_ghost_sync_h__
