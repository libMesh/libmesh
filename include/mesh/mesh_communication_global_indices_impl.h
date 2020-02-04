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

#ifndef LIBMESH_MESH_COMMUNICATION_GLOBAL_INDICES_IMPL_H
#define LIBMESH_MESH_COMMUNICATION_GLOBAL_INDICES_IMPL_H

// Local includes
#include "libmesh/mesh_communication.h"

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_hilbert.h"
#include "libmesh/parallel_sort.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/node_range.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"

// C/C++ includes
#ifdef LIBMESH_HAVE_LIBHILBERT
#  include "hilbert.h"
#endif


#ifdef LIBMESH_HAVE_LIBHILBERT
namespace { // anonymous namespace for helper functions

using namespace libMesh;

// Utility function to map (x,y,z) in [bbox.min, bbox.max]^3 into
// [0,max_inttype]^3 for computing Hilbert keys
template <typename RealType>
void get_hilbert_coords (const PointTempl<RealType> & p,
                         const libMesh::BoundingBoxTempl<RealType> & bbox,
                         CFixBitVec icoords[3])
{
  static const Hilbert::inttype max_inttype = static_cast<Hilbert::inttype>(-1);

  const long double // put (x,y,z) in [0,1]^3 (don't divide by 0)
    x = static_cast<long double>
        ((bbox.first(0) == bbox.second(0)) ? 0. :
         (p(0)-bbox.first(0))/(bbox.second(0)-bbox.first(0))),

#if LIBMESH_DIM > 1
    y = static_cast<long double>
        ((bbox.first(1) == bbox.second(1)) ? 0. :
         (p(1)-bbox.first(1))/(bbox.second(1)-bbox.first(1))),
#else
    y = 0.,
#endif

#if LIBMESH_DIM > 2
    z = static_cast<long double>
        ((bbox.first(2) == bbox.second(2)) ? 0. :
         (p(2)-bbox.first(2))/(bbox.second(2)-bbox.first(2)));
#else
  z = 0.;
#endif

  // (icoords) in [0,max_inttype]^3
  icoords[0] = static_cast<Hilbert::inttype>(x*max_inttype);
  icoords[1] = static_cast<Hilbert::inttype>(y*max_inttype);
  icoords[2] = static_cast<Hilbert::inttype>(z*max_inttype);
}



template <typename RealType>
Parallel::DofObjectKey
get_dofobject_key (const ElemTempl<RealType> * e,
                   const libMesh::BoundingBoxTempl<RealType> & bbox)
{
  static const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

  Hilbert::HilbertIndices index;
  CFixBitVec icoords[3];
  Hilbert::BitVecType bv;
  get_hilbert_coords (e->centroid(), bbox, icoords);
  Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, bv);
  index = bv;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  return std::make_pair(index, e->unique_id());
#else
  return index;
#endif
}



// Compute the hilbert index
template <typename RealType>
Parallel::DofObjectKey
get_dofobject_key (const NodeTempl<RealType> * n,
                   const libMesh::BoundingBoxTempl<RealType> & bbox)
{
  static const unsigned int sizeof_inttype = sizeof(Hilbert::inttype);

  Hilbert::HilbertIndices index;
  CFixBitVec icoords[3];
  Hilbert::BitVecType bv;
  get_hilbert_coords (*n, bbox, icoords);
  Hilbert::coordsToIndex (icoords, 8*sizeof_inttype, 3, bv);
  index = bv;

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  return std::make_pair(index, n->unique_id());
#else
  return index;
#endif
}

// Helper class for threaded Hilbert key computation
template <typename RealType>
class ComputeHilbertKeys
{
public:
  ComputeHilbertKeys (const libMesh::BoundingBoxTempl<RealType> & bbox,
                      std::vector<Parallel::DofObjectKey> & keys) :
    _bbox(bbox),
    _keys(keys)
  {}

  // computes the hilbert index for a node
  void operator() (const ConstNodeRangeTempl<RealType> & range) const
  {
    std::size_t pos = range.first_idx();
    for (const auto & node : range)
      {
        libmesh_assert(node);
        libmesh_assert_less (pos, _keys.size());
        _keys[pos++] = get_dofobject_key (node, _bbox);
      }
  }

  // computes the hilbert index for an element
  void operator() (const ConstElemRangeTempl<RealType> & range) const
  {
    std::size_t pos = range.first_idx();
    for (const auto & elem : range)
      {
        libmesh_assert(elem);
        libmesh_assert_less (pos, _keys.size());
        _keys[pos++] = get_dofobject_key (elem, _bbox);
      }
  }

private:
  const libMesh::BoundingBoxTempl<RealType> & _bbox;
  std::vector<Parallel::DofObjectKey> & _keys;
};
}
#endif


namespace libMesh
{

// ------------------------------------------------------------
// MeshCommunication class members
#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
template <typename RealType>
void MeshCommunication::assign_global_indices (MeshBaseTempl<RealType> & mesh) const
{
  LOG_SCOPE ("assign_global_indices()", "MeshCommunication");

  // This method determines partition-agnostic global indices
  // for nodes and elements.

  // Algorithm:
  // (1) compute the Hilbert key for each local node/element
  // (2) perform a parallel sort of the Hilbert key
  // (3) get the min/max value on each processor
  // (4) determine the position in the global ranking for
  //     each local object

  const Parallel::Communicator & communicator (mesh.comm());

  // Global bounding box.  We choose the nodal bounding box for
  // backwards compatibility; the element bounding box may be looser
  // on curved elements.
  BoundingBox bbox =
    MeshTools::create_nodal_bounding_box (mesh);

  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  std::vector<Parallel::DofObjectKey>
    node_keys, elem_keys;

  {
    // Nodes first
    {
      ConstNodeRangeTempl<RealType> nr (mesh.local_nodes_begin(),
                         mesh.local_nodes_end());
      node_keys.resize (nr.size());
      Threads::parallel_for (nr, ComputeHilbertKeys<RealType> (bbox, node_keys));

      // // It's O(N^2) to check that these keys don't duplicate before the
      // // sort...
      // typename MeshBaseTempl<RealType>::const_node_iterator nodei = mesh.local_nodes_begin();
      // for (std::size_t i = 0; i != node_keys.size(); ++i, ++nodei)
      //   {
      //     typename MeshBaseTempl<RealType>::const_node_iterator nodej = mesh.local_nodes_begin();
      //     for (std::size_t j = 0; j != i; ++j, ++nodej)
      //       {
      //         if (node_keys[i] == node_keys[j])
      //           {
      //             CFixBitVec icoords[3], jcoords[3];
      //             get_hilbert_coords(**nodej, bbox, jcoords);
      //             libMesh::err << "node " << (*nodej)->id() << ", " << static_cast<PointTempl<RealType> &>(**nodej) << " has HilbertIndices " << node_keys[j] << std::endl;
      //             get_hilbert_coords(**nodei, bbox, icoords);
      //             libMesh::err << "node " << (*nodei)->id() << ", " << static_cast<PointTempl<RealType> &>(**nodei) << " has HilbertIndices " << node_keys[i] << std::endl;
      //             libmesh_error_msg("Error: nodes with duplicate Hilbert keys!");
      //           }
      //       }
      //   }
    }

    // Elements next
    {
      ConstElemRangeTempl<RealType> er (mesh.local_elements_begin(),
                         mesh.local_elements_end());
      elem_keys.resize (er.size());
      Threads::parallel_for (er, ComputeHilbertKeys<RealType> (bbox, elem_keys));

      // // For elements, the keys can be (and in the case of TRI, are
      // // expected to be) duplicates, but only if the elements are at
      // // different levels
      // typename MeshBaseTempl<RealType>::const_element_iterator elemi = mesh.local_elements_begin();
      // for (std::size_t i = 0; i != elem_keys.size(); ++i, ++elemi)
      //   {
      //     typename MeshBaseTempl<RealType>::const_element_iterator elemj = mesh.local_elements_begin();
      //     for (std::size_t j = 0; j != i; ++j, ++elemj)
      //       {
      //         if ((elem_keys[i] == elem_keys[j]) &&
      //             ((*elemi)->level() == (*elemj)->level()))
      //           {
      //             libMesh::err << "level " << (*elemj)->level()
      //                          << " elem\n" << (**elemj)
      //                          << " centroid " << (*elemj)->centroid()
      //                          << " has HilbertIndices " << elem_keys[j]
      //                          << " or " << get_dofobject_key((*elemj), bbox)
      //                          << std::endl;
      //             libMesh::err << "level " << (*elemi)->level()
      //                          << " elem\n" << (**elemi)
      //                          << " centroid " << (*elemi)->centroid()
      //                          << " has HilbertIndices " << elem_keys[i]
      //                          << " or " << get_dofobject_key((*elemi), bbox)
      //                          << std::endl;
      //             libmesh_error_msg("Error: level " << (*elemi)->level() << " elements with duplicate Hilbert keys!");
      //           }
      //       }
      //  }
    }
  } // done computing Hilbert keys



  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  Parallel::Sort<Parallel::DofObjectKey> node_sorter (communicator,
                                                      node_keys);
  node_sorter.sort(); /* done with node_keys */ //node_keys.clear();

  const std::vector<Parallel::DofObjectKey> & my_node_bin =
    node_sorter.bin();

  Parallel::Sort<Parallel::DofObjectKey> elem_sorter (communicator,
                                                      elem_keys);
  elem_sorter.sort(); /* done with elem_keys */ //elem_keys.clear();

  const std::vector<Parallel::DofObjectKey> & my_elem_bin =
    elem_sorter.bin();



  //-------------------------------------------------------------
  // (3) get the max value on each processor
  std::vector<Parallel::DofObjectKey>
    node_upper_bounds(communicator.size()),
    elem_upper_bounds(communicator.size());

  { // limit scope of temporaries
    std::vector<Parallel::DofObjectKey> recvbuf(2*communicator.size());
    std::vector<unsigned short int> /* do not use a vector of bools here since it is not always so! */
      empty_nodes (communicator.size()),
      empty_elem  (communicator.size());
    std::vector<Parallel::DofObjectKey> my_max(2);

    communicator.allgather (static_cast<unsigned short int>(my_node_bin.empty()), empty_nodes);
    communicator.allgather (static_cast<unsigned short int>(my_elem_bin.empty()),  empty_elem);

    if (!my_node_bin.empty()) my_max[0] = my_node_bin.back();
    if (!my_elem_bin.empty()) my_max[1] = my_elem_bin.back();

    communicator.allgather (my_max, /* identical_buffer_sizes = */ true);

    // Be careful here.  The *_upper_bounds will be used to find the processor
    // a given object belongs to.  So, if a processor contains no objects (possible!)
    // then copy the bound from the lower processor id.
    for (auto p : IntRange<processor_id_type>(0, communicator.size()))
      {
        node_upper_bounds[p] = my_max[2*p+0];
        elem_upper_bounds[p] = my_max[2*p+1];

        if (p > 0) // default hilbert index value is the OK upper bound for processor 0.
          {
            if (empty_nodes[p]) node_upper_bounds[p] = node_upper_bounds[p-1];
            if (empty_elem[p])  elem_upper_bounds[p] = elem_upper_bounds[p-1];
          }
      }
  }



  //-------------------------------------------------------------
  // (4) determine the position in the global ranking for
  //     each local object
  {
    //----------------------------------------------
    // Nodes first -- all nodes, not just local ones
    {
      // Request sets to send to each processor
      std::map<processor_id_type, std::vector<Parallel::DofObjectKey>>
        requested_ids;
      // Results to gather from each processor - kept in a map so we
      // do only one loop over nodes after all receives are done.
      std::map<dof_id_type, std::vector<dof_id_type>>
        filled_request;

      // build up list of requests
      for (const auto & node : mesh.node_ptr_range())
        {
          libmesh_assert(node);
          const Parallel::DofObjectKey hi =
            get_dofobject_key (node, bbox);
          const processor_id_type pid =
            cast_int<processor_id_type>
            (std::distance (node_upper_bounds.begin(),
                            std::lower_bound(node_upper_bounds.begin(),
                                             node_upper_bounds.end(),
                                             hi)));

          libmesh_assert_less (pid, communicator.size());

          requested_ids[pid].push_back(hi);
        }

      // The number of objects in my_node_bin on each processor
      std::vector<dof_id_type> node_bin_sizes(communicator.size());
      communicator.allgather (static_cast<dof_id_type>(my_node_bin.size()), node_bin_sizes);

      // The offset of my first global index
      dof_id_type my_offset = 0;
      for (auto pid : IntRange<processor_id_type>(0, communicator.rank()))
        my_offset += node_bin_sizes[pid];

      auto gather_functor =
        [
#ifndef NDEBUG
         & node_upper_bounds,
         & communicator,
#endif
         & my_node_bin,
         my_offset
        ]
        (processor_id_type,
         const std::vector<Parallel::DofObjectKey> & keys,
         std::vector<dof_id_type> & global_ids)
        {
          // Fill the requests
          const std::size_t keys_size = keys.size();
          global_ids.reserve(keys_size);
          for (std::size_t idx=0; idx != keys_size; idx++)
            {
              const Parallel::DofObjectKey & hi = keys[idx];
              libmesh_assert_less_equal (hi, node_upper_bounds[communicator.rank()]);

              // find the requested index in my node bin
              std::vector<Parallel::DofObjectKey>::const_iterator pos =
                std::lower_bound (my_node_bin.begin(), my_node_bin.end(), hi);
              libmesh_assert (pos != my_node_bin.end());
              libmesh_assert_equal_to (*pos, hi);

              // Finally, assign the global index based off the position of the index
              // in my array, properly offset.
              global_ids.push_back(cast_int<dof_id_type>(std::distance(my_node_bin.begin(), pos) + my_offset));
            }
        };

      auto action_functor =
        [&filled_request]
        (processor_id_type pid,
         const std::vector<Parallel::DofObjectKey> &,
         const std::vector<dof_id_type> & global_ids)
        {
          filled_request[pid] = global_ids;
        };

      // Trade requests with other processors
      const dof_id_type * ex = nullptr;
      Parallel::pull_parallel_vector_data
        (communicator, requested_ids, gather_functor, action_functor, ex);

      // We now have all the filled requests, so we can loop through our
      // nodes once and assign the global index to each one.
      {
        std::map<dof_id_type, std::vector<dof_id_type>::const_iterator>
          next_obj_on_proc;
        for (auto & p : filled_request)
          next_obj_on_proc[p.first] = p.second.begin();

        for (auto & node : mesh.node_ptr_range())
          {
            libmesh_assert(node);
            const Parallel::DofObjectKey hi =
              get_dofobject_key (node, bbox);
            const processor_id_type pid =
              cast_int<processor_id_type>
              (std::distance (node_upper_bounds.begin(),
                              std::lower_bound(node_upper_bounds.begin(),
                                               node_upper_bounds.end(),
                                               hi)));

            libmesh_assert_less (pid, communicator.size());
            libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

            const dof_id_type global_index = *next_obj_on_proc[pid];
            libmesh_assert_less (global_index, mesh.n_nodes());
            node->set_id() = global_index;

            ++next_obj_on_proc[pid];
          }
      }
    }

    //---------------------------------------------------
    // elements next -- all elements, not just local ones
    {
      // Request sets to send to each processor
      std::map<processor_id_type, std::vector<Parallel::DofObjectKey>>
        requested_ids;
      // Results to gather from each processor - kept in a map so we
      // do only one loop over elements after all receives are done.
      std::map<dof_id_type, std::vector<dof_id_type>>
        filled_request;

      for (const auto & elem : mesh.element_ptr_range())
        {
          libmesh_assert(elem);
          const Parallel::DofObjectKey hi =
            get_dofobject_key (elem, bbox);
          const processor_id_type pid =
            cast_int<processor_id_type>
            (std::distance (elem_upper_bounds.begin(),
                            std::lower_bound(elem_upper_bounds.begin(),
                                             elem_upper_bounds.end(),
                                             hi)));

          libmesh_assert_less (pid, communicator.size());

          requested_ids[pid].push_back(hi);
        }

      // The number of objects in my_elem_bin on each processor
      std::vector<dof_id_type> elem_bin_sizes(communicator.size());
      communicator.allgather (static_cast<dof_id_type>(my_elem_bin.size()), elem_bin_sizes);

      // The offset of my first global index
      dof_id_type my_offset = 0;
      for (auto pid : IntRange<processor_id_type>(0, communicator.rank()))
        my_offset += elem_bin_sizes[pid];

      auto gather_functor =
        [
#ifndef NDEBUG
         & elem_upper_bounds,
         & communicator,
#endif
         & my_elem_bin,
         my_offset
        ]
        (processor_id_type,
         const std::vector<Parallel::DofObjectKey> & keys,
         std::vector<dof_id_type> & global_ids)
        {
          // Fill the requests
          const std::size_t keys_size = keys.size();
          global_ids.reserve(keys_size);
          for (std::size_t idx=0; idx != keys_size; idx++)
            {
              const Parallel::DofObjectKey & hi = keys[idx];
              libmesh_assert_less_equal (hi, elem_upper_bounds[communicator.rank()]);

              // find the requested index in my elem bin
              std::vector<Parallel::DofObjectKey>::const_iterator pos =
                std::lower_bound (my_elem_bin.begin(), my_elem_bin.end(), hi);
              libmesh_assert (pos != my_elem_bin.end());
              libmesh_assert_equal_to (*pos, hi);

              // Finally, assign the global index based off the position of the index
              // in my array, properly offset.
              global_ids.push_back (cast_int<dof_id_type>(std::distance(my_elem_bin.begin(), pos) + my_offset));
            }
        };

      auto action_functor =
        [&filled_request]
        (processor_id_type pid,
         const std::vector<Parallel::DofObjectKey> &,
         const std::vector<dof_id_type> & global_ids)
        {
          filled_request[pid] = global_ids;
        };

      // Trade requests with other processors
      const dof_id_type * ex = nullptr;
      Parallel::pull_parallel_vector_data
        (communicator, requested_ids, gather_functor, action_functor, ex);

      // We now have all the filled requests, so we can loop through our
      // elements once and assign the global index to each one.
      {
        std::vector<std::vector<dof_id_type>::const_iterator>
          next_obj_on_proc; next_obj_on_proc.reserve(communicator.size());
        for (auto pid : IntRange<processor_id_type>(0, communicator.size()))
          next_obj_on_proc.push_back(filled_request[pid].begin());

        for (auto & elem : mesh.element_ptr_range())
          {
            libmesh_assert(elem);
            const Parallel::DofObjectKey hi =
              get_dofobject_key (elem, bbox);
            const processor_id_type pid =
              cast_int<processor_id_type>
              (std::distance (elem_upper_bounds.begin(),
                              std::lower_bound(elem_upper_bounds.begin(),
                                               elem_upper_bounds.end(),
                                               hi)));

            libmesh_assert_less (pid, communicator.size());
            libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

            const dof_id_type global_index = *next_obj_on_proc[pid];
            libmesh_assert_less (global_index, mesh.n_elem());
            elem->set_id() = global_index;

            ++next_obj_on_proc[pid];
          }
      }
    }
  }
}
#else // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI
template <typename RealType>
void MeshCommunication::assign_global_indices (MeshBaseTempl<RealType> &) const
{
}
#endif // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI


#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
template <typename RealType>
void MeshCommunication::check_for_duplicate_global_indices (MeshBaseTempl<RealType> & mesh) const
{
  LOG_SCOPE ("check_for_duplicate_global_indices()", "MeshCommunication");

  // Global bounding box.  We choose the nodal bounding box for
  // backwards compatibility; the element bounding box may be looser
  // on curved elements.
  BoundingBox bbox =
    MeshTools::create_nodal_bounding_box (mesh);


  std::vector<Parallel::DofObjectKey>
    node_keys, elem_keys;

  {
    // Nodes first
    {
      ConstNodeRangeTempl<RealType> nr (mesh.local_nodes_begin(),
                         mesh.local_nodes_end());
      node_keys.resize (nr.size());
      Threads::parallel_for (nr, ComputeHilbertKeys<RealType> (bbox, node_keys));

      // It's O(N^2) to check that these keys don't duplicate before the
      // sort...
      typename MeshBaseTempl<RealType>::const_node_iterator nodei = mesh.local_nodes_begin();
      for (std::size_t i = 0; i != node_keys.size(); ++i, ++nodei)
        {
          typename MeshBaseTempl<RealType>::const_node_iterator nodej = mesh.local_nodes_begin();
          for (std::size_t j = 0; j != i; ++j, ++nodej)
            {
              if (node_keys[i] == node_keys[j])
                {
                  CFixBitVec icoords[3], jcoords[3];
                  get_hilbert_coords(**nodej, bbox, jcoords);
                  libMesh::err <<
                    "node " << (*nodej)->id() << ", " <<
                    *(PointTempl<RealType> *)(*nodej) << " has HilbertIndices " <<
                    node_keys[j] << std::endl;
                  get_hilbert_coords(**nodei, bbox, icoords);
                  libMesh::err <<
                    "node " << (*nodei)->id() << ", " <<
                    *(PointTempl<RealType> *)(*nodei) << " has HilbertIndices " <<
                    node_keys[i] << std::endl;
                  libmesh_error_msg("Error: nodes with duplicate Hilbert keys!");
                }
            }
        }
    }

    // Elements next
    {
      ConstElemRangeTempl<RealType> er (mesh.local_elements_begin(),
                         mesh.local_elements_end());
      elem_keys.resize (er.size());
      Threads::parallel_for (er, ComputeHilbertKeys<RealType> (bbox, elem_keys));

      // For elements, the keys can be (and in the case of TRI, are
      // expected to be) duplicates, but only if the elements are at
      // different levels
      typename MeshBaseTempl<RealType>::const_element_iterator elemi = mesh.local_elements_begin();
      for (std::size_t i = 0; i != elem_keys.size(); ++i, ++elemi)
        {
          typename MeshBaseTempl<RealType>::const_element_iterator elemj = mesh.local_elements_begin();
          for (std::size_t j = 0; j != i; ++j, ++elemj)
            {
              if ((elem_keys[i] == elem_keys[j]) &&
                  ((*elemi)->level() == (*elemj)->level()))
                {
                  libMesh::err <<
                    "level " << (*elemj)->level() << " elem\n" <<
                    (**elemj) << " centroid " <<
                    (*elemj)->centroid() << " has HilbertIndices " <<
                    elem_keys[j] << " or " <<
                    get_dofobject_key((*elemj), bbox) <<
                    std::endl;
                  libMesh::err <<
                    "level " << (*elemi)->level() << " elem\n" <<
                    (**elemi) << " centroid " <<
                    (*elemi)->centroid() << " has HilbertIndices " <<
                    elem_keys[i] << " or " <<
                    get_dofobject_key((*elemi), bbox) <<
                    std::endl;
                  libmesh_error_msg("Error: level " << (*elemi)->level() << " elements with duplicate Hilbert keys!");
                }
            }
        }
    }
  } // done checking Hilbert keys
}
#else // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI
template <typename RealType>
void MeshCommunication::check_for_duplicate_global_indices (MeshBaseTempl<RealType> &) const
{
}
#endif // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI

#if defined(LIBMESH_HAVE_LIBHILBERT) && defined(LIBMESH_HAVE_MPI)
template <typename RealType, typename ForwardIterator>
void MeshCommunication::find_local_indices (const libMesh::BoundingBoxTempl<RealType> & bbox,
                                            const ForwardIterator & begin,
                                            const ForwardIterator & end,
                                            std::unordered_map<dof_id_type, dof_id_type> & index_map) const
{
  LOG_SCOPE ("find_local_indices()", "MeshCommunication");

  // This method determines id-agnostic local indices
  // for nodes and elements by sorting Hilbert keys.

  index_map.clear();

  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  // These aren't trivial to compute, and we will need them again.
  // But the binsort will sort the input vector, trashing the order
  // that we'd like to rely on.  So, two vectors...
  std::map<Parallel::DofObjectKey, dof_id_type> hilbert_keys;
  {
    LOG_SCOPE("local_hilbert_indices", "MeshCommunication");
    for (ForwardIterator it=begin; it!=end; ++it)
      {
        const Parallel::DofObjectKey hi(get_dofobject_key ((*it), bbox));
        hilbert_keys.emplace(hi, (*it)->id());
      }
  }

  {
    dof_id_type cnt = 0;
    for (auto key_val : hilbert_keys)
      index_map[key_val.second] = cnt++;
  }
}


template <typename RealType, typename ForwardIterator>
void MeshCommunication::find_global_indices (const Parallel::Communicator & communicator,
                                             const libMesh::BoundingBoxTempl<RealType> & bbox,
                                             const ForwardIterator & begin,
                                             const ForwardIterator & end,
                                             std::vector<dof_id_type> & index_map) const
{
  LOG_SCOPE ("find_global_indices()", "MeshCommunication");

  // This method determines partition-agnostic global indices
  // for nodes and elements.

  // Algorithm:
  // (1) compute the Hilbert key for each local node/element
  // (2) perform a parallel sort of the Hilbert key
  // (3) get the min/max value on each processor
  // (4) determine the position in the global ranking for
  //     each local object
  index_map.clear();
  std::size_t n_objects = std::distance (begin, end);
  index_map.reserve(n_objects);

  //-------------------------------------------------------------
  // (1) compute Hilbert keys
  // These aren't trivial to compute, and we will need them again.
  // But the binsort will sort the input vector, trashing the order
  // that we'd like to rely on.  So, two vectors...
  std::vector<Parallel::DofObjectKey>
    sorted_hilbert_keys,
    hilbert_keys;
  sorted_hilbert_keys.reserve(n_objects);
  hilbert_keys.reserve(n_objects);
  {
    LOG_SCOPE("compute_hilbert_indices()", "MeshCommunication");
    for (ForwardIterator it=begin; it!=end; ++it)
      {
        const Parallel::DofObjectKey hi(get_dofobject_key (*it, bbox));
        hilbert_keys.push_back(hi);

        if ((*it)->processor_id() == communicator.rank())
          sorted_hilbert_keys.push_back(hi);

        // someone needs to take care of unpartitioned objects!
        if ((communicator.rank() == 0) &&
            ((*it)->processor_id() == DofObject::invalid_processor_id))
          sorted_hilbert_keys.push_back(hi);
      }
  }

  //-------------------------------------------------------------
  // (2) parallel sort the Hilbert keys
  START_LOG ("parallel_sort()", "MeshCommunication");
  Parallel::Sort<Parallel::DofObjectKey> sorter (communicator,
                                                 sorted_hilbert_keys);
  sorter.sort();
  STOP_LOG ("parallel_sort()", "MeshCommunication");
  const std::vector<Parallel::DofObjectKey> & my_bin = sorter.bin();

  // The number of objects in my_bin on each processor
  std::vector<unsigned int> bin_sizes(communicator.size());
  communicator.allgather (static_cast<unsigned int>(my_bin.size()), bin_sizes);

  // The offset of my first global index
  unsigned int my_offset = 0;
  for (auto pid : IntRange<processor_id_type>(0, communicator.rank()))
    my_offset += bin_sizes[pid];

  //-------------------------------------------------------------
  // (3) get the max value on each processor
  std::vector<Parallel::DofObjectKey>
    upper_bounds(1);

  if (!my_bin.empty())
    upper_bounds[0] = my_bin.back();

  communicator.allgather (upper_bounds, /* identical_buffer_sizes = */ true);

  // Be careful here.  The *_upper_bounds will be used to find the processor
  // a given object belongs to.  So, if a processor contains no objects (possible!)
  // then copy the bound from the lower processor id.
  for (auto p : IntRange<processor_id_type>(1, communicator.size()))
    if (!bin_sizes[p]) upper_bounds[p] = upper_bounds[p-1];


  //-------------------------------------------------------------
  // (4) determine the position in the global ranking for
  //     each local object
  {
    //----------------------------------------------
    // all objects, not just local ones

    // Request sets to send to each processor
    std::map<processor_id_type, std::vector<Parallel::DofObjectKey>>
      requested_ids;
    // Results to gather from each processor
    std::map<processor_id_type, std::vector<dof_id_type>>
      filled_request;

    // build up list of requests
    std::vector<Parallel::DofObjectKey>::const_iterator hi =
      hilbert_keys.begin();

    for (ForwardIterator it = begin; it != end; ++it)
      {
        libmesh_assert (hi != hilbert_keys.end());

        std::vector<Parallel::DofObjectKey>::iterator lb =
          std::lower_bound(upper_bounds.begin(), upper_bounds.end(),
                           *hi);

        const processor_id_type pid =
          cast_int<processor_id_type>
          (std::distance (upper_bounds.begin(), lb));

        libmesh_assert_less (pid, communicator.size());

        requested_ids[pid].push_back(*hi);

        ++hi;
        // go ahead and put pid in index_map, that way we
        // don't have to repeat the std::lower_bound()
        index_map.push_back(pid);
      }

  auto gather_functor =
    [
#ifndef NDEBUG
     & upper_bounds,
     & communicator,
#endif
     & bbox,
     & my_bin,
       my_offset
    ]
    (processor_id_type, const std::vector<Parallel::DofObjectKey> & keys,
     std::vector<dof_id_type> & global_ids)
    {
      // Ignore unused lambda capture warnings in devel mode
      libmesh_ignore(bbox);

      // Fill the requests
      const std::size_t keys_size = keys.size();
      global_ids.clear();
      global_ids.reserve(keys_size);
      for (std::size_t idx=0; idx != keys_size; idx++)
        {
          const Parallel::DofObjectKey & hilbert_indices = keys[idx];
          libmesh_assert_less_equal (hilbert_indices, upper_bounds[communicator.rank()]);

          // find the requested index in my node bin
          std::vector<Parallel::DofObjectKey>::const_iterator pos =
            std::lower_bound (my_bin.begin(), my_bin.end(), hilbert_indices);
          libmesh_assert (pos != my_bin.end());
#ifdef DEBUG
          // If we could not find the requested Hilbert index in
          // my_bin, something went terribly wrong, possibly the
          // Mesh was displaced differently on different processors,
          // and therefore the Hilbert indices don't agree.
          if (*pos != hilbert_indices)
            {
              // The input will be hilbert_indices.  We convert it
              // to BitVecType using the operator= provided by the
              // BitVecType class. BitVecType is a CBigBitVec!
              Hilbert::BitVecType input;
#ifdef LIBMESH_ENABLE_UNIQUE_ID
              input = hilbert_indices.first;
#else
              input = hilbert_indices;
#endif

              // Get output in a vector of CBigBitVec
              std::vector<CBigBitVec> output(3);

              // Call the indexToCoords function
              Hilbert::indexToCoords(output.data(), 8*sizeof(Hilbert::inttype), 3, input);

              // The entries in the output racks are integers in the
              // range [0, Hilbert::inttype::max] which can be
              // converted to floating point values in [0,1] and
              // finally to actual values using the bounding box.
              const Real max_int_as_real =
                static_cast<Real>(std::numeric_limits<Hilbert::inttype>::max());

              // Get the points in [0,1]^3.  The zeroth rack of each entry in
              // 'output' maps to the normalized x, y, and z locations,
              // respectively.
              Point p_hat(static_cast<Real>(output[0].racks()[0]) / max_int_as_real,
                          static_cast<Real>(output[1].racks()[0]) / max_int_as_real,
                          static_cast<Real>(output[2].racks()[0]) / max_int_as_real);

              // Convert the points from [0,1]^3 to their actual (x,y,z) locations
              Real
                xmin = bbox.first(0),
                xmax = bbox.second(0),
                ymin = bbox.first(1),
                ymax = bbox.second(1),
                zmin = bbox.first(2),
                zmax = bbox.second(2);

              // Convert the points from [0,1]^3 to their actual (x,y,z) locations
              Point p(xmin + (xmax-xmin)*p_hat(0),
                      ymin + (ymax-ymin)*p_hat(1),
                      zmin + (zmax-zmin)*p_hat(2));

              libmesh_error_msg("Could not find hilbert indices: "
                                << hilbert_indices
                                << " corresponding to point " << p);
            }
#endif

          // Finally, assign the global index based off the position of the index
          // in my array, properly offset.
          global_ids.push_back (cast_int<dof_id_type>(std::distance(my_bin.begin(), pos) + my_offset));
        }
    };

  auto action_functor =
    [&filled_request]
    (processor_id_type pid,
     const std::vector<Parallel::DofObjectKey> &,
     const std::vector<dof_id_type> & global_ids)
    {
      filled_request[pid] = global_ids;
    };

  const dof_id_type * ex = nullptr;
  Parallel::pull_parallel_vector_data
    (communicator, requested_ids, gather_functor, action_functor, ex);

    // We now have all the filled requests, so we can loop through our
    // nodes once and assign the global index to each one.
    {
      std::vector<std::vector<dof_id_type>::const_iterator>
        next_obj_on_proc; next_obj_on_proc.reserve(communicator.size());
      for (auto pid : IntRange<processor_id_type>(0, communicator.size()))
        next_obj_on_proc.push_back(filled_request[pid].begin());

      unsigned int cnt=0;
      for (ForwardIterator it = begin; it != end; ++it, cnt++)
        {
          const processor_id_type pid = cast_int<processor_id_type>
            (index_map[cnt]);

          libmesh_assert_less (pid, communicator.size());
          libmesh_assert (next_obj_on_proc[pid] != filled_request[pid].end());

          const dof_id_type global_index = *next_obj_on_proc[pid];
          index_map[cnt] = global_index;

          ++next_obj_on_proc[pid];
        }
    }
  }

  libmesh_assert_equal_to(index_map.size(), n_objects);
}
#else // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI
template <typename RealType, typename ForwardIterator>
void MeshCommunication::find_local_indices (const libMesh::BoundingBoxTempl<RealType> &,
                                            const ForwardIterator & begin,
                                            const ForwardIterator & end,
                                            std::unordered_map<dof_id_type, dof_id_type> & index_map) const
{
  index_map.clear();
  dof_id_type index = 0;
  for (ForwardIterator it=begin; it!=end; ++it)
    index_map[(*it)->id()] = (index++);
}

template <typename RealType, typename ForwardIterator>
void MeshCommunication::find_global_indices (const Parallel::Communicator &,
                                             const libMesh::BoundingBoxTempl<RealType> &,
                                             const ForwardIterator & begin,
                                             const ForwardIterator & end,
                                             std::vector<dof_id_type> & index_map) const
{
  index_map.clear();
  index_map.reserve(std::distance (begin, end));
  unsigned int index = 0;
  for (ForwardIterator it=begin; it!=end; ++it)
    index_map.push_back(index++);
}
#endif // LIBMESH_HAVE_LIBHILBERT, LIBMESH_HAVE_MPI

#define MESH_GLOBAL_INDICES_INSTANTIATIONS(RealType) \
  template void MeshCommunication::find_local_indices(                                             \
      const libMesh::BoundingBoxTempl<RealType> &,                                                 \
      const typename MeshBaseTempl<RealType>::const_element_iterator &,                            \
      const typename MeshBaseTempl<RealType>::const_element_iterator &,                            \
      std::unordered_map<dof_id_type, dof_id_type> &) const;                                       \
  template void MeshCommunication::find_global_indices(                                            \
      const Parallel::Communicator & communicator,                                                 \
      const libMesh::BoundingBoxTempl<RealType> &,                                                 \
      const typename MeshBaseTempl<RealType>::const_element_iterator &,                            \
      const typename MeshBaseTempl<RealType>::const_element_iterator &,                            \
      std::vector<dof_id_type> &) const;                                                           \
  template void MeshCommunication::find_global_indices(                                            \
      const Parallel::Communicator & communicator,                                                 \
      const libMesh::BoundingBoxTempl<RealType> &,                                                 \
      const typename MeshBaseTempl<RealType>::element_iterator &,                                  \
      const typename MeshBaseTempl<RealType>::element_iterator &,                                  \
      std::vector<dof_id_type> &) const;                                                           \
  template void MeshCommunication::find_global_indices(                                            \
      const Parallel::Communicator & communicator,                                                 \
      const libMesh::BoundingBoxTempl<RealType> &,                                                 \
      const typename MeshBaseTempl<RealType>::const_node_iterator &,                               \
      const typename MeshBaseTempl<RealType>::const_node_iterator &,                               \
      std::vector<dof_id_type> &) const;                                                           \
  template void MeshCommunication::find_global_indices(                                            \
      const Parallel::Communicator & communicator,                                                 \
      const libMesh::BoundingBoxTempl<RealType> &,                                                 \
      const typename MeshBaseTempl<RealType>::node_iterator &,                                     \
      const typename MeshBaseTempl<RealType>::node_iterator &,                                     \
      std::vector<dof_id_type> &) const;                                                           \
  template void MeshCommunication::assign_global_indices(MeshBaseTempl<RealType> &) const;         \
  template void MeshCommunication::check_for_duplicate_global_indices(MeshBaseTempl<RealType> &)   \
      const


} // namespace libMesh

#endif // LIBMESH_MESH_COMMUNICATION_GLOBAL_INDICES_IMPL_H
