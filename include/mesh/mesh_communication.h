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



#ifndef LIBMESH_MESH_COMMUNICATION_H
#define LIBMESH_MESH_COMMUNICATION_H

// Local Includes
#include "libmesh/compare_elems_by_level.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_tools.h"

// C++ Includes
#include <unordered_map>

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
template <typename> class DistributedMeshTempl;

// This is for backwards compatibility, but if your code relies on
// forward declarations in our headers then fix it.
template <typename> class ParallelMeshTempl;

/**
 * This is the \p MeshCommunication class.  It handles all the details
 * of communicating mesh information from one processor to another.  All
 * parallelization of the \p Mesh data structures is done via this class.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 */
class MeshCommunication
{
public:

  /**
   * Constructor.
   */
  MeshCommunication () {}

  /**
   * Destructor.
   */
  ~MeshCommunication () {}

  /**
   * Clears all data structures and resets to a pristine state.
   */
  void clear ();

  //   /**
  //    * Finds all the processors that may contain
  //    * elements that neighbor my elements.  This list
  //    * is guaranteed to include all processors that border
  //    * any of my elements, but may include additional ones as
  //    * well.  This method computes bounding boxes for the
  //    * elements on each processor and checks for overlaps.
  //    */
  //   void find_neighboring_processors(const MeshBase &);

  /**
   * This method takes a mesh (which is assumed to reside on
   * processor 0) and broadcasts it to all the other processors.
   * It also broadcasts any boundary information the mesh has
   * associated with it.
   */
  template <typename RealType>
  void broadcast (MeshBaseTempl<RealType> &) const;

  /**
   * This method takes a parallel distributed mesh and redistributes
   * the elements.  Specifically, any elements stored on a given
   * processor are sent to the processor which "owns" them.  Similarly,
   * any elements assigned to the current processor but stored on
   * another are received. Once this step is completed any required ghost
   * elements are updated.  The final result is that each processor stores
   * only the elements it actually owns and any ghost elements required
   * to satisfy data dependencies. This method can be invoked after a
   * partitioning step to affect the new partitioning.
   *
   * Redistribution can also be done with newly coarsened elements'
   * neighbors only.
   */
  template <typename RealType>
  void redistribute (DistributedMeshTempl<RealType> & mesh,
                     bool newly_coarsened_only = false) const;

  /**
   *
   */
  template <typename RealType>
  void gather_neighboring_elements (DistributedMeshTempl<RealType> &) const;

  /**
   * Examine a just-coarsened mesh, and for any newly-coarsened elements,
   * send the associated ghosted elements to the processor which needs them.
   */
  template <typename RealType>
  void send_coarse_ghosts (MeshBaseTempl<RealType> &) const;

  /**
   * This method takes an input \p DistributedMesh which may be
   * distributed among all the processors.  Each processor then
   * sends its local nodes and elements to processor \p root_id.
   * The end result is that a previously distributed \p DistributedMesh
   * will be serialized on processor \p root_id.  Since this method is
   * collective it must be called by all processors. For the special
   * case of \p root_id equal to \p DofObject::invalid_processor_id
   * this function performs an allgather.
   */
  template <typename RealType>
  void gather (const processor_id_type root_id, DistributedMeshTempl<RealType> &) const;

  /**
   * This method takes an input \p DistributedMesh which may be
   * distributed among all the processors.  Each processor then
   * sends its local nodes and elements to the other processors.
   * The end result is that a previously distributed \p DistributedMesh
   * will be serialized on each processor.  Since this method is
   * collective it must be called by all processors.
   */
  template <typename RealType>
  void allgather (DistributedMeshTempl<RealType> & mesh) const
  { MeshCommunication::gather(DofObject::invalid_processor_id, mesh); }

  /**
   * This method takes an input \p DistributedMesh which may be
   * distributed among all the processors.  Each processor
   * deletes all elements which are neither local elements nor "ghost"
   * elements which touch local elements, and deletes all nodes which
   * are not contained in local or ghost elements.
   * The end result is that a previously serial \p DistributedMesh
   * will be distributed between processors.  Since this method is
   * collective it must be called by all processors.
   *
   * The std::set is a list of extra elements that you _don't_ want
   * to delete.  These will be left on the current processor along with
   * local elements and ghosted neighbors.
   */
  template <typename RealType>
  void delete_remote_elements (DistributedMeshTempl<RealType> &, const std::set<ElemTempl<RealType> *> &) const;

  /**
   * This method assigns globally unique, partition-agnostic
   * indices to the nodes and elements in the mesh.  The approach
   * is to compute the Hilbert space-filling curve key and use its
   * value to assign an index in [0,N_global). Since the Hilbert key
   * is unique for each spatial location, two objects occupying the
   * same location will be assigned the same global id.  Thus, this
   * method can also be useful for identifying duplicate nodes
   * which may occur during parallel refinement.
   */
  template <typename RealType>
  void assign_global_indices (MeshBaseTempl<RealType> &) const;

  /**
   * Throw an error if we have any index clashes in the numbering used by
   * assign_global_indices.
   */
  template <typename RealType>
  void check_for_duplicate_global_indices (MeshBaseTempl<RealType> & ) const;

  /**
   * This method determines a locally unique, contiguous
   * index for each object in the input range.
   */
  template <typename RealType, typename ForwardIterator>
  void find_local_indices (const libMesh::BoundingBoxTempl<RealType> &,
                           const ForwardIterator &,
                           const ForwardIterator &,
                           std::unordered_map<dof_id_type, dof_id_type> &) const;

  /**
   * This method determines a globally unique, partition-agnostic
   * index for each object in the input range.
   */
  template <typename RealType, typename ForwardIterator>
  void find_global_indices (const Parallel::Communicator & communicator,
                            const libMesh::BoundingBoxTempl<RealType> &,
                            const ForwardIterator &,
                            const ForwardIterator &,
                            std::vector<dof_id_type> &) const;

  /**
   * Copy ids of ghost elements from their local processors.
   */
  template <typename RealType>
  void make_elems_parallel_consistent (MeshBaseTempl<RealType> &);

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Copy p levels of ghost elements from their local processors.
   */
  template <typename RealType>
  void make_p_levels_parallel_consistent (MeshBaseTempl<RealType> &);
#endif // LIBMESH_ENABLE_AMR

  /**
   * Assuming all ids on local nodes are globally unique, and
   * assuming all processor ids are parallel consistent, this function makes
   * all other ids parallel consistent.
   */
  template <typename RealType>
  void make_node_ids_parallel_consistent (MeshBaseTempl<RealType> &);

  /**
   * Assuming all unique_ids on local nodes are globally unique, and
   * assuming all processor ids are parallel consistent, this function makes
   * all ghost unique_ids parallel consistent.
   */
  template <typename RealType>
  void make_node_unique_ids_parallel_consistent (MeshBaseTempl<RealType> &);

  /**
   * Assuming all processor ids on nodes touching local elements
   * are parallel consistent, this function makes all other processor ids
   * parallel consistent as well.
   */
  template <typename RealType>
  void make_node_proc_ids_parallel_consistent (MeshBaseTempl<RealType> &);

  /**
   * Assuming all processor ids on nodes touching local elements
   * are parallel consistent, this function makes processor ids
   * on new nodes on other processors parallel consistent as well.
   */
  template <typename RealType>
  void make_new_node_proc_ids_parallel_consistent (MeshBaseTempl<RealType> &);

  /**
   * Copy processor_ids and ids on ghost nodes from their
   * local processors.  This is useful for code which wants to add
   * nodes to a distributed mesh.
   */
  template <typename RealType>
  void make_nodes_parallel_consistent (MeshBaseTempl<RealType> &);

  /**
   * Copy processor_ids and ids on new nodes from their local
   * processors.
   */
  template <typename RealType>
  void make_new_nodes_parallel_consistent (MeshBaseTempl<RealType> &);
};


// Related utilities

// Ask a mesh's ghosting functors to insert into a set all elements
// that are either on or connected to processor id \p pid.  Ask only
// for elements in the range from \p elem_it before \p elem_end;
// typically this may be mesh.active_pid_elements_*(pid)
template <typename RealType>
void query_ghosting_functors(const MeshBaseTempl<RealType> & mesh,
                             processor_id_type pid,
                             typename MeshBaseTempl<RealType>::const_element_iterator elem_it,
                             typename MeshBaseTempl<RealType>::const_element_iterator elem_end,
                             std::set<const ElemTempl<RealType> *, CompareElemIdsByLevel> & connected_elements);

// Take a set of elements and insert all immediate
// children of elements in the given range
template <typename RealType>
void connect_children(const MeshBaseTempl<RealType> & mesh,
                      typename MeshBaseTempl<RealType>::const_element_iterator elem_it,
                      typename MeshBaseTempl<RealType>::const_element_iterator elem_end,
                      std::set<const ElemTempl<RealType> *, CompareElemIdsByLevel> & connected_elements);

// Take a set of elements and insert all elements' ancestors and
// subactive descendants as well.
template <typename RealType>
void connect_families(std::set<const ElemTempl<RealType> *, CompareElemIdsByLevel> & connected_elements);

// Take a set of elements and create a set of connected nodes.
template <typename RealType>
void reconnect_nodes (const std::set<const ElemTempl<RealType> *, CompareElemIdsByLevel> & connected_elements,
                      std::set<const NodeTempl<RealType> *> & connected_nodes);



//--------------------------------------------------------------
// MeshCommunication inline members


} // namespace libMesh

#endif // LIBMESH_MESH_COMMUNICATION_H
