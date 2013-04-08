// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_tools.h"

// C++ Includes   -----------------------------------

namespace libMesh
{

// Forward declarations
class MeshBase;
class ParallelMesh;
class BoundaryInfo;
template <typename T> class LocationMap;



/**
 * This is the \p MeshCommunication class.  It handles all the details
 * of communicating mesh information from one processor to another.  All
 * parallelization of the \p Mesh data structures is done via this class.
 *
 * @author Benjamin S. Kirk, 2003
 */


// ------------------------------------------------------------
// MeshCommunication class definition
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
   * Clears all data structures and returns to a pristine state.
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
  void broadcast (MeshBase &) const;

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
   */
  void redistribute (ParallelMesh &) const;

  /**
   *
   */
  void gather_neighboring_elements (ParallelMesh &) const;

  /**
   * This method takes an input \p ParallelMesh which may be
   * distributed among all the processors.  Each processor then
   * sends its local nodes and elements to the other processors.
   * The end result is that a previously distributed \p ParallelMesh
   * will be serialized on each processor.  Since this method is
   * collective it must be called by all processors.
   */
  void allgather (ParallelMesh &) const;

  /**
   * This method takes an input \p ParallelMesh which may be
   * distributed among all the processors.  Each processor
   * deletes all elements which are neither local elements nor "ghost"
   * elements which touch local elements, and deletes all nodes which
   * are not contained in local or ghost elements.
   * The end result is that a previously serial \p ParallelMesh
   * will be distributed between processors.  Since this method is
   * collective it must be called by all processors.
   *
   * The std::set is a list of extra elements that you _don't_ want
   * to delete.  These will be left on the current processor along with
   * local elements and ghosted neighbors.
   */
  void delete_remote_elements (ParallelMesh &, const std::set<Elem *> &) const;

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
  void assign_global_indices (MeshBase &) const;


  /**
   * This method determines a globally unique, partition-agnostic
   * index for each object in the input range.
   */
  template <typename ForwardIterator>
  void find_global_indices (const Parallel::Communicator &communicator,
			    const MeshTools::BoundingBox &,
			    const ForwardIterator &,
			    const ForwardIterator &,
			    std::vector<dof_id_type> &) const;

  /**
   * Copy ids of ghost elements from their local processors.
   */
  void make_elems_parallel_consistent (MeshBase &);

  /**
   * Assuming all ids on local nodes are globally unique, and
   * assuming all processor ids are parallel consistent, this function makes
   * all other ids parallel consistent.
   */
  void make_node_ids_parallel_consistent (MeshBase &,
                                          LocationMap<Node> &);

  /**
   * Assuming all processor ids on nodes touching local elements
   * are parallel consistent, this function makes all other processor ids
   * parallel consistent as well.
   */
  void make_node_proc_ids_parallel_consistent (MeshBase &,
                                               LocationMap<Node> &);

  /**
   * Copy processor_ids and ids on ghost nodes from their
   * local processors.  This is an internal function of MeshRefinement
   * which turns out to be useful for other code which wants to add
   * nodes to a distributed mesh.
   */
  void make_nodes_parallel_consistent (MeshBase &,
                                       LocationMap<Node> &);
};



//--------------------------------------------------------------
// MeshCommunication inline members


} // namespace libMesh

#endif // LIBMESH_MESH_COMMUNICATION_H
