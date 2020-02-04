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



#ifndef LIBMESH_MESH_TOOLS_H
#define LIBMESH_MESH_TOOLS_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/bounding_box.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh_base.h"

#ifdef LIBMESH_FORWARD_DECLARE_ENUMS
namespace libMesh
{
enum ElemType : int;
}
#else
#include "libmesh/enum_elem_type.h"
#endif

// C++ Includes
#include <limits>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace libMesh
{

// forward declarations
class Sphere;
template <typename> class ElemTempl;
template <typename> class NodeTempl;
template <typename> class NodeTempl;
template <typename> class PointTempl;

/**
 * Utility functions for operations on a \p Mesh object.  Here is where
 * useful functions for interfacing with a \p Mesh should be defined.
 * In general this namespace should be used to prevent the \p Mesh class
 * from becoming too cluttered.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
namespace MeshTools
{

/**
 * Backwards compatibility with forward declarations.
 *
 * \deprecated Use libMesh::BoundingBox instead.
 */
#ifdef LIBMESH_ENABLE_DEPRECATED
class BoundingBox : public libMesh::BoundingBoxTempl<Real>
{
public:
  BoundingBox (const PointTempl<Real> & new_min,
               const PointTempl<Real> & new_max) :
    libMesh::BoundingBoxTempl<Real>(new_min, new_max) {
    libmesh_deprecated(); // Switch to libMesh::BoundingBox
  }

  BoundingBox (const std::pair<PointTempl<Real>, PointTempl<Real>> & bbox) :
    libMesh::BoundingBoxTempl<Real>(bbox) {
    libmesh_deprecated(); // Switch to libMesh::BoundingBox
  }

  BoundingBox () {
    libmesh_deprecated(); // Switch to libMesh::BoundingBox
  }
};
#endif


/**
 * \returns The sum over all the elements of the number
 * of nodes per element.
 *
 * This can be useful for partitioning hybrid meshes.  A feasible load
 * balancing scheme is to keep the weight per processor as uniform as
 * possible.
 */
template <typename RealType>
dof_id_type total_weight (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The sum over all the elements on processor \p pid
 * of nodes per element.
 *
 * This can be useful for partitioning hybrid meshes.  A feasible load
 * balancing scheme is to keep the weight per processor as uniform as
 * possible.
 */
template <typename RealType>
dof_id_type weight (const MeshBaseTempl<RealType> & mesh,
                    const processor_id_type pid);

template <typename RealType>
inline
dof_id_type weight (const MeshBaseTempl<RealType> & mesh)
{ return MeshTools::weight (mesh, mesh.processor_id()); }

/**
 * After calling this function the input vector \p nodes_to_elem_map
 * will contain the node to element connectivity.  That is to say
 * \p nodes_to_elem_map[i][j] is the global number of \f$ j^{th} \f$
 * element connected to node \p i.
 */
template <typename RealType>
void build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                              std::vector<std::vector<dof_id_type>> & nodes_to_elem_map);

/**
 * The same, except element pointers are returned instead of indices.
 */
template <typename RealType>
void build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                              std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map);

/**
 * After calling this function the input map \p nodes_to_elem_map
 * will contain the node to element connectivity.  That is to say
 * \p nodes_to_elem_map[i][j] is the global number of \f$ j^{th} \f$
 * element connected to node \p i.
 */
template <typename RealType>
void build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                              std::unordered_map<dof_id_type, std::vector<dof_id_type>> & nodes_to_elem_map);

/**
 * The same, except element pointers are returned instead of indices.
 */
template <typename RealType>
void build_nodes_to_elem_map (const MeshBaseTempl<RealType> & mesh,
                              std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map);


//   /**
//    * Calling this function on a 2D mesh will convert all the elements
//    * to triangles.  \p QUAD4s will be converted to \p TRI3s, \p QUAD8s
//    * and \p QUAD9s will be converted to \p TRI6s.
//    */
//   void all_tri (MeshBaseTempl<RealType> & mesh);

#ifdef LIBMESH_ENABLE_DEPRECATED
/**
 * Fills the vector "on_boundary" with flags that tell whether each node
 * is on the domain boundary (true)) or not (false).
 */
template <typename RealType>
void find_boundary_nodes (const MeshBaseTempl<RealType> & mesh,
                          std::vector<bool> & on_boundary);
#endif

/**
 * Returns a std::set containing Node IDs for all of the boundary nodes
 */
template <typename RealType>
std::unordered_set<dof_id_type> find_boundary_nodes(const MeshBaseTempl<RealType> & mesh);

/**
 * Returns a std::set containing Node IDs for all of the block boundary nodes
 *
 * A "block boundary node" is a node that is connected to elemenents from 2
 * or more blockse
 */
template <typename RealType>
std::unordered_set<dof_id_type> find_block_boundary_nodes(const MeshBaseTempl<RealType> & mesh);

/**
 * \returns Two points defining a cartesian box that bounds the
 * mesh.  The first entry in the pair is the minimum, the second
 * is the maximum.
 *
 * \deprecated Use create_bounding_box() instead.
 */
#ifdef LIBMESH_ENABLE_DEPRECATED
BoundingBox
bounding_box (const MeshBaseTempl<Real> & mesh);
#endif

/**
 * The same functionality as the deprecated MeshTools::bounding_box().
 *
 * \returns The non-deprecated libMesh::BoundingBox type.
 */
template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
create_bounding_box (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns A bounding sphere for \p mesh instead of a bounding box.
 */
template <typename RealType>
Sphere
bounding_sphere (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns Two points defining a cartesian box that bounds the
 * nodes of the mesh.
 *
 * In the case of curved elements, this box might *not* bound the
 * elements of the mesh.
 */
template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
create_nodal_bounding_box (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns Two points defining a cartesian box that bounds the
 * elements belonging to the local processor.
 *
 * Unlike the other bounding box creation functions, this does *not*
 * need to be run in parallel, because this is the only function we
 * can guarantee can be resolved with only local information.
 */
template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
create_local_bounding_box (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns Two points defining a cartesian box that bounds the
 * elements belonging to processor pid.
 *
 * \deprecated Use create_processor_bounding_box() instead.
 */
#ifdef LIBMESH_ENABLE_DEPRECATED
BoundingBox
processor_bounding_box (const MeshBaseTempl<Real> & mesh,
                        const processor_id_type pid);
#endif

/**
 * The same functionality as the deprecated MeshTools::processor_bounding_box().
 *
 * \returns The non-deprecated libMesh::BoundingBox type.
 */
template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
create_processor_bounding_box (const MeshBaseTempl<RealType> & mesh,
                               const processor_id_type pid);

/**
 * \returns A processor bounding sphere instead of a processor bounding box.
 */
template <typename RealType>
Sphere
processor_bounding_sphere (const MeshBaseTempl<RealType> & mesh,
                           const processor_id_type pid);

/**
 * \returns Two points defining a Cartesian box that bounds the
 * elements belonging to subdomain sid.
 *
 * \deprecated Use create_subdomain_bounding_box() instead.
 */
#ifdef LIBMESH_ENABLE_DEPRECATED
BoundingBox
subdomain_bounding_box (const MeshBaseTempl<Real> & mesh,
                        const subdomain_id_type sid);
#endif


/**
 * The same functionality as the deprecated MeshTools::subdomain_bounding_box().
 *
 * \returns The non-deprecated libMesh::BoundingBox type.
 */
template <typename RealType>
libMesh::BoundingBoxTempl<RealType>
create_subdomain_bounding_box (const MeshBaseTempl<RealType> & mesh,
                               const subdomain_id_type sid);

/**
 * \returns A subdomain bounding sphere instead of a subdomain bounding box.
 */
template <typename RealType>
Sphere
subdomain_bounding_sphere (const MeshBaseTempl<RealType> & mesh,
                           const subdomain_id_type sid);


/**
 * Fills in a vector of all element types in the mesh.  Implemented
 * in terms of element_iterators.
 */
template <typename RealType>
void elem_types (const MeshBaseTempl<RealType> & mesh,
                 std::vector<ElemType> & et);

/**
 * \returns The number of elements of type \p type.
 *
 * Implemented in terms of type_element_iterators.
 */
template <typename RealType>
dof_id_type n_elem_of_type (const MeshBaseTempl<RealType> & mesh,
                            const ElemType type);

/**
 * \returns The number of active elements of type \p type.
 *
 * Implemented in terms of active_type_element_iterators.
 */
template <typename RealType>
dof_id_type n_active_elem_of_type (const MeshBaseTempl<RealType> & mesh,
                                   const ElemType type);

/**
 * \returns The number of elements of type \p type at the specified
 * refinement level.
 *
 * \todo Replace all of the n_xxx_elem() functions like this with
 * a single function which takes a range of iterators and computes the
 * std::distance between them.
 */
template <typename RealType>
dof_id_type n_non_subactive_elem_of_type_at_level(const MeshBaseTempl<RealType> & mesh,
                                                  const ElemType type,
                                                  const unsigned int level);

/**
 * \returns The number of levels of refinement in the mesh.
 *
 * Implemented by looping over all the local elements and
 * unpartitioned elements and finding the maximum level, then summing
 * in parallel.
 */
template <typename RealType>
unsigned int n_levels(const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The number of levels of refinement in the local mesh.
 *
 * Implemented by looping over all the local elements and finding the
 * maximum level.
 */
template <typename RealType>
unsigned int n_local_levels(const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The number of levels of refinement in the active mesh.
 *
 * Implemented by looping over all the active local elements and finding
 * the maximum level, then taking the max in parallel.
 */
template <typename RealType>
unsigned int n_active_levels(const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The number of levels of refinement in the active local mesh.
 *
 * Implemented by looping over all the active local elements and finding
 * the maximum level.
 */
template <typename RealType>
unsigned int n_active_local_levels(const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The number of p-levels of refinement in the mesh.
 *
 * Implemented by looping over all the local elements and finding the
 * maximum p-level, then summing in parallel.
 */
template <typename RealType>
unsigned int n_p_levels (const MeshBaseTempl<RealType> & mesh);

/**
 * \returns The number of levels of refinement in the mesh, even if that
 * mesh is not currently properly distributed or properly serialized.
 *
 * Implemented by looping over all elements and finding the maximum
 * level, then summing in parallel.  This is much slower than
 * n_levels() but will return correct values even when the mesh is in
 * an inconsistent parallel state.
 */
template <typename RealType>
unsigned int paranoid_n_levels(const MeshBaseTempl<RealType> & mesh);

/**
 * Builds a set of node IDs for nodes which belong to non-subactive
 * elements.  Non-subactive elements are those which are either active
 * or inactive.  This is useful for determining which nodes should be
 * written to a data file, and is used by the XDA mesh writing methods.
 */
template <typename RealType>
void get_not_subactive_node_ids(const MeshBaseTempl<RealType> & mesh,
                                std::set<dof_id_type> & not_subactive_node_ids);

/**
 * Count up the number of elements of a specific type
 * (as defined by an iterator range).
 */
template <typename IteratorType>
dof_id_type n_elem (const IteratorType & begin,
                    const IteratorType & end);


/**
 * Count up the number of nodes of a specific type
 * (as defined by an iterator range).
 */
template <typename IteratorType>
dof_id_type n_nodes (const IteratorType & begin,
                     const IteratorType & end);


/**
 * Find the maximum h-refinement level in a mesh.
 */
template <typename RealType>
unsigned int max_level (const MeshBaseTempl<RealType> & mesh);

/**
 * Given a mesh and a node in the mesh, the vector will be filled with
 * every node directly attached to the given one.
 */
template <typename RealType>
void find_nodal_neighbors(const MeshBaseTempl<RealType> & mesh,
                          const NodeTempl<RealType> & n,
                          const std::vector<std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map,
                          std::vector<const NodeTempl<RealType> *> & neighbors);

/**
 * Given a mesh and a node in the mesh, the vector will be filled with
 * every node directly attached to the given one.
 */
template <typename RealType>
void find_nodal_neighbors(const MeshBaseTempl<RealType> & mesh,
                          const NodeTempl<RealType> & n,
                          const std::unordered_map<dof_id_type, std::vector<const ElemTempl<RealType> *>> & nodes_to_elem_map,
                          std::vector<const NodeTempl<RealType> *> & neighbors);

/**
 * Given a mesh hanging_nodes will be filled with an associative array keyed off the
 * global id of all the hanging nodes in the mesh.  It will hold an array of the
 * parents of the node (meaning the two nodes to either side of it that make up
 * the side the hanging node is on.
 */
template <typename RealType>
void find_hanging_nodes_and_parents(const MeshBaseTempl<RealType> & mesh,
                                    std::map<dof_id_type, std::vector<dof_id_type>> & hanging_nodes);

/**
 * Changes the processor ids on each node so be the same as the id of the
 * lowest element touching that node.
 *
 * This corrects "orphaned" processor ids that may occur from element
 * coarsening.
 *
 * On a distributed mesh, this function must be called in parallel
 * to sync everyone's corrected processor ids on ghost nodes.
 */
template <typename RealType>
void correct_node_proc_ids(MeshBaseTempl<RealType> &);


#ifdef DEBUG
/**
 * A function for verifying that an element has been cut off
 * from the rest of the mesh
 */
template <typename RealType>
void libmesh_assert_no_links_to_elem(const MeshBaseTempl<RealType> & mesh,
                                     const ElemTempl<RealType> * bad_elem);

/**
 * A function for testing that all DofObjects within a mesh
 * have the same n_systems count
 */
template <typename RealType>
void libmesh_assert_equal_n_systems (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for testing that all non-recently-created DofObjects
 * within a mesh have old_dof_object data.  This is not expected to
 * be true at all points within a simulation code.
 */
template <typename RealType>
void libmesh_assert_old_dof_objects (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for walking across the mesh to try and ferret out
 * invalidated or misassigned pointers
 */
template <typename RealType>
void libmesh_assert_valid_node_pointers (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that active local elements' neighbors
 * are never remote elements
 */
template <typename RealType>
void libmesh_assert_valid_remote_elems (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that ids and processor assignment of elements
 * are correctly sorted (monotone increasing)
 */
template <typename RealType>
void libmesh_assert_valid_elem_ids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that ids of elements are correctly
 * sorted for AMR (parents have lower ids than children)
 */
template <typename RealType>
void libmesh_assert_valid_amr_elem_ids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that any interior_parent pointers on
 * elements are consistent with AMR (parents' interior_parents are
 * interior_parents' parents)
 */
template <typename RealType>
void libmesh_assert_valid_amr_interior_parents (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that all nodes are connected to at least
 * one element.
 *
 * This will fail in the most general case.  When DistributedMesh and
 * NodeConstraints are enabled, we expect the possibility that a
 * processor will be given remote nodes to satisfy node constraints
 * without also being given the remote elements connected to those
 * nodes.
 */
template <typename RealType>
void libmesh_assert_connected_nodes (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that boundary condition ids match
 * across processors.
 */
template <typename RealType>
void libmesh_assert_valid_boundary_ids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that degree of freedom indexing matches
 * across processors.
 *
 * Verify a particular system by specifying that system's number, or
 * verify all systems at once by leaving \p sysnum unspecified.
 */
template <typename RealType>
void libmesh_assert_valid_dof_ids (const MeshBaseTempl<RealType> & mesh,
                                   unsigned int sysnum = libMesh::invalid_uint);

/**
 * A function for verifying that degree of freedom indexes are
 * contiguous on each processors, as is required by libMesh numeric
 * classes.
 *
 * Verify a particular system by specifying that system's number.
 */
template <typename RealType>
void libmesh_assert_contiguous_dof_ids (const MeshBaseTempl<RealType> & mesh,
                                        unsigned int sysnum);

#ifdef LIBMESH_ENABLE_UNIQUE_ID
/**
 * A function for verifying that unique ids match across processors.
 *
 * FIXME: we ought to check for uniqueness too.
 */
template <typename RealType>
void libmesh_assert_valid_unique_ids (const MeshBaseTempl<RealType> & mesh);
#endif

/**
 * A function for verifying that distribution of dof objects is
 * parallel consistent (every processor can see every node or element
 * it owns)
 */
template <typename RealType>
void libmesh_assert_consistent_distributed(const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that distribution of nodes is parallel
 * consistent (every processor can see every node it owns) even before
 * node ids have been made consistent
 */
template <typename RealType>
void libmesh_assert_consistent_distributed_nodes(const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that processor assignment is parallel
 * consistent (every processor agrees on the processor id of each node
 * it can see) even on nodes which have not yet recieved consistent
 * DofObject::id(), using element topology to identify matching nodes.
 */
template <typename RealType>
void libmesh_assert_parallel_consistent_new_node_procids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that processor assignment is parallel
 * consistent (every processor agrees on the processor id of each dof
 * object it can see)
 */
template <typename DofObjectSubclass, typename RealType>
void libmesh_assert_parallel_consistent_procids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that processor assignment is
 * topologically consistent on nodes (each node part of an active
 * element on its processor) or elements (each parent has the
 * processor id of one of its children).
 */
template <typename DofObjectSubclass, typename RealType>
void libmesh_assert_topology_consistent_procids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that processor assignment is
 * both parallel and topologically consistent.
 */
template <typename DofObjectSubclass, typename RealType>
void libmesh_assert_valid_procids (const MeshBaseTempl<RealType> & mesh) {
  libmesh_assert_parallel_consistent_procids<DofObjectSubclass>(mesh);
  libmesh_assert_topology_consistent_procids<DofObjectSubclass>(mesh);
}

/**
 * A function for verifying that processor assignment of nodes matches
 * the heuristic specified in Node::choose_processor_id()
 */
template <typename RealType>
void libmesh_assert_canonical_node_procids (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that refinement flags on elements
 * are consistent between processors
 */
template <typename RealType>
void libmesh_assert_valid_refinement_flags (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that elements on this processor have
 * valid descendants and consistent active flags.
 */
template <typename RealType>
void libmesh_assert_valid_refinement_tree (const MeshBaseTempl<RealType> & mesh);

/**
 * A function for verifying that neighbor connectivity is correct (each
 * element is a neighbor of or descendant of a neighbor of its neighbors)
 * and consistent (each neighbor link goes to either the same neighbor
 * or to a RemoteElem on each processor)
 *
 * If assert_valid_remote_elems is set to false, then no error will be
 * thrown for neighbor links where a remote_elem should exist but a nullptr
 * exists instead.
 */
template <typename RealType>
void libmesh_assert_valid_neighbors (const MeshBaseTempl<RealType> & mesh,
                                     bool assert_valid_remote_elems=true);
#endif

// There is no reason for users to call functions in the MeshTools::Private namespace.
namespace Private {
/**
 * There is no reason for a user to ever call this function.
 *
 * This function determines partition-agnostic global indices for all
 * nodes and elements in the mesh.
 *
 * \note After this function is called, the mesh will likely be in an
 * inconsistent state, i.e. \p mesh.nodes(i)->id() != i in the nodes
 * container.  Direct node/element access via the \p mesh.node(n) or
 * \p mesh.elem(e) functions will likely fail. The original numbering
 * can (and should) be restored with a subsequent call to \p
 * fix_node_and_element_numbering().
 *
 */
template <typename RealType>
void globally_renumber_nodes_and_elements (MeshBaseTempl<RealType> &);
} // end namespace Private

} // end namespace MeshTools

} // namespace libMesh


#endif // LIBMESH_MESH_TOOLS_H
