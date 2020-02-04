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



#ifndef LIBMESH_BOUNDARY_INFO_H
#define LIBMESH_BOUNDARY_INFO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <cstddef>
#include <map>
#include <set>
#include <vector>
#include <tuple>

namespace libMesh
{


// Forward declarations
template <typename> class ElemTempl;
template <typename> class NodeTempl;
template <typename> class MeshBaseTempl;
template <typename> class UnstructuredMeshTempl;


/**
 * The \p BoundaryInfo class contains information relevant to boundary
 * conditions including storing faces, edges, and nodes on the
 * boundary, along with ids that can be used to identify the type of
 * boundary each entity is part of. It can also build a mesh that just
 * includes boundary elements/faces.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Used by the Mesh to keep track of boundary nodes and elements.
 */
template <typename RealType = Real>
class BoundaryInfoTempl : public ParallelObject
{
public:
  typedef BoundaryInfoTempl<RealType> BoundaryInfo;
  typedef ElemTempl<RealType> Elem;
  typedef NodeTempl<RealType> Node;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef UnstructuredMeshTempl<RealType> UnstructuredMesh;

protected:
  friend class MeshBaseTempl<RealType>;

  /**
   * Constructor.  Takes a reference to the mesh.
   * The BoundaryInfo class is only used internally
   * by the Mesh class.  A user should never instantiate
   * this class.  Therefore the constructor is protected.
   */
  BoundaryInfoTempl (MeshBase & m);

public:
  /**
   * Actual copying operation.
   *
   * \note The copy will have a reference to the same Mesh as the original.
   */
  BoundaryInfo & operator=(const BoundaryInfo & other_boundary_info);


  /**
   * Destructor.  Not much to do.
   */
  ~BoundaryInfoTempl ();

  /**
   * Clears the underlying data structures and restores the object to
   * a pristine state with no data stored.
   */
  void clear ();

  /**
   * Clears and regenerates the cached sets of ids.
   * This is necessary after use of remove_*() functions, which remove
   * individual id associations (an O(1) process) without checking to
   * see whether that is the last association with the id (an O(N)
   * process.
   */
  void regenerate_id_sets ();


  /**
   * Generates \p boundary_mesh data structures corresponding to the
   * \p mesh data structures.  Allows the \p boundary_mesh to be used
   * like any other mesh, except with interior_parent() values defined
   * for algorithms which couple boundary and interior mesh
   * information.  Any pre-existing \p boundary_mesh data is cleared.
   */
  void sync (UnstructuredMesh & boundary_mesh);

  /**
   * Generates \p boundary_mesh data structures corresponding to the
   * \p mesh data structures.  Allows the \p boundary_mesh to be used
   * like any other mesh, except with interior_parent() values defined
   * for algorithms which couple boundary and interior mesh
   * information.  Any pre-existing \p boundary_mesh data is cleared.
   * Only boundary elements with the specified ids are extracted.
   * Boundary IDs for the nodes on \p requested_boundary_ids
   * will also be copied over to \p boundary_mesh. We do not
   * currently copy edge boundary IDs over to \p boundary_mesh.
   */
  void sync (const std::set<boundary_id_type> & requested_boundary_ids,
             UnstructuredMesh & boundary_mesh);

  /**
   * Like the other sync() implementations, but specifically intended
   * for building "boundary" meshes from internal sidesets. In the
   * case of an internal sideset, each side may belong to 2
   * higher-dimensional parent elements, and typically we do not want
   * to add the same side to the boundary mesh twice. The
   * std::set<subdomain_id_type> passed into this function specifies
   * which subdomain the sides in question should relative to, so that
   * they are only added once.
   */
  void sync (const std::set<boundary_id_type> & requested_boundary_ids,
             UnstructuredMesh & boundary_mesh,
             const std::set<subdomain_id_type> & subdomains_relative_to);

  /**
   * Suppose we have used sync to create \p boundary_mesh. Then each
   * element in \p boundary_mesh will have interior_parent defined.
   * This method gets extra data for us:
   *  - \p node_id_map stores a map from the node ids on the interior mesh
   *    to the corresponding node ids of \p boundary_mesh.
   *  - \p side_id_map stores a map from the element ids of the boundary mesh
   *    to the side index of the interior_parent that the boundary element
   *    corresponds to.
   * \p tolerance is used to identify when we have matching elements.
   */
  void get_side_and_node_maps (UnstructuredMesh & boundary_mesh,
                               std::map<dof_id_type, dof_id_type> & node_id_map,
                               std::map<dof_id_type, unsigned char> & side_id_map,
                               Real tolerance=1.e-6);

  /**
   * Generates \p elements along the boundary of our _mesh, which
   * use pre-existing nodes on the boundary_mesh, and which have
   * interior_parent values properly defined.
   *
   * The \p boundary_mesh may be the *same* as the interior mesh; this
   * generates a mesh with elements of mixed dimension.
   *
   * Only boundary elements with the specified ids are created.
   */
  void add_elements (const std::set<boundary_id_type> & requested_boundary_ids,
                     UnstructuredMesh & boundary_mesh);

  /**
   * Same as the add_elements() function above, but takes a set of
   * subdomains for which the sides must be relative to. This is
   * necessary to avoid double-adding sides of internal sidesets to
   * the BoundaryMesh.
   */
  void add_elements(const std::set<boundary_id_type> & requested_boundary_ids,
                    UnstructuredMesh & boundary_mesh,
                    const std::set<subdomain_id_type> & subdomains_relative_to);

  /**
   * Add \p Node \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const Node * node,
                 const boundary_id_type id);

  /**
   * Add node number \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const dof_id_type node,
                 const boundary_id_type id);

  /**
   * Add \p Node \p node with boundary ids \p ids to the boundary
   * information data structure.
   */
  void add_node (const Node * node,
                 const std::vector<boundary_id_type> & ids);

  /**
   * Clears all the boundary information from all of the nodes in the mesh
   */
  void clear_boundary_node_ids();

  /**
   * Add edge \p edge of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const dof_id_type elem,
                 const unsigned short int edge,
                 const boundary_id_type id);

  /**
   * Add edge \p edge of element \p elem with boundary id \p id
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const Elem * elem,
                 const unsigned short int edge,
                 const boundary_id_type id);

  /**
   * Add edge \p edge of element \p elem with boundary ids \p ids
   * to the boundary information data structure.
   * Edge-based boundary IDs should only be used in 3D.
   */
  void add_edge (const Elem * elem,
                 const unsigned short int edge,
                 const std::vector<boundary_id_type> & ids);

  /**
   * Add shell face \p shellface of element number \p elem with boundary id \p id
   * to the boundary information data structure. This is only relevant for shell
   * elements.
   */
  void add_shellface (const dof_id_type elem,
                      const unsigned short int shellface,
                      const boundary_id_type id);

  /**
   * Add shell face \p shellface of element \p elem with boundary id \p id
   * to the boundary information data structure. This is only relevant for shell
   * elements.
   */
  void add_shellface (const Elem * elem,
                      const unsigned short int shellface,
                      const boundary_id_type id);

  /**
   * Add shell face \p shellface of element \p elem with boundary ids \p ids
   * to the boundary information data structure. This is only relevant for shell
   * elements.
   */
  void add_shellface (const Elem * elem,
                      const unsigned short int shellface,
                      const std::vector<boundary_id_type> & ids);

  /**
   * Add side \p side of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const dof_id_type elem,
                 const unsigned short int side,
                 const boundary_id_type id);

  /**
   * Add side \p side of element \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const Elem * elem,
                 const unsigned short int side,
                 const boundary_id_type id);

  /**
   * Add side \p side of element \p elem with boundary ids \p ids
   * to the boundary information data structure.
   */
  void add_side (const Elem * elem,
                 const unsigned short int side,
                 const std::vector<boundary_id_type> & ids);

  /**
   * Removes the boundary conditions associated with node \p node,
   * if any exist.
   */
  void remove (const Node * node);

  /**
   * Removes the boundary conditions associated with element \p elem,
   * if any exist.
   */
  void remove (const Elem * elem);

  /**
   * Removes boundary id \p id from node \p node, if it exists.
   */
  void remove_node (const Node * node,
                    const boundary_id_type id);

  /**
   * Removes all boundary conditions associated with edge \p edge of
   * element \p elem, if any exist.
   */
  void remove_edge (const Elem * elem,
                    const unsigned short int edge);

  /**
   * Removes the boundary id \p id from edge \p edge of element \p
   * elem, if it exists.
   */
  void remove_edge (const Elem * elem,
                    const unsigned short int edge,
                    const boundary_id_type id);

  /**
   * Removes all boundary conditions associated with shell face
   * \p shellface of element \p elem, if any exist.
   */
  void remove_shellface (const Elem * elem,
                         const unsigned short int shellface);

  /**
   * Removes all boundary conditions associated with shell face
   * \p shellface of element \p elem, if any exist.
   */
  void remove_shellface (const Elem * elem,
                         const unsigned short int shellface,
                         const boundary_id_type id);

  /**
   * Removes all boundary conditions associated with side \p side of
   * element \p elem, if any exist.
   */
  void remove_side (const Elem * elem,
                    const unsigned short int side);

  /**
   * Removes the boundary id \p id from side \p side of element \p
   * elem, if it exists.
   */
  void remove_side (const Elem * elem,
                    const unsigned short int side,
                    const boundary_id_type id);

  /**
   * Removes all entities (nodes, sides, edges, shellfaces) with boundary
   * id \p id from their respective containers and erases any record of
   * \p id's existence from the BoundaryInfo object.  That is, after
   * calling remove_id(), \p id will no longer be in the sets returned by
   * get_boundary_ids(), get_side_boundary_ids(), etc., and will not
   * be in the bc_id_list vector returned by build_side_list(), etc.
   */
  void remove_id (boundary_id_type id);

  /**
   * \returns The number of user-specified boundary ids on the
   * semilocal part of the mesh.
   *
   * \note DistributedMesh users may need to compare boundary_ids sets
   * via inter-processor communication.
   */
  std::size_t n_boundary_ids () const { return _boundary_ids.size(); }

  /**
   * \returns \p true if \p node is associated with boundary \p id.
   */
  bool has_boundary_id (const Node * const node,
                        const boundary_id_type id) const;

  /**
   * \returns The boundary ids associated with \p Node \p node.
   *
   * \deprecated Instead, use the version of this function that fills
   * a std::vector.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> boundary_ids (const Node * node) const;
#endif

  /**
   * Fills a user-provided std::vector with the boundary ids associated
   * with \p Node \p node.
   *
   * This is the non-deprecated version of the function.
   */
  void boundary_ids (const Node * node,
                     std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The number of boundary ids associated with \p Node \p node.
   */
  unsigned int n_boundary_ids (const Node * node) const;

  /**
   * \returns The number of boundary ids associated with the \p edge
   * edge of element \p elem.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   */
  unsigned int n_edge_boundary_ids (const Elem * const elem,
                                    const unsigned short int edge) const;

  /**
   * \returns The list of boundary ids associated with the \p edge edge of
   * element \p elem.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   *
   * \deprecated Instead, use the version of this function that fills
   * a std::vector.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> edge_boundary_ids (const Elem * const elem,
                                                   const unsigned short int edge) const;
#endif

  /**
   * \returns The list of boundary ids associated with the \p edge edge of
   * element \p elem.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   *
   * This is the non-deprecated version of the function.
   */
  void edge_boundary_ids (const Elem * const elem,
                          const unsigned short int edge,
                          std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The list of raw boundary ids associated with the \p edge
   * edge of element \p elem.
   *
   * These ids are "raw" because they exclude ids which are implicit,
   * such as a child's inheritance of its ancestors' boundary id.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   *
   * \deprecated Instead, use the version of this function that fills
   * a std::vector.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> raw_edge_boundary_ids (const Elem * const elem,
                                                       const unsigned short int edge) const;
#endif

  /**
   * \returns The list of raw boundary ids associated with the \p edge
   * edge of element \p elem.
   *
   * These ids are "raw" because they exclude ids which are implicit,
   * such as a child's inheritance of its ancestors' boundary id.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   *
   * This is the non-deprecated version of the function.
   */
  void raw_edge_boundary_ids (const Elem * const elem,
                              const unsigned short int edge,
                              std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The number of boundary ids associated with the specified
   * shell face of element \p elem.
   *
   * \note This is only relevant for shell elements.
   */
  unsigned int n_shellface_boundary_ids (const Elem * const elem,
                                         const unsigned short int shellface) const;

  /**
   * \returns The list of boundary ids associated with the specified shell face
   * of element \p elem.
   *
   * \note This is only relevant for shell elements.
   */
  void shellface_boundary_ids (const Elem * const elem,
                               const unsigned short int shellface,
                               std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The list of raw boundary ids associated with the
   * specified shell face of element \p elem.
   *
   * These ids are "raw" because they exclude ids which are implicit,
   * such as a child's inheritance of its ancestors' boundary id.
   *
   * \note This is only relevant for shell elements.
   */
  void raw_shellface_boundary_ids (const Elem * const elem,
                                   const unsigned short int shellface,
                                   std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns \p true if side \p side of Elem \p elem is associated
   * with \p id.
   */
  bool has_boundary_id (const Elem * const elem,
                        const unsigned short int side,
                        const boundary_id_type id) const;

  /**
   * \returns The boundary id associated with the \p side side of
   * element \p elem, or \p invalid_id if the \p side does not have an
   * associated boundary id.
   *
   * \note Only one id per side is allowed, however multiple sides per
   * element are allowed.
   *
   * \deprecated Asking for just one boundary id means your code isn't
   * safe to use on meshes with overlapping boundary ids.  Try using
   * BoundaryInfo::boundary_ids() or BoundaryInfo::has_boundary_id()
   * instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  boundary_id_type boundary_id (const Elem * const elem,
                                const unsigned short int side) const;
#endif

  /**
   * \returns The number of boundary ids associated with the \p side
   * side of element \p elem.
   */
  unsigned int n_boundary_ids (const Elem * const elem,
                               const unsigned short int side) const;

  /**
   * \returns The list of boundary ids associated with the \p side side of
   * element \p elem.
   *
   * \deprecated Instead, use the version of this function that fills
   * a std::vector.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> boundary_ids (const Elem * const elem,
                                              const unsigned short int side) const;
#endif

  /**
   * \returns The list of boundary ids associated with the \p side side of
   * element \p elem.
   *
   * This is the non-deprecated version of the function.
   */
  void boundary_ids (const Elem * const elem,
                     const unsigned short int side,
                     std::vector<boundary_id_type> & vec_to_fill) const;

  /**
   * \returns The list of raw boundary ids associated with the \p side
   * side of element \p elem.
   *
   * These ids are "raw" because they exclude ids which are implicit,
   * such as a child's inheritance of its ancestors' boundary id.
   *
   * \deprecated Instead, use the version of this function that fills
   * a std::vector.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  std::vector<boundary_id_type> raw_boundary_ids (const Elem * const elem,
                                                  const unsigned short int side) const;
#endif

  /**
   * \returns The list of raw boundary ids associated with the \p side
   * side of element \p elem.
   *
   * These ids are "raw" because they exclude ids which are implicit,
   * such as a child's inheritance of its ancestors' boundary id.
   *
   * This is the non-deprecated version of the function.
   */
  void raw_boundary_ids (const Elem * const elem,
                         const unsigned short int side,
                         std::vector<boundary_id_type> & vec_to_fill) const;

  /*
   * Copy boundary ids associated with old_elem (but not its nodes)
   * from old_boundary_info (which may be this) into this boundary
   * info, associating them with new_elem.
   */
  void copy_boundary_ids (const BoundaryInfo & old_boundary_info,
                          const Elem * const old_elem,
                          const Elem * const new_elem);

  /**
   * \returns A side of element \p elem whose associated boundary id is
   * \p boundary_id if such a side exists, and \p invalid_uint otherwise.
   *
   * \note If multiple sides of \p elem have the same id, only the lowest numbered
   * such side is returned.
   */
  unsigned int side_with_boundary_id(const Elem * const elem,
                                     const boundary_id_type boundary_id) const;

  /**
   * Builds the list of unique node boundary ids.
   *
   * On a ReplicatedMesh this will be all ids; on a DistributedMesh
   * only ids on semilocal nodes will be included.
   */
  void build_node_boundary_ids(std::vector<boundary_id_type> & b_ids) const;

  /**
   * Builds the list of unique side boundary ids.
   *
   * On a ReplicatedMesh this will be all ids; on a DistributedMesh
   * only ids on sides of semilocal elements will be included.
   */
  void build_side_boundary_ids(std::vector<boundary_id_type> & b_ids) const;

  /**
   * Builds the list of unique shellface boundary ids.
   *
   * On a ReplicatedMesh this will be all ids; on a DistributedMesh
   * only ids on shellfaces of semilocal elements will be included.
   */
  void build_shellface_boundary_ids(std::vector<boundary_id_type> & b_ids) const;

  /**
   * \returns The number of element-side-based boundary conditions.
   *
   * This will be the correct global count even on a distributed mesh.
   */
  std::size_t n_boundary_conds () const;

  /**
   * \returns The number of edge-based boundary conditions.
   * Edge-based boundary IDs should only be used in 3D.
   *
   * This will be the correct global count even on a distributed mesh.
   */
  std::size_t n_edge_conds () const;

  /**
   * \returns The number of shellface-based boundary conditions.
   * This is only relevant on shell elements.
   *
   * This will be the correct global count even on a distributed mesh.
   */
  std::size_t n_shellface_conds () const;

  /**
   * \returns The number of node-based boundary conditions.
   *
   * This will be the correct global count even on a distributed mesh.
   */
  std::size_t n_nodeset_conds () const;

  /**
   * Creates a list of nodes and ids for those nodes.
   *
   * On a ReplicatedMesh this will include all nodes; on a
   * DistributedMesh only semilocal nodes will be included.
   *
   * \deprecated Use the version of build_node_list() below that
   * returns a std::vector of tuples instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void build_node_list (std::vector<dof_id_type> &      node_id_list,
                        std::vector<boundary_id_type> & bc_id_list) const;
#endif

  /**
   * As above, but the library creates and fills in a vector of
   * (node-id, bc-id) pairs and returns it to the user, taking
   * advantage of guaranteed RVO. Note: we could use std::pairs for
   * this, but for consistency with the other build_XYZ_list
   * functions, we're using tuples.
   */
  std::vector<std::tuple<dof_id_type, boundary_id_type>>
  build_node_list() const;

  /**
   * Adds nodes with boundary ids based on the side's boundary
   * ids they are connected to.
   */
  void build_node_list_from_side_list();

  /**
   * Adds sides to a sideset if every node on that side are in the same
   * sideset
   */
  void build_side_list_from_node_list();

  /**
   * Creates a list of element numbers, sides, and ids for those sides.
   *
   * On a ReplicatedMesh this will include all sides; on a
   * DistributedMesh only sides of semilocal elements will be
   * included.
   *
   * \deprecated Use the version of build_side_list() below that
   * returns a std::vector of tuples instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void build_side_list (std::vector<dof_id_type> &        element_id_list,
                        std::vector<unsigned short int> & side_list,
                        std::vector<boundary_id_type> &   bc_id_list) const;
#endif

  /**
   * As above, but the library creates and fills in a vector of
   * (elem-id, side-id, bc-id) triplets and returns it to the user,
   * taking advantage of guaranteed RVO.
   */
  typedef std::tuple<dof_id_type, unsigned short int, boundary_id_type> BCTuple;
  std::vector<BCTuple> build_side_list() const;

  /**
   * Creates a list of active element numbers, sides, and ids for those sides.
   *
   * On a ReplicatedMesh this will include all sides; on a
   * DistributedMesh only sides of semilocal elements will be
   * included.
   *
   * \deprecated Use the version of build_active_side_list() below
   * that returns a std::vector of tuples instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void build_active_side_list (std::vector<dof_id_type> &        element_id_list,
                               std::vector<unsigned short int> & side_list,
                               std::vector<boundary_id_type> &   bc_id_list) const;
#endif

  /**
   * As above, but the library creates and fills in a vector of
   * (elem-id, side-id, bc-id) triplets and returns it to the user,
   * taking advantage of guaranteed RVO.
   */
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
  build_active_side_list () const;

  /**
   * Creates a list of element numbers, edges, and boundary ids for those edges.
   *
   * On a ReplicatedMesh this will include all edges; on a
   * DistributedMesh only edges of semilocal elements will be
   * included.
   *
   * \deprecated Use the version of build_edge_list() below
   * that returns a std::vector of tuples instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void build_edge_list (std::vector<dof_id_type> &        element_id_list,
                        std::vector<unsigned short int> & edge_list,
                        std::vector<boundary_id_type> &   bc_id_list) const;
#endif

  /**
   * As above, but the library creates and fills in a vector of
   * (elem-id, side-id, bc-id) triplets and returns it to the user,
   * taking advantage of guaranteed RVO.
   */
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
  build_edge_list() const;

  /**
   * Creates a list of element numbers, shellfaces, and boundary ids for those shellfaces.
   *
   * On a ReplicatedMesh this will include all shellfaces; on a
   * DistributedMesh only shellfaces of semilocal elements will be
   * included.
   *
   * \deprecated Use the version of build_shellface_list() below
   * that returns a std::vector of tuples instead.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  void build_shellface_list (std::vector<dof_id_type> &        element_id_list,
                             std::vector<unsigned short int> & shellface_list,
                             std::vector<boundary_id_type> &   bc_id_list) const;
#endif

  /**
   * As above, but the library creates and fills in a vector of
   * (elem-id, side-id, bc-id) triplets and returns it to the user,
   * taking advantage of guaranteed RVO.
   */
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>>
  build_shellface_list() const;

  /**
   * \returns A set of the boundary ids which exist on semilocal parts
   * of the mesh.
   *
   * DistributedMesh-compatible code may need a set_union or other
   * manipulations to work with sets of boundary ids which include ids
   * on remote parts of the mesh.
   */
  const std::set<boundary_id_type> & get_boundary_ids () const
  { return _boundary_ids; }

  /**
   * \returns A reference to the set of the boundary IDs specified on
   * sides of semilocal mesh elements.
   */
  const std::set<boundary_id_type> & get_side_boundary_ids () const
  { return _side_boundary_ids; }

  /**
   * \returns A reference to the set of all boundary IDs specified on
   * edges of semilocal mesh elements.
   *
   * \note Edge-based boundary IDs should only be used in 3D.
   */
  const std::set<boundary_id_type> & get_edge_boundary_ids () const
  { return _edge_boundary_ids; }

  /**
   * \returns A reference to the set of all boundary IDs specified on
   * shell faces.
   *
   * \note This is only relevant on shell elements.
   */
  const std::set<boundary_id_type> & get_shellface_boundary_ids () const
  { return _shellface_boundary_ids; }

  /**
   * \returns A reference to the set of all boundary IDs specified on
   * semilocal mesh nodes.
   */
  const std::set<boundary_id_type> & get_node_boundary_ids () const
  { return _node_boundary_ids; }


  /**
   * Prints the boundary information data structure.
   */
  void print_info (std::ostream & out=libMesh::out) const;

  /**
   * Prints a summary of the boundary information.
   */
  void print_summary (std::ostream & out=libMesh::out) const;

  /**
   * \returns A reference for getting an optional name for a sideset.
   */
  const std::string & get_sideset_name(boundary_id_type id) const;

  /**
   * \returns A writable reference for setting an optional
   * name for a sideset.
   */
  std::string & sideset_name(boundary_id_type id);

  /**
   * \returns A reference for getting an optional name for a nodeset.
   */
  const std::string & get_nodeset_name(boundary_id_type id) const;

  /**
   * \returns A writable reference for setting an optional
   * name for a nodeset.
   */
  std::string & nodeset_name(boundary_id_type id);

  /**
   * \returns A const reference to an optional edgeset name.
   */
  const std::string & get_edgeset_name(boundary_id_type id) const;

  /**
   * \returns A writable reference to an optional edgeset name.
   */
  std::string & edgeset_name(boundary_id_type id);

  /**
   * \returns The id of the named boundary if it exists, \p invalid_id
   * otherwise.
   */
  boundary_id_type get_id_by_name(const std::string & name) const;

  /**
   * \returns A writable reference to the sideset name map.
   */
  std::map<boundary_id_type, std::string> & set_sideset_name_map ()
  { return _ss_id_to_name; }
  const std::map<boundary_id_type, std::string> & get_sideset_name_map () const
  { return _ss_id_to_name; }

  /**
   * \returns A writable reference to the nodeset name map.
   */
  std::map<boundary_id_type, std::string> & set_nodeset_name_map ()
  { return _ns_id_to_name; }
  const std::map<boundary_id_type, std::string> & get_nodeset_name_map () const
  { return _ns_id_to_name; }

  /**
   * \returns Writable/const reference to the edgeset name map.
   */
  std::map<boundary_id_type, std::string> & set_edgeset_name_map ()
  { return _es_id_to_name; }
  const std::map<boundary_id_type, std::string> & get_edgeset_name_map () const
  { return _es_id_to_name; }

  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const boundary_id_type invalid_id;


private:

  /**
   * Helper method for finding consistent maps of interior to boundary
   * dof_object ids.  Either node_id_map or side_id_map can be nullptr,
   * in which case it will not be filled.
   */
  void _find_id_maps (const std::set<boundary_id_type> & requested_boundary_ids,
                      dof_id_type first_free_node_id,
                      std::map<dof_id_type, dof_id_type> * node_id_map,
                      dof_id_type first_free_elem_id,
                      std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> * side_id_map,
                      const std::set<subdomain_id_type> & subdomains_relative_to);

  /**
   * The Mesh this boundary info pertains to.
   */
  MeshBase & _mesh;

  /**
   * Data structure that maps nodes in the mesh
   * to boundary ids.
   */
  std::multimap<const Node *,
                boundary_id_type> _boundary_node_id;

  /**
   * Data structure that maps edges of elements
   * to boundary ids. This is only relevant in 3D.
   */
  std::multimap<const Elem *,
                std::pair<unsigned short int, boundary_id_type>>
  _boundary_edge_id;

  /**
   * Data structure that maps faces of shell elements
   * to boundary ids. This is only relevant for shell elements.
   */
  std::multimap<const Elem *,
                std::pair<unsigned short int, boundary_id_type>>
  _boundary_shellface_id;

  /**
   * Data structure that maps sides of elements
   * to boundary ids.
   */
  std::multimap<const Elem *,
                std::pair<unsigned short int, boundary_id_type>>
  _boundary_side_id;

  /**
   * A collection of user-specified boundary ids for sides, edges, nodes,
   * and shell faces.
   * See _side_boundary_ids, _edge_boundary_ids, _node_boundary_ids, and
   * _shellface_boundary_ids for sets containing IDs for only sides, edges,
   * nodes, and shell faces, respectively.
   */
  std::set<boundary_id_type> _boundary_ids;

  /**
   * Set of user-specified boundary IDs for sides *only*.
   *
   * \note \p _boundary_ids is the union of this set, \p
   * _edge_boundary_ids, \p _node_boundary_ids, and \p
   * _shellface_boundary_ids.
   */
  std::set<boundary_id_type> _side_boundary_ids;

  /**
   * Set of user-specified boundary IDs for edges *only*.
   * This is only relevant in 3D.
   *
   * \note \p _boundary_ids is the union of this set, \p _side_boundary_ids,
   * \p _node_boundary_ids, and \p _shellface_boundary_ids.
   */
  std::set<boundary_id_type> _edge_boundary_ids;

  /**
   * Set of user-specified boundary IDs for nodes *only*.
   *
   * \note \p _boundary_ids is the union of this set, \p
   * _edge_boundary_ids, \p _side_boundary_ids, and \p
   * _shellface_boundary_ids.
   */
  std::set<boundary_id_type> _node_boundary_ids;

  /**
   * Set of user-specified boundary IDs for shellfaces *only*.
   * This is only relevant for shell elements.
   *
   * \note \p _boundary_ids is the union of this set, \p
   * _side_boundary_ids, \p _edge_boundary_ids, and \p
   * _node_boundary_ids.
   */
  std::set<boundary_id_type> _shellface_boundary_ids;

  /**
   * This structure maintains the mapping of named side sets
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<boundary_id_type, std::string> _ss_id_to_name;

  /**
   * This structure maintains the mapping of named node sets
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<boundary_id_type, std::string> _ns_id_to_name;

  /**
   * This structure maintains the mapping of named edge sets
   * for file formats that support named blocks.  Currently
   * this is only implemented for ExodusII
   */
  std::map<boundary_id_type, std::string> _es_id_to_name;
};

//------------------------------------------------------
// BoundaryInfo static member initializations
template <typename RealType>
const boundary_id_type BoundaryInfoTempl<RealType>::invalid_id = -123;

typedef BoundaryInfoTempl<Real> BoundaryInfo;

} // namespace libMesh

#endif // LIBMESH_BOUNDARY_INFO_H
