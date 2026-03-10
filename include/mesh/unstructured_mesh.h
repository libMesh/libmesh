// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_UNSTRUCTURED_MESH_H
#define LIBMESH_UNSTRUCTURED_MESH_H

// Local Includes
#include "libmesh/mesh_base.h"

// C++ Includes
#include <cstddef>

namespace libMesh
{

/**
 * The \p UnstructuredMesh class is derived from the \p MeshBase class.  The
 * user will typically want to instantiate and use the
 * Mesh class in her applications, which is currently a simple
 * derived class of UnstructuredMesh.
 * In order to use the adaptive mesh refinement capabilities
 * of the library, first instantiate a MeshRefinement object
 * with a reference to this class.  Then call the appropriate
 * refinement functions from that object.  To interact with the
 * boundary, instantiate a BoundaryMesh with a reference to
 * this class, and then use that object's functionality.
 *
 * \author Roy Stogner
 * \date 2007
 * \brief Base class for Replicated and Distributed meshes.
 */
class UnstructuredMesh : public MeshBase
{
public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  UnstructuredMesh (const Parallel::Communicator & comm_in,
                    unsigned char dim=1);

  /**
   * UnstructuredMesh uses a defaulted copy constructor.
   */
  UnstructuredMesh(const UnstructuredMesh &) = default;

  /**
   * UnstructuredMesh constructor from arbitrary (e.g. Cartesian
   * someday?) meshes
   */
  UnstructuredMesh(const MeshBase &);

  /**
   * Move-constructor deleted in MeshBase.
   */
  UnstructuredMesh(UnstructuredMesh &&) = delete;

  /**
   * Copy assignment is not allowed.
   */
  UnstructuredMesh & operator= (const UnstructuredMesh &) = delete;

  /**
   * Move assignment is allowed, by subclasses who handle
   * post_dofobject_moves()
   */
  UnstructuredMesh & operator= (UnstructuredMesh && other_mesh) = default;

  virtual MeshBase & assign(MeshBase && other_mesh) override = 0;

  /**
   * Destructor.
   */
  virtual ~UnstructuredMesh();

  /**
   * Reads the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.  This is now the only
   * way to read a mesh.  The \p UnstructuredMesh then initializes its data
   * structures and is ready for use.
   *
   * The skip_renumber_nodes_and_elements argument is now deprecated -
   * to disallow renumbering, set \p MeshBase::allow_renumbering(false).
   *
   * Set skip_find_neighbors=true to skip the find-neighbors operation
   * during prepare_for_use. This operation isn't always necessary
   * and it can be time-consuming, which is why we provide an option to
   * skip it.
   */
  virtual void read (const std::string & name,
                     void * mesh_data=nullptr,
                     bool skip_renumber_nodes_and_elements=false,
                     bool skip_find_neighbors=false) override;
  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   */
  virtual void write (const std::string & name) const override;

  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write (const std::string & name,
              const std::vector<Number> & values,
              const std::vector<std::string> & variable_names) const;

  /**
   * Converts a mesh with higher-order
   * elements into a mesh with linear elements.  For
   * example, a mesh consisting of \p Tet10 will be converted
   * to a mesh with \p Tet4 etc.
   */
  virtual void all_first_order () override;

  /**
   * Converts a (conforming, non-refined) mesh with linear elements
   * into a mesh with second-order elements.  For example, a mesh
   * consisting of \p Tet4 will be converted to a mesh with \p Tet10
   * etc.
   *
   * \note For some elements like \p Hex8 there exist two higher order
   * equivalents, \p Hex20 and \p Hex27.  When \p full_ordered is \p
   * true (default), then \p Hex27 is built.  Otherwise, \p Hex20 is
   * built.  The same holds obviously for \p Quad4, \p Prism6, etc.
   */
  virtual void all_second_order_range (const SimpleRange<element_iterator> & range,
                                       const bool full_ordered=true) override;

  /**
   * Converts a (conforming, non-refined) mesh with linear elements
   * into a mesh with "complete" order elements, i.e. elements which
   * can store degrees of freedom on any vertex, edge, or face.  For
   * example, a mesh consisting of \p Tet4 or \p Tet10 will be
   * converted to a mesh with \p Tet14 etc.
   */
  virtual void all_complete_order_range (const SimpleRange<element_iterator> & range) override;

  /**
   * Generates a new mesh containing all the elements which
   * are assigned to processor \p pid.  This mesh is written
   * to the pid_mesh reference which you must create and pass
   * to the function.
   */
  void create_pid_mesh (UnstructuredMesh & pid_mesh,
                        const processor_id_type pid) const;

  /**
   * Constructs a mesh called "new_mesh" from the current mesh by
   * iterating over the elements between it and it_end and adding
   * them to the new mesh.
   */
  void create_submesh (UnstructuredMesh & new_mesh,
                       const const_element_iterator & it,
                       const const_element_iterator & it_end) const;

  /**
   * Stitch \p other_mesh to this mesh so that this mesh is the union of the two meshes.
   * \p this_mesh_boundary and \p other_mesh_boundary are used to specify a dim-1 dimensional
   * surface on which we seek to merge any "overlapping" nodes, where we use the parameter
   * \p tol as a relative tolerance (relative to the smallest edge length on the surfaces
   * being stitched) to determine whether or not nodes are overlapping.
   * If \p clear_stitched_boundary_ids==true, this function clears boundary_info IDs in this
   * mesh associated \p this_mesh_boundary and \p other_mesh_boundary.
   * If \p use_binary_search is true, we use an optimized "sort then binary search" algorithm
   * for finding matching nodes. Otherwise we use a N^2 algorithm (which can be more reliable
   * at dealing with slightly misaligned meshes).
   * If \p enforce_all_nodes_match_on_boundaries is true, we throw an error if the number of
   * nodes on the specified boundaries don't match the number of nodes that were merged.
   * This is a helpful error check in some cases. If this is true, it overrides the value of
   * \p merge_boundary_nodes_all_or_nothing.
   * If \p skip_find_neighbors is true, a faster stitching method is used, where the lists of
   * neighbors for each elements are copied as well and patched, without calling the time-consuming
   * find_neighbors() function. This option is now hard-coded to true.
   * If \p merge_boundary_nodes_all_or_nothing is true, instead of throwing an error
   * like \p enforce_all_nodes_match_on_boundaries, the meshes are combined anyway but coincident
   * nodes are not merged into single nodes. This is useful in cases where you are not sure if the
   * boundaries are fully conforming beforehand and you want to handle the non-conforming cases
   * differently.
   *
   * Note that the element IDs for elements in the stitched mesh corresponding to "this" mesh
   * will be unchanged. The IDs for elements corresponding to \p other_mesh will be incremented
   * by this->max_elem_id().
   *
   * There is no simple a priori relationship between node IDs in "this" mesh
   * and other_mesh and node IDs in the stitched mesh because the number of nodes (and hence
   * the node IDs) in the stitched mesh depend on how many nodes are stitched.
   *
   * If \p remap_subdomain_ids is true then we assume that some
   * subdomain ids might have been autogenerated, so we remap them as
   * necessary, treating subdomain names as the important thing for
   * consistency; if we have missing names and cannot infer a
   * consistent resolution to an id conflict then we exit with an
   * error.  If \p remap_subdomain_ids is false then we revert to the
   * older libMesh behavior: leave all subdomain ids alone and woe
   * unto you if you weren't keeping track of them.
   *
   * If \p prepare_after_stitching is true then we prepare the newly
   * stitched mesh for use immediately after stitching.  If the mesh
   * does not need to be completely prepared for use yet (e.g. because
   * it will undergo further stitching etc. before repartitioning for
   * load balancing is desired), this can be set to false to
   * potentially improve performance.
   *
   * \returns the count of how many nodes were merged between the two meshes.
   * This can be zero in the case of no matching nodes or if
   * \p merge_boundary_nodes_all_or_nothing was active and relevant.
   */
  std::size_t stitch_meshes (const MeshBase & other_mesh,
                             boundary_id_type this_mesh_boundary,
                             boundary_id_type other_mesh_boundary,
                             Real tol=TOLERANCE,
                             bool clear_stitched_boundary_ids=false,
                             bool verbose=true,
                             bool use_binary_search=true,
                             bool enforce_all_nodes_match_on_boundaries=false,
                             bool merge_boundary_nodes_all_or_nothing=false,
                             bool remap_subdomain_ids=false,
                             bool prepare_after_stitching=true);

  /**
   * Similar to stitch_meshes, except that we stitch two adjacent surfaces within this mesh.
   */
  std::size_t stitch_surfaces (boundary_id_type boundary_id_1,
                               boundary_id_type boundary_id_2,
                               Real tol=TOLERANCE,
                               bool clear_stitched_boundary_ids=false,
                               bool verbose=true,
                               bool use_binary_search=true,
                               bool enforce_all_nodes_match_on_boundaries=false,
                               bool merge_boundary_nodes_all_or_nothing=false,
                               bool prepare_after_stitching=true);

  /**
   * Deep copy of nodes and elements from another mesh object (used by
   * subclass copy constructors and by mesh merging operations)
   *
   * This will not copy most "high level" data in the mesh; that is
   * done separately by constructors.  An exception is that, if the
   * \p other_mesh has element or node extra_integer data, any names
   * for that data which do not already exist on \p this mesh are
   * added so that all such data can be copied.
   *
   * If an \p id_remapping map is provided, then element subdomain ids
   * in \p other_mesh will be converted using it before adding them to
   * \p this mesh.
   *
   * For backwards compatibility, this does some limited mesh
   * preparation after the copy: everything except for renumbering,
   * remote element removal, and partitioning.  To skip just the
   * step of that preparation which finds new neighbor_ptr links
   * between elements, set \p skip_find_neighbors.  To skip all of
   * that preparation, set \p skip_preparation.  If preparation is
   * skipped, it is the users responsibility to set the flags
   * indicating what preparation may still be necessary before using
   * the mesh later.
   */
  virtual void copy_nodes_and_elements (const MeshBase & other_mesh,
                                        const bool skip_find_neighbors = false,
                                        dof_id_type element_id_offset = 0,
                                        dof_id_type node_id_offset = 0,
                                        unique_id_type unique_id_offset = 0,
                                        std::unordered_map<subdomain_id_type, subdomain_id_type> *
                                          id_remapping = nullptr,
                                        const bool skip_preparation = false);

  /**
   * Move node and elements from other_mesh to this mesh.
   */
  virtual void move_nodes_and_elements(MeshBase && other_mesh) = 0;


  /**
   * Other functions from MeshBase requiring re-definition.
   */
  virtual void find_neighbors (const bool reset_remote_elements = false,
                               const bool reset_current_list    = true,
                               const bool assert_valid          = true,
                               const bool check_non_remote      = true) override;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract () override;
#endif // #ifdef LIBMESH_ENABLE_AMR

private:

  /**
   * Helper function for stitch_meshes and stitch_surfaces
   * that does the mesh stitching.
   */
  std::size_t stitching_helper (const MeshBase * other_mesh,
                                boundary_id_type boundary_id_1,
                                boundary_id_type boundary_id_2,
                                Real tol,
                                bool clear_stitched_boundary_ids,
                                bool verbose,
                                bool use_binary_search,
                                bool enforce_all_nodes_match_on_boundaries,
                                bool skip_find_neighbors,
                                bool merge_boundary_nodes_all_or_nothing,
                                bool remap_subdomain_ids,
                                bool prepare_after_stitching);
};


} // namespace libMesh

#endif // LIBMESH_UNSTRUCTURED_MESH_H
