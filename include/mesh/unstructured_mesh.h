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
 * In order to use the adaptive mesh refinment capabilities
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

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  UnstructuredMesh (unsigned char dim=1);
#endif

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
                     void * mesh_data=libmesh_nullptr,
                     bool skip_renumber_nodes_and_elements=false,
                     bool skip_find_neighbors=false) libmesh_override;
  /**
   * Write the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension.
   */
  virtual void write (const std::string & name) libmesh_override;

  /**
   * Write to the file specified by \p name.  Attempts to figure out the
   * proper method by the file extension. Also writes data.
   */
  void write (const std::string & name,
              const std::vector<Number> & values,
              const std::vector<std::string> & variable_names);

  /**
   * Converts a mesh with higher-order
   * elements into a mesh with linear elements.  For
   * example, a mesh consisting of \p Tet10 will be converted
   * to a mesh with \p Tet4 etc.
   */
  virtual void all_first_order () libmesh_override;

  /**
   * Converts a (conforming, non-refined) mesh with linear
   * elements into a mesh with second-order elements.  For
   * example, a mesh consisting of \p Tet4 will be converted
   * to a mesh with \p Tet10 etc.  Note that for some elements
   * like \p Hex8 there exist @e two higher order equivalents,
   * \p Hex20 and \p Hex27.  When \p full_ordered is \p true
   * (default), then \p Hex27 is built.  Otherwise, \p Hex20
   * is built.  The same holds obviously for \p Quad4, \p Prism6
   * ...
   */
  virtual void all_second_order (const bool full_ordered=true) libmesh_override;

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
                       const_element_iterator & it,
                       const const_element_iterator & it_end) const;


  /**
   * Deep copy of another unstructured mesh class (used by subclass
   * copy constructors)
   */
  virtual void copy_nodes_and_elements(const UnstructuredMesh & other_mesh, const bool skip_find_neighbors=false);


  /**
   * Other functions from MeshBase requiring re-definition.
   */
  virtual void find_neighbors (const bool reset_remote_elements = false,
                               const bool reset_current_list    = true) libmesh_override;

#ifdef LIBMESH_ENABLE_AMR
  /**
   * Delete subactive (i.e. children of coarsened) elements.
   * This removes all elements descended from currently active
   * elements in the mesh.
   */
  virtual bool contract () libmesh_override;
#endif // #ifdef LIBMESH_ENABLE_AMR

};


} // namespace libMesh

#endif // LIBMESH_UNSTRUCTURED_MESH_H
