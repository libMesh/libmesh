// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_TET_INTERFACE_H
#define LIBMESH_MESH_TET_INTERFACE_H

#include "libmesh/libmesh_config.h"

// Local includes
#include "libmesh/bounding_box.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/point.h" // used for specifying holes

// C++ includes
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace libMesh
{
// Forward Declarations
class UnstructuredMesh;
class Elem;

/**
 * Class \p MeshTetInterface provides an abstract interface for
 * tetrahedralization of meshes by subclasses.
 *
 * \author Roy H. Stogner
 * \date 2024
 */
class MeshTetInterface
{
public:

  /**
   * Constructor. Takes a reference to the mesh.
   */
  explicit
  MeshTetInterface (UnstructuredMesh & mesh);

  /**
   * Default destructor in base class.
   */
  virtual ~MeshTetInterface();

  /**
   * Sets and/or gets the desired tetrahedron volume.  Set to zero to
   * disable volume constraint.
   */
  Real & desired_volume() {return _desired_volume;}

  /**
   * Sets/gets flag which tells whether to do two steps of Laplace
   * mesh smoothing after generating the grid.  False by default (for
   * compatibility with old TetGenMeshInterface behavior).
   */
  bool & smooth_after_generating() {return _smooth_after_generating;}

  /**
   * Sets and/or gets the desired element type.  This should be a Tet
   * type.
   */
  ElemType & elem_type() {return _elem_type;}

  /**
   * Attaches a vector of Mesh pointers defining holes which will be
   * meshed around.  We use unique_ptr here because we expect that we
   * may need to modify these meshes internally.
   */
  void attach_hole_list(std::unique_ptr<std::vector<std::unique_ptr<UnstructuredMesh>>> holes);

  /**
   * This is the main public interface for this function.
   */
  virtual void triangulate () = 0;

protected:

  /**
   * Remove volume elements from the given mesh, after converting
   * their outer boundary faces to surface elements.
   *
   * Returns the bounding box of the mesh; this is useful for
   * detecting misplaced holes later.
   */
  static BoundingBox volume_to_surface_mesh (UnstructuredMesh & mesh);

  /**
   * This function checks the integrity of the current set of
   * elements in the Mesh to see if they comprise a hull,
   * that is:
   * - If they are all TRI3 elements
   * - They all have non-nullptr neighbors
   *
   * \returns
   * - 0 if the mesh forms a valid hull
   * - 1 if a non-TRI3 element is found
   * - 2 if an element with a nullptr-neighbor is found
   */
  unsigned check_hull_integrity();

  /**
   * This function prints an informative message and
   * crashes based on the output of the check_hull_integrity()
   * function.  It is a separate function so that you
   * can check hull integrity without crashing if you desire.
   */
  void process_hull_integrity_result(unsigned result);

  /**
   * Delete original convex hull elements from the Mesh
   * after performing a Delaunay tetrahedralization.
   */
  void delete_2D_hull_elements();

  /**
   * The desired volume for the elements in the resulting mesh.
   * Unlimited (indicated by 0) by default
   */
  Real _desired_volume;

  /**
   * Flag which tells whether we should smooth the mesh after
   * it is generated.  False by default.
   */
  bool _smooth_after_generating;

  /**
   * The exact type of tetrahedra we intend to construct
   */
  ElemType _elem_type;

  /**
   * Local reference to the mesh we are working with.
   */
  UnstructuredMesh & _mesh;

  /**
   * A pointer to a vector of meshes each defining a hole.  If this is
   * nullptr, there are no holes!
   */
  std::unique_ptr<std::vector<std::unique_ptr<UnstructuredMesh>>> _holes;
};

} // namespace libMesh

#endif // LIBMESH_MESH_TET_INTERFACE_H
