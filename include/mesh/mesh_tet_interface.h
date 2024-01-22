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
#include "libmesh/enum_elem_type.h"
#include "libmesh/point.h" // used for specifying holes

// C++ includes
#include <cstddef>
#include <map>
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
   * Empty destructor.
   */
  ~MeshTetInterface() = default;

  /**
   * Sets and/or gets the desired element type.
   */
  ElemType & elem_type() {return _elem_type;}

  /**
   * This is the main public interface for this function.
   */
  virtual void triangulate () = 0;

protected:
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
   * The exact type of tetrahedra we intend to construct
   */
  ElemType _elem_type;

  /**
   * Local reference to the mesh we are working with.
   */
  UnstructuredMesh & _mesh;
};

} // namespace libMesh

#endif // LIBMESH_MESH_TET_INTERFACE_H
