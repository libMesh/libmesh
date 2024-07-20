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


#ifndef LIBMESH_MESH_NETGEN_INTERFACE_H
#define LIBMESH_MESH_NETGEN_INTERFACE_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_NETGEN

// Local includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_tet_interface.h"
#include "libmesh/point.h" // used for specifying holes

namespace nglib {
#include "netgen/nglib/nglib.h"
}

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
 * Class \p NetGenMeshInterface provides an interface for
 * tetrahedralization of meshes using the NetGen library.  For
 * information about TetGen cf.
 * <a href="https://github.com/NGSolve/netgen/">NetGen github
 * repository</a>.
 *
 * \author Roy H. Stogner
 * \date 2024
 */
class NetGenMeshInterface : public MeshTetInterface
{
public:

  /**
   * Constructor. Takes a reference to the mesh containing the
   * triangulated surface which is to be tetrahedralized.
   */
  explicit
  NetGenMeshInterface (UnstructuredMesh & mesh);

  /**
   * Empty destructor.
   */
  virtual ~NetGenMeshInterface() override = default;

  /**
   * Method invokes NetGen library to compute a tetrahedralization
   */
  virtual void triangulate () override;

protected:

  /**
   * Tetgen only operates on serial meshes.
   */
  MeshSerializer _serializer;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_NETGEN

#endif // LIBMESH_MESH_NETGEN_INTERFACE_H
