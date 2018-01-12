// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_H
#define LIBMESH_MESH_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PARMESH
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/distributed_mesh.h"
namespace libMesh {
typedef DistributedMesh DefaultMesh;
}
#else
#include "libmesh/replicated_mesh.h"
namespace libMesh {
typedef ReplicatedMesh DefaultMesh;
}
#endif

namespace libMesh
{

// Forward declarations don't like typedefs...
// typedef ReplicatedMesh Mesh;

/**
 * The \p Mesh class is a thin wrapper, around the \p ReplicatedMesh
 * class by default.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Manages a collection of Nodes and Elems.
 */
class Mesh : public DefaultMesh
{
public:

  /**
   * Constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  Mesh (const Parallel::Communicator & comm_in,
        unsigned char dim=1)
    : DefaultMesh(comm_in,dim) {}

#ifndef LIBMESH_DISABLE_COMMWORLD
  /**
   * Deprecated constructor.  Takes \p dim, the dimension of the mesh.
   * The mesh dimension can be changed (and may automatically be
   * changed by mesh generation/loading) later.
   */
  explicit
  Mesh (unsigned char dim=1)
    : DefaultMesh(dim) {}
#endif

  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  Mesh (const UnstructuredMesh & other_mesh) : DefaultMesh(other_mesh) {}

  /**
   * Destructor.
   */
  ~Mesh() {}
};



} // namespace libMesh



#endif // LIBMESH_MESH_H
