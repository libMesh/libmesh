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



#ifndef LIBMESH_SERIAL_MESH_H
#define LIBMESH_SERIAL_MESH_H

#include "libmesh/replicated_mesh.h"


namespace libMesh
{

// For backwards compatibility with our original name, users can still
// refer to a ReplicatedMesh as a SerialMesh.
//
// This has to be a shim class rather than a typedef so that forward
// declarations still work.
class SerialMesh : public ReplicatedMesh
{
public:
  explicit
  SerialMesh (const Parallel::Communicator & comm_in,
              unsigned char dim=1)
    : ReplicatedMesh(comm_in,dim) {}

#ifndef LIBMESH_DISABLE_COMMWORLD
  explicit
  SerialMesh (unsigned char dim=1)
    : ReplicatedMesh(dim)
  {
    libmesh_deprecated();
  }
#endif

  SerialMesh (const UnstructuredMesh & other_mesh) : ReplicatedMesh(other_mesh) {}

  virtual UniquePtr<MeshBase> clone () const libmesh_override
  { return UniquePtr<MeshBase>(new SerialMesh(*this)); }

  ~SerialMesh() {}
};

} // namespace libMesh


#endif // LIBMESH_SERIAL_MESH_H
