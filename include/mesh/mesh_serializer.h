// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_SERIALIZER_H
#define LIBMESH_MESH_SERIALIZER_H

// Local includes

// C++ includes


namespace libMesh
{
// Forward declarations
class MeshBase;

/**
 * Temporarily serialize a ParallelMesh for output; a distributed
 * mesh is allgathered by the MeshSerializer constructor if
 * need_serial is true, then remote elements are deleted again by the
 * destructor.
 */
class MeshSerializer
{
public:
  MeshSerializer(MeshBase & mesh, bool need_serial = true);

  ~MeshSerializer();

private:
  MeshBase & _mesh;
  bool reparallelize;
};

} // namespace libMesh

#endif // LIBMESH_MESH_SERIALIZER_H
