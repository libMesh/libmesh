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


#ifndef LIBMESH_MESH_SERIALIZER_H
#define LIBMESH_MESH_SERIALIZER_H

// Local includes
#include "libmesh/mesh_base.h"
#include "libmesh/parallel_only.h"

namespace libMesh
{
// Forward declarations
template <typename> class MeshBaseTempl;

/**
 * Temporarily serialize a DistributedMesh for output; a distributed
 * mesh is allgathered by the MeshSerializer constructor if
 * need_serial is true, then remote elements are deleted again by the
 * destructor.
 *
 * \author Roy Stogner
 * \date 2011
 * \brief Temporarily serializes a DistributedMesh for output.
 */
template <typename RealType = Real>
class MeshSerializerTempl
{
  typedef MeshBaseTempl<RealType> MeshBase;

public:
  MeshSerializerTempl(MeshBase & mesh, bool need_serial = true, bool serial_only_needed_on_proc_0 = false);

  ~MeshSerializerTempl();

private:
  MeshBase & _mesh;
  bool reparallelize;
};

template <typename RealType>
MeshSerializerTempl<RealType>::MeshSerializerTempl(
  MeshBase & mesh, bool need_serial, bool serial_only_needed_on_proc_0) :
    _mesh(mesh),
    reparallelize(false)
{
  libmesh_parallel_only(mesh.comm());
  if (need_serial && !_mesh.is_serial() && !serial_only_needed_on_proc_0) {
    reparallelize = true;
    _mesh.allgather();
  }
  else if (need_serial && !_mesh.is_serial() && serial_only_needed_on_proc_0) {
    // Note: NOT reparallelizing on purpose.
    // Just waste a bit of space on processor 0 to speed things up
    _mesh.gather_to_zero();
  }
}



template <typename RealType>
MeshSerializerTempl<RealType>::~MeshSerializerTempl()
{
  if (reparallelize)
    _mesh.delete_remote_elements();
}

typedef MeshSerializerTempl<Real> MeshSerializer;

} // namespace libMesh

#endif // LIBMESH_MESH_SERIALIZER_H
