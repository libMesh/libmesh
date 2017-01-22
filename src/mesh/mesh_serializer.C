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


// Local includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h" // parallel_only() macro

namespace libMesh
{

MeshSerializer::MeshSerializer(MeshBase & mesh, bool need_serial, bool serial_only_needed_on_proc_0) :
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



MeshSerializer::~MeshSerializer()
{
  if (reparallelize)
    _mesh.delete_remote_elements();
}

} // namespace libMesh
