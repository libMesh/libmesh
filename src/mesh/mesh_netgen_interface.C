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

#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_HAVE_NETGEN


// C++ includes
#include <sstream>

// Local includes
#include "libmesh/mesh_netgen_interface.h"

#include "libmesh/boundary_info.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_tri3.h"
#include "libmesh/unstructured_mesh.h"

namespace {

// RAII for exception safety
class WrappedNgMesh
{
public:
  WrappedNgMesh()
  {
    _ngmesh = nglib::Ng_NewMesh();
  }

  ~WrappedNgMesh()
  {
    nglib::Ng_DeleteMesh(_ngmesh);
  }

  operator nglib::Ng_Mesh* () {
    return _ngmesh;
  }

private:
  nglib::Ng_Mesh * _ngmesh;
};

}

namespace libMesh
{

//----------------------------------------------------------------------
// NetGenMeshInterface class members
NetGenMeshInterface::NetGenMeshInterface (UnstructuredMesh & mesh) :
  MeshTetInterface(mesh),
  _desired_volume(0),
  _smooth_after_generating(true),
  _serializer(mesh)
{
}



void NetGenMeshInterface::triangulate ()
{
  this->check_hull_integrity();

  WrappedNgMesh ngmesh;
}



} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_NETGEN
