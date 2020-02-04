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



#ifndef LIBMESH_FRO_IO_H
#define LIBMESH_FRO_IO_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_output.h"

// C++ includes

namespace libMesh
{

// Forward declarations
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;

/**
 * This class implements writing meshes in the .fro format
 * used by the MIT ACDL.  Valid only for triangular meshes.
 *
 * \author Benjamin S. Kirk
 * \date 2007
 */
class FroIO : public MeshOutput<MeshBase>
{
public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  explicit
  FroIO (const MeshBase &);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string &) override;
};



// ------------------------------------------------------------
// FroIO inline members
inline
FroIO::FroIO (const MeshBase & mesh_in) :
  MeshOutput<MeshBase> (mesh_in)
{
}


} // namespace libMesh


#endif // LIBMESH_FRO_IO_H
