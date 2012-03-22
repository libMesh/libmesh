// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __fro_io_h__
#define __fro_io_h__

// Local includes
#include "libmesh_common.h"
#include "mesh_output.h"

namespace libMesh
{

// Forward declarations
class MeshBase;



/**
 * This class implements writing meshes in the .fro format
 * used by the MIT ACDL.  Valid only for triangular meshes.
 *
 * @author Benjamin S. Kirk, 2007
 */

// ------------------------------------------------------------
// FroIO class definition
class FroIO : public MeshOutput<MeshBase>
{
 public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  FroIO (const MeshBase&);

  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& );


 private:

};



// ------------------------------------------------------------
// FroIO inline members
inline
FroIO::FroIO (const MeshBase& mesh) :
  MeshOutput<MeshBase> (mesh)
{
}


} // namespace libMesh


#endif // #define __fro_io_h__
