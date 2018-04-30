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



#ifndef LIBMESH_BOUNDARY_MESH_H
#define LIBMESH_BOUNDARY_MESH_H

// Local Includes
#include "libmesh/mesh.h"

namespace libMesh
{

/**
 * The \p BoundaryMesh is a \p Mesh in its own right, but it
 * contains a description of the boundary of some other mesh.
 * This is useful for writing the boundary of a domain for inspecting
 * boundary conditions and other things.
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief Specialized mesh used for keeping track of boundary elements.
 */
class BoundaryMesh : public Mesh
{
public:
  /**
   * Constructor. Initializes dimension and processor id.
   */
  explicit
  BoundaryMesh (const Parallel::Communicator & comm_in,
                unsigned char dim=1);

  /**
   * Destructor.
   */
  ~BoundaryMesh();
};

} // namespace libMesh


#endif // LIBMESH_BOUNDARY_MESH_H
