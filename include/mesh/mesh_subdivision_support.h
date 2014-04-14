// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_SUBDIVISION_SUPPORT_H
#define LIBMESH_MESH_SUBDIVISION_SUPPORT_H



// Local Includes -----------------------------------
#include "libmesh/libmesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/elem.h"

// C++ Includes   -----------------------------------

namespace libMesh
{

// ------------------------------------------------------------
// MeshTools::Subdivision namespace
namespace MeshTools
{
/**
 * Utility functions for subdivision surface operations on a \p Mesh.
 */
namespace Subdivision
{
/**
 * Determines the 1-ring of element \p elem, and writes it to the
 * \p nodes vector. This is necessary because subdivision
 * elements have a larger local support than conventionally
 * interpolated elements. The 1-ring may, for instance, look like
 * this:
 * \verbatim
 *    N+4 - N+1 - N+2
 *    / \   / \   / \
 *   /   \ /   \ /   \
 * N+5 -- N --- 1 -- N+3
 *   \   / \ e / \   /
 *    \ /   \ /   \ /
 *    N-1--- 0 --- 2
 *      \   /|\   /
 *       \ / | \ /
 *        5--4--3
 * \endverbatim
 */
void find_one_ring(const Tri3Subdivision* elem, std::vector<Node*>& nodes);

/**
 * Turns a triangulated \p mesh into a subdivision mesh. This
 * function normally needn't be called by the user, because it is
 * invoked by \p prepare_subdivision_mesh.
 */
void all_subdivision(MeshBase& mesh);

/**
 * Prepares the \p mesh for use with subdivision elements. The
 * \p ghosted flag determines how boundaries are treated. If \p false,
 * a new layer of "ghost" elements is appended along the domain
 * boundaries. If \p true, the outermost element layer is taken as
 * ghosts, i.e. no new elements are added.
 */
void prepare_subdivision_mesh(MeshBase& mesh, bool ghosted = false);

/**
 * Adds a new layer of "ghost" elements along the domain boundaries.
 * This function normally needn't be called by the user, because it
 * is invoked by \p prepare_subdivision_mesh.
 */
void add_boundary_ghosts(MeshBase& mesh);

/**
 * Flags the outermost element layer along the domain boundaries as
 * "ghost" elements. This function normally needn't be called by the
 * user, because it is invoked by \p prepare_subdivision_mesh.
 */
void tag_boundary_ghosts(MeshBase& mesh);

/**
 * A lookup table for the increment modulo 3 operation, for iterating
 * through the three nodes per element in positive direction.
 */
static const unsigned int next[3] = {1,2,0};

/**
 * A lookup table for the decrement modulo 3 operation, for iterating
 * through the three nodes per element in negative direction.
 */
static const unsigned int prev[3] = {2,0,1};

} // end namespace MeshTools::Subdivision
} // end namespace MeshTools

} // namespace libMesh


#endif // LIBMESH_MESH_SUBDIVISION_SUPPORT_H
