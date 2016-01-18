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


#ifndef LIBMESH_MESH_TRIANGLE_WRAPPER_H
#define LIBMESH_MESH_TRIANGLE_WRAPPER_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local Includes
#include "libmesh/libmesh_common.h" // Real
#include "libmesh/enum_elem_type.h" // For ElemType declaration below

// C++ includes
#include <cstddef>

namespace libMesh
{
// Forward declarations
class UnstructuredMesh;

// Make sure Triangle uses our "Real" as its "REAL"
typedef Real REAL;

/**
 * A special namespace for wrapping the standard Triangle API,
 * as well as some helper functions for initializing/destroying
 * the structs triangle uses to communicate.
 */
namespace TriangleWrapper
{
extern "C"
{
#include "triangle.h"
}

enum IO_Type {
  INPUT  = 0,
  OUTPUT = 1,
  BOTH   = 2};

/**
 * Initializes the fields of t to NULL/0 as necessary.
 * This is helpful for preventing the access of uninitialized
 * memory when working with C, which has no constructors or
 * destructors.
 */
void init(triangulateio & t);

/**
 * Frees any memory which has been dynamically allocated by
 * Triangle.  Note the following facts:
 * 1) Triangle does not free any memory itself
 * 2) It is always safe to call free on a NULL pointer.
 *
 * However, triangle *does* shallow-copy (for example)
 * the holelist pointer from the input to output struct **without**
 * performing a deep copy of the holelist itself.  Therefore, double-free
 * will occur without additional care!
 */
void destroy(triangulateio & t, IO_Type);

/**
 * Copies triangulation data computed by triange from a triangulateio object
 * to a LibMesh mesh.   This routine is used internally by the
 * MeshTools::Generation::build_delaunay_square(...) and
 * MeshTools::Generation::build_delaunay_square_with_hole(...) routines.
 */
void copy_tri_to_mesh(const triangulateio & triangle_data_input,
                      UnstructuredMesh & mesh_output,
                      const ElemType type);
} // namespace TriangleWrapper
} // namespace libMesh


#endif // LIBMESH_HAVE_TRIANGLE

#endif // LIBMESH_MESH_TRIANGLE_WRAPPER_H
