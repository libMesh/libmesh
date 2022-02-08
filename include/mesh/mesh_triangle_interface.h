// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_MESH_TRIANGLE_INTERFACE_H
#define LIBMESH_MESH_TRIANGLE_INTERFACE_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_TRIANGLE

// Local Includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/triangulator_interface.h"

namespace libMesh
{

/**
 * A C++ interface between LibMesh and the Triangle library written by
 * J.R. Shewchuk.
 *
 * \author John W. Peterson
 * \date 2011
 */
class TriangleInterface : public TriangulatorInterface
{
public:
  /**
   * The constructor.  A reference to the mesh containing the points
   * which are to be triangulated must be provided.  Unless otherwise
   * specified, a convex hull will be computed for the set of input points
   * and the convex hull will be meshed.
   */
  explicit
  TriangleInterface(UnstructuredMesh & mesh);

  /**
   * Empty destructor.
   */
  ~TriangleInterface() = default;

  /**
   * Internally, this calls Triangle's triangulate routine.
   */
  virtual void triangulate() override;

  /**
   * Sets and/or gets additional flags to be passed to triangle
   */
  std::string & extra_flags() {return _extra_flags;}

private:
  /**
   * Additional flags to be passed to triangle
   */
  std::string _extra_flags;

  /**
   * Triangle only operates on serial meshes.
   */
  MeshSerializer _serializer;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE

#endif // ifndef LIBMESH_MESH_TRIANGLE_INTERFACE_H
