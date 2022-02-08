// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_POLY2TRI_TRIANGULATOR_H
#define LIBMESH_POLY2TRI_TRIANGULATOR_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_POLY2TRI

// Local Includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/triangulator_interface.h"

namespace libMesh
{

/**
 * A C++ interface between LibMesh and the poly2tri library, with
 * custom code for Steiner point insertion.
 *
 * \author Roy H. Stogner
 * \date 2022
 */
class Poly2TriTriangulator : public TriangulatorInterface
{
public:
  /**
   * The constructor.  A reference to the mesh containing the points
   * which are to be triangulated must be provided.  Unless otherwise
   * specified, a convex hull will be computed for the set of input points
   * and the convex hull will be meshed.
   */
  explicit
  Poly2TriTriangulator(UnstructuredMesh & mesh);

  /**
   * Empty destructor.
   */
  ~Poly2TriTriangulator() = default;

  /**
   * Internally, this calls the poly2tri triangulation code in a loop,
   * inserting our owner Steiner points as necessary to promote mesh
   * quality.
   */
  virtual void triangulate() override;

private:

  /**
   * We only operate on serialized meshes.
   */
  MeshSerializer _serializer;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE

#endif // ifndef LIBMESH_POLY2TRI_TRIANGULATOR_H
