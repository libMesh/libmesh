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

#ifndef LIBMESH_INF_ELEM_BUILDER_H
#define LIBMESH_INF_ELEM_BUILDER_H


#include "libmesh/libmesh_config.h"


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/id_types.h"
#include "libmesh/point.h"

// C++ includes
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace libMesh
{

// Forward Declarations
class MeshBase;
class Node;

/**
 * This class is used to build infinite elements on
 * top of an existing mesh.  It only makes sense to
 * use this if LIBMESH_ENABLE_INFINITE_ELEMENTS is true.
 *
 * \author Daniel Dreyer
 * \author John W. Peterson
 * \date 2004
 */
class InfElemBuilder
{
public:
  /**
   * Constructor.
   */
  explicit
  InfElemBuilder(MeshBase & mesh) : _mesh(mesh) {}

  /**
   * Useful typedef
   */
  typedef std::pair<bool, double> InfElemOriginValue;

  /**
   * Build infinite elements atop a volume-based mesh,
   * determine origin automatically.  Also returns the
   * origin as a \p const \p Point to make it more obvious that
   * the origin should not change after the infinite elements
   * have been built.  When symmetry planes are present, use
   * the version with optional symmetry switches.
   * The flag \p be_verbose enables some diagnostic output.
   */
  const Point build_inf_elem (const bool be_verbose = false);

  /**
   * @returns the origin of the infinite elements.
   * Builds infinite elements atop a volume-based mesh.
   * Finds all faces on the outer boundary and build infinite elements
   * on them.  Using the \p InfElemOriginValue the user can
   * prescribe only selected origin coordinates.  The remaining
   * coordinates are computed from the center of the bounding box
   * of the mesh.
   *
   * During the search for faces on which infinite elements are built,
   * @e interior faces that are not on symmetry planes are found, too.
   * When an (optional) pointer to \p inner_boundary_nodes is provided,
   * then this vector will be filled with the nodes that lie on the
   * inner boundary.
   *
   * Faces which lie in at least one symmetry plane are skipped.
   * The three optional booleans \p x_sym, \p y_sym,
   * \p z_sym indicate symmetry planes (through the origin, obviously)
   * perpendicular to the \p x, \p y and \p z direction,
   * respectively.
   * The flag \p be_verbose enables some diagnostic output.
   */
  const Point build_inf_elem (const InfElemOriginValue & origin_x,
                              const InfElemOriginValue & origin_y,
                              const InfElemOriginValue & origin_z,
                              const bool x_sym = false,
                              const bool y_sym = false,
                              const bool z_sym = false,
                              const bool be_verbose = false,
                              std::vector<const Node *> * inner_boundary_nodes = libmesh_nullptr);



private:
  /**
   * Build infinite elements atop a volume-based mesh.
   * Actual implementation.
   */
  void build_inf_elem (const Point & origin,
                       const bool x_sym = false,
                       const bool y_sym = false,
                       const bool z_sym = false,
                       const bool be_verbose = false,
                       std::set<std::pair<dof_id_type,
                       unsigned int> > * inner_faces = libmesh_nullptr);
  /**
   * Reference to the mesh we're building infinite
   * elements for.
   */
  MeshBase & _mesh;
};


} // namespace libMesh

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
#endif // LIBMESH_INF_ELEM_BUILDER_H
