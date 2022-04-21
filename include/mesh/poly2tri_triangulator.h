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
#include "libmesh/dof_object.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/triangulator_interface.h"

namespace libMesh
{

// Forward Declarations
class BoundaryInfo;
class Elem;

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
   * which are to be triangulated must be provided.  The first
   * \p n_boundary_nodes are expected to form a closed loop around the
   * mesh domain; any subsequent nodes are expected to be interior
   * nodes or in the middle of (internal hole or external) boundary
   * segments.
   *
   * If \p n_boundary_nodes is not supplied or is \p invalid_id then
   * all mesh points are expected to be boundary polyline points.
   */
  explicit
  Poly2TriTriangulator(UnstructuredMesh & mesh,
                       dof_id_type n_boundary_nodes =
                         DofObject::invalid_id);

  /**
   * Empty destructor.  Defaulted in the .C so we can forward declare
   * unique_ptr contents.
   */
  virtual ~Poly2TriTriangulator();

  /**
   * Internally, this calls the poly2tri triangulation code in a loop,
   * inserting our owner Steiner points as necessary to promote mesh
   * quality.
   */
  virtual void triangulate() override;

  /**
   * Set a function giving desired triangle area as a function of
   * position.  Set this to nullptr to disable position-dependent area
   * constraint (falling back on desired_area()).
   */
  virtual void set_desired_area_function (FunctionBase<Real> * desired) override;

  /**
   * Get the function giving desired triangle area as a function of
   * position, or \p nullptr if no such function has been set.
   */
  virtual FunctionBase<Real> * get_desired_area_function () override;

  /**
   * Set whether or not the triangulation is allowed to refine the
   * mesh boundary when refining the interior.  This is true by
   * default, but may be set to false to make the mesh boundary more
   * predictable (and so easier to stitch to other meshes) later.
   *
   * This may not be implemented in all subclasses.
   */
  virtual void set_refine_boundary_allowed (bool refine_bdy_allowed)
  { _refine_bdy_allowed = refine_bdy_allowed; }

  /**
   * Get whether or not the triangulation is allowed to refine the
   * mesh boundary when refining the interior.  True by default.
   */
  virtual bool refine_boundary_allowed ()
  { return _refine_bdy_allowed; }


protected:
  /**
   * Is refining this element's boundary side allowed?
   */
  bool is_refine_boundary_allowed(const BoundaryInfo & boundary_info,
                                  const Elem & elem,
                                  unsigned int side);

  /**
   * Triangulate the current mesh and hole points.
   */
  void triangulate_current_points();

  /**
   * Add Steiner points as new mesh nodes, as necessary to refine an
   * existing trangulation.  Returns true iff new points were added.
   */
  bool insert_refinement_points();

  /**
   * Returns true if the given element ought to be refined according
   * to current criteria.
   */
  bool should_refine_elem(Elem & elem);

private:

  /**
   * We might have to replace the user-provided holes with refined
   * versions
   */
  std::map<const Hole *, std::unique_ptr<ArbitraryHole>> replaced_holes;

  /**
   * We only operate on serialized meshes.
   */
  MeshSerializer _serializer;

  /**
   * Keep track of how many mesh nodes are boundary nodes.
   */
  dof_id_type _n_boundary_nodes;

  /**
   * Location-dependent area requirements
   */
  std::unique_ptr<FunctionBase<Real>> _desired_area_func;

  /**
   * Whether to allow boundary refinement
   */
  bool _refine_bdy_allowed;
};

} // namespace libMesh

#endif // LIBMESH_HAVE_TRIANGLE

#endif // ifndef LIBMESH_POLY2TRI_TRIANGULATOR_H
