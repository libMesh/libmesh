
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



#ifndef LIBMESH_FEM_PHYSICS_H
#define LIBMESH_FEM_PHYSICS_H

// Local Includes
#include "libmesh/libmesh.h" // for libMesh::invalid_uint
#include "libmesh/diff_physics.h"

// C++ includes

namespace libMesh
{

/**
 * This class provides a specific system class.  It aims
 * to generalize any system, linear or nonlinear, which
 * provides both a residual and a Jacobian.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * \author Roy H. Stogner
 * \date 2012
 */
class FEMPhysics : public virtual DifferentiablePhysics
{
public:

  /**
   * Constructor.
   */
  FEMPhysics () :
    DifferentiablePhysics()
  {}

  /**
   * Destructor.
   */
  virtual ~FEMPhysics () {}

  /**
   * Adds a pseudo-convection contribution on \p elem to
   * elem_residual, if the nodes of \p elem are being translated by a
   * moving mesh.
   *
   * This function assumes that the user's time derivative equations
   * (except for any equations involving unknown mesh xyz coordinates
   * themselves) are expressed in an Eulerian frame of reference, and
   * that the user is satisfied with an unstabilized convection term.
   * Lagrangian equations will probably require overriding
   * eulerian_residual() with a blank function; ALE or stabilized
   * formulations will require reimplementing eulerian_residual()
   * entirely.
   */
  virtual bool eulerian_residual (bool request_jacobian,
                                  DiffContext & context) libmesh_override;

  /**
   * Subtracts a mass vector contribution on \p elem from
   * elem_residual.
   *
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Many problems can use the reimplementation in
   * FEMPhysics::mass_residual which subtracts (du/dt,v) for each
   * transient variable u; users with more complicated transient
   * problems will need to reimplement this themselves.
   */
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext &) libmesh_override;
};



} // namespace libMesh


#endif // LIBMESH_FEM_PHYSICS_H
