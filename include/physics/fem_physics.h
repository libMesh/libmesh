
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
 * @author Roy H. Stogner 2012
 */

// ------------------------------------------------------------
// Finite Element Method Physics class definition

class FEMPhysics : public virtual DifferentiablePhysics
{
public:

  /**
   * Constructor.
   */
  FEMPhysics () :
    DifferentiablePhysics(),
    _mesh_sys                         (NULL),
    _mesh_x_var                       (libMesh::invalid_uint),
    _mesh_y_var                       (libMesh::invalid_uint),
    _mesh_z_var                       (libMesh::invalid_uint)
    {}

  /**
   * Destructor.
   */
  virtual ~FEMPhysics () {}

  /**
   * Tells the FEMPhysics that system \p sys contains the
   * isoparametric Lagrangian variables which correspond to the
   * coordinates of mesh nodes, in problems where the mesh itself is
   * expected to move in time.
   *
   * The system with mesh coordinate data (which may be \p this system
   * itself, for fully coupled moving mesh problems) is currently
   * assumed to have new (end of time step) mesh coordinates stored in
   * solution, old (beginning of time step) mesh coordinates stored in
   * _old_nonlinear_solution, and constant velocity motion during each
   * time step.
   *
   * Activating this function ensures that local (but not neighbor!) element
   * geometry is correctly repositioned when evaluating element residuals.
   *
   * Currently \p sys must be \p *this for a tightly coupled moving
   * mesh problem or NULL to stop mesh movement; loosely coupled
   * moving mesh problems are not implemented.
   *
   * This code is experimental.  "Trust but verify, and not in that
   * order"
   */
  virtual void set_mesh_system(System* sys);

  /**
   * Returns a reference to the system with variables corresponding to
   * mesh nodal coordinates, or NULL if the mesh is fixed.
   * Useful for ALE calculations.
   */
  const System* get_mesh_system() const;

  /**
   * Tells the FEMPhysics that variable \p var from the mesh system
   * should be used to update the x coordinate of mesh nodes, in problems where
   * the mesh itself is expected to move in time.
   *
   * The system with mesh coordinate data (which may be this system itself, for
   * fully coupled moving mesh problems) is currently assumed to have new (end
   * of time step) mesh coordinates stored in solution, old (beginning of time
   * step) mesh coordinates stored in _old_nonlinear_solution, and constant
   * velocity motion during each time step.
   *
   * Activating this function ensures that local (but not neighbor!) element
   * geometry is correctly repositioned when evaluating element residuals.
   */
  virtual void set_mesh_x_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh x coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_x_var() const;

  /**
   * Tells the FEMPhysics that variable \p var from the mesh system
   * should be used to update the y coordinate of mesh nodes.
   */
  virtual void set_mesh_y_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh y coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_y_var() const;

  /**
   * Tells the FEMPhysics that variable \p var from the mesh system
   * should be used to update the z coordinate of mesh nodes.
   */
  virtual void set_mesh_z_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh z coordinate. Useful for ALE calculations.
   */
  unsigned int get_mesh_z_var() const;

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
                                  DiffContext &context);


  /**
   * Adds a mass vector contribution on \p elem to elem_residual.
   * If this method receives request_jacobian = true, then it
   * should compute elem_jacobian and return true if possible.  If
   * elem_jacobian has not been computed then the method should
   * return false.
   *
   * Most problems can use the reimplementation in
   * FEMPhysics::mass_residual; few users will need to reimplement
   * this themselves.
   */
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext &);

protected:

  /**
   * System from which to acquire moving mesh information
   */
  System *_mesh_sys;

  /**
   * Variables from which to acquire moving mesh information
   */
  unsigned int _mesh_x_var, _mesh_y_var, _mesh_z_var;
};



// ------------------------------------------------------------
// FEMPhysics inline methods


inline
void FEMPhysics::set_mesh_system(System* sys)
{
  // For now we assume that we're doing fully coupled mesh motion
//  if (sys && sys != this)
//    libmesh_not_implemented();

  // For the foreseeable future we'll assume that we keep these
  // Systems in the same EquationSystems
  // libmesh_assert_equal_to (&this->get_equation_systems(),
  //                          &sys->get_equation_systems());

  // And for the immediate future this code may not even work
  libmesh_experimental();

  _mesh_sys = sys;
}



inline
void FEMPhysics::set_mesh_x_var (unsigned int var)
{
  _mesh_x_var = var;
}



inline
void FEMPhysics::set_mesh_y_var (unsigned int var)
{
  _mesh_y_var = var;
}



inline
void FEMPhysics::set_mesh_z_var (unsigned int var)
{
  _mesh_z_var = var;
}



inline
const System* FEMPhysics::get_mesh_system() const
{
  return _mesh_sys;
}

inline
unsigned int FEMPhysics::get_mesh_x_var() const
{
  return _mesh_x_var;
}

inline
unsigned int FEMPhysics::get_mesh_y_var() const
{
  return _mesh_y_var;
}

inline
unsigned int FEMPhysics::get_mesh_z_var() const
{
  return _mesh_z_var;
}



} // namespace libMesh


#endif // LIBMESH_FEM_PHYSICS_H
