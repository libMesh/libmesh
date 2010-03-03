
// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __fem_system_h__
#define __fem_system_h__

// C++ includes

// Local Includes
#include "diff_system.h"
#include "fem_context.h"

// Forward Declarations


/**
 * This class provides a specific system class.  It aims
 * at nonlinear implicit systems, requiring only a
 * cell residual calculation from the user.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the class \p ExplicitSystem.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2006
 */

// ------------------------------------------------------------
// FEMSystem class definition

class FEMSystem : public DifferentiableSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FEMSystem (EquationSystems& es,
	         const std::string& name,
	         const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~FEMSystem ();

  /**
   * The type of system.
   */
  typedef FEMSystem sys_type;

  /**
   * The type of the parent.
   */
  typedef DifferentiableSystem Parent;
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

  /**
   * Prepares \p matrix or \p rhs for matrix assembly.
   * Users may reimplement this to add pre- or post-assembly
   * code before or after calling FEMSystem::assembly()
   */
  virtual void assembly (bool get_residual, bool get_jacobian);

  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   *
   * For moving mesh systems, this also translates the mesh to the
   * solution position.
   */
  virtual void solve ();

  /**
   * Tells the FEMSystem that variable \p var is evolving with
   * respect to time.  In general, the user's init() function
   * should call time_evolving() for any variables which
   * behave like du/dt = F(u), and should not call time_evolving()
   * for any variables which behave like 0 = G(u).
   *
   * Most derived systems will not have to reimplment this function; however
   * any system which reimplements mass_residual() may have to reimplement
   * time_evolving() to prepare data structures.
   */
  virtual void time_evolving (unsigned int var);

  /**
   * Tells the FEMSystem that system \p sys contains the
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
   */
  const System* get_mesh_system() const;

  /**
   * Tells the FEMSystem that variable \p var from the mesh system
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
   * mesh x coordinate
   */
  unsigned int get_mesh_x_var() const;

  /**
   * Tells the FEMSystem that variable \p var from the mesh system
   * should be used to update the y coordinate of mesh nodes.
   */
  virtual void set_mesh_y_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh y coordinate
   */
  unsigned int get_mesh_y_var() const;

  /**
   * Tells the FEMSystem that variable \p var from the mesh system
   * should be used to update the z coordinate of mesh nodes.
   */
  virtual void set_mesh_z_var(unsigned int var);

  /**
   * Returns the variable number corresponding to the
   * mesh z coordinate
   */
  unsigned int get_mesh_z_var() const;

  /**
   * Tells the FEMSystem to set the degree of freedom coefficients
   * which should correspond to mesh nodal coordinates.
   */
  void mesh_position_get();

  /**
   * Tells the FEMSystem to set the mesh nodal coordinates
   * which should correspond to degree of freedom coefficients.
   */
  void mesh_position_set();

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
   * Most problems can use the FEMSystem::mass_residual implementation,
   * which calculates the residual (u, phi_i) and jacobian (phi_i, phi_j);
   * few users will need to reimplement this themselves.  Using a custom
   * mass matrix (e.g. for divergence-free elements or mass lumping)
   * requires reimplementing mass_residual().
   */
  virtual bool mass_residual (bool request_jacobian,
                              DiffContext &context);

  /**
   * Builds a FEMContext object with enough information to do
   * evaluations on each element.
   *
   * For most problems, the default FEMSystem implementation is correct; users
   * who subclass FEMContext will need to also reimplement this method to build
   * it.
   */
  virtual AutoPtr<DiffContext> build_context();

  /*
   * Prepares the result of a build_context() call for use.
   * 
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular physics requires.
   */
  virtual void init_context(DiffContext &);
 
  /**
   * Runs a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess ();

  /**
   * Runs a qoi assembly loop over all elements, and if
   * \p assemble_qoi_sides is true over all sides.
   *
   * Users may have to override this function for quantities of
   * interest that are not expressible as a sum of element qois.
   */
  virtual void assemble_qoi
    (const QoISet& indices = QoISet());

  /**
   * Runs a qoi derivative assembly loop over all elements, and if
   * \p assemble_qoi_sides is true over all sides.
   *
   * Users may have to override this function for quantities of
   * interest that are not expressible as a sum of element qois.
   */
  virtual void assemble_qoi_derivative
    (const QoISet& indices = QoISet());

  /**
   * If fe_reinit_during_postprocess is true (it is true by default), FE
   * objects will be reinit()ed with their default quadrature rules.  If false,
   * FE objects will need to be reinit()ed by the user or will be in an
   * undefined state.
   */
  bool fe_reinit_during_postprocess;

  /**
   * By default, when calling the user-defined residual functions, the
   * FEMSystem will first set up an appropriate
   * FEType::default_quadrature_rule() object for performing the integration.
   * This rule will integrate elements of order up to 2*p+1 exactly (where p is
   * the sum of the base FEType and local p refinement levels), but if
   * additional (or reduced) quadrature accuracy is desired then this
   * extra_quadrature_order (default 0) will be added.
   */
  int extra_quadrature_order;

  /**
   * If calculating numeric jacobians is required, the FEMSystem
   * will perturb each solution vector entry by numerical_jacobian_h
   * when calculating finite differences.
   */ 
  Real numerical_jacobian_h;

  /**
   * If verify_analytic_jacobian is equal to zero (as it is by
   * default), no numeric jacobians will be calculated unless
   * an overloaded element_time_derivative(), element_constraint(),
   * side_time_derivative(), or side_constraint() function cannot
   * provide an analytic jacobian upon request.
   * 
   * If verify_analytic_jacobian is equal to the positive value tol,
   * then any time a full analytic element jacobian can be calculated
   * it will be tested against a numerical jacobian on the same element,
   * and the program will abort if the relative error (in matrix l1 norms)
   * exceeds tol.
   */
  Real verify_analytic_jacobians;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Syntax sugar to make numerical_jacobian() declaration easier.
   */
  typedef bool (TimeSolver::*TimeSolverResPtr)(bool, DiffContext&);

  /**
   * Uses the results of multiple \p res calls
   * to numerically differentiate the corresponding jacobian.
   */
  void numerical_jacobian (TimeSolverResPtr res,
                           FEMContext &context);

  /**
   * Uses the results of multiple element_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element.
   */
  void numerical_elem_jacobian (FEMContext &context);

  /**
   * Uses the results of multiple side_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element's side.
   */
  void numerical_side_jacobian (FEMContext &context);

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
// FEMSystem inline methods



inline
void FEMSystem::set_mesh_system(System* sys)
{
  // For now we assume that we're doing fully coupled mesh motion
  if (sys && sys != this)
    libmesh_not_implemented();

  // For the foreseeable future we'll assume that we keep these
  // Systems in the same EquationSystems
  libmesh_assert(&this->get_equation_systems() ==
                 &sys->get_equation_systems());

  _mesh_sys = sys;
}



inline
const System* FEMSystem::get_mesh_system() const
{
  return _mesh_sys;
}



inline
void FEMSystem::set_mesh_x_var (unsigned int var)
{
  _mesh_x_var = var;
}



inline
unsigned int FEMSystem::get_mesh_x_var() const
{
  return _mesh_x_var;
}



inline
void FEMSystem::set_mesh_y_var (unsigned int var)
{
  _mesh_y_var = var;
}



inline
unsigned int FEMSystem::get_mesh_y_var() const
{
  return _mesh_y_var;
}



inline
void FEMSystem::set_mesh_z_var (unsigned int var)
{
  _mesh_z_var = var;
}



inline
unsigned int FEMSystem::get_mesh_z_var() const
{
  return _mesh_z_var;
}


#endif
