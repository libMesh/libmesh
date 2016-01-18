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



#ifndef LIBMESH_FEM_SYSTEM_H
#define LIBMESH_FEM_SYSTEM_H

// Local Includes
#include "libmesh/diff_system.h"
#include "libmesh/fem_physics.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward Declarations
class DiffContext;
class FEMContext;


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
 * \author Roy H. Stogner
 * \date 2006
 */
class FEMSystem : public DifferentiableSystem,
                  public FEMPhysics
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  FEMSystem (EquationSystems & es,
             const std::string & name,
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
   * Prepares \p matrix or \p rhs for matrix assembly.
   * Users may reimplement this to add pre- or post-assembly
   * code before or after calling FEMSystem::assembly()
   */
  virtual void assembly (bool get_residual,
                         bool get_jacobian,
                         bool apply_heterogeneous_constraints = false) libmesh_override;

  /**
   * Invokes the solver associated with the system.  For steady state
   * solvers, this will find a root x where F(x) = 0.  For transient
   * solvers, this will integrate dx/dt = F(x).
   *
   * For moving mesh systems, this also translates the mesh to the
   * solution position.
   */
  virtual void solve () libmesh_override;

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
   * Builds a FEMContext object with enough information to do
   * evaluations on each element.
   *
   * For most problems, the default FEMSystem implementation is correct; users
   * who subclass FEMContext will need to also reimplement this method to build
   * it.
   */
  virtual UniquePtr<DiffContext> build_context() libmesh_override;

  /*
   * Prepares the result of a build_context() call for use.
   *
   * Most FEMSystem-based problems will need to reimplement this in order to
   * call FE::get_*() as their particular physics requires.
   */
  virtual void init_context(DiffContext &) libmesh_override;

  /**
   * Runs a postprocessing loop over all elements, and if
   * \p postprocess_sides is true over all sides.
   */
  virtual void postprocess () libmesh_override;

  /**
   * Runs a qoi assembly loop over all elements, and if
   * \p assemble_qoi_sides is true over all sides.
   *
   * Users may have to override this function if they have any
   * quantities of interest that are not expressible as a sum of
   * element qois.
   */
  virtual void assemble_qoi (const QoISet & indices = QoISet()) libmesh_override;

  /**
   * Runs a qoi derivative assembly loop over all elements, and if
   * \p assemble_qoi_sides is true over all sides.
   *
   * Users may have to override this function for quantities of
   * interest that are not expressible as a sum of element qois.
   */
  virtual void assemble_qoi_derivative (const QoISet & qoi_indices = QoISet(),
                                        bool include_liftfunc = true,
                                        bool apply_constraints = true) libmesh_override;

  /**
   * If fe_reinit_during_postprocess is true (it is true by default), FE
   * objects will be reinit()ed with their default quadrature rules.  If false,
   * FE objects will need to be reinit()ed by the user or will be in an
   * undefined state.
   */
  bool fe_reinit_during_postprocess;

  /**
   * If calculating numeric jacobians is required, the FEMSystem
   * will perturb each solution vector entry by numerical_jacobian_h
   * when calculating finite differences.  This defaults to the
   * libMesh TOLERANCE but can be set manually.
   *
   * For ALE terms, the FEMSystem will perturb each mesh point in an
   * element by numerical_jacobian_h * Elem::hmin()
   */
  Real numerical_jacobian_h;

  /**
   * If numerical_jacobian_h_for_var(var_num) is changed from its
   * default value (numerical_jacobian_h), the FEMSystem will perturb
   * solution vector entries for variable var_num by that amount when
   * calculating finite differences with respect to that variable.
   *
   * This is useful in multiphysics problems which have not been
   * normalized.
   */
  Real numerical_jacobian_h_for_var(unsigned int var_num) const;

  void set_numerical_jacobian_h_for_var(unsigned int var_num, Real new_h);

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

  /**
   * Syntax sugar to make numerical_jacobian() declaration easier.
   */
  typedef bool (TimeSolver::*TimeSolverResPtr)(bool, DiffContext &);

  /**
   * Uses the results of multiple \p res calls
   * to numerically differentiate the corresponding jacobian.
   */
  void numerical_jacobian (TimeSolverResPtr res,
                           FEMContext & context) const;

  /**
   * Uses the results of multiple element_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element.
   */
  void numerical_elem_jacobian (FEMContext & context) const;

  /**
   * Uses the results of multiple side_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on an element's side.
   */
  void numerical_side_jacobian (FEMContext & context) const;

  /**
   * Uses the results of multiple side_residual() calls
   * to numerically differentiate the corresponding jacobian
   * on nonlocal DoFs.
   */
  void numerical_nonlocal_jacobian (FEMContext & context) const;

protected:
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data () libmesh_override;

private:
  std::vector<Real> _numerical_jacobian_h_for_var;
};

// --------------------------------------------------------------
// FEMSystem inline methods
inline
Real
FEMSystem::numerical_jacobian_h_for_var(unsigned int var_num) const
{
  if ((var_num >= _numerical_jacobian_h_for_var.size()) ||
      _numerical_jacobian_h_for_var[var_num] == Real(0))
    return numerical_jacobian_h;

  return _numerical_jacobian_h_for_var[var_num];
}

inline
void FEMSystem::set_numerical_jacobian_h_for_var(unsigned int var_num,
                                                 Real new_h)
{
  if (_numerical_jacobian_h_for_var.size() <= var_num)
    _numerical_jacobian_h_for_var.resize(var_num+1,Real(0));

  libmesh_assert_greater(new_h, 0);

  _numerical_jacobian_h_for_var[var_num] = new_h;
}

} // namespace libMesh


#endif // LIBMESH_FEM_SYSTEM_H
