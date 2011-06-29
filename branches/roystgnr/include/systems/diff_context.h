
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



#ifndef __diff_context_h__
#define __diff_context_h__

// C++ includes

// Local Includes
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"

namespace libMesh
{

// Forward declarations
class System;

/**
 * This class provides all data required for a physics package
 * (e.g. a DifferentiableSystem subclass) to perform local element
 * residual and jacobian integrations.
 *
 * This class is part of the new DifferentiableSystem framework,
 * which is still experimental.  Users of this framework should
 * beware of bugs and future API changes.
 *
 * @author Roy H. Stogner 2009
 */

// ------------------------------------------------------------
// DifferentiableSystem class definition

class DiffContext
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  DiffContext (const System &);

  /**
   * Destructor.
   */
  virtual ~DiffContext ();

  /**
   * If \p postprocess_sides is true (it is false by default), the
   * postprocessing loop will loop over all sides as well as all elements.
   */
  bool postprocess_sides;

  /**
   * Gives derived classes the opportunity to reinitialize data (FE objects in
   * FEMSystem, for example) needed for an interior integration at a new point
   * within a timestep
   */
  virtual void elem_reinit(Real) {}

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for a side integration at a new point within a timestep
   */
  virtual void elem_side_reinit(Real) {}

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for an edge integration at a new point within a timestep
   */
  virtual void elem_edge_reinit(Real) {}

  /**
   * For time-dependent problems, this is the time t for which the current
   * nonlinear_solution is defined.
   * FIXME - this needs to be tweaked mid-timestep by all transient solvers!
   */
  Real time;

  /**
   * This is the time stored in the System class at the time this context
   * was created, i.e. the time at the beginning of the current timestep.
   * This value gets set in the constructor and unlike DiffContext::time,
   * is not tweaked mid-timestep by transient solvers: it remains equal
   * to the value it was assigned at construction.
   */
  const Real system_time;

  /**
   * Element by element components of nonlinear_solution
   * as adjusted by a time_solver
   */
  DenseVector<Number> elem_solution;
  std::vector<DenseSubVector<Number> *> elem_subsolutions;

  /**
   * Element by element components of nonlinear_solution
   * at a fixed point in a timestep, for optional use by e.g.
   * stabilized methods
   */
  DenseVector<Number> elem_fixed_solution;
  std::vector<DenseSubVector<Number> *> elem_fixed_subsolutions;

  /**
   * The derivative of elem_solution with respect to the nonlinear solution,
   * for use by systems constructing jacobians with elem_fixed_solution
   * based methods
   */
  Real elem_solution_derivative;

  /**
   * The derivative of elem_fixed_solution with respect to the nonlinear
   * solution, for use by systems constructing jacobians with
   * elem_fixed_solution based methods
   */
  Real fixed_solution_derivative;

  /**
   * Element residual vector
   */
  DenseVector<Number> elem_residual;

  /**
   * Element jacobian: derivatives of elem_residual with respect to
   * elem_solution
   */
  DenseMatrix<Number> elem_jacobian;

  /**
   * Element quantity of interest contributions
   */
  std::vector<Number> elem_qoi;

  /**
   * Element quantity of interest derivative contributions
   */
  std::vector<DenseVector<Number> > elem_qoi_derivative;
  std::vector<std::vector<DenseSubVector<Number> *> > elem_qoi_subderivatives;

  /**
   * Element residual subvectors and Jacobian submatrices
   */
  std::vector<DenseSubVector<Number> *> elem_subresiduals;
  std::vector<std::vector<DenseSubMatrix<Number> *> > elem_subjacobians;

  /** 
   * Global Degree of freedom index lists
   */
  std::vector<unsigned int> dof_indices;
  std::vector<std::vector<unsigned int> > dof_indices_var;

  /**
   * Points the _deltat member of this class at a timestep value
   * stored in the creating System, for example DiffSystem::deltat
   */
  void set_deltat_pointer(Real* dt);

  /**
   * Returns the value currently pointed to by this class's _deltat
   * member
   */
  Real get_deltat_value();
  
private:
  /**
   * Default NULL, can optionally be used to point to a timestep value
   * in the System-derived class responsible for creating this DiffContext.
   *
   * In DiffSystem's build_context() function, is assigned to point to
   * the deltat member of that class.
   *
   * Accessible via public get_deltat()/set_deltat() methods of this class.
   *
   * Always test for NULL before using!
   */
  Real* _deltat;
};

} // namespace libMesh


#endif
