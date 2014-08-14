// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_DIFF_CONTEXT_H
#define LIBMESH_DIFF_CONTEXT_H

// Local Includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dense_vector.h"
#include "libmesh/id_types.h"

// C++ includes
#include <cstddef>
#include <map>
#include <vector>

namespace libMesh
{

// Forward declarations
template <typename T> class NumericVector;
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
  explicit
  DiffContext (const System &);

  /**
   * Destructor.
   */
  virtual ~DiffContext ();

  /**
   * Gives derived classes the opportunity to reinitialize data (FE objects in
   * FEMSystem, for example) needed for an interior integration at a new point
   * within a timestep
   */
  virtual void elem_reinit(Real /* theta */) {}

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for a side integration at a new point within a timestep
   */
  virtual void elem_side_reinit(Real /* theta */) {}

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for an edge integration at a new point within a timestep
   */
  virtual void elem_edge_reinit(Real /* theta */) {}

  /**
   * Gives derived classes the opportunity to reinitialize data needed
   * for nonlocal calculations at a new point within a timestep
   */
  virtual void nonlocal_reinit(Real /* theta */) {}

  /**
   * Number of variables in solution.
   */
  unsigned int n_vars() const
  { return cast_int<unsigned int>(dof_indices_var.size()); }

  /**
   * Accessor for associated system.
   */
  const System& get_system() const
  { return _system; }

  /**
   * Accessor for element solution.
   */
  const DenseVector<Number>& get_elem_solution() const
  { return elem_solution; }

  /**
   * Non-const accessor for element solution.
   */
  DenseVector<Number>& get_elem_solution()
  { return elem_solution; }

  /**
   * Accessor for element solution of a particular variable corresponding
   * to the variable index argument.
   */
  const DenseSubVector<Number>& get_elem_solution( unsigned int var ) const
  { return *(elem_subsolutions[var]); }

  /**
   * Accessor for element fixed solution.
   */
  const DenseVector<Number>& get_elem_fixed_solution() const
  { return elem_fixed_solution; }

  /**
   * Non-const accessor for element fixed solution.
   */
  DenseVector<Number>& get_elem_fixed_solution()
  { return elem_fixed_solution; }

  /**
   * Accessor for element fixed solution of a particular variable corresponding
   * to the variable index argument.
   */
  const DenseSubVector<Number>& get_elem_fixed_solution( unsigned int var ) const
  { return *(elem_fixed_subsolutions[var]); }

  /**
   * Const accessor for element residual.
   */
  const DenseVector<Number>& get_elem_residual() const
  { return elem_residual; }

  /**
   * Non-const accessor for element residual.
   */
  DenseVector<Number>& get_elem_residual()
  { return elem_residual; }

  /**
   * Const accessor for element residual of a particular variable corresponding
   * to the variable index argument.
   */
  const DenseSubVector<Number>& get_elem_residual( unsigned int var ) const
  { return *(elem_subresiduals[var]); }

  /**
   * Non-const accessor for element residual of a particular variable corresponding
   * to the variable index argument.
   */
  DenseSubVector<Number>& get_elem_residual( unsigned int var )
  { return *(elem_subresiduals[var]); }

  /**
   * Const accessor for element Jacobian.
   */
  const DenseMatrix<Number>& get_elem_jacobian() const
  { return elem_jacobian; }

  /**
   * Non-const accessor for element Jacobian.
   */
  DenseMatrix<Number>& get_elem_jacobian()
  { return elem_jacobian; }

  /**
   * Const accessor for element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  const DenseSubMatrix<Number>& get_elem_jacobian( unsigned int var1, unsigned int var2 ) const
  { return *(elem_subjacobians[var1][var2]); }

  /**
   * Non-const accessor for element Jacobian of particular variables corresponding
   * to the variable index arguments.
   */
  DenseSubMatrix<Number>& get_elem_jacobian( unsigned int var1, unsigned int var2 )
  { return *(elem_subjacobians[var1][var2]); }

  /**
   * Const accessor for QoI vector.
   */
  const std::vector<Number>& get_qois() const
  { return elem_qoi; }

  /**
   * Non-const accessor for QoI vector.
   */
  std::vector<Number>& get_qois()
  { return elem_qoi; }

  /**
   * Const accessor for QoI derivatives.
   */
  const std::vector<DenseVector<Number> > & get_qoi_derivatives() const
  { return elem_qoi_derivative; }

  /**
   * Non-const accessor for QoI derivatives.
   */
  std::vector<DenseVector<Number> > & get_qoi_derivatives()
  { return elem_qoi_derivative; }

  /**
   * Const accessor for QoI derivative of a particular qoi and variable corresponding
   * to the index arguments.
   */
  const DenseSubVector<Number>& get_qoi_derivatives( unsigned int qoi, unsigned int var ) const
  { return *(elem_qoi_subderivatives[qoi][var]); }

  /**
   * Non-const accessor for QoI derivative of a particular qoi and variable corresponding
   * to the index arguments.
   */
  DenseSubVector<Number>& get_qoi_derivatives( unsigned int qoi, unsigned int var )
  { return *(elem_qoi_subderivatives[qoi][var]); }

  /**
   * Accessor for element dof indices
   */
  const std::vector<dof_id_type>& get_dof_indices() const
  { return dof_indices; }

  /**
   * Non-const accessor for element dof indices
   */
  std::vector<dof_id_type>& get_dof_indices()
  { return dof_indices; }

  /**
   * Accessor for element dof indices of a particular variable corresponding
   * to the index argument.
   */
  const std::vector<dof_id_type>& get_dof_indices( unsigned int var ) const
  { return dof_indices_var[var]; }

  /**
   * Accessor for the time variable stored in the system class.
   */
  Real get_system_time() const
  { return system_time; }

  /**
   * Accessor for the time for which the current nonlinear_solution is defined.
   */
  Real get_time() const
  { return time; }

  /**
   * Set the time for which the current nonlinear_solution is defined.
   */
  void set_time( Real time_in )
  { time = time_in; }

  Real get_elem_solution_derivative() const
  { return elem_solution_derivative; }

  Real get_fixed_solution_derivative() const
  { return fixed_solution_derivative; }

  /**
   * Accessor for querying whether we need to do a primal
   * or adjoint solve
   */
  bool is_adjoint() const
  { return _is_adjoint; }

  /**
   * Accessor for setting whether we need to do a primal
   * or adjoint solve
   */
  bool& is_adjoint()
  { return _is_adjoint; }

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
   * Points the _deltat member of this class at a timestep value
   * stored in the creating System, for example DiffSystem::deltat
   */
  void set_deltat_pointer(Real* dt);

  /**
   * Returns the value currently pointed to by this class's _deltat
   * member
   */
  Real get_deltat_value();

  /**
   * Adds a vector to the map of localized vectors. We can later evaluate interior_values,
   * interior_gradients and side_values for these fields these vectors represent.
   */
  void add_localized_vector (NumericVector<Number> & _localized_vector, const System & _sys);

  /**
   * Typedef for the localized_vectors iterator
   */
  typedef std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > >::iterator localized_vectors_iterator;

  /**
   * Return a reference to DenseVector localization of _localized_vector
   * contained in the localized_vectors map
   */
  DenseVector<Number> & get_localized_vector (const NumericVector<Number> & _localized_vector);

  /**
   * const accessible version of get_localized_vector function
   */
  const DenseVector<Number> & get_localized_vector (const NumericVector<Number> & _localized_vector) const;

  /**
   * Return a reference to DenseSubVector localization of _localized_vector at variable _var
   * contained in the localized_vectors map
   */
  DenseSubVector<Number> & get_localized_subvector (const NumericVector<Number> & _localized_vector, unsigned int _var);

  /**
   * const accessible version of get_localized_subvector function
   */
  const DenseSubVector<Number> & get_localized_subvector (const NumericVector<Number> & _localized_vector, unsigned int _var) const;

protected:

  /**
   * Contains pointers to vectors the user has asked to be localized, keyed with
   * pairs of element localized versions of that vector and per variable views
   */

  std::map<const NumericVector<Number>*, std::pair<DenseVector<Number>, std::vector<DenseSubVector<Number>*> > > localized_vectors;

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
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type> > dof_indices_var;

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

  /**
   * A reference to the system this context is constructed with
   */
  const System& _system;

  /**
   * Is this context to be used for a primal or adjoint solve?
   */
  bool _is_adjoint;

};

} // namespace libMesh


#endif // LIBMESH_DIFF_CONTEXT_H
