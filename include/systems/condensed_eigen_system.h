// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_CONDENSED_EIGEN_SYSTEM_H
#define LIBMESH_CONDENSED_EIGEN_SYSTEM_H

#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// Local Includes
#include "libmesh/eigen_system.h"
#include "libmesh/sparse_matrix.h"

namespace libMesh
{

/**
 * This class extends EigenSystem to allow a simple way of solving
 * (standard or generalized) eigenvalue problems in the case where
 * we want to remove certain degrees of freedom from the system.
 * This is useful, for example, in the case that one wants to solve
 * eigenvalue problems with Dirichlet boundary conditions.
 *
 * \author David Knezevic
 * \date 2011
 * \brief Extends EigenSystem to allow certain DOFs to be condensed out.
 */
class CondensedEigenSystem : public EigenSystem
{
public:

  /**
   * Constructor.
   */
  CondensedEigenSystem (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in);

  /**
   * Special functions.
   * - This class has the same restrictions/defaults as its base class.
   * - The destructor is defaulted out-of-line.
   */
  CondensedEigenSystem (const CondensedEigenSystem &) = delete;
  CondensedEigenSystem & operator= (const CondensedEigenSystem &) = delete;
  CondensedEigenSystem (CondensedEigenSystem &&) = default;
  CondensedEigenSystem & operator= (CondensedEigenSystem &&) = delete;
  virtual ~CondensedEigenSystem ();

  /**
   * The type of system.
   */
  typedef CondensedEigenSystem sys_type;

  /**
   * The type of the parent
   */
  typedef EigenSystem Parent;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Loop over the dofs on each processor to initialize the list
   * of non-condensed dofs. These are the dofs in the system that
   * are not contained in \p global_dirichlet_dofs_set and are not
   * subject to constraints due to adaptive mesh hanging nodes,
   * periodic boundary conditions, or Dirichlet boundary conditions.
   *
   * Most users will not need to use the \p global_condensed_dofs_set
   * argument; simply call initialize_condensed_dofs() after any time
   * the EquationSystems (and therefore its constraint equations) gets
   * initialized or reinitialized.
   */
  void initialize_condensed_dofs(const std::set<dof_id_type> &
                                 global_condensed_dofs_set =
                                 std::set<dof_id_type>());

  /**
   * \returns The global number of non-condensed dofs in the system.
   */
  dof_id_type n_global_non_condensed_dofs() const;

  /**
   * Override to solve the condensed eigenproblem with
   * the dofs in local_non_condensed_dofs_vector
   * stripped out of the system matrices on each processor.
   */
  virtual void solve() override;

  virtual void clear () override;

  /**
   * Override \p get_eigenpair() to retrieve the eigenpair for
   * the condensed eigensolve. We only set the non-condensed
   * entries of the solution vector (the condensed
   * entries are set to zero by default).
   */
  virtual std::pair<Real, Real> get_eigenpair(dof_id_type i) override;

  bool has_condensed_matrix_A() const { return _condensed_matrix_A.get(); }
  bool has_condensed_matrix_B() const { return _condensed_matrix_B.get(); }
  bool has_condensed_precond_matrix() const { return _condensed_precond_matrix.get(); }

  /**
   * \returns The system matrix used for standard eigenvalue problems
   */
  SparseMatrix<Number> & get_condensed_matrix_A();

  /**
   * \returns The second system matrix used for generalized eigenvalue problems.
   */
  SparseMatrix<Number> & get_condensed_matrix_B();

  /**
   * \returns The condensed preconditioning matrix
   */
  SparseMatrix<Number> & get_condensed_precond_matrix();

  /**
   * Copy a logically sub-vector into a super-vector. The sub-vector local size should correspond to
   * the size of our \p local_non_condensed_dofs_vector, e.g. the sub-vector should represent
   * non-condensed degree of freedom. The \p super should contain both condensed and non-condensed
   * dofs
   */
  void copy_sub_to_super(const NumericVector<Number> & sub, NumericVector<Number> & super);

  /**
   * Copy a logically super-vector into a sub-vector. The \p super vector should represent all dofs,
   * both condensed and non-condensed. The \p sub should contain only non-condensed dofs, e.g. its
   * local size should match the size of our \p local_non_condensed_dofs_vector
   */
  void copy_super_to_sub(NumericVector<Number> & super, NumericVector<Number> & sub);

  /**
   * Copy a logically super-matrix into a sub-matrix. The \p super matrix should represent all dofs,
   * both condensed and non-condensed. The \p sub should contain only non-condensed dofs, e.g. its
   * local row ownership size should match the size of our \p local_non_condensed_dofs_vector
   */
  void copy_super_to_sub(const SparseMatrix<Number> & super, SparseMatrix<Number> & sub);

  /**
   * Instructs not to create the condensed submatrices from the global matrices right before the
   * solve. This API should be used when assembly occurs during callbacks during the solve, e.g. for
   * a nonlinear eigen problem, solved with a Newton-based solver to determine the fundamental mode,
   * in which function and matrix evalations are liable to change with every nonlinear iteration.
   * Moreover, calling \p create_submatrix wipes away callbacks associated with the condensed matrix
   */
  void dont_create_submatrices_in_solve() { _create_submatrices_in_solve = false; }

#ifdef LIBMESH_ENABLE_DEPRECATED
  /**
   * The (condensed) system matrix for standard eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_condensed_matrix_A() instead.
   */
  SparseMatrix<Number> * condensed_matrix_A;

  /**
   * A second (condensed) system matrix for generalized eigenvalue problems.
   *
   * Public access to this member variable will be deprecated in
   * the future! Use get_condensed_matrix_B() instead.
   */
  SparseMatrix<Number> * condensed_matrix_B;
#endif

  /**
   * Vector storing the local dof indices that will not be condensed.
   * All dofs that are not in this vector will be eliminated from
   * the system when we perform a solve.
   */
  std::vector<dof_id_type> local_non_condensed_dofs_vector;

  /**
   * Initializes the condensed matrices. This API should be used if the condensed matrices will be
   * assembled during the solve as opposed to pre-solve
   */
  void initialize_condensed_matrices();

  /**
   * @returns Whether there are condensed degrees of freedom
   */
  bool have_condensed_dofs() const
  { libmesh_assert(_condensed_dofs_initialized); return _have_condensed_dofs; }

  virtual void reinit() override;

protected:
  virtual void add_matrices () override;

private:
  virtual bool condense_constrained_dofs() const override final { return true; }

#ifdef LIBMESH_ENABLE_DEPRECATED
  void set_raw_pointers();
#endif

  /**
   * A private flag to indicate whether the condensed dofs
   * have been initialized.
   */
  bool _condensed_dofs_initialized;

  /**
   * Whether there are any condensed degrees of freedom
   */
  bool _have_condensed_dofs;

  /**
   * Denotes whether to create the condensed submatrices from the global matrices in the solve
   */
  bool _create_submatrices_in_solve;

  /**
   * The (condensed) system matrix for standard eigenvalue problems.
   */
  std::unique_ptr<SparseMatrix<Number>> _condensed_matrix_A;

  /**
   * A second (condensed) system matrix for generalized eigenvalue problems.
   */
  std::unique_ptr<SparseMatrix<Number>> _condensed_matrix_B;

  /**
   * The condensed preconditioning matrix
   */
  std::unique_ptr<SparseMatrix<Number>> _condensed_precond_matrix;
};


} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC

#endif // LIBMESH_CONDENSED_EIGEN_SYSTEM_H
