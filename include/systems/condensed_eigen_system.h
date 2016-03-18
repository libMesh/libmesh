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


#ifndef LIBMESH_CONDENSED_EIGEN_SYSTEM_H
#define LIBMESH_CONDENSED_EIGEN_SYSTEM_H

#include "libmesh/libmesh_config.h"

// Currently, the EigenSystem should only be available
// if SLEPc support is enabled.
#if defined(LIBMESH_HAVE_SLEPC)

// Local Includes
#include "libmesh/eigen_system.h"
#include "libmesh/sparse_matrix.h"

// C++ includes

namespace libMesh
{

/**
 * This class extends EigenSystem to allow a simple way of solving
 * (standard or generalized) eigenvalue problems in the case where
 * we want to remove certain degrees of freedom from the system.
 * This is useful, for example, in the case that one wants to solve
 * eigenvalue problems with Dirichlet boundary conditions.
 */
class CondensedEigenSystem : public EigenSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  CondensedEigenSystem (EquationSystems & es,
                        const std::string & name_in,
                        const unsigned int number_in);

  /**
   * The type of system.
   */
  typedef CondensedEigenSystem sys_type;

  /**
   * The type of the parent
   */
  typedef EigenSystem Parent;

  /**
   * @returns a clever pointer to the system.
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
   * @return the global number of non-condensed dofs in the system.
   */
  dof_id_type n_global_non_condensed_dofs() const;

  /**
   * Override to solve the condensed eigenproblem with
   * the dofs in local_non_condensed_dofs_vector
   * stripped out of the system matrices on each processor.
   */
  virtual void solve() libmesh_override;

  /**
   * Overload get_eigenpair to retrieve the eigenpair for
   * the condensed eigensolve. We only set the non-condensed
   * entries of the solution vector (the condensed
   * entries are set to zero by default).
   */
  virtual std::pair<Real, Real> get_eigenpair(dof_id_type i) libmesh_override;

  /**
   * The (condensed) system matrix for standard eigenvalue problems.
   */
  UniquePtr< SparseMatrix<Number> > condensed_matrix_A;

  /**
   * A second (condensed) system matrix for generalized eigenvalue problems.
   */
  UniquePtr< SparseMatrix<Number> > condensed_matrix_B;

  /**
   * Vector storing the local dof indices that will not be condensed.
   * All dofs that are not in this vector will be eliminated from
   * the system when we perform a solve.
   */
  std::vector<dof_id_type> local_non_condensed_dofs_vector;

private:

  /**
   * A private flag to indicate whether the condensed dofs
   * have been initialized.
   */
  bool condensed_dofs_initialized;
};


} // namespace libMesh


#endif // LIBMESH_HAVE_SLEPC

#endif // LIBMESH_CONDENSED_EIGEN_SYSTEM_H
