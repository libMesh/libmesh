// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_STATIC_CONDENSATION_H
#define LIBMESH_STATIC_CONDENSATION_H

#include "libmesh/libmesh_config.h"

// subvectors currently only work with petsc

#ifdef LIBMESH_HAVE_PETSC

#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/preconditioner.h"

#include <map>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace libMesh
{
class MeshBase;
class DofMap;
class System;
class Elem;
template <typename>
class DenseMatrix;
template <typename>
class LinearSolver;
template <typename>
class SparseMatrix;
template <typename>
class NumericVector;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
typedef Eigen::MatrixXcd EigenMatrix;
typedef Eigen::VectorXcd EigenVector;
#else
typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::VectorXd EigenVector;
#endif

class StaticCondensation : public Preconditioner<Number>
{
public:
  StaticCondensation(const MeshBase & mesh, System & sys, const DofMap & dof_map);
  virtual ~StaticCondensation();

  /**
   * Build the element global to local index maps and size the element matrices
   */
  virtual void init() override;

  virtual void setup() override;

  virtual void apply(const NumericVector<Number> & full_rhs,
                     NumericVector<Number> & full_sol) override;

  virtual void clear() override;

  void add_matrix(const Elem & elem,
                  const unsigned int i_var,
                  const unsigned int j_var,
                  const DenseMatrix<Number> & k);

  const SparseMatrix<Number> & get_condensed_mat() const;

  /**
   * Add \p vars to the list of variables not to condense. This can be useful when some variable's
   * equation is discretized with a DG method or if including the variable in the condensed block
   * diagonal would result in it being singular
   */
  void dont_condense_vars(const std::set<unsigned int> & vars);

  /**
   * @returns our list of variables for whom we do not condense out any dofs
   */
  const std::set<unsigned int> & uncondensed_vars() const { return _uncondensed_vars; }

private:
  static void set_local_vectors(const NumericVector<Number> & global_vector,
                                const std::vector<dof_id_type> & elem_dof_indices,
                                std::vector<Number> & elem_dof_values_vec,
                                EigenVector & elem_dof_values);

  void forward_elimination(const NumericVector<Number> & full_rhs);
  void backwards_substitution(const NumericVector<Number> & full_rhs,
                              NumericVector<Number> & full_sol);

  static auto
  total_and_condensed_from_scalar_dofs_functor(std::vector<dof_id_type> & elem_condensed_dofs);

  auto total_and_condensed_from_field_dofs_functor(std::vector<dof_id_type> & elem_condensed_dofs);

  struct LocalData
  {
    /// condensed-condensed matrix entries
    EigenMatrix Acc;
    /// condensed-uncondensed matrix entries
    EigenMatrix Acu;
    /// uncondensed-condensed matrix entries
    EigenMatrix Auc;
    /// uncondensed-uncondensed matrix entries
    EigenMatrix Auu;

    // Acc LU decompositions
    typename std::remove_const<decltype(Acc.partialPivLu())>::type AccFactor;

    struct VarData
    {
      unsigned int uncondensed_dofs_offset;
      unsigned int condensed_dofs_offset;
      unsigned int num_uncondensed_dofs;
      unsigned int num_condensed_dofs;
    };

    /// The uncondensed degrees of freedom with global numbering corresponding to the the \emph reduced
    /// system. Note that initially this will actually hold the indices corresponding to the fully
    /// sized problem, but we will swap it out by the time we are done initializing
    std::vector<dof_id_type> reduced_space_indices;

    std::map<unsigned int, VarData> var_to_data;
  };

  std::map<dof_id_type, LocalData> _elem_to_local_data;
  std::vector<dof_id_type> _local_uncondensed_dofs;

  const MeshBase & _mesh;
  const System & _sys;
  const DofMap & _dof_map;
  std::unique_ptr<SparseMatrix<Number>> _reduced_sys_mat;
  std::unique_ptr<NumericVector<Number>> _reduced_sol;
  std::unique_ptr<NumericVector<Number>> _reduced_rhs;
  std::unique_ptr<LinearSolver<Number>> _reduced_solver;

  /// Variables for which we will keep all dofs
  std::set<unsigned int> _uncondensed_vars;
};

inline const SparseMatrix<Number> &
StaticCondensation::get_condensed_mat() const
{
  libmesh_assert(_reduced_sys_mat);
  return *_reduced_sys_mat;
}

inline void
StaticCondensation::dont_condense_vars(const std::set<unsigned int> & vars)
{
  _uncondensed_vars.insert(vars.begin(), vars.end());
}

}

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_STATIC_CONDENSATION_H
