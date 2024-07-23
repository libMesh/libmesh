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
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"

#include <unordered_map>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace libMesh
{
class MeshBase;
class System;
class DofMap;
class Elem;
template <typename>
class LinearSolver;
template <typename>
class NumericVector;
template <typename>
class Preconditioner;
class StaticCondensationPreconditioner;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
typedef Eigen::MatrixXcd EigenMatrix;
typedef Eigen::VectorXcd EigenVector;
#else
typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::VectorXd EigenVector;
#endif

class StaticCondensation : public SparseMatrix<Number>
{
public:
  StaticCondensation(const MeshBase & mesh, const System & system, const DofMap & dof_map);
  virtual ~StaticCondensation();

  //
  // SparseMatrix overrides
  //

  virtual SparseMatrix<Number> & operator=(const SparseMatrix<Number> &) override;

  virtual SolverPackage solver_package() override;

  virtual void init(const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type = 30,
                    const numeric_index_type = 10,
                    const numeric_index_type = 1) override;

  virtual void init(ParallelType) override { this->init(); }

  virtual void clear() override;

  virtual void zero() override;

  virtual std::unique_ptr<SparseMatrix<Number>> zero_clone() const override;

  virtual std::unique_ptr<SparseMatrix<Number>> clone() const override;

  virtual void close() override;

  virtual numeric_index_type m() const override;

  virtual numeric_index_type n() const override { return this->m(); }

  virtual numeric_index_type row_start() const override;

  virtual numeric_index_type row_stop() const override;

  virtual numeric_index_type col_start() const override { return this->row_start(); }

  virtual numeric_index_type col_stop() const override { return this->row_stop(); }

  virtual void
  set(const numeric_index_type i, const numeric_index_type j, const Number value) override;

  virtual void
  add(const numeric_index_type i, const numeric_index_type j, const Number value) override;

  virtual void add_matrix(const DenseMatrix<Number> & dm,
                          const std::vector<numeric_index_type> & rows,
                          const std::vector<numeric_index_type> & cols) override;

  virtual void add_matrix(const DenseMatrix<Number> & dm,
                          const std::vector<numeric_index_type> & dof_indices) override;

  virtual void add(const Number a, const SparseMatrix<Number> & X) override;

  virtual Number operator()(const numeric_index_type i, const numeric_index_type j) const override;

  virtual Real l1_norm() const override;

  virtual Real linfty_norm() const override;

  virtual bool closed() const override;

  virtual void print_personal(std::ostream & os = libMesh::out) const override;

  virtual void get_diagonal(NumericVector<Number> & dest) const override;

  virtual void get_transpose(SparseMatrix<Number> & dest) const override;

  virtual void get_row(numeric_index_type i,
                       std::vector<numeric_index_type> & indices,
                       std::vector<Number> & values) const override;

  //
  // Unique methods
  //

  /**
   * Build the element global to local index maps and size the element matrices
   */
  void init();

  /**
   * A no-op to be consistent with shimming from the \p StaticCondenstionPreconditioner. "setup"
   * should take place at the end of matrix assembly, e.g. during a call to \p close()
   */
  void setup();

  /**
   * Perform our three stages, \p forward_elimination(), a \p solve() on the condensed system, and
   * finally a \p backwards_substitution()
   */
  void apply(const NumericVector<Number> & full_rhs, NumericVector<Number> & full_sol);

  const SparseMatrix<Number> & get_condensed_mat() const;

  /**
   * Add \p vars to the list of variables not to condense. This can be useful when some variable's
   * equation is discretized with a DG method or if including the variable in the condensed block
   * diagonal would result in it being singular
   */
  void dont_condense_vars(const std::unordered_set<unsigned int> & vars);

  /**
   * @returns our list of variables for whom we do not condense out any dofs
   */
  const std::unordered_set<unsigned int> & uncondensed_vars() const { return _uncondensed_vars; }

  /**
   * Set the current element. This enables fast lookups of local indices from global indices
   */
  void set_current_elem(const Elem & elem);

  /**
   * Get the preconditioning wrapper
   */
  StaticCondensationPreconditioner & get_preconditioner() { return *_scp; }

private:
  static void set_local_vectors(const NumericVector<Number> & global_vector,
                                const std::vector<dof_id_type> & elem_dof_indices,
                                std::vector<Number> & elem_dof_values_vec,
                                EigenVector & elem_dof_values);

  void forward_elimination(const NumericVector<Number> & full_rhs);
  void backwards_substitution(const NumericVector<Number> & full_rhs,
                              NumericVector<Number> & full_sol);

  static void total_dofs_from_scalar_dofs(std::vector<dof_id_type> & dofs,
                                          const std::vector<dof_id_type> & scalar_dofs);

  static void condensed_dofs_from_scalar_dofs(std::vector<dof_id_type> & /*condensed_dofs*/,
                                              const std::vector<dof_id_type> & /*scalar_dofs*/);

  static void uncondensed_dofs_from_scalar_dofs(std::vector<dof_id_type> & uncondensed_dofs,
                                                const std::vector<dof_id_type> & scalar_dofs);

  static void total_dofs_from_field_dof(std::vector<dof_id_type> & dofs,
                                        const Elem & /*elem*/,
                                        const unsigned int /*node_num*/,
                                        const unsigned int /*var_num*/,
                                        const dof_id_type field_dof);

  void condensed_dofs_from_field_dof(std::vector<dof_id_type> & condensed_dofs,
                                     const Elem & elem,
                                     const unsigned int node_num,
                                     const unsigned int var_num,
                                     const dof_id_type field_dof) const;

  void uncondensed_dofs_from_field_dof(std::vector<dof_id_type> & uncondensed_dofs,
                                       const Elem & elem,
                                       const unsigned int node_num,
                                       const unsigned int var_num,
                                       const dof_id_type field_dof) const;

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

    /// The uncondensed degrees of freedom with global numbering corresponding to the the \emph reduced
    /// system. Note that initially this will actually hold the indices corresponding to the fully
    /// sized problem, but we will swap it out by the time we are done initializing
    std::vector<dof_id_type> reduced_space_indices;

    std::unordered_map<dof_id_type, dof_id_type> uncondensed_global_to_local_map;
    std::unordered_map<dof_id_type, dof_id_type> condensed_global_to_local_map;
  };

  std::unordered_map<dof_id_type, LocalData> _elem_to_local_data;
  std::vector<dof_id_type> _local_uncondensed_dofs;

  const MeshBase & _mesh;
  const System & _system;
  const DofMap & _dof_map;
  std::unique_ptr<SparseMatrix<Number>> _reduced_sys_mat;
  std::unique_ptr<NumericVector<Number>> _reduced_sol;
  std::unique_ptr<NumericVector<Number>> _reduced_rhs;
  std::unique_ptr<LinearSolver<Number>> _reduced_solver;
  std::unique_ptr<NumericVector<Number>> _ghosted_full_sol;

  /// Variables for which we will keep all dofs
  std::unordered_set<unsigned int> _uncondensed_vars;

  /// The current element ID. This is one half of a key, along with the global index, that maps to a
  /// local index
  dof_id_type _current_elem_id;

  /// Helper data member for adding individual matrix elements
  DenseMatrix<Number> _size_one_mat;

  /// Preconditioner object which will call back to us for the preconditioning action
  std::unique_ptr<StaticCondensationPreconditioner> _scp;
};

inline const SparseMatrix<Number> &
StaticCondensation::get_condensed_mat() const
{
  libmesh_assert(_reduced_sys_mat);
  return *_reduced_sys_mat;
}

inline void
StaticCondensation::dont_condense_vars(const std::unordered_set<unsigned int> & vars)
{
  _uncondensed_vars.insert(vars.begin(), vars.end());
}

}

#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_STATIC_CONDENSATION_H
