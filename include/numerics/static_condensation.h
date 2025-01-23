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

#ifndef LIBMESH_STATIC_CONDENSATION_H
#define LIBMESH_STATIC_CONDENSATION_H

#include "libmesh/libmesh_config.h"

// subvectors currently only work with petsc and we rely on Eigen for local LU factorizations
#if defined(LIBMESH_HAVE_EIGEN) && defined(LIBMESH_HAVE_PETSC)

#include "libmesh/petsc_matrix_shell_matrix.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/dense_matrix.h"

#include <unordered_map>
#include <memory>
#include <vector>

// Warnings about stack protection
#include "libmesh/ignore_warnings.h"
#include <Eigen/Dense>
#include "libmesh/restore_warnings.h"

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

class StaticCondensation : public PetscMatrixShellMatrix<Number>
{
public:
  StaticCondensation(const MeshBase & mesh, const System & system, const DofMap & dof_map);
  virtual ~StaticCondensation();

  //
  // SparseMatrix overrides
  //

  virtual SparseMatrix<Number> & operator=(const SparseMatrix<Number> &) override;

  virtual SolverPackage solver_package() override;

  virtual void init(const numeric_index_type m,
                    const numeric_index_type n,
                    const numeric_index_type m_l,
                    const numeric_index_type n_l,
                    const numeric_index_type nnz = 30,
                    const numeric_index_type noz = 10,
                    const numeric_index_type blocksize = 1) override;

  virtual void init(const ParallelType type) override;

  virtual bool initialized() const override;

  virtual void clear() noexcept override;

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
  /**
   * Retrieves the degree of freedom values from \p global_vector corresponding to \p
   * elem_dof_indices, filling both \p elem_dof_values_vec and \p elem_dof_values
   */
  static void set_local_vectors(const NumericVector<Number> & global_vector,
                                const std::vector<dof_id_type> & elem_dof_indices,
                                std::vector<Number> & elem_dof_values_vec,
                                EigenVector & elem_dof_values);

  /**
   * Takes an incoming "full" RHS from the solver, where full means the union of uncondensed and
   * condennsed dofs, and condenses it down into the condensed \p _reduced_rhs data member using
   * element Schur complements
   */
  void forward_elimination(const NumericVector<Number> & full_rhs);

  /**
   * After performing the solve with the \p _reduced_rhs and Schur complement matrix (\p
   * _reduced_sys_mat) to determine the \p _reduced_sol, we use the uncondensed \p full_rhs along
   * with the \p _reduced_sol to back substitute into the \p full_sol, which is the final output
   * data of the static condensation preconditioner application
   */
  void backwards_substitution(const NumericVector<Number> & full_rhs,
                              NumericVector<Number> & full_sol);

  /**
   * Data stored on a per-element basis used to compute element Schur complements and their
   * applications to vectors
   */
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

    /// A map from the global degree of freedom number for the full system (condensed + uncondensed)
    /// to an element local number. If this map is queried with a condensed dof, nothing will be
    /// found. The size of this container will be the number of uncondensed degrees of freedom whose
    /// basis functions are nonzero on the element
    std::unordered_map<dof_id_type, dof_id_type> uncondensed_global_to_local_map;
    /// A map from the global degree of freedom number for the full system (condensed + uncondensed)
    /// to an element local number. If this map is queried with an uncondensed dof, nothing will be
    /// found. The size of this container will be the number of condensed degrees of freedom whose
    /// basis functions are nonzero on the element
    std::unordered_map<dof_id_type, dof_id_type> condensed_global_to_local_map;
  };

  /// A map from element ID to Schur complement data
  std::unordered_map<dof_id_type, LocalData> _elem_to_local_data;

  /// All the uncondensed degrees of freedom (numbered in the "full" uncondensed + condensed
  /// space). This data member is used for creating subvectors corresponding to only uncondensed
  /// dofs
  std::vector<dof_id_type> _local_uncondensed_dofs;

  const MeshBase & _mesh;
  const System & _system;
  const DofMap & _dof_map;

  /// global sparse matrix for the uncondensed degrees of freedom
  std::unique_ptr<SparseMatrix<Number>> _reduced_sys_mat;
  /// solution for the uncondensed degrees of freedom
  std::unique_ptr<NumericVector<Number>> _reduced_sol;
  /// RHS corresponding to the uncondensed degrees of freedom. This is constructed by applying the
  /// element Schur complements to an incoming RHS which contains both condensed and uncondensed
  /// degrees of freedom
  std::unique_ptr<NumericVector<Number>> _reduced_rhs;
  /// The solver for the uncondensed degrees of freedom
  std::unique_ptr<LinearSolver<Number>> _reduced_solver;
  /// This is a ghosted representation of the full (uncondensed + condensed) solution. Note that
  // this is, in general, *not equal* to the system solution, e.g. this may correspond to the
  // solution for the Newton *update*
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

  /// Whether our object has been initialized
  bool _sc_is_initialized;

  /// Whether we have cached values via add_XXX()
  bool _have_cached_values;
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

#else

#include "libmesh/sparse_matrix.h"
#include <unordered_set>

namespace libMesh
{
class MeshBase;
class System;
class DofMap;
class StaticCondensationPreconditioner;

class StaticCondensation : public SparseMatrix<Number>
{
public:
  StaticCondensation(const MeshBase &, const System &, const DofMap & dof_map);

  const std::unordered_set<unsigned int> & uncondensed_vars() const { libmesh_not_implemented(); }
  StaticCondensationPreconditioner & get_preconditioner() { libmesh_not_implemented(); }
  virtual SparseMatrix<Number> & operator=(const SparseMatrix<Number> &) override
  {
    libmesh_not_implemented();
  }
  virtual SolverPackage solver_package() override { libmesh_not_implemented(); }
  virtual void init(const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type,
                    const numeric_index_type = 30,
                    const numeric_index_type = 10,
                    const numeric_index_type = 1) override
  {
    libmesh_not_implemented();
  }
  virtual void init(ParallelType) override { libmesh_not_implemented(); }
  virtual void clear() override { libmesh_not_implemented(); }
  virtual void zero() override { libmesh_not_implemented(); }
  virtual std::unique_ptr<SparseMatrix<Number>> zero_clone() const override
  {
    libmesh_not_implemented();
  }
  virtual std::unique_ptr<SparseMatrix<Number>> clone() const override
  {
    libmesh_not_implemented();
  }
  virtual void close() override { libmesh_not_implemented(); }
  virtual numeric_index_type m() const override { libmesh_not_implemented(); }
  virtual numeric_index_type n() const override { libmesh_not_implemented(); }
  virtual numeric_index_type row_start() const override { libmesh_not_implemented(); }
  virtual numeric_index_type row_stop() const override { libmesh_not_implemented(); }
  virtual numeric_index_type col_start() const override { libmesh_not_implemented(); }
  virtual numeric_index_type col_stop() const override { libmesh_not_implemented(); }
  virtual void set(const numeric_index_type, const numeric_index_type, const Number) override
  {
    libmesh_not_implemented();
  }
  virtual void add(const numeric_index_type, const numeric_index_type, const Number) override
  {
    libmesh_not_implemented();
  }
  virtual void add_matrix(const DenseMatrix<Number> &,
                          const std::vector<numeric_index_type> &,
                          const std::vector<numeric_index_type> &) override
  {
    libmesh_not_implemented();
  }
  virtual void add_matrix(const DenseMatrix<Number> &,
                          const std::vector<numeric_index_type> &) override
  {
    libmesh_not_implemented();
  }
  virtual void add(const Number, const SparseMatrix<Number> &) override
  {
    libmesh_not_implemented();
  }
  virtual Number operator()(const numeric_index_type, const numeric_index_type) const override
  {
    libmesh_not_implemented();
  }
  virtual Real l1_norm() const override { libmesh_not_implemented(); }
  virtual Real linfty_norm() const override { libmesh_not_implemented(); }
  virtual bool closed() const override { libmesh_not_implemented(); }
  virtual void print_personal(std::ostream & = libMesh::out) const override
  {
    libmesh_not_implemented();
  }
  virtual void get_diagonal(NumericVector<Number> &) const override { libmesh_not_implemented(); }
  virtual void get_transpose(SparseMatrix<Number> &) const override { libmesh_not_implemented(); }
  virtual void get_row(numeric_index_type,
                       std::vector<numeric_index_type> &,
                       std::vector<Number> &) const override
  {
    libmesh_not_implemented();
  }
  void init() { libmesh_not_implemented(); }
  void setup() { libmesh_not_implemented(); }
  void apply(const NumericVector<Number> &, NumericVector<Number> &) { libmesh_not_implemented(); }
};
}

#endif // LIBMESH_HAVE_EIGEN && LIBMESH_HAVE_PETSC
#endif // LIBMESH_STATIC_CONDENSATION_H
