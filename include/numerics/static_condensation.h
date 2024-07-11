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

#include <unordered_map>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace libMesh
{
class MeshBase;
class ImplicitSystem;
class DofMap;
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
  StaticCondensation(const MeshBase & mesh, ImplicitSystem & system, const DofMap & dof_map);
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
  void dont_condense_vars(const std::unordered_set<unsigned int> & vars);

  /**
   * @returns our list of variables for whom we do not condense out any dofs
   */
  const std::unordered_set<unsigned int> & uncondensed_vars() const { return _uncondensed_vars; }

private:
  static void set_local_vectors(const NumericVector<Number> & global_vector,
                                const std::vector<dof_id_type> & elem_dof_indices,
                                std::vector<Number> & elem_dof_values_vec,
                                EigenVector & elem_dof_values);

  void forward_elimination(NumericVector<Number> & rhs);
  void backwards_substitution(const NumericVector<Number> & rhs, NumericVector<Number> & sol);

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

    struct VarData
    {
      unsigned int uncondensed_dofs_offset;
      unsigned int condensed_dofs_offset;
      unsigned int num_uncondensed_dofs;
      unsigned int num_condensed_dofs;
    };

    std::unordered_map<unsigned int, VarData> var_to_data;
  };

  std::unordered_map<dof_id_type, LocalData> _elem_to_local_data;
  std::vector<dof_id_type> _local_uncondensed_dofs;

  const MeshBase & _mesh;
  ImplicitSystem & _system;
  const DofMap & _dof_map;
  std::unique_ptr<SparseMatrix<Number>> _reduced_sys_mat;
  std::unique_ptr<NumericVector<Number>> _reduced_sol;
  std::unique_ptr<NumericVector<Number>> _reduced_rhs;
  std::unique_ptr<LinearSolver<Number>> _reduced_solver;
  std::unique_ptr<NumericVector<Number>> _eliminated_rhs;
  std::unique_ptr<NumericVector<Number>> _ghosted_full_sol;

  /// Variables for which we will keep all dofs
  std::unordered_set<unsigned int> _uncondensed_vars;
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
