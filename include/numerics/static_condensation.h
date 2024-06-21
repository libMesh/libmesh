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

#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
// These next three includes are not strictly necessary but unique_ptr will error with incomplete
// type if user code does not have these includes, which could be confusing for the user because
// they're only making calls on the StaticCondensation class
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"
#include <unordered_map>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace libMesh
{
class MeshBase;
class DofMap;
class Elem;
template <typename>
class DenseMatrix;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
typedef Eigen::MatrixXcd EigenMatrix;
typedef Eigen::VectorXcd EigenVector;
#else
typedef Eigen::MatrixXd EigenMatrix;
typedef Eigen::VectorXd EigenVector;
#endif

class StaticCondensation
{
public:
  StaticCondensation(const MeshBase & mesh, const DofMap & dof_map);

  /**
   * Build the element global to local index maps and size the element matrices
   */
  void init();

  void add_matrix(const Elem & elem,
                  const unsigned int i_var,
                  const unsigned int j_var,
                  const DenseMatrix<Number> & k);

  void solve(const NumericVector<Number> & full_rhs, NumericVector<Number> & full_sol);

private:
  static void set_local_vectors(const NumericVector<Number> & global_vector,
                                const std::vector<dof_id_type> & elem_dof_indices,
                                std::vector<Number> & elem_dof_values_vec,
                                EigenVector & elem_dof_values);

  void assemble_reduced_mat();
  void forward_elimination(const NumericVector<Number> & full_rhs);
  void backwards_substitution(const NumericVector<Number> & full_rhs,
                              NumericVector<Number> & full_sol);

  static auto computeElemDofsScalar(std::vector<dof_id_type> & elem_interior_dofs,
                                    std::vector<dof_id_type> & elem_trace_dofs) -> decltype(auto);

  static auto computeElemDofsField(std::vector<dof_id_type> & elem_interior_dofs,
                                   std::vector<dof_id_type> & elem_trace_dofs) -> decltype(auto);

  struct LocalData
  {
    /// interior-interior matrix entries
    EigenMatrix Aii;
    /// interior-boundary matrix entries
    EigenMatrix Aib;
    /// boundary-interior matrix entries
    EigenMatrix Abi;
    /// boundary-boundary matrix entries
    EigenMatrix Abb;

    // Aii LU decompositions
    typename std::remove_const<decltype(Aii.partialPivLu())>::type AiiFactor;

    struct VarData
    {
      unsigned int boundary_dofs_offset;
      unsigned int interior_dofs_offset;
      unsigned int num_boundary_dofs;
      unsigned int num_interior_dofs;
    };

    /// The trace degrees of freedom with global numbering corresponding to the the \emph reduced
    /// system. Note that initially this will actually hold the indices corresponding to the fully
    /// sized problem, but we will swap it out by the time we are done initializing
    std::vector<dof_id_type> reduced_space_indices;

    std::unordered_map<unsigned int, VarData> var_to_data;
  };

  std::unordered_map<dof_id_type, LocalData> _elem_to_local_data;

  const MeshBase & _mesh;
  const DofMap & _dof_map;
  std::unique_ptr<SparseMatrix<Number>> _reduced_sys_mat;
  std::unique_ptr<NumericVector<Number>> _reduced_sol;
  std::unique_ptr<NumericVector<Number>> _reduced_rhs;
  std::unique_ptr<LinearSolver<Number>> _reduced_solver;
};
}

#endif // LIBMESH_STATIC_CONDENSATION_H
