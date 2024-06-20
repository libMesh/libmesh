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
template <typename>
class NumericVector;
template <typename>
class SparseMatrix;

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

  void assemble_reduced_mat();
  void solve(const NumericVector<Number> & full_rhs);

private:
  void computeElemDofsScalar(const Elem & elem,
                             const std::vector<dof_id_type> & scalar_dof_indices,
                             std::vector<dof_id_type> & elem_dof_indices,
                             std::vector<dof_id_type> & elem_interior_dofs,
                             std::vector<dof_id_type> & elem_trace_dofs) const;

  void computeElemDofsField(const Elem & elem,
                            const unsigned int node_num,
                            const dof_id_type field_dof,
                            std::vector<dof_id_type> & elem_dof_indices,
                            std::vector<dof_id_type> & elem_interior_dofs,
                            std::vector<dof_id_type> & elem_trace_dofs) const;

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
};
}

#endif // LIBMESH_STATIC_CONDENSATION_H
