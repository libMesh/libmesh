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
#include <unordered_map>

namespace libMesh
{
class MeshBase;
class DofMap;
class Elem;
template <typename>
class DenseMatrix;

template <typename T>
class StaticCondensation
{
public:
  StaticCondensation(const MeshBase & mesh, const DofMap & dof_map);

  /**
   * Build the element global to local index maps
   */
  void build_idx_maps();

  void add_matrix(const Elem & elem,
                  const unsigned int i_var,
                  const unsigned int j_var,
                  const DenseMatrix<T> & k);

private:
  struct LocalData
  {
    /// interior-interior matrix entries
    DenseMatrix<T> Aii;
    /// interior-boundary matrix entries
    DenseMatrix<T> Aib;
    /// boundary-interior matrix entries
    DenseMatrix<T> Abi;
    /// boundary-boundary matrix entries
    DenseMatrix<T> Abb;

    struct VarData
    {
      unsigned int boundary_dofs_offset;
      unsigned int interior_dofs_offset;
      unsigned int num_boundary_dofs;
      unsigned int num_interior_dofs;
    };

    std::unordered_map<unsigned int, VarData> var_to_data;
  };

  std::unordered_map<dof_id_type, LocalData> _elem_to_local_data;

  const MeshBase & _mesh;
  const DofMap & _dof_map;
};
}

#endif // LIBMESH_STATIC_CONDENSATION_H
