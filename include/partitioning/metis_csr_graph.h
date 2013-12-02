// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_METIS_CSR_GRAPH_H
#define LIBMESH_METIS_CSR_GRAPH_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"

// C++ Includes   -----------------------------------
#include <vector>
#include <numeric>



namespace libMesh
{

  /**
   * This utility class provides a convenient implementation for
   * building the compressed-row-storage graph required for the METIS/ParMETIS
   * graph partitioning schemes.
   */
  class METIS_CSR_Graph
  {
  public:
    std::vector<int> offsets, vals;

    void prep_n_nonzeros(const libMesh::dof_id_type row, const libMesh::dof_id_type n_nonzeros)
    {
      libmesh_assert_less (row+1, offsets.size());
      offsets[row+1] = n_nonzeros;
    }



    libMesh::dof_id_type n_nonzeros (const libMesh::dof_id_type row) const
    {
      libmesh_assert_less (row+1, offsets.size());
      return (offsets[row+1] - offsets[row]);
    }


    void prepare_for_use()
    {
      std::partial_sum (offsets.begin(), offsets.end(), offsets.begin());
      libmesh_assert (!offsets.empty());
      vals.resize(offsets.back());

      if (vals.empty())
	vals.push_back(0);
    }



    int& operator()(const libMesh::dof_id_type row, const libMesh::dof_id_type nonzero)
    {
      libmesh_assert_greater (vals.size(), offsets[row]+nonzero);

      return vals[offsets[row]+nonzero];
    }



    const int& operator()(const libMesh::dof_id_type row, const libMesh::dof_id_type nonzero) const
    {
      libmesh_assert_greater (vals.size(), offsets[row]+nonzero);

      return vals[offsets[row]+nonzero];
    }

  };

} // namespace libMesh


#endif // LIBMESH_METIS_CSR_GRAPH_H
