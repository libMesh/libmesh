// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local Includes -----------------------------------
#include "libmesh/dof_map.h"
#include "libmesh/sparsity_pattern.h"



namespace libMesh
{
  namespace SparsityPattern
  {


    //-------------------------------------------------------
    // we need to implement these constructors here so that
    // a full DofMap definition is available.
    Build::Build (const MeshBase &mesh_in,
		  const DofMap &dof_map_in,
		  const CouplingMatrix *dof_coupling_in,
		  const bool implicit_neighbor_dofs_in,
		  const bool need_full_sparsity_pattern_in) :
      ParallelObject(dof_map_in),
      mesh(mesh_in),
      dof_map(dof_map_in),
      dof_coupling(dof_coupling_in),
      implicit_neighbor_dofs(implicit_neighbor_dofs_in),
      need_full_sparsity_pattern(need_full_sparsity_pattern_in),
      sparsity_pattern(),
      nonlocal_pattern(),
      n_nz(),
      n_oz()
    {}



    Build::Build (Build &other, Threads::split) :
      ParallelObject(other),
      mesh(other.mesh),
      dof_map(other.dof_map),
      dof_coupling(other.dof_coupling),
      implicit_neighbor_dofs(other.implicit_neighbor_dofs),
      need_full_sparsity_pattern(other.need_full_sparsity_pattern),
      sparsity_pattern(),
      nonlocal_pattern(),
      n_nz(),
      n_oz()
    {}



  } // namespace SparsityPattern
} // namespace libMesh
