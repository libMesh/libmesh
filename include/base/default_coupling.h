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



#ifndef LIBMESH_DEFAULT_COUPLING_H
#define LIBMESH_DEFAULT_COUPLING_H

// Local Includes -----------------------------------
#include "libmesh/ghosting_functor.h"

namespace libMesh
{

/**
 * This class implements the default geometry ghosting in libMesh:
 * point neighbors and interior_parent elements are ghosted.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
class DefaultCoupling : public GhostingFunctor
{
public:

  /**
   * Constructor.
   */
  DefaultCoupling(const MeshBase & mesh,
                  const CouplingMatrix * dof_coupling,
                  bool couple_neighbor_dofs) :
    _mesh(mesh), _dof_coupling(dof_coupling),
    _couple_neighbor_dofs(couple_neighbor_dofs) {}

  /**
   * For the specified range of active elements, find the elements
   * which will be coupled to them in the sparsity pattern.
   *
   * This will be only the elements themselves by default, but will
   * include side neighbors if an all-discontinuous-variable system is
   * detected and/or if the user specified --implicit_neighbor_dofs on
   * the command line or used set_implicit_neighbor_dofs() in their
   * code.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements);

private:

  const MeshBase & _mesh;
  const CouplingMatrix * _dof_coupling;
  bool _couple_neighbor_dofs;
};

} // namespace libMesh

#endif // LIBMESH_DEFAULT_COUPLING_H
