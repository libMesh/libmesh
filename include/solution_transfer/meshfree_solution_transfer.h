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



#ifndef MESHFREE_SOLUTION_TRANSFER_H
#define MESHFREE_SOLUTION_TRANSFER_H

#include "libmesh/solution_transfer.h"
#include "libmesh/meshfree_interpolation_function.h"

#include <string>

namespace libMesh
{

/**
 * Implementation of a SolutionTransfer object that utilizes the
 * MeshfreeInterpolation system to interpolate one solution to
 * another.
 *
 * \author Derek Gaston
 * \date 2013
 * \brief A SolutionTransfer object that does "mesh free" interpolation.
 */
class MeshfreeSolutionTransfer : public SolutionTransfer
{
public:
  MeshfreeSolutionTransfer(const libMesh::Parallel::Communicator & comm_in) :
    SolutionTransfer(comm_in)
  {}

  virtual ~MeshfreeSolutionTransfer() = default;

  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) override;
};

} // namespace libMesh

#endif // #define MESHFREE_SOLUTION_TRANSFER_H
