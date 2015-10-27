// The libMesh Finite Element Library.
// Copyright (C) 2002-2015 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef MESHFUNCTIONSOLUTIONTRANSFER_H
#define MESHFUNCTIONSOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"

#include <string>

namespace libMesh
{

/**
 * Implementation of a SolutionTransfer object that only works for
 * transferring the solution using a MeshFunction
 *
 * Note: A serialization of the "from" solution vector will be
 * performed!  This can be slow in parallel and take a lot of memory!
 */
class MeshFunctionSolutionTransfer : public SolutionTransfer
{
public:
  MeshFunctionSolutionTransfer(const libMesh::Parallel::Communicator
                               &comm LIBMESH_CAN_DEFAULT_TO_COMMWORLD);
  virtual ~MeshFunctionSolutionTransfer();

  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) libmesh_override;
};

} // namespace libMesh

#endif
