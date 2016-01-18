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



#ifndef MESHFREESOLUTIONTRANSFER_H
#define MESHFREESOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"

#include <string>

namespace libMesh
{

/**
 * Implementation of a SolutionTransfer object that utilizes the
 * MeshfreeInterpolation system to interpolate one solution to
 * another.
 */
class MeshfreeSolutionTransfer : public SolutionTransfer
{
public:
  MeshfreeSolutionTransfer(const libMesh::Parallel::Communicator & comm_in
                           LIBMESH_CAN_DEFAULT_TO_COMMWORLD) :
    SolutionTransfer(comm_in)
  {}

  virtual ~MeshfreeSolutionTransfer() {}

  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) libmesh_override;
};

} // namespace libMesh

#endif // #define MESHFREESOLUTIONTRANSFER_H
