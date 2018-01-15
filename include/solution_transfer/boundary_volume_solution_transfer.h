// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef BOUNDARY_VOLUME_SOLUTION_TRANSFER
#define BOUNDARY_VOLUME_SOLUTION_TRANSFER

#include "libmesh/solution_transfer.h"


namespace libMesh {

/**
 * SolutionTransfer derived class which is specifically for
 * transferring solutions back and forth between a VolumeMesh and its
 * associated BoundaryMesh.  That is, when calling transfer() below,
 * we must have either:
 *
 * .) Transfer from the boundary of a VolumeMesh to its associated
 *    BoundaryMesh, or
 * .) Transfer from a BoundaryMesh to its corresponding VolumeMesh.
 *
 * Because of this assumption, we can use the interior_parent() to
 * directly define the mapping between solutions.  The transfer()
 * function internally determines which direction the transfer is
 * going, and calls the appropriate algorithm internally.
 *
 * \author Xikai Jiang
 * \author John W. Peterson
 * \date 2016
 */
class BoundaryVolumeSolutionTransfer : public SolutionTransfer
{
public:
  BoundaryVolumeSolutionTransfer (const Parallel::Communicator & comm_in LIBMESH_CAN_DEFAULT_TO_COMMWORLD) :
    SolutionTransfer(comm_in)
  {}

  virtual ~BoundaryVolumeSolutionTransfer() {}

  /**
   * Transfer values from a Variable in a System associated with a
   * volume mesh to a Variable in a System associated with the
   * corresponding BoundaryMesh, or vice-versa.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) libmesh_override;

private:
  /**
   * Transfer values *from* volume mesh *to* boundary mesh.  Called
   * when needed from transfer().
   */
  void transfer_volume_boundary(const Variable & from_var, const Variable & to_var);

  /**
   * Transfer values *from* boundary mesh *to* volume mesh.  Called
   * when needed from transfer().
   */
  void transfer_boundary_volume(const Variable & from_var, const Variable & to_var);
};

} // namespace libMesh

#endif
