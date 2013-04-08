#ifndef MESHFUNCTIONSOLUTIONTRANSFER_H
#define MESHFUNCTIONSOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"

#include <string>

namespace libMesh {

/**
 * Implementation of a SolutionTransfer object that only works for transferring the solution using a MeshFunction
 *
 * Note: A serialization of the "from" solution vector will be performed!  This can be slow in parallel and take a lot of memory!
 */
class MeshFunctionSolutionTransfer : public SolutionTransfer
{
public:
  MeshFunctionSolutionTransfer(const libMesh::Parallel::Communicator &comm = libMesh::CommWorld);
  virtual ~MeshFunctionSolutionTransfer();

  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var);
};

} // namespace libMesh

#endif
