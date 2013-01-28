#ifndef DIRECTSOLUTIONTRANSFER_H
#define DIRECTSOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"

#include <string>

namespace libMesh {

/**
 * Implementation of a SolutionTransfer object that only works for transferring the solution but only in the case of:
 *
 * 1.  The Systems must have EXACTLY the same mesh.
 * 2.  The two variables involved must have EXACTLY the same finite element family and order.
 */
class DirectSolutionTransfer : public SolutionTransfer
{
public:
  DirectSolutionTransfer();
  virtual ~DirectSolutionTransfer();
  
  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var);
};

} // namespace libMesh

#endif
