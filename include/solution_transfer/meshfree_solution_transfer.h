#ifndef MESHFREESOLUTIONTRANSFER_H
#define MESHFREESOLUTIONTRANSFER_H

#include "libmesh/solution_transfer.h"

#include <string>

namespace libMesh {

/**
 * Implementation of a SolutionTransfer object that utilizes the MeshfreeInterpolation system to interpolate one solution to another.
 */
class MeshfreeSolutionTransfer : public SolutionTransfer
{
public:
  MeshfreeSolutionTransfer() {}
  virtual ~MeshfreeSolutionTransfer() {}
  
  /**
   * Transfer the values of a variable to another.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var);
};

} // namespace libMesh

#endif
