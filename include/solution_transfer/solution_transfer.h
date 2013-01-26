#ifndef SOLUTIONTRANSFER_H
#define SOLUTIONTRANSFER_H

#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"

#include <string>
#include <map>

namespace libMesh {

/**
 * Base class for objects that allow transfering variable values between different systems with different meshes.
 */
class SolutionTransfer
{
public:
  SolutionTransfer() {};
  virtual ~SolutionTransfer() {};

  /**
   * Transfer the values of a variable to another.
   *
   * This is meant for transferring values from one EquationSystems to another
   * even in the case of having different meshes.
   */
  virtual void transfer(const Variable & from_var, const Variable & to_var) = 0;
};

} // namespace libMesh

#endif
