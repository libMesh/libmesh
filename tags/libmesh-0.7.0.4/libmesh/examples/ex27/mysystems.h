#include "femparameters.h"

#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

FEMSystem &build_system(EquationSystems &es, GetPot &, FEMParameters &param)
{
  LaplaceSystem &system = es.add_system<LaplaceSystem> ("LaplaceSystem");

  // Use the prescribed FE type
  system.fe_family() = param.fe_family[0];
  system.fe_order() = param.fe_order[0];

  // Use analytical jacobians?
  system.analytic_jacobians() = param.analytic_jacobians;

  return system;
}



