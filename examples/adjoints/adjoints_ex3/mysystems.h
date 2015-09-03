#include "coupled_system.h"
#include "femparameters.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

FEMSystem &build_system(EquationSystems &es, GetPot &, FEMParameters& /* param */)
{
  CoupledSystem &system = es.add_system<CoupledSystem> ("CoupledSystem");

  return system;
}

Number exact_value(const Point& /* p */,       // xyz location
                   const Parameters& /* param */,  // EquationSystem parameters
                   const std::string&, // sys_name
                   const std::string&) // unknown_name
{
  std::cout<<"Warning ! No exact value function specified ! Returning 0." << std::endl;

  return Number(0.);
}

Gradient exact_grad(const Point& /* p */,       // xyz location
                    const Parameters& /* param */,  // EquationSystems parameters
                    const std::string&, // sys_name
                    const std::string&) // unknown_name
{
  std::cout<<"Warning ! No exact value function specified ! Returning 0." << std::endl;

  return Gradient(0.);
}
