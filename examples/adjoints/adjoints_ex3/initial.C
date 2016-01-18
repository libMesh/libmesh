#include "initial.h"
#include "femparameters.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void read_initial_parameters()
{
}

void finish_initialization()
{
}



// Initial conditions
Number initial_value(const Point & /* p */,
                     const Parameters & /* param */,
                     const std::string &,
                     const std::string &)
{

  return Number(1.);

}



Gradient initial_grad(const Point & /* p */,
                      const Parameters & /* param */,
                      const std::string &,
                      const std::string &)
{
  return Gradient(0.);
}
