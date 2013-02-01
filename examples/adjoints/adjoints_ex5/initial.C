#include "initial.h"

using namespace libMesh;

void read_initial_parameters()
{
}

void finish_initialization()
{
}



// Initial conditions
Number initial_value(const Point& p,
                     const Parameters&,
                     const std::string&,
                     const std::string&)
{
  Real x = p(0), y = p(1);
  
  return sin(M_PI * x) * sin(M_PI * y);
}



Gradient initial_grad(const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  Real x = p(0), y = p(1);

  return Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
  M_PI*sin(M_PI * x) * cos(M_PI * y));  
}
