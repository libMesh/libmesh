#include "adjoint_initial.h"

using namespace libMesh;

void adjoint_read_initial_parameters()
{
}

void adjoint_finish_initialization()
{
}



// Initial conditions
Number adjoint_initial_value(const Point& p,
                     const Parameters&,
                     const std::string&,
                     const std::string&)
{
  Real x = p(0), y = p(1);

  Number val = 0.;
 
  val = sin(M_PI * x) * sin(M_PI * y);
 
  return val;
}



Gradient adjoint_initial_grad(const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  Real x = p(0), y = p(1);

  return Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
  M_PI*sin(M_PI * x) * cos(M_PI * y));    
}
