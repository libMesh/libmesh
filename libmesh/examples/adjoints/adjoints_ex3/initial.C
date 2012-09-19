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
Number initial_value(const Point& p,
                     const Parameters& param,
                     const std::string&,
                     const std::string&)
{
  const Real x1 = p(0);
  const Real x2 = p(1);
  //  const Real pie = libMesh::pi;

  //Real alpha = param.get<Real>("alpha");
  // Real q = param.get<Real>("q");
  
  // Real c1 = 1./(2*alpha*pie*pie);
  // Real c2 = (q/alpha + 1./(2*alpha*pie))/(pie * cosh(pie));
  
  return 1.;
      
}



Gradient initial_grad(const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  Real x = p(0), y = p(1);
/*
  if ((y < x) && (y < (1.0-x)))
    return Gradient(0.0,1.0);
  if ((y > x) && (y < (1.0-x)))
    return Gradient(-1.0,0.0);
  if ((y > x) && (y > (1.0-x)))
    return Gradient(0.0,-1.0);
  return Gradient(1.0,0.0);
*/
  return 0.;
}
