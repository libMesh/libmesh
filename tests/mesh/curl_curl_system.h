// DiffSystem framework files
#include "libmesh/fem_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dirichlet_boundaries.h"

using namespace libMesh;

#ifndef CURL_CURL_SYSTEM_H
#define CURL_CURL_SYSTEM_H

class CurlCurlExactSolution
{
public:
  CurlCurlExactSolution(){}

  ~CurlCurlExactSolution(){}

  RealGradient operator() (Real x, Real y)
  {
    const Real ux =  cos(pi*x)*sin(pi*y);
    const Real uy = -sin(pi*x)*cos(pi*y);

    return RealGradient(ux, uy);
  }

  RealTensor grad(Real x, Real y)
  {
    const Real dux_dx = -pi*sin(pi*x)*sin(pi*y);
    const Real dux_dy = pi*cos(pi*x)*cos(pi*y);
    const Real duy_dx = -pi*cos(pi*x)*cos(pi*y);
    const Real duy_dy = pi*sin(pi*x)*sin(pi*y);

    return RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy);
  }

  RealGradient curl(Real x, Real y)
  {
    const Real dux_dy =  pi*cos(pi*x)*cos(pi*y);
    const Real duy_dx = -pi*cos(pi*x)*cos(pi*y);

    return RealGradient(Real(0), Real(0), duy_dx - dux_dy);
  }

  RealGradient forcing(Real x, Real y)
  {
    const Real fx =  (2*pi*pi + 1)*cos(pi*x)*sin(pi*y);
    const Real fy = -(2*pi*pi + 1)*sin(pi*x)*cos(pi*y);

    return RealGradient(fx, fy);
  }
};

/**
 * FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
 * but we must specify element residuals.
 */
class CurlCurlSystem : public FEMSystem
{
public:
  // Constructor
  CurlCurlSystem(EquationSystems & es,
                 const std::string & name,
                 const unsigned int number);

  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context(DiffContext & context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext & context);

  virtual bool side_time_derivative(bool request_jacobian,
                                    DiffContext & context);

protected:
  // Indices for each variable;
  unsigned int u_var;
  unsigned int v_var;

  void init_dirichlet_bcs();

  // Returns the value of a forcing function at point p.
  RealGradient forcing(const Point & p);

  // Returns the value of the exact solution for this model problem.
  RealGradient exact_solution(const Point & p);

  CurlCurlExactSolution soln;
};

#endif // CURL_CURL_SYSTEM_H
