#ifndef __fem_parameters_h__
#define __fem_parameters_h__

#include <limits>
#include <string>

#include "libmesh_common.h"
#include "getpot.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class FEMParameters
{
public:
    FEMParameters() :  
      domainfile("lshaped.xda"),
      coarserefinements(0),            
      solver_quiet(true),
      require_residual_reduction(true),
      min_step_length(1e-5),
      max_linear_iterations(200000), max_nonlinear_iterations(20),
      relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
      initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
      linear_tolerance_multiplier(1.e-3),
      nelem_target(30000), global_tolerance(0.0),
      refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
      refine_uniformly(false),
      max_adaptivesteps(1),      
      indicator_type("kelly"),
      fe_family(1, "LAGRANGE"), fe_order(1, 1),	
      analytic_jacobians(true), verify_analytic_jacobians(0.0),
      print_solution_norms(false), print_solutions(false),
      print_residual_norms(false), print_residuals(false),
      print_jacobian_norms(false), print_jacobians(false) {}

    void read(GetPot &input);
   
    unsigned int dimension;  
    Real elementorder;    
    std::string domainfile;
    unsigned int coarserefinements;
            
    bool solver_quiet, require_residual_reduction;
    Real min_step_length;
    unsigned int max_linear_iterations, max_nonlinear_iterations;
    Real relative_step_tolerance, relative_residual_tolerance,
	 initial_linear_tolerance, minimum_linear_tolerance,
	 linear_tolerance_multiplier;

    unsigned int nelem_target;
    Real global_tolerance;
    Real refine_fraction, coarsen_fraction, coarsen_threshold;
    bool refine_uniformly;
    unsigned int max_adaptivesteps;
    
    std::string indicator_type;    

    std::vector<std::string> fe_family;
    std::vector<unsigned int> fe_order;

    bool analytic_jacobians;
    Real verify_analytic_jacobians;

    bool print_solution_norms, print_solutions,
         print_residual_norms, print_residuals,
         print_jacobian_norms, print_jacobians;
};

#endif // __fem_parameters_h__
