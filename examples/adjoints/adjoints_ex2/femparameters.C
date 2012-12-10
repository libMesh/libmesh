#include "femparameters.h"

#define GETPOT_INPUT(A) { A = input(#A, A);\
  std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); };
#define GETPOT_INT_INPUT(A) { A = input(#A, (int)A);\
  std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); };
#define GETPOT_REGISTER(A) { \
  std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); };


void FEMParameters::read(GetPot &input)
{
    std::vector<std::string> variable_names;
   
    
    GETPOT_INT_INPUT(coarserefinements);  
    GETPOT_INPUT(domainfile);

    GETPOT_INPUT(solver_quiet);
    GETPOT_INPUT(require_residual_reduction);
    GETPOT_INPUT(min_step_length);
    GETPOT_INT_INPUT(max_linear_iterations);
    GETPOT_INT_INPUT(max_nonlinear_iterations);
    GETPOT_INPUT(relative_step_tolerance);
    GETPOT_INPUT(relative_residual_tolerance);
    GETPOT_INPUT(initial_linear_tolerance);
    GETPOT_INPUT(minimum_linear_tolerance);
    GETPOT_INPUT(linear_tolerance_multiplier);
    GETPOT_INT_INPUT(nelem_target);
    GETPOT_INPUT(global_tolerance);
    GETPOT_INPUT(refine_fraction);
    GETPOT_INPUT(coarsen_fraction);
    GETPOT_INPUT(coarsen_threshold);
    GETPOT_INT_INPUT(max_adaptivesteps);   
    GETPOT_INPUT(refine_uniformly);
    GETPOT_INPUT(indicator_type);

    GETPOT_REGISTER(fe_family);
    const unsigned int n_fe_family =
      std::max(1u, input.vector_variable_size("fe_family"));
    fe_family.resize(n_fe_family, "LAGRANGE");
    for (unsigned int i=0; i != n_fe_family; ++i)
      fe_family[i]              = input("fe_family", fe_family[i].c_str(), i);
    GETPOT_REGISTER(fe_order);
    const unsigned int n_fe_order =
      input.vector_variable_size("fe_order");
    fe_order.resize(n_fe_order, 1);
    for (unsigned int i=0; i != n_fe_order; ++i)
      fe_order[i]               = input("fe_order", (int)fe_order[i], i);

    GETPOT_INPUT(analytic_jacobians);
    GETPOT_INPUT(verify_analytic_jacobians);
    GETPOT_INPUT(print_solution_norms);
    GETPOT_INPUT(print_solutions);
    GETPOT_INPUT(print_residual_norms);
    GETPOT_INPUT(print_residuals);
    GETPOT_INPUT(print_jacobian_norms);
    GETPOT_INPUT(print_jacobians);

  std::vector<std::string> bad_variables =
    input.unidentified_arguments(variable_names);

  if (libMesh::processor_id() == 0 && !bad_variables.empty())
    {
      std::cerr << "ERROR: Unrecognized variables:" << std::endl;
      for (unsigned int i = 0; i != bad_variables.size(); ++i)
        std::cerr << bad_variables[i] << std::endl;
      std::cerr << "not found among recognized variables." << std::endl;
      for (unsigned int i = 0; i != variable_names.size(); ++i)
        std::cerr << variable_names[i] << std::endl;
      libmesh_error();
    }
 }
