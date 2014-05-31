
#include "femparameters.h"

using namespace libMesh;

#define GETPOT_INPUT(A) { A = input(#A, A);                             \
    std::string stringval = input(#A, std::string());                   \
    variable_names.push_back(std::string(#A "=") + stringval); };
#define GETPOT_INT_INPUT(A) { A = input(#A, (int)A);                    \
    std::string stringval = input(#A, std::string());                   \
    variable_names.push_back(std::string(#A "=") + stringval); };
#define GETPOT_REGISTER(A) {                                            \
    std::string stringval = input(#A, std::string());                   \
    variable_names.push_back(std::string(#A "=") + stringval); };


void FEMParameters::read(GetPot &input)
{
  std::vector<std::string> variable_names;

  GETPOT_INT_INPUT(initial_timestep);
  GETPOT_INT_INPUT(n_timesteps);
  GETPOT_INPUT(transient);
  GETPOT_INT_INPUT(deltat_reductions);
  GETPOT_INPUT(timesolver_core);
  GETPOT_INPUT(end_time);
  GETPOT_INPUT(deltat);
  GETPOT_INPUT(timesolver_theta);
  GETPOT_INPUT(timesolver_maxgrowth);
  GETPOT_INPUT(timesolver_upper_tolerance);
  GETPOT_INPUT(timesolver_tolerance);
  GETPOT_INPUT(steadystate_tolerance);

  GETPOT_REGISTER(timesolver_norm);
  const unsigned int n_timesolver_norm = input.vector_variable_size("timesolver_norm");
  timesolver_norm.resize(n_timesolver_norm, L2);
  for (unsigned int i=0; i != n_timesolver_norm; ++i)
    {
      int current_norm = 0; // L2
      if (timesolver_norm[i] == H1)
        current_norm = 1;
      if (timesolver_norm[i] == H2)
        current_norm = 2;
      current_norm = input("timesolver_norm", current_norm, i);
      if (current_norm == 0)
        timesolver_norm[i] = L2;
      else if (current_norm == 1)
        timesolver_norm[i] = H1;
      else if (current_norm == 2)
        timesolver_norm[i] = H2;
      else
        timesolver_norm[i] = DISCRETE_L2;
    }

  GETPOT_INT_INPUT(dimension);
  GETPOT_INPUT(domaintype);
  GETPOT_INPUT(domainfile);
  GETPOT_INPUT(elementtype);
  GETPOT_INPUT(fine_mesh_file_primal);
  GETPOT_INPUT(fine_mesh_soln_primal);
  GETPOT_INPUT(fine_mesh_file_adjoint);
  GETPOT_INPUT(fine_mesh_soln_adjoint);
  GETPOT_INPUT(elementorder);
  GETPOT_INPUT(domain_xmin);
  GETPOT_INPUT(domain_ymin);
  GETPOT_INPUT(domain_zmin);
  GETPOT_INPUT(domain_edge_width);
  GETPOT_INPUT(domain_edge_length);
  GETPOT_INPUT(domain_edge_height);
  GETPOT_INT_INPUT(coarsegridx);
  GETPOT_INT_INPUT(coarsegridy);
  GETPOT_INT_INPUT(coarsegridz);
  GETPOT_INT_INPUT(coarserefinements);
  GETPOT_INT_INPUT(coarsecoarsenings);
  GETPOT_INT_INPUT(extrarefinements);
  GETPOT_INPUT(use_petsc_snes);
  GETPOT_INPUT(time_solver_quiet);
  GETPOT_INPUT(solver_quiet);
  GETPOT_INPUT(reuse_preconditioner);
  GETPOT_INPUT(require_residual_reduction);
  GETPOT_INPUT(min_step_length);
  GETPOT_INT_INPUT(max_linear_iterations);
  GETPOT_INT_INPUT(max_nonlinear_iterations);
  GETPOT_INPUT(relative_step_tolerance);
  GETPOT_INPUT(absolute_residual_tolerance);
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
  GETPOT_INT_INPUT(initial_adaptivesteps);
  GETPOT_INT_INPUT(initial_sobolev_order);
  GETPOT_INT_INPUT(initial_extra_quadrature);
  GETPOT_INPUT(refine_uniformly);
  GETPOT_INPUT(indicator_type);
  GETPOT_INPUT(adjoint_residual_type);
  GETPOT_INPUT(patch_reuse);
  GETPOT_INT_INPUT(sobolev_order);
  GETPOT_INPUT(alternate_with_uniform_steps);
  GETPOT_INPUT(alternate_step_number);
  GETPOT_INPUT(component_wise_error);
  GETPOT_INPUT(compute_sensitivities);
  GETPOT_INPUT(compare_to_fine_solution_primal);
  GETPOT_INPUT(compare_to_fine_solution_adjoint);
  GETPOT_INPUT(do_forward_sensitivity);
  GETPOT_INPUT(do_adjoint_sensitivity);
  GETPOT_INPUT(postprocess_adjoint);
  GETPOT_INT_INPUT(write_interval);
  GETPOT_INPUT(output_xda);
  GETPOT_INPUT(output_xdr);

  GETPOT_INPUT(output_gz);
#ifndef LIBMESH_HAVE_GZSTREAM
  output_gz                   = false;
#endif

  GETPOT_INPUT(output_bz2);
#ifndef LIBMESH_HAVE_BZ2
  output_bz2                  = false;
#endif

  GETPOT_INPUT(output_gmv);
  GETPOT_INPUT(write_gmv_error);

#ifndef LIBMESH_HAVE_GMV
  output_gmv                  = false;
  write_gmv_error             = false;
#endif

  GETPOT_INPUT(output_tecplot);
  GETPOT_INPUT(write_tecplot_error);
#ifndef LIBMESH_HAVE_TECPLOT_API
  output_tecplot              = false;
  write_tecplot_error         = false;
#endif

  GETPOT_INPUT(run_simulation);
  GETPOT_INPUT(run_postprocess);

  run_simulation              = input("run_simulation", run_simulation);
  run_postprocess             = input("run_postprocess", run_postprocess);

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

  GETPOT_INPUT(extra_quadrature_order);
  GETPOT_INPUT(analytic_jacobians);
  GETPOT_INPUT(verify_analytic_jacobians);
  GETPOT_INPUT(print_solution_norms);
  GETPOT_INPUT(print_solutions);
  GETPOT_INPUT(print_residual_norms);
  GETPOT_INPUT(print_residuals);
  GETPOT_INPUT(print_jacobian_norms);
  GETPOT_INPUT(print_jacobians);


  // std::vector<std::string> bad_variables =
  //   input.unidentified_arguments(variable_names);

  // if (libMesh::processor_id() == 0 && !bad_variables.empty())
  //   {
  //     libMesh::err << "ERROR: Unrecognized variable(s):" << std::endl;
  //     for (unsigned int i = 0; i != bad_variables.size(); ++i)
  //       libMesh::err << bad_variables[i] << std::endl;
  //     libMesh::err << "\nRecognized variables are:" << std::endl;
  //     for (unsigned int i = 0; i != variable_names.size(); ++i)
  //       libMesh::err << variable_names[i] << std::endl;
  //     libmesh_error_msg("Fix/remove the unrecognized variables and try again.");
  //   }
}
