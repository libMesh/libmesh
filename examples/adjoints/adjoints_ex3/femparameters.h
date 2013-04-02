
#ifndef __fem_parameters_h__
#define __fem_parameters_h__

#include <limits>

#include "libmesh/libmesh_common.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/getpot.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class FEMParameters
{
public:
    FEMParameters() :
      initial_timestep(0), n_timesteps(100),
      transient(true),
      deltat_reductions(0),
      timesolver_core("euler"),
      end_time(std::numeric_limits<Real>::max()),
      deltat(0.0001), timesolver_theta(0.5),
      timesolver_maxgrowth(0.), timesolver_tolerance(0.),
      timesolver_upper_tolerance(0.),
      steadystate_tolerance(0.),
      timesolver_norm(0,L2),
      dimension(2),
	domaintype("square"), domainfile("mesh.xda"), elementtype("quad"),
	fine_mesh_file_primal("fine_mesh.xda"), fine_mesh_soln_primal("fine_mesh_soln.xda"),
	fine_mesh_file_adjoint("fine_mesh.xda"), fine_mesh_soln_adjoint("fine_mesh_soln.xda"),
      elementorder(2),
      domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
      domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
      coarsegridx(1), coarsegridy(1), coarsegridz(1),
	coarserefinements(0), coarsecoarsenings(0), extrarefinements(0),
      use_petsc_snes(false),
      time_solver_quiet(true), solver_quiet(true),
	reuse_preconditioner(false),
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
      initial_adaptivesteps(0), initial_sobolev_order(1),
	initial_extra_quadrature(0),
	indicator_type("kelly"), patch_reuse(true), sobolev_order(1), adjoint_residual_type("patchpatch"), alternate_with_uniform_steps("false"), alternate_step_number(2),
	component_wise_error("false"), compare_to_fine_solution_primal(false), compare_to_fine_solution_adjoint(false), compute_sensitivities(false), do_forward_sensitivity(true), do_adjoint_sensitivity(false),
	postprocess_adjoint(false),
      write_interval(10),
      write_gmv_error(false), write_tecplot_error(false),
      output_xda(false), output_xdr(false),
      output_bz2(true), output_gz(true),
      output_gmv(false), output_tecplot(false),
      run_simulation(true), run_postprocess(false),
      fe_family(1, "LAGRANGE"), fe_order(1, 1),
      extra_quadrature_order(0),
      analytic_jacobians(true), verify_analytic_jacobians(0.0),
      print_solution_norms(false), print_solutions(false),
      print_residual_norms(false), print_residuals(false),
      print_jacobian_norms(false), print_jacobians(false) {}

    void read(GetPot &input);



    unsigned int initial_timestep, n_timesteps;
    bool transient;
    unsigned int deltat_reductions;
    std::string timesolver_core;
    Real end_time, deltat, timesolver_theta,
	 timesolver_maxgrowth, timesolver_tolerance,
	 timesolver_upper_tolerance, steadystate_tolerance;
    std::vector<FEMNormType> timesolver_norm;

    unsigned int dimension;
    std::string domaintype, domainfile, elementtype, fine_mesh_file_primal, fine_mesh_soln_primal, fine_mesh_file_adjoint, fine_mesh_soln_adjoint;
    Real elementorder;
    Real domain_xmin, domain_ymin, domain_zmin;
    Real domain_edge_width, domain_edge_length, domain_edge_height;
    unsigned int coarsegridx, coarsegridy, coarsegridz;
    unsigned int coarserefinements, coarsecoarsenings, extrarefinements;

    bool use_petsc_snes;
    bool time_solver_quiet, solver_quiet, reuse_preconditioner, require_residual_reduction;
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
    unsigned int initial_adaptivesteps;
    unsigned int initial_sobolev_order;
    unsigned int initial_extra_quadrature;

    std::string indicator_type;
    bool patch_reuse;
    unsigned int sobolev_order;
    std::string adjoint_residual_type;
    bool alternate_with_uniform_steps;
    unsigned int alternate_step_number;
    bool component_wise_error;
    bool compare_to_fine_solution_primal, compare_to_fine_solution_adjoint;

    bool compute_sensitivities;
    bool do_forward_sensitivity;
    bool do_adjoint_sensitivity;

    bool postprocess_adjoint;

    unsigned int write_interval;
    bool write_gmv_error, write_tecplot_error, output_xda, output_xdr,
	 output_bz2, output_gz, output_gmv, output_tecplot;
    bool run_simulation, run_postprocess;

    std::vector<std::string> fe_family;
    std::vector<unsigned int> fe_order;
    int extra_quadrature_order;

    bool analytic_jacobians;
    Real verify_analytic_jacobians;

    bool print_solution_norms, print_solutions,
         print_residual_norms, print_residuals,
         print_jacobian_norms, print_jacobians;
};

#endif // __fem_parameters_h__
