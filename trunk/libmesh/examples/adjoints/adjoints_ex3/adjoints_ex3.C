// C++ includes
#include <iostream>
#include <sys/time.h>

// Libmesh includes
#include "libmesh/equation_systems.h"

#include "libmesh/twostep_time_solver.h"
#include "libmesh/euler_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/exact_error_estimator.h"
#include "libmesh/steady_solver.h"

#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_diff_solver.h"

#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/system.h"
#include "libmesh/system_norm.h"

#include "libmesh/adjoint_residual_error_estimator.h"
#include "libmesh/const_fem_function.h"
#include "libmesh/error_vector.h"
#include "libmesh/exact_solution.h"
#include "libmesh/fem_function_base.h"
#include "libmesh/fourth_error_estimators.h"
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/refinement_selector.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/qoi_set.h"
#include "libmesh/weighted_patch_recovery_error_estimator.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "libmesh/o_string_stream.h"

// Local includes
#include "coupled_system.h"
#include "domain.h"
#include "initial.h"
#include "femparameters.h"
#include "mysystems.h"
#include "output.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Number output files

std::string numbered_filename(unsigned int t_step, // The timestep count
                              unsigned int a_step, // The adaptive step count
                              std::string solution_type, // primal or adjoint solve
                              std::string type,
                              std::string extension,
                              FEMParameters &param)
{
  OStringStream file_name;
  file_name << solution_type << ".out." << type << '.';
  OSSRealzeroright(file_name, 3, 0, t_step);
  file_name << '.' ;
  OSSRealzeroright(file_name, 2, 0, a_step);
  file_name << '.' << extension;

  if (param.output_bz2)
    file_name << ".bz2";
  else if (param.output_gz)
    file_name << ".gz";
  return file_name.str();
}

// Write tecplot, gmv and xda/xdr output

void write_output(EquationSystems &es,                  
                  unsigned int t_step, // The timestep count
                  unsigned int a_step, // The adaptive step count
                  std::string solution_type, // primal or adjoint solve
                  FEMParameters &param)
{
  MeshBase &mesh = es.get_mesh();

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      OStringStream file_name_gmv;
      file_name_gmv << solution_type << ".out.gmv.";
      OSSRealzeroright(file_name_gmv,3,0,t_step);
      file_name_gmv << '.' ;
      OSSRealzeroright(file_name_gmv,2,0,a_step);

      GMVIO(mesh).write_equation_systems
        (file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.output_tecplot)
    {
      OStringStream file_name_tecplot;
      file_name_tecplot << solution_type << ".out.";
      OSSRealzeroright(file_name_tecplot,3,0,t_step);
      file_name_tecplot << '.' ;
      OSSRealzeroright(file_name_tecplot,2,0,a_step);
      file_name_tecplot << ".plt";

      TecplotIO(mesh).write_equation_systems
        (file_name_tecplot.str(), es);
    }
#endif

  if (param.output_xda || param.output_xdr)
    mesh.renumber_nodes_and_elements();
  if (param.output_xda)
    {
      mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param));
      es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xda", param),
               libMeshEnums::WRITE, EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
  if (param.output_xdr)
    {
      mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param));
      es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param),
               libMeshEnums::WRITE, EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
}

void write_output_headers(FEMParameters &param)
{
  // Only one processor needs to take care of headers.
  if (libMesh::processor_id() != 0)
    return;

  start_output(param.initial_timestep, "out_clocktime.m", "vector_clocktime");

  if (param.run_simulation)
    {
      start_output(param.initial_timestep, "out_activemesh.m", "table_activemesh");
      start_output(param.initial_timestep, "out_solvesteps.m", "table_solvesteps");

      if (param.timesolver_tolerance)
        {
          start_output(param.initial_timestep, "out_time.m", "vector_time");
          start_output(param.initial_timestep, "out_timesteps.m", "vector_timesteps");
        }
      if (param.steadystate_tolerance)
        start_output(param.initial_timestep, "out_changerate.m", "vector_changerate");
    }
}

void write_output_solvedata(EquationSystems &es,
                            unsigned int a_step,
                            unsigned int newton_steps,
                            unsigned int krylov_steps,
                            unsigned int tv_sec,
                            unsigned int tv_usec)
{
  MeshBase &mesh = es.get_mesh();
  unsigned int n_active_elem = mesh.n_active_elem();
  unsigned int n_active_dofs = es.n_active_dofs();

  if (libMesh::processor_id() == 0)
    {
      // Write out the number of elements/dofs used
      std::ofstream activemesh ("out_activemesh.m",
        std::ios_base::app | std::ios_base::out);
      activemesh.precision(17);
      activemesh << (a_step + 1) << ' '
                 << n_active_elem << ' '
                 << n_active_dofs << std::endl;

      // Write out the number of solver steps used
      std::ofstream solvesteps ("out_solvesteps.m",
        std::ios_base::app | std::ios_base::out);
      solvesteps.precision(17);
      solvesteps << newton_steps << ' '
                 << krylov_steps << std::endl;

      // Write out the clock time
      std::ofstream clocktime ("out_clocktime.m",
        std::ios_base::app | std::ios_base::out);
      clocktime.precision(17);
      clocktime << tv_sec << '.' << tv_usec << std::endl;
    }
}

void write_output_footers(FEMParameters &param)
{
  // Write footers on output .m files
  if (libMesh::processor_id() == 0)
    {
      std::ofstream clocktime ("out_clocktime.m",
        std::ios_base::app | std::ios_base::out);
      clocktime << "];" << std::endl;

      if (param.run_simulation)
        {
          std::ofstream activemesh ("out_activemesh.m",
            std::ios_base::app | std::ios_base::out);
          activemesh << "];" << std::endl;

          std::ofstream solvesteps ("out_solvesteps.m",
            std::ios_base::app | std::ios_base::out);
          solvesteps << "];" << std::endl;

          if (param.timesolver_tolerance)
            {
              std::ofstream times ("out_time.m",
                std::ios_base::app | std::ios_base::out);
              times << "];" << std::endl;
              std::ofstream timesteps ("out_timesteps.m",
                std::ios_base::app | std::ios_base::out);
              timesteps << "];" << std::endl;
            }
          if (param.steadystate_tolerance)
            {
              std::ofstream changerate ("out_changerate.m",
                std::ios_base::app | std::ios_base::out);
              changerate << "];" << std::endl;
            }
        }

      // We're done, so write out a file (for e.g. Make to check)
      std::ofstream complete ("complete");
      complete << "complete" << std::endl;
    }
}

#if defined(LIBMESH_HAVE_GMV) || defined(LIBMESH_HAVE_TECPLOT_API)
void write_error(EquationSystems &es,
                 ErrorVector &error,
                 unsigned int t_number,
                 unsigned int a_number,
                 FEMParameters &param,
                 std::string error_type)
#else
void write_error(EquationSystems &,
                 ErrorVector &,
                 unsigned int,
                 unsigned int,
                 FEMParameters &,
                 std::string)
#endif
{
#ifdef LIBMESH_HAVE_GMV
  if (param.write_gmv_error)
    {
      OStringStream error_gmv;
      error_gmv << "error.gmv.";
      OSSRealzeroright(error_gmv,3,0, a_number);
      error_gmv << ".";
      OSSRealzeroright(error_gmv,2,0, t_number);
      error_gmv << error_type;
      error.plot_error(error_gmv.str(), es.get_mesh());
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.write_tecplot_error)
    {
      OStringStream error_tecplot;
      error_tecplot << "error.plt.";
      OSSRealzeroright(error_tecplot,3,0, a_number);
      error_tecplot << ".";
      OSSRealzeroright(error_tecplot,2,0, t_number);
      error_tecplot << error_type;
      error.plot_error(error_tecplot.str(), es.get_mesh());
    }
#endif
}

void read_output(EquationSystems &es,
                 unsigned int t_step,
                 unsigned int a_step,
                 std::string solution_type, 
                 FEMParameters &param)
{
  MeshBase &mesh = es.get_mesh();

  std::string file_name_mesh, file_name_soln;
  // Look for ASCII files first
  if (param.output_xda)
    {
      file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param);
      file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xda", param);
    }
  else if (param.output_xdr)
    {
      file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param);
      file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param);
    }

  // Read in the mesh
  mesh.read(file_name_mesh);

  // And the stored solution
  es.read(file_name_soln, libMeshEnums::READ,
          EquationSystems::READ_HEADER |
          EquationSystems::READ_DATA |
          EquationSystems::READ_ADDITIONAL_DATA);

  // Put systems in a consistent state
  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).update();

  // Figure out the current time
  Real current_time = 0., current_timestep = 0.;

  if (param.timesolver_tolerance)
    {
      std::ifstream times ("out_time.m");
      std::ifstream timesteps ("out_timesteps.m");
      if (times.is_open() && timesteps.is_open())
        {
          // Read headers
          const unsigned int headersize = 25;
          char header[headersize];
          timesteps.getline (header, headersize);
          if (strcmp(header, "vector_timesteps = [") != 0)
            {
              std::cout << "Bad header in out_timesteps.m:" << std::endl
                        << header
                        << std::endl;
              libmesh_error();
            }

          times.getline (header, headersize);
          if (strcmp(header, "vector_time = [") != 0)
            {
              std::cout << "Bad header in out_time.m:" << std::endl
                        << header
                        << std::endl;
              libmesh_error();
            }

          // Read each timestep
          for (unsigned int i = 0; i != t_step; ++i)
            {
              if (!times.good())
                libmesh_error();
              times >> current_time;
              timesteps >> current_timestep;
            }
          // Remember to increment the last timestep; out_times.m
          // lists each *start* time
          current_time += current_timestep;
        }
      else
        libmesh_error();
    }
  else
    current_time = t_step * param.deltat;

  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).time = current_time;
}

void set_system_parameters(FEMSystem &system, FEMParameters &param)
{
  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;

  // More desperate debugging options
  system.print_solution_norms = param.print_solution_norms;
  system.print_solutions      = param.print_solutions;
  system.print_residual_norms = param.print_residual_norms;
  system.print_residuals      = param.print_residuals;
  system.print_jacobian_norms = param.print_jacobian_norms;
  system.print_jacobians      = param.print_jacobians;

  // Solve this as a time-dependent or steady system
  if (param.transient)
    {
      UnsteadySolver *innersolver;
      if (param.timesolver_core == "euler2")
        {
          Euler2Solver *euler2solver =
            new Euler2Solver(system);

          euler2solver->theta = param.timesolver_theta;
          innersolver = euler2solver;
        }
      else if (param.timesolver_core == "euler")
        {
          EulerSolver *eulersolver =
            new EulerSolver(system);

          eulersolver->theta = param.timesolver_theta;
          innersolver = eulersolver;
        }
      else
        {
          std::cerr << "Don't recognize core TimeSolver type: "
                    << param.timesolver_core << std::endl;
          libmesh_error();
        }

      if (param.timesolver_tolerance)
        {
          TwostepTimeSolver *timesolver =
            new TwostepTimeSolver(system);

          timesolver->max_growth       = param.timesolver_maxgrowth;
          timesolver->target_tolerance = param.timesolver_tolerance;
          timesolver->upper_tolerance  = param.timesolver_upper_tolerance;
          timesolver->component_norm   = SystemNorm(param.timesolver_norm);

          timesolver->core_time_solver =
            AutoPtr<UnsteadySolver>(innersolver);
          system.time_solver =
            AutoPtr<UnsteadySolver>(timesolver);
        }
      else
        system.time_solver =
          AutoPtr<TimeSolver>(innersolver);
    }
  else
    system.time_solver =
      AutoPtr<TimeSolver>(new SteadySolver(system));

  system.time_solver->reduce_deltat_on_diffsolver_failure =
                                        param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

  // Create any periodic boundary conditions
  for (unsigned int i=0; i != param.periodic_boundaries.size(); ++i)
    system.get_dof_map().add_periodic_boundary(param.periodic_boundaries[i]);

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);
#else
      libmesh_error();
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);

      solver->quiet                       = param.solver_quiet;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
      solver->relative_residual_tolerance = param.relative_residual_tolerance;
      solver->require_residual_reduction  = param.require_residual_reduction;
      solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      if (system.time_solver->reduce_deltat_on_diffsolver_failure)
        {
          solver->continue_after_max_iterations = true;
          solver->continue_after_backtrack_failure = true;
        }

      // And the linear solver options
      solver->max_linear_iterations       = param.max_linear_iterations;
      solver->initial_linear_tolerance    = param.initial_linear_tolerance;
      solver->minimum_linear_tolerance    = param.minimum_linear_tolerance;      
    }
}

AutoPtr<MeshRefinement> build_mesh_refinement(MeshBase &mesh,
                                              FEMParameters &param)
{
  AutoPtr<MeshRefinement> mesh_refinement(new MeshRefinement(mesh));
  mesh_refinement->coarsen_by_parents() = true;
  mesh_refinement->absolute_global_tolerance() = param.global_tolerance;
  mesh_refinement->nelem_target()      = param.nelem_target;  
  mesh_refinement->refine_fraction()   = param.refine_fraction;
  mesh_refinement->coarsen_fraction()  = param.coarsen_fraction;  
  mesh_refinement->coarsen_threshold() = param.coarsen_threshold;

  return mesh_refinement;
}

AutoPtr<ErrorEstimator> build_error_estimator(FEMParameters &param)
{
  AutoPtr<ErrorEstimator> error_estimator;

  error_estimator.reset(new KellyErrorEstimator);

  return error_estimator;
}


AutoPtr<ErrorEstimator>
build_error_estimator_component_wise
  (FEMParameters &param, 
   unsigned int a_number,
   std::vector<std::vector<Number> > &term_weights,
   std::vector<libMeshEnums::FEMNormType> &primal_error_norm_type,
   std::vector<libMeshEnums::FEMNormType> &dual_error_norm_type,
   unsigned int gmv_flag)
{  
  AutoPtr<ErrorEstimator> error_estimator;

  AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;

  error_estimator.reset (adjoint_residual_estimator);

  // We solve the adjoint problem beforehand reusing the preconditioner from the forward solve
  adjoint_residual_estimator->adjoint_already_solved = true;

  if(param.write_gmv_error)
     {
       OStringStream error_gmv;
       error_gmv << "error.gmv.";
       OSSRealzeroright(error_gmv,3,0, a_number);
       error_gmv <<".";
       OSSRealzeroright(error_gmv,1,0, gmv_flag);      

       adjoint_residual_estimator->error_plot_suffix = error_gmv.str();
     }

  PatchRecoveryErrorEstimator *p1 =
    new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->primal_error_estimator().reset(p1);

  PatchRecoveryErrorEstimator *p2 =
    new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->dual_error_estimator().reset(p2);  

  p1->set_patch_reuse(param.patch_reuse);
  p2->set_patch_reuse(param.patch_reuse);

  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(0, primal_error_norm_type[0]);      
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(1, primal_error_norm_type[1]);
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(2, primal_error_norm_type[2]);
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(3, primal_error_norm_type[3]);

  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(0, dual_error_norm_type[0]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(1, dual_error_norm_type[1]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(2, dual_error_norm_type[2]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(3, dual_error_norm_type[3]);  

  // Now set the right weights for each term in the error estimate
  adjoint_residual_estimator->error_norm.set_weight(0, term_weights[0][0]);
  adjoint_residual_estimator->error_norm.set_weight(1, term_weights[1][1]);
  adjoint_residual_estimator->error_norm.set_weight(2, term_weights[2][2]);
  adjoint_residual_estimator->error_norm.set_weight(3, term_weights[3][3]);

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 1, term_weights[0][1]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 2, term_weights[0][2]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 3, term_weights[0][3]);

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 0, term_weights[1][0]);  
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 2, term_weights[1][2]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 3, term_weights[1][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 0, term_weights[2][0]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 1, term_weights[2][1]);  
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 3, term_weights[2][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 0, term_weights[3][0]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 1, term_weights[3][1]);           
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 2, term_weights[3][2]);

  return error_estimator;
}

AutoPtr<ErrorEstimator>
build_weighted_error_estimator_component_wise
  (FEMParameters &param,
   unsigned int a_number,
   std::vector<std::vector<Number> > &term_weights,
   std::vector<libMeshEnums::FEMNormType> &primal_error_norm_type,
   std::vector<libMeshEnums::FEMNormType> &dual_error_norm_type,
   std::vector<FEMFunctionBase<Number>*> coupled_system_weight_functions,
   unsigned int gmv_flag)
{  
  AutoPtr<ErrorEstimator> error_estimator;

  AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;

  error_estimator.reset (adjoint_residual_estimator);

  // We solve the adjoint problem beforehand reusing the preconditioner from the forward solve
  adjoint_residual_estimator->adjoint_already_solved = true;

  if(param.write_gmv_error)
     {
       OStringStream error_gmv;
       error_gmv << "error.gmv.";
       OSSRealzeroright(error_gmv,3,0, a_number);
       error_gmv <<".";
       OSSRealzeroright(error_gmv,1,0, gmv_flag);      

       adjoint_residual_estimator->error_plot_suffix = error_gmv.str();
     }

  WeightedPatchRecoveryErrorEstimator *p1 =
    new WeightedPatchRecoveryErrorEstimator;
  adjoint_residual_estimator->primal_error_estimator().reset(p1);

  PatchRecoveryErrorEstimator *p2 =
    new PatchRecoveryErrorEstimator;
  adjoint_residual_estimator->dual_error_estimator().reset(p2);  

  p1->set_patch_reuse(param.patch_reuse);
  p2->set_patch_reuse(param.patch_reuse);

  p1->weight_functions.clear();

  p1->weight_functions.push_back(coupled_system_weight_functions[0]);
  p1->weight_functions.push_back(coupled_system_weight_functions[1]);
  p1->weight_functions.push_back(coupled_system_weight_functions[2]);
  p1->weight_functions.push_back(coupled_system_weight_functions[3]);  

  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(0, primal_error_norm_type[0]);      
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(1, primal_error_norm_type[1]);
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(2, primal_error_norm_type[2]);
  adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(3, primal_error_norm_type[3]);  

  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(0, dual_error_norm_type[0]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(1, dual_error_norm_type[1]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(2, dual_error_norm_type[2]);
  adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(3, dual_error_norm_type[3]);  

  // Now set the right weights for each term in the error estimate
  adjoint_residual_estimator->error_norm.set_weight(0, term_weights[0][0]);
  adjoint_residual_estimator->error_norm.set_weight(1, term_weights[1][1]);
  adjoint_residual_estimator->error_norm.set_weight(2, term_weights[2][2]);
  adjoint_residual_estimator->error_norm.set_weight(3, term_weights[3][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 1, term_weights[0][1]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 2, term_weights[0][2]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(0, 3, term_weights[0][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 0, term_weights[1][0]);  
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 2, term_weights[1][2]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(1, 3, term_weights[1][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 0, term_weights[2][0]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 1, term_weights[2][1]);  
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(2, 3, term_weights[2][3]);  

  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 0, term_weights[3][0]);
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 1, term_weights[3][1]);           
  adjoint_residual_estimator->error_norm.set_off_diagonal_weight(3, 2, term_weights[3][2]);

  return error_estimator;
}

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  std::cout << "Started " << argv[0] << std::endl;

  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("general.in");
    if (!i)
      {
        std::cerr << '[' << libMesh::processor_id() 
                  << "] Can't find general.in; exiting early." 
                  << std::endl; 
        libmesh_error();
      }
  }

  GetPot infile("general.in");

  // Read in parameters from the input file
  FEMParameters param;
  param.read(infile);

  // Create a mesh.
  Mesh mesh (param.dimension);

  // And an object to refine it
  AutoPtr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  std::cout << "Building mesh" << std::endl;

  if (!param.initial_timestep && param.run_simulation)
    build_domain(mesh, param);

  // Uniformly refine it
  //mesh_refinement->uniformly_refine(1);

  std::cout << "Building system" << std::endl;

  FEMSystem &system = build_system(equation_systems, infile, param);

  set_system_parameters(system, param);

  std::cout << "Initializing systems" << std::endl;

  equation_systems.init ();

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  std::cout<<"Starting adaptive loop"<<std::endl<<std::endl;

  // Count the number of steps used
  unsigned int newton_steps = 0;
  unsigned int krylov_steps = 0;
  unsigned int a_step = 0;

  // Get the linear solver object to set the preconditioner reuse flag
  LinearSolver<Number> *linear_solver = system.get_linear_solver();

  // Adaptively solve the timestep
  for (; a_step <= param.max_adaptivesteps; ++a_step)
    {                                    
      // We have to ensure that we are refining to either an error
      // tolerance or to a target number of elements, and, not to a
      // non-zero tolerance and a target number of elements
      // simultaneously
      if(param.global_tolerance > 0 && param.nelem_target > 0)
        {
          std::cout<<"We cant refine to both a non-zero tolerance and a target number of elements, EXITING adaptive loop. "<<std::endl<<std::endl;
          break;
        }        

      // We dont want to reuse the preconditioner for the forward solve
      linear_solver->reuse_preconditioner(false);

      std::cout<< "Adaptive step " << a_step << std::endl;

      std::cout<< "Solving the forward problem" <<std::endl;

      std::cout << "We have " << mesh.n_active_elem()
                    << " active elements and " << equation_systems.n_active_dofs()
                << " active dofs." << std::endl << std::endl;

      system.solve();

      write_output(equation_systems, 0, a_step, "primal", param);

      // Our QoI is integrated on the side
      system.postprocess_sides = true;

      // Postprocess to get the averaged corner velocities and the QoI
      system.postprocess();

      Number QoI_0_computed = (dynamic_cast<CoupledSystem&>(system)).get_QoI_value();
      std::cout<< "The boundary QoI is " << std::setprecision(17) << QoI_0_computed << std::endl << std::endl;

      // You dont need to compute error estimates and refine at the last
      // adaptive step, only before that
      if(a_step!=param.max_adaptivesteps)
        {
          NumericVector<Number> &primal_solution = (*system.solution);

          // Make sure we get the contributions to the adjoint RHS from the sides
          system.assemble_qoi_sides = true;

          // We are about to solve the adjoint system, but before we
          // do this we see the same preconditioner flag to reuse the
          // preconditioner from the forward solver                
          linear_solver->reuse_preconditioner(param.reuse_preconditioner);

          // Solve the adjoint system
          std::cout<< "Solving the adjoint problem" <<std::endl;
          system.adjoint_solve();

          NumericVector<Number> &dual_solution = system.get_adjoint_solution();
          primal_solution.swap(dual_solution);

          write_output(equation_systems, 0, a_step, "adjoint", param);

          // Swap back
          primal_solution.swap(dual_solution);

          // We need the values of the parameters beta and Pe from the
          // system for the adjoint error estimate
          Number Pe = (dynamic_cast<CoupledSystem&>(system)).get_Pe();

          ErrorVector error;

          // The total error is the sum of the error_non_pressure +
          // error_with_pressure + ...
          // error_estimator_convection_diffusion_x +
          // error_estimator_convection_diffusion_y

          // We construct the non-pressure and pressure contributions
          // to the error separately
          ErrorVector error_non_pressure;

          // First building norm_type_vector_non_pressure and
          // weights_matrix_non_pressure for the non-pressure term
          // error contributions
          std::vector<libMeshEnums::FEMNormType>
            primal_norm_type_vector_non_pressure;                        
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          primal_norm_type_vector_non_pressure.push_back(L2);
          primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);

          std::vector<libMeshEnums::FEMNormType>
            dual_norm_type_vector_non_pressure;                        
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
          dual_norm_type_vector_non_pressure.push_back(L2);
          dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);

          std::vector<std::vector<Number> >
            weights_matrix_non_pressure(system.n_vars(),
              std::vector<Number>(system.n_vars(), 0.0));
          weights_matrix_non_pressure[0][0] = 1.;
          weights_matrix_non_pressure[0][1] = 0.;
          weights_matrix_non_pressure[0][2] = 0.;
          weights_matrix_non_pressure[0][3] = 0.;                  

          weights_matrix_non_pressure[1][0] = 0.;
          weights_matrix_non_pressure[1][1] = 1.;
          weights_matrix_non_pressure[1][2] = 0.;
          weights_matrix_non_pressure[1][3] = 0.;                  

          weights_matrix_non_pressure[2][0] = 0.;
          weights_matrix_non_pressure[2][1] = 0.;
          weights_matrix_non_pressure[2][2] = 0.;
          weights_matrix_non_pressure[2][3] = 0.;


          weights_matrix_non_pressure[3][0] = 0.;
          weights_matrix_non_pressure[3][1] = 0.;
          weights_matrix_non_pressure[3][2] = 0.;
          weights_matrix_non_pressure[3][3] = 1./Pe;

          // Build the error estimator to estimate the contributions
          // to the QoI error from the non pressure term
          AutoPtr<ErrorEstimator> error_estimator_non_pressure = 
            build_error_estimator_component_wise
              (param, a_step, weights_matrix_non_pressure,
               primal_norm_type_vector_non_pressure,
               dual_norm_type_vector_non_pressure, 10);

          // Estimate the contributions to the QoI error from the non
          // pressure terms
          error_estimator_non_pressure->estimate_error(system, error_non_pressure);

          write_error(equation_systems, error_non_pressure, 0, a_step, param, "_non_pressure");

          // Now for the pressure contribs
          ErrorVector error_with_pressure;

          // Now building norm_type_vector_with_pressure and
          // weights_matrix_with_pressure for the pressure term error
          // contributions
          std::vector<libMeshEnums::FEMNormType>
            primal_norm_type_vector_with_pressure;                        
          primal_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
          primal_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
          primal_norm_type_vector_with_pressure.push_back(L2);
          primal_norm_type_vector_with_pressure.push_back(L2);

          std::vector<libMeshEnums::FEMNormType>
            dual_norm_type_vector_with_pressure;                        
          dual_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
          dual_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
          dual_norm_type_vector_with_pressure.push_back(L2);
          dual_norm_type_vector_with_pressure.push_back(L2);

          std::vector<std::vector<Number> >
            weights_matrix_with_pressure
              (system.n_vars(),
               std::vector<Number>(system.n_vars(), 0.0));
          weights_matrix_with_pressure[0][0] = 0.;
          weights_matrix_with_pressure[0][1] = 0.;
          weights_matrix_with_pressure[0][2] = 1.;
          weights_matrix_with_pressure[0][3] = 0.;                  

          weights_matrix_with_pressure[1][0] = 0.;
          weights_matrix_with_pressure[1][1] = 0.;
          weights_matrix_with_pressure[1][2] = 1.;
          weights_matrix_with_pressure[1][3] = 0.;                  

          weights_matrix_with_pressure[2][0] = 1.;
          weights_matrix_with_pressure[2][1] = 1.;
          weights_matrix_with_pressure[2][2] = 0.;
          weights_matrix_with_pressure[2][3] = 0.;                  

          weights_matrix_with_pressure[3][0] = 0.;
          weights_matrix_with_pressure[3][1] = 0.;
          weights_matrix_with_pressure[3][2] = 0.;
          weights_matrix_with_pressure[3][3] = 0.;                  

          // Build the error estimator to estimate the contributions
          // to the QoI error from the pressure term
          AutoPtr<ErrorEstimator> error_estimator_with_pressure =
            build_error_estimator_component_wise
              (param, a_step, weights_matrix_with_pressure,
               primal_norm_type_vector_with_pressure,
               dual_norm_type_vector_with_pressure, 20);

          // Estimate the contributions to the QoI error from the pressure terms  
          error_estimator_with_pressure->estimate_error(system, error_with_pressure);

          write_error(equation_systems, error_with_pressure, 0, a_step, param, "_with_pressure");

          // Get the error contribution from the convection diffusion
          // equation weighted with the dual concentration
          ErrorVector error_convection_diffusion_x;

          // Now building norm_type_vector_with_pressure and
          // weights_matrix_with_pressure for the pressure term error
          // contributions
          std::vector<libMeshEnums::FEMNormType>
            primal_norm_type_vector_convection_diffusion_x;                        
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(L2);
          primal_norm_type_vector_convection_diffusion_x.push_back(H1_X_SEMINORM);

          std::vector<libMeshEnums::FEMNormType>
            dual_norm_type_vector_convection_diffusion_x;                        
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);
          // Note that we need the error of the dual concentration in L2
          dual_norm_type_vector_convection_diffusion_x.push_back(L2);

          std::vector<std::vector<Number> >
            weights_matrix_convection_diffusion_x
              (system.n_vars(),
               std::vector<Number>(system.n_vars(), 0.0));
          weights_matrix_convection_diffusion_x[0][0] = 0.;
          weights_matrix_convection_diffusion_x[0][1] = 0.;
          weights_matrix_convection_diffusion_x[0][2] = 0.;
          weights_matrix_convection_diffusion_x[0][3] = 0.;                  

          weights_matrix_convection_diffusion_x[1][0] = 0.;
          weights_matrix_convection_diffusion_x[1][1] = 0.;
          weights_matrix_convection_diffusion_x[1][2] = 0.;
          weights_matrix_convection_diffusion_x[1][3] = 0.;                  

          weights_matrix_convection_diffusion_x[2][0] = 0.;
          weights_matrix_convection_diffusion_x[2][1] = 0.;
          weights_matrix_convection_diffusion_x[2][2] = 0.;
          weights_matrix_convection_diffusion_x[2][3] = 0.;                  

          weights_matrix_convection_diffusion_x[3][0] = 0.;
          weights_matrix_convection_diffusion_x[3][1] = 0.;
          weights_matrix_convection_diffusion_x[3][2] = 0.;
          weights_matrix_convection_diffusion_x[3][3] = 1.;                  

          ConstFEMFunction<Number> identity(1);

          // Declare objects of class CoupledFEMFunctionsx
          CoupledFEMFunctionsx convdiffx(system);                  

          // Make a vector of pointers to these objects
          std::vector<FEMFunctionBase<Number> *> coupled_system_weight_functions_x;
          coupled_system_weight_functions_x.push_back(&identity);
          coupled_system_weight_functions_x.push_back(&identity);
          coupled_system_weight_functions_x.push_back(&identity);                  
          coupled_system_weight_functions_x.push_back(&convdiffx);                  

          // Build the error estimator to estimate the contributions
          // to the QoI error from the pressure term
          AutoPtr<ErrorEstimator> error_estimator_convection_diffusion_x =
          build_weighted_error_estimator_component_wise
            (param, a_step, weights_matrix_convection_diffusion_x,
             primal_norm_type_vector_convection_diffusion_x,
             dual_norm_type_vector_convection_diffusion_x,
             coupled_system_weight_functions_x, 30);

          // Estimate the contributions to the QoI error from the
          // convection diffusion terms  
          error_estimator_convection_diffusion_x->estimate_error
            (system, error_convection_diffusion_x);

          write_error(equation_systems, error_convection_diffusion_x,
                      0, a_step, param, "_convection_diffusion_x");

          ErrorVector error_convection_diffusion_y;

          // Now building norm_type_vector_with_pressure and
          // weights_matrix_with_pressure for the pressure term error
          // contributions
          std::vector<libMeshEnums::FEMNormType>
            primal_norm_type_vector_convection_diffusion_y;                        
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(L2);
          primal_norm_type_vector_convection_diffusion_y.push_back(H1_Y_SEMINORM);

          std::vector<libMeshEnums::FEMNormType>
            dual_norm_type_vector_convection_diffusion_y;                        
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);
          // Note that we need the error of the dual concentration in L2
          dual_norm_type_vector_convection_diffusion_y.push_back(L2);

          std::vector<std::vector<Number> >
            weights_matrix_convection_diffusion_y
              (system.n_vars(), std::vector<Number>(system.n_vars(), 0.0));
          weights_matrix_convection_diffusion_y[0][0] = 0.;
          weights_matrix_convection_diffusion_y[0][1] = 0.;
          weights_matrix_convection_diffusion_y[0][2] = 0.;
          weights_matrix_convection_diffusion_y[0][3] = 0.;                  

          weights_matrix_convection_diffusion_y[1][0] = 0.;
          weights_matrix_convection_diffusion_y[1][1] = 0.;
          weights_matrix_convection_diffusion_y[1][2] = 0.;
          weights_matrix_convection_diffusion_y[1][3] = 0.;                  

          weights_matrix_convection_diffusion_y[2][0] = 0.;
          weights_matrix_convection_diffusion_y[2][1] = 0.;
          weights_matrix_convection_diffusion_y[2][2] = 0.;
          weights_matrix_convection_diffusion_y[2][3] = 0.;                  

          weights_matrix_convection_diffusion_y[3][0] = 0.;
          weights_matrix_convection_diffusion_y[3][1] = 0.;
          weights_matrix_convection_diffusion_y[3][2] = 0.;
          weights_matrix_convection_diffusion_y[3][3] = 1.;                  

          // Declare objects of class CoupledFEMFunctionsy                  
          CoupledFEMFunctionsy convdiffy(system);                  

          // Make a vector of pointers to these objects
          std::vector<FEMFunctionBase<Number> *> coupled_system_weight_functions_y;
          coupled_system_weight_functions_y.push_back(&identity);
          coupled_system_weight_functions_y.push_back(&identity);
          coupled_system_weight_functions_y.push_back(&identity);
          coupled_system_weight_functions_y.push_back(&convdiffy);                  

          // Build the error estimator to estimate the contributions
          // to the QoI error from the pressure term
          AutoPtr<ErrorEstimator> error_estimator_convection_diffusion_y =
            build_weighted_error_estimator_component_wise
              (param, a_step, weights_matrix_convection_diffusion_y,
               primal_norm_type_vector_convection_diffusion_y,
               dual_norm_type_vector_convection_diffusion_y,
               coupled_system_weight_functions_y, 40);

          // Estimate the contributions to the QoI error from the
          // convection diffusion terms  
          error_estimator_convection_diffusion_y->estimate_error(system, error_convection_diffusion_y);

          write_error(equation_systems, error_convection_diffusion_y, 0, a_step, param, "_convection_diffusion_y");

          libmesh_assert_equal_to(error_non_pressure.size(), error_with_pressure.size()); 
          libmesh_assert_equal_to(error_convection_diffusion_x.size(), error_non_pressure.size());

          if(param.indicator_type == "adjoint_residual")
            {
              error.resize(error_non_pressure.size());

              // Now combine the contribs from the pressure and non
              // pressure terms to get the complete estimate
              for(unsigned int i = 0; i < error.size(); i++)
                {                          
                  error[i] = error_non_pressure[i] + error_with_pressure[i] + error_convection_diffusion_x[i] + error_convection_diffusion_y[i];
                }
            }
          else
            {
              std::cout<<"Using Kelly Estimator"<<std::endl;
              // Build the Kelly error estimator
              AutoPtr<ErrorEstimator> error_estimator =  build_error_estimator(param);

              // Estimate the error
              error_estimator->estimate_error(system, error);
            }

          write_error(equation_systems, error, 0, a_step, param, "_total");

          // We have to refine either based on reaching an error
          // tolerance or a number of elements target, which should be
          // verified above

          // Otherwise we flag elements by error tolerance or nelem
          // target

          if(param.refine_uniformly)
            {
              mesh_refinement->uniformly_refine(1);
            }
          else if(param.global_tolerance >= 0. && param.nelem_target == 0.) // Refine based on reaching an error tolerance
            {                              
              mesh_refinement->flag_elements_by_error_tolerance (error);              

              // Carry out the adaptive mesh refinement/coarsening
              mesh_refinement->refine_and_coarsen_elements();
            }
          else if(param.alternate_with_uniform_steps) // Do alternate adaptive and uniform steps
            {
              //Refine uniformly every param.alternate_step_number, For eg. if param.alternate_step_number = 5, every 5th step will be uniform
              if((a_step !=0 && a_step >= 10) && ((a_step % param.alternate_step_number == 0) || (a_step % param.alternate_step_number == 1)) )
                {
                  std::cout<<"Not refining adaptively ! Doing a uniform refinement step !"<<std::endl<<std::endl;
                  mesh_refinement->uniformly_refine(1);

                }
              else // Refine adaptively on the other steps
                {
                  std::cout<<"Refining adaptively based on n elem target"<<std::endl<<std::endl;
                  mesh_refinement->flag_elements_by_nelem_target (error);              
                  // Carry out the adaptive mesh refinement/coarsening
                  mesh_refinement->refine_and_coarsen_elements();
                }
            }
          else // Refine based on reaching a target number of elements
            {                                
              //If we have reached the desired number of elements, we
              //dont do any more adaptive steps
              if (mesh.n_active_elem() >= param.nelem_target)
                {
                  std::cout<<"We reached the target number of elements."<<std::endl <<std::endl;
                  break;
                }

              mesh_refinement->flag_elements_by_nelem_target (error);

              // Carry out the adaptive mesh refinement/coarsening
              mesh_refinement->refine_and_coarsen_elements();
            }

          equation_systems.reinit();

        }

    }


  write_output_footers(param);

  // All done.  
  return 0;    

}

