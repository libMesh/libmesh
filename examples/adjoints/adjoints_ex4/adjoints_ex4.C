// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Adjoints Example 4 - Laplace Equation in the L-Shaped Domain with AdjointRefinementErrorEstimator</h1>
// \author Vikram Garg
// \date 2012
//
// This example solves the Laplace equation on the classic "L-shaped"
// domain with adaptive mesh refinement. The exact
// solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). The kelly and
// adjoint residual error estimators are used to develop error indicators and
// guide mesh adaptation. Since we use the adjoint capabilities of libMesh in
// this example, we use the DiffSystem framework. This file (adjoints_ex1.C)
// contains the declaration of mesh and equation system objects, L-shaped.C
// contains the assembly of the system, element_qoi_derivative.C and
// side_qoi_derivative.C contain the RHS for the adjoint systems.
// Postprocessing to compute the QoIs is done in element_postprocess.C and
// side_postprocess.C.

// The initial mesh contains three QUAD9 elements which represent the
// standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
// i.e.
// Element 0: [-1,0]x[ 0,1]
// Element 1: [ 0,1]x[ 0,1]
// Element 2: [-1,0]x[-1,0]
// The mesh is provided in the standard libMesh ASCII format file
// named "lshaped.xda".  In addition, an input file named "general.in"
// is provided which allows the user to set several parameters for
// the solution so that the problem can be re-run without a
// re-compile.  The solution technique employed is to have a
// refinement loop with a linear (forward and adjoint) solve inside followed by a
// refinement of the grid and projection of the solution to the new grid
// In the final loop iteration, there is no additional
// refinement after the solve.  In the input file "general.in", the variable
// "max_adaptivesteps" controls the number of refinement steps, and
// "refine_fraction" / "coarsen_fraction" determine the number of
// elements which will be refined / coarsened at each step.

// C++ includes
#include <iostream>
#include <iomanip>

// General libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/steady_solver.h"
#include "libmesh/system_norm.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// Adjoint Related includes
#include "libmesh/qoi_set.h"
#include "libmesh/adjoint_refinement_estimator.h"

// libMesh I/O includes
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exodusII_io.h"

// Local includes
#include "femparameters.h"
#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Local function declarations

// Number output files, the files are give a prefix of primal or adjoint_i depending on
// whether the output is the primal solution or the dual solution for the ith QoI

// Optionally write different types of output files.
void write_output(EquationSystems & es,
                  unsigned int a_step,       // The adaptive step count
                  std::string solution_type, // primal or adjoint solve
                  FEMParameters & param)
{
  // Ignore parameters when there are no output formats available.
  libmesh_ignore(es);
  libmesh_ignore(a_step);
  libmesh_ignore(solution_type);
  libmesh_ignore(param);

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      MeshBase & mesh = es.get_mesh();

      std::ostringstream file_name_gmv;
      file_name_gmv << solution_type
                    << ".out.gmv."
                    << std::setw(2)
                    << std::setfill('0')
                    << std::right
                    << a_step;

      GMVIO(mesh).write_equation_systems
        (file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_EXODUS_API
  if (param.output_exodus)
    {
      MeshBase & mesh = es.get_mesh();

      // We write out one file per adaptive step. The files are named in
      // the following way:
      // foo.e
      // foo.e-s002
      // foo.e-s003
      // ...
      // so that, if you open the first one with Paraview, it actually
      // opens the entire sequence of adapted files.
      std::ostringstream file_name_exodus;

      file_name_exodus << solution_type << ".e";
      if (a_step > 0)
        file_name_exodus << "-s"
                         << std::setw(3)
                         << std::setfill('0')
                         << std::right
                         << a_step + 1;

      // We write each adaptive step as a pseudo "time" step, where the
      // time simply matches the (1-based) adaptive step we are on.
      ExodusII_IO(mesh).write_timestep(file_name_exodus.str(),
                                       es,
                                       1,
                                       /*time=*/a_step + 1);
    }
#endif
}

// Set the parameters for the nonlinear and linear solvers to be used during the simulation

void set_system_parameters(LaplaceSystem & system,
                           FEMParameters & param)
{
  // Use analytical jacobians?
  system.analytic_jacobians() = param.analytic_jacobians;

  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;

  // Use the prescribed FE type
  system.fe_family() = param.fe_family[0];
  system.fe_order() = param.fe_order[0];

  // More desperate debugging options
  system.print_solution_norms = param.print_solution_norms;
  system.print_solutions      = param.print_solutions;
  system.print_residual_norms = param.print_residual_norms;
  system.print_residuals      = param.print_residuals;
  system.print_jacobian_norms = param.print_jacobian_norms;
  system.print_jacobians      = param.print_jacobians;

  // No transient time solver
  system.time_solver = libmesh_make_unique<SteadySolver>(system);

  // Nonlinear solver options
  {
    NewtonSolver * solver = new NewtonSolver(system);
    system.time_solver->diff_solver() = std::unique_ptr<DiffSolver>(solver);

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

// Build the mesh refinement object and set parameters for refining/coarsening etc

#ifdef LIBMESH_ENABLE_AMR

std::unique_ptr<MeshRefinement> build_mesh_refinement(MeshBase & mesh,
                                                      FEMParameters & param)
{
  MeshRefinement * mesh_refinement = new MeshRefinement(mesh);
  mesh_refinement->coarsen_by_parents() = true;
  mesh_refinement->absolute_global_tolerance() = param.global_tolerance;
  mesh_refinement->nelem_target()      = param.nelem_target;
  mesh_refinement->refine_fraction()   = param.refine_fraction;
  mesh_refinement->coarsen_fraction()  = param.coarsen_fraction;
  mesh_refinement->coarsen_threshold() = param.coarsen_threshold;

  return std::unique_ptr<MeshRefinement>(mesh_refinement);
}

// This is where declare the adjoint refined error estimator. This estimator builds an error bound
// for Q(u) - Q(u_h), by solving the adjoint problem on a finer Finite Element space. For more details
// see the description of the Adjoint Refinement Error Estimator in adjoint_refinement_error_estimator.C
std::unique_ptr<AdjointRefinementEstimator> build_adjoint_refinement_error_estimator(QoISet & qois)
{
  libMesh::out << "Computing the error estimate using the Adjoint Refinement Error Estimator" << std::endl << std::endl;

  AdjointRefinementEstimator *adjoint_refinement_estimator = new AdjointRefinementEstimator;

  adjoint_refinement_estimator->qoi_set() = qois;

  // We enrich the FE space for the dual problem by doing 2 uniform h refinements
  adjoint_refinement_estimator->number_h_refinements = 2;

  return std::unique_ptr<AdjointRefinementEstimator>(adjoint_refinement_estimator);
}

#endif // LIBMESH_ENABLE_AMR


// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // This doesn't converge with Eigen BICGSTAB for some reason...
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc");

  libMesh::out << "Started " << argv[0] << std::endl;

  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("general.in");
    if (!i)
      libmesh_error_msg('[' << init.comm().rank() << "] Can't find general.in; exiting early.");
  }
  GetPot infile("general.in");

  // Read in parameters from the input file
  FEMParameters param(init.comm());
  param.read(infile);

  // Skip this default-2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // And an object to refine it
  std::unique_ptr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  libMesh::out << "Reading in and building the mesh" << std::endl;

  // Read in the mesh
  mesh.read(param.domainfile.c_str());
  // Make all the elements of the mesh second order so we can compute
  // with a higher order basis
  mesh.all_second_order();

  // Create a mesh refinement object to do the initial uniform refinements
  // on the coarse grid read in from lshaped.xda
  MeshRefinement initial_uniform_refinements(mesh);
  initial_uniform_refinements.uniformly_refine(param.coarserefinements);

  libMesh::out << "Building system" << std::endl;

  // Build the FEMSystem
  LaplaceSystem & system = equation_systems.add_system<LaplaceSystem> ("LaplaceSystem");

  // Set its parameters
  set_system_parameters(system, param);

  libMesh::out << "Initializing systems" << std::endl;

  equation_systems.init ();

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();
  LinearSolver<Number> *linear_solver = system.get_linear_solver();

  {
    // Adaptively solve the timestep
    unsigned int a_step = 0;
    for (; a_step != param.max_adaptivesteps; ++a_step)
      {
        // We can't adapt to both a tolerance and a
        // target mesh size
        if (param.global_tolerance != 0.)
          libmesh_assert_equal_to (param.nelem_target, 0);
        // If we aren't adapting to a tolerance we need a
        // target mesh size
        else
          libmesh_assert_greater (param.nelem_target, 0);

        linear_solver->reuse_preconditioner(false);

        // Solve the forward problem
        system.solve();

        // Write out the computed primal solution
        write_output(equation_systems, a_step, "primal", param);

        // Get a pointer to the primal solution vector
        NumericVector<Number> & primal_solution = *system.solution;

        // Declare a QoISet object, we need this object to set weights for our QoI error contributions
        QoISet qois;

        // Declare a qoi_indices vector, each index will correspond to a QoI
        std::vector<unsigned int> qoi_indices;
        qoi_indices.push_back(0);
        qoi_indices.push_back(1);
        qois.add_indices(qoi_indices);

        // Set weights for each index, these will weight the contribution of each QoI in the final error
        // estimate to be used for flagging elements for refinement
        qois.set_weight(0, 0.5);
        qois.set_weight(1, 0.5);

        // Make sure we get the contributions to the adjoint RHS from the sides
        system.assemble_qoi_sides = true;

        // We are about to solve the adjoint system, but before we do this we see the same preconditioner
        // flag to reuse the preconditioner from the forward solver
        linear_solver->reuse_preconditioner(param.reuse_preconditioner);

        // Solve the adjoint system. This takes the transpose of the stiffness matrix and then
        // solves the resulting system
        system.adjoint_solve();

        // Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unnecessarily in the error estimator
        system.set_adjoint_already_solved(true);

        // Get a pointer to the solution vector of the adjoint problem for QoI 0
        NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

        // Swap the primal and dual solutions so we can write out the adjoint solution
        primal_solution.swap(dual_solution_0);
        write_output(equation_systems, a_step, "adjoint_0", param);

        // Swap back
        primal_solution.swap(dual_solution_0);

        // Get a pointer to the solution vector of the adjoint problem for QoI 0
        NumericVector<Number> & dual_solution_1 = system.get_adjoint_solution(1);

        // Swap again
        primal_solution.swap(dual_solution_1);
        write_output(equation_systems, a_step, "adjoint_1", param);

        // Swap back again
        primal_solution.swap(dual_solution_1);

        libMesh::out << "Adaptive step "
                     << a_step
                     << ", we have "
                     << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs."
                     << std::endl;

        // Postprocess, compute the approximate QoIs and write them out to the console
        libMesh::out << "Postprocessing: " << std::endl;
        system.postprocess_sides = true;
        system.postprocess();
        Number QoI_0_computed = system.get_QoI_value("computed", 0);
        Number QoI_0_exact = system.get_QoI_value("exact", 0);
        Number QoI_1_computed = system.get_QoI_value("computed", 1);
        Number QoI_1_exact = system.get_QoI_value("exact", 1);

        libMesh::out << "The relative error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(QoI_0_computed - QoI_0_exact) / std::abs(QoI_0_exact)
                     << std::endl;

        libMesh::out << "The relative error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(QoI_1_computed - QoI_1_exact) / std::abs(QoI_1_exact)
                     << std::endl
                     << std::endl;

        // We will declare an error vector for passing to the adjoint refinement error estimator
        ErrorVector QoI_elementwise_error;

        // Build an adjoint refinement error estimator object
        std::unique_ptr<AdjointRefinementEstimator> adjoint_refinement_error_estimator =
          build_adjoint_refinement_error_estimator(qois);

        // Estimate the error in each element using the Adjoint Refinement estimator
        adjoint_refinement_error_estimator->estimate_error(system, QoI_elementwise_error);

        // Print out the computed error estimate, note that we access the global error estimates
        // using an accessor function, right now sum(QoI_elementwise_error) != global_QoI_error_estimate
        libMesh::out << "The computed relative error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_exact)
                     << std::endl;

        libMesh::out << "The computed relative error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_exact)
                     << std::endl
                     << std::endl;

        // Also print out effectivity indices (estimated error/true error)
        libMesh::out << "The effectivity index for the computed error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact)
                     << std::endl;

        libMesh::out << "The effectivity index for the computed error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact)
                     << std::endl
                     << std::endl;

        // For refinement purposes we need to sort by error
        // *magnitudes*, but AdjointRefinement gives us signed errors.
        if (!param.refine_uniformly)
          for (std::size_t i=0; i<QoI_elementwise_error.size(); i++)
            if (QoI_elementwise_error[i] != 0.)
              QoI_elementwise_error[i] = std::abs(QoI_elementwise_error[i]);

        // We have to refine either based on reaching an error tolerance or
        // a number of elements target, which should be verified above
        // Otherwise we flag elements by error tolerance or nelem target

        // Uniform refinement
        if (param.refine_uniformly)
          {
            mesh_refinement->uniformly_refine(1);
          }
        // Adaptively refine based on reaching an error tolerance
        else if (param.global_tolerance >= 0. && param.nelem_target == 0.)
          {
            mesh_refinement->flag_elements_by_error_tolerance (QoI_elementwise_error);

            mesh_refinement->refine_and_coarsen_elements();
          }
        // Adaptively refine based on reaching a target number of elements
        else
          {
            if (mesh.n_active_elem() >= param.nelem_target)
              {
                libMesh::out << "We reached the target number of elements." << std::endl << std::endl;
                break;
              }

            mesh_refinement->flag_elements_by_nelem_target (QoI_elementwise_error);

            mesh_refinement->refine_and_coarsen_elements();
          }

        // Dont forget to reinit the system after each adaptive refinement !
        equation_systems.reinit();

        libMesh::out << "Refined mesh to "
                     << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs."
                     << std::endl;
      }

    // Do one last solve if necessary
    if (a_step == param.max_adaptivesteps)
      {
        linear_solver->reuse_preconditioner(false);
        system.solve();

        write_output(equation_systems, a_step, "primal", param);

        NumericVector<Number> & primal_solution = *system.solution;

        QoISet qois;
        std::vector<unsigned int> qoi_indices;

        qoi_indices.push_back(0);
        qoi_indices.push_back(1);
        qois.add_indices(qoi_indices);

        qois.set_weight(0, 0.5);
        qois.set_weight(1, 0.5);

        system.assemble_qoi_sides = true;
        linear_solver->reuse_preconditioner(param.reuse_preconditioner);
        system.adjoint_solve();

        // Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unnecessarily in the error estimator
        system.set_adjoint_already_solved(true);

        NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

        primal_solution.swap(dual_solution_0);
        write_output(equation_systems, a_step, "adjoint_0", param);

        primal_solution.swap(dual_solution_0);

        NumericVector<Number> & dual_solution_1 = system.get_adjoint_solution(1);

        primal_solution.swap(dual_solution_1);
        write_output(equation_systems, a_step, "adjoint_1", param);

        primal_solution.swap(dual_solution_1);

        libMesh::out << "Adaptive step "
                     << a_step
                     << ", we have "
                     << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs."
                     << std::endl;

        libMesh::out << "Postprocessing: " << std::endl;
        system.postprocess_sides = true;
        system.postprocess();

        Number QoI_0_computed = system.get_QoI_value("computed", 0);
        Number QoI_0_exact = system.get_QoI_value("exact", 0);
        Number QoI_1_computed = system.get_QoI_value("computed", 1);
        Number QoI_1_exact = system.get_QoI_value("exact", 1);

        libMesh::out << "The relative error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(QoI_0_computed - QoI_0_exact) / std::abs(QoI_0_exact)
                     << std::endl;

        libMesh::out << "The relative error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(QoI_1_computed - QoI_1_exact) / std::abs(QoI_1_exact)
                     << std::endl
                     << std::endl;

        // We will declare an error vector for passing to the adjoint refinement error estimator
        // Right now, only the first entry of this vector will be filled (with the global QoI error estimate)
        // Later, each entry of the vector will contain elementwise error that the user can sum to get the total error
        ErrorVector QoI_elementwise_error;

        // Build an adjoint refinement error estimator object
        std::unique_ptr<AdjointRefinementEstimator> adjoint_refinement_error_estimator =
          build_adjoint_refinement_error_estimator(qois);

        // Estimate the error in each element using the Adjoint Refinement estimator
        adjoint_refinement_error_estimator->estimate_error(system, QoI_elementwise_error);

        // Print out the computed error estimate, note that we access the global error estimates
        // using an accessor function, right now sum(QoI_elementwise_error) != global_QoI_error_estimate
        libMesh::out << "The computed relative error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_exact)
                     << std::endl;

        libMesh::out << "The computed relative error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_exact)
                     << std::endl
                     << std::endl;

        // Also print out effectivity indices (estimated error/true error)
        libMesh::out << "The effectivity index for the computed error in QoI 0 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact)
                     << std::endl;

        libMesh::out << "The effectivity index for the computed error in QoI 1 is "
                     << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact)
                     << std::endl
                     << std::endl;

        // Hard coded assert to ensure that the actual numbers we are getting are what they should be

        // The effectivity index isn't exactly reproducible at single precision
        // libmesh_assert_less(std::abs(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact) - 0.84010976704434637), 1.e-5);
        // libmesh_assert_less(std::abs(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact) - 0.48294428289950514), 1.e-5);

        // But the effectivity indices should always be sane
        libmesh_assert_less(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact), 2.5);
        libmesh_assert_greater(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact), .4);
        libmesh_assert_less(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact), 2.5);
        libmesh_assert_greater(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact), .4);

        // And the computed errors should still be low
        libmesh_assert_less(std::abs(QoI_0_computed - QoI_0_exact), 2e-4);
        libmesh_assert_less(std::abs(QoI_1_computed - QoI_1_exact), 2e-4);
      }
  }

  libMesh::err << '[' << mesh.processor_id()
               << "] Completing output."
               << std::endl;

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
