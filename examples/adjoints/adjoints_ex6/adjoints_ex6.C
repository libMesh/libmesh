// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Adjoints Example 6 - Poisson Equation in square domain with AdjointRefinementErrorEstimator</h1>
// \author Vikram Garg
// \date 2017
//
// This example solves the Poisson equation, whose solution displays a
// sharp layer, with QoI based adjoint error estimation and adaptive
// mesh refinement. The exact QoI value is in poisson.in. This example
// also illustrates the use of the adjoint Dirichlet boundary
// condition capability, necessary for handling flux QoIs. We access
// the adjoint capabilities of libMesh via the DiffSystem
// framework. This file (adjoints_ex6.C) contains the declaration of
// mesh and equation system objects, poissonsystem.C contains the
// assembly of the system. Postprocessing to compute the QoI is done
// in element_postprocess.C.  There is no need for element and side
// qoi derivative functions, since the adjoint RHS is supplied by the
// adjoint dirichlet boundary condition.

// WARNING: Adjoint-dirichlet based weighted flux quantities of
// interest are computed internally by FEMSystem as R^h(u^h, L^h)
// where R^h is the physics residual, u^h is the solution, and L^h the
// lift function at the current mesh.  This is appropriate for
// unstabilized weak residuals.  Non-FEMSystem users and stabilized
// method users will need to override the default flux QoI behavior.

// An input file named "general.in" is provided which allows the user
// to set several parameters for the solution so that the problem can
// be re-run without a re-compile.  The solution technique employed is
// to have a refinement loop with a linear (forward and adjoint) solve
// inside followed by a refinement of the grid and projection of the
// solution to the new grid In the final loop iteration, there is no
// additional refinement after the solve.  In the input file
// "general.in", the variable "max_adaptivesteps" controls the number
// of refinement steps, and "refine_fraction" / "coarsen_fraction"
// determine the number of elements which will be refined / coarsened
// at each step.

// C++ includes
#include <iostream>
#include <iomanip>

// General libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/newton_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/steady_solver.h"
#include "libmesh/system_norm.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/enum_solver_package.h"

// Adjoint Related includes
#include "libmesh/qoi_set.h"
#include "libmesh/adjoint_refinement_estimator.h"

// libMesh I/O includes
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"

// Local includes
#include "femparameters.h"
#include "poisson.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Local function declarations

// Number output files, the files are give a prefix of primal or adjoint_i depending on
// whether the output is the primal solution or the dual solution for the ith QoI

// Write gmv output
void write_output(EquationSystems & es,
                  unsigned int a_step,       // The adaptive step count
                  std::string solution_type) // primal or adjoint solve
{
#ifdef LIBMESH_HAVE_GMV
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
#endif
}

// Set the parameters for the nonlinear and linear solvers to be used during the simulation
void set_system_parameters(PoissonSystem & system, FEMParameters & param)
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


#ifdef LIBMESH_ENABLE_AMR

// Build the mesh refinement object and set parameters for refining/coarsening etc
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


// This is where declare the adjoint refined error estimator. This
// estimator builds an error bound for Q(u) - Q(u_h), by solving the
// adjoint problem on a finer Finite Element space. For more details
// see the description of the Adjoint Refinement Error Estimator in
// adjoint_refinement_error_estimator.C
std::unique_ptr<AdjointRefinementEstimator> build_adjoint_refinement_error_estimator(QoISet & qois)
{
  libMesh::out << "Computing the error estimate using the Adjoint Refinement Error Estimator\n" << std::endl;

  AdjointRefinementEstimator * adjoint_refinement_estimator = new AdjointRefinementEstimator;

  adjoint_refinement_estimator->qoi_set() = qois;

  // We enrich the FE space for the dual problem by doing 2 uniform h refinements
  adjoint_refinement_estimator->number_h_refinements = 2;

  return std::unique_ptr<AdjointRefinementEstimator>(adjoint_refinement_estimator);
}

#endif // LIBMESH_ENABLE_AMR


// The main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // This doesn't converge with Eigen BICGSTAB for some reason...
  libmesh_example_requires((libMesh::default_solver_package() != EIGEN_SOLVERS) &&
                           (libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE),
                           "--enable-petsc or --enable-trilinos");

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
  MeshTools::Generation::build_square
    (mesh, 2, 2,
     0.0, 1.0,
     0.0, 1.0,
     QUAD4);

  // Make all the elements of the mesh second order so we can compute
  // with a higher order basis
  mesh.all_second_order();

  // Create a mesh refinement object to do the initial uniform refinements
  // on the coarse grid read in from lshaped.xda
  MeshRefinement initial_uniform_refinements(mesh);
  initial_uniform_refinements.uniformly_refine(param.coarserefinements);

  libMesh::out << "Building system" << std::endl;

  // Build the FEMSystem
  PoissonSystem & system = equation_systems.add_system<PoissonSystem> ("PoissonSystem");

  // Set its parameters
  set_system_parameters(system, param);

  libMesh::out << "Initializing systems" << std::endl;

  equation_systems.init ();

  // Add an adjoint_solution0 vector to the system
  system.add_vector("adjoint_solution0", false, GHOSTED);

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  // Get a pointer to the linear solver object to be able to reuse preconditioner
  LinearSolver<Number> * linear_solver = system.get_linear_solver();

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

        // Dont reuse preconditioners before the primal solve
        linear_solver->reuse_preconditioner(false);

        // Solve the forward problem
        system.solve();

        // Write out the computed primal solution
        write_output(equation_systems, a_step, "primal");

        // Declare a QoISet object, we need this object to set weights for our QoI error contributions
        QoISet qois;

        // Declare a qoi_indices vector, each index will correspond to a QoI
        std::vector<unsigned int> qoi_indices;
        qoi_indices.push_back(0);
        qois.add_indices(qoi_indices);

        // Set weights for each index, these will weight the contribution of each QoI in the final error
        // estimate to be used for flagging elements for refinement
        qois.set_weight(0, 1.0);

        // Make sure we get the contributions to the adjoint RHS from the sides
        system.assemble_qoi_sides = true;

        // We are about to solve the adjoint system, but before we do this we see the same preconditioner
        // flag to reuse the preconditioner from the forward solver
        linear_solver->reuse_preconditioner(param.reuse_preconditioner);

        // Solve the adjoint system. This takes the transpose of the stiffness matrix and then
        // solves the resulting system
        system.adjoint_solve();

        //Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unnecessarily in the error estimator
        system.set_adjoint_already_solved(true);

        // Get a pointer to the primal solution vector
        NumericVector<Number> & primal_solution = *system.solution;

        //Get a pointer to the solution vector of the adjoint problem for QoI 0
        NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

        //Swap the primal and dual solutions so we can write out the adjoint solution
        primal_solution.swap(dual_solution_0);
        write_output(equation_systems, a_step, "adjoint_0");

        //Swap back
        primal_solution.swap(dual_solution_0);

        libMesh::out << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs." << std::endl ;

        // Postprocess, compute the approximate QoIs and write them out to the console
        libMesh::out << "Postprocessing: " << std::endl;
        system.postprocess_sides = true;
        system.postprocess();

        Number QoI_0_computed = system.get_QoI_value("computed", 0);
        Number QoI_0_exact = system.get_QoI_value("exact", 0);

        libMesh::out << "The computed QoI 0 is " << std::setprecision(17)
                     << QoI_0_computed << std::endl;
        libMesh::out << "The relative error in QoI 0 is " << std::setprecision(17)
                     << std::abs(QoI_0_computed - QoI_0_exact) << std::endl; // / std::abs(QoI_0_exact)

        // We will declare an error vector for passing to the adjoint refinement error estimator
        ErrorVector QoI_elementwise_error;

        // Build an adjoint refinement error estimator object
        std::unique_ptr<AdjointRefinementEstimator> adjoint_refinement_error_estimator =
          build_adjoint_refinement_error_estimator(qois);

        // Estimate the error in each element using the Adjoint Refinement estimator
        adjoint_refinement_error_estimator->estimate_error(system, QoI_elementwise_error);

        // Print out the computed error estimate, note that we access the global error estimates
        // using an accessor function, right now sum(QoI_elementwise_error) != global_QoI_error_estimate
        libMesh::out << "The computed relative error in QoI 0 is " << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) << std::endl; // / std::abs(QoI_0_exact)

        // Also print out effectivity indices (estimated error/true error)
        libMesh::out << "The effectivity index for the computed error in QoI 0 is " << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) /
          std::abs(QoI_0_computed - QoI_0_exact) << std::endl;

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
                libMesh::out << "We reached the target number of elements.\n" << std::endl;
                break;
              }

            mesh_refinement->flag_elements_by_nelem_target (QoI_elementwise_error);

            mesh_refinement->refine_and_coarsen_elements();
          }

        // Dont forget to reinit the system after each adaptive refinement!
        equation_systems.reinit();

        libMesh::out << "Refined mesh to "
                     << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs." << std::endl;
      }

    // On the last adaptive step, dont refine elements and check regressions via asserts
    if (a_step == param.max_adaptivesteps)
      {
        linear_solver->reuse_preconditioner(false);
        system.solve();

        write_output(equation_systems, a_step, "primal");

        NumericVector<Number> & primal_solution = *system.solution;

        QoISet qois;
        std::vector<unsigned int> qoi_indices;

        qoi_indices.push_back(0);
        qois.add_indices(qoi_indices);

        qois.set_weight(0, 1.0);

        system.assemble_qoi_sides = true;

        linear_solver->reuse_preconditioner(param.reuse_preconditioner);
        system.adjoint_solve();
        system.set_adjoint_already_solved(true);

        NumericVector<Number> & dual_solution_0 = system.get_adjoint_solution(0);

        primal_solution.swap(dual_solution_0);
        write_output(equation_systems, a_step, "adjoint_0");

        primal_solution.swap(dual_solution_0);

        libMesh::out << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
                     << " active elements and "
                     << equation_systems.n_active_dofs()
                     << " active dofs." << std::endl ;

        libMesh::out << "Postprocessing: " << std::endl;
        system.postprocess_sides = true;
        system.postprocess();

        Number QoI_0_computed = system.get_QoI_value("computed", 0);
        Number QoI_0_exact = system.get_QoI_value("exact", 0);

        libMesh::out << "The computed QoI 0 is " << std::setprecision(17)
                     << QoI_0_computed << std::endl;
        libMesh::out << "The relative error in QoI 0 is " << std::setprecision(17)
                     << std::abs(QoI_0_computed - QoI_0_exact) << std::endl; // / std::abs(QoI_0_exact)


        // We will declare an error vector for passing to the adjoint
        // refinement error estimator Right now, only the first entry
        // of this vector will be filled (with the global QoI error
        // estimate) Later, each entry of the vector will contain
        // elementwise error that the user can sum to get the total
        // error
        ErrorVector QoI_elementwise_error;

        // Build an adjoint refinement error estimator object
        std::unique_ptr<AdjointRefinementEstimator> adjoint_refinement_error_estimator =
          build_adjoint_refinement_error_estimator(qois);

        // Estimate the error in each element using the Adjoint Refinement estimator
        adjoint_refinement_error_estimator->estimate_error(system, QoI_elementwise_error);

        // Print out the computed error estimate, note that we access
        // the global error estimates using an accessor function,
        // right now sum(QoI_elementwise_error) != global_QoI_error_estimate
        libMesh::out << "The computed relative error in QoI 0 is " << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) << std::endl; // / std::abs(QoI_0_exact)

        // Also print out effectivity indices (estimated error/true error)
        libMesh::out << "The effectivity index for the computed error in QoI 0 is " << std::setprecision(17)
                     << std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) /
          std::abs(QoI_0_computed - QoI_0_exact) << std::endl;

        // Hard coded assert to ensure that the actual numbers we are getting are what they should be

        // The effectivity index isn't exactly reproducible at single precision
        // libmesh_assert_less(std::abs(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact) - 0.84010976704434637), 1.e-5);
        // libmesh_assert_less(std::abs(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact) - 0.48294428289950514), 1.e-5);

        // But the effectivity indices should always be sane
        // libmesh_assert_less(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact), 2.5);
        // libmesh_assert_greater(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(0)) / std::abs(QoI_0_computed - QoI_0_exact), .4);
        // libmesh_assert_less(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact), 2.5);
        // libmesh_assert_greater(std::abs(adjoint_refinement_error_estimator->get_global_QoI_error_estimate(1)) / std::abs(QoI_1_computed - QoI_1_exact), .4);

        // And the computed errors should still be low
        // libmesh_assert_less(std::abs(QoI_0_computed - QoI_0_exact), 2e-4);
        // libmesh_assert_less(std::abs(QoI_1_computed - QoI_1_exact), 2e-4);
      }
  }

  libMesh::err << '[' << mesh.processor_id()
               << "] Completing output." << std::endl;

#endif // #ifndef LIBMESH_ENABLE_AMR

  // All done.
  return 0;
}
