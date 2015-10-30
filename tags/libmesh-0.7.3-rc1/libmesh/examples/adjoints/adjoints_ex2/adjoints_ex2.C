/* The Next Great Finite Element Library. */
/* Copyright (C) 2004  Benjamin S. Kirk, John W. Peterson */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  U$ */

// <h1>Example 27 - Laplace Equation in the L-Shaped Domain with
// Adjoint based sensitivity</h1>
 //
 // This example solves the Laplace equation on the classic "L-shaped"
 // domain with adaptive mesh refinement. The exact solution is
 // u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). We scale this
 // exact solution by a combination of parameters, (\alpha_{1} + 2
 // *\alpha_{2}) to get u = (\alpha_{1} + 2 *\alpha_{2}) * r^{2/3} *
 // \sin ( (2/3) * \theta), which again satisfies the Laplace
 // Equation. We define the Quantity of Interest in element_qoi.C, and
 // compute the sensitivity of the QoI to \alpha_{1} and \alpha_{2}
 // using the adjoint sensitivity method. Since we use the adjoint
 // capabilities of libMesh in this example, we use the DiffSystem
 // framework. This file (main.C) contains the declaration of mesh and
 // equation system objects, L-shaped.C contains the assembly of the
 // system, element_qoi_derivative.C contains
 // the RHS for the adjoint system. Postprocessing to compute the the
 // QoIs is done in element_qoi.C

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

// General libMesh includes
#include "equation_systems.h"
#include "error_vector.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "newton_solver.h"
#include "numeric_vector.h"
#include "steady_solver.h"
#include "system_norm.h"

// Sensitivity Calculation related includes
#include "parameter_vector.h"
#include "sensitivity_data.h"

// Error Estimator includes
#include "kelly_error_estimator.h"
#include "patch_recovery_error_estimator.h"

// Adjoint Related includes
#include "adjoint_residual_error_estimator.h"
#include "qoi_set.h"

// libMesh I/O includes
#include "getpot.h"
#include "gmv_io.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// Local includes
#include "femparameters.h"
#include "L-shaped.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Local function declarations

// Number output files, the files are give a prefix of primal or adjoint_i depending on
// whether the output is the primal solution or the dual solution for the ith QoI

// Write gmv output

void write_output(EquationSystems &es,		                   
		  unsigned int a_step, // The adaptive step count
		  std::string solution_type) // primal or adjoint solve
{
  MeshBase &mesh = es.get_mesh();

#ifdef LIBMESH_HAVE_GMV
  
  OStringStream file_name_gmv;
  file_name_gmv << solution_type << ".out.gmv.";      
  OSSRealzeroright(file_name_gmv,2,0,a_step);
  
  GMVIO(mesh).write_equation_systems
    (file_name_gmv.str(), es);
    
#endif    
}

// Set the parameters for the nonlinear and linear solvers to be used during the simulation

void set_system_parameters(LaplaceSystem &system, FEMParameters &param)
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
  system.time_solver =
      AutoPtr<TimeSolver>(new SteadySolver(system));

  // Nonlinear solver options
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

// Build the mesh refinement object and set parameters for refining/coarsening etc

#ifdef LIBMESH_ENABLE_AMR

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

#endif // LIBMESH_ENABLE_AMR

// This is where we declare the error estimators to be built and used for
// mesh refinement. The adjoint residual estimator needs two estimators.
// One for the forward component of the estimate and one for the adjoint
// weighting factor. Here we use the Patch Recovery indicator to estimate both the
// forward and adjoint weights. The H1 seminorm component of the error is used
// as dictated by the weak form the Laplace equation.

AutoPtr<ErrorEstimator> build_error_estimator(FEMParameters &param)
{
  AutoPtr<ErrorEstimator> error_estimator;

  if (param.indicator_type == "kelly")
    {
      std::cout<<"Using Kelly Error Estimator"<<std::endl;

      error_estimator.reset(new KellyErrorEstimator);
    }
  else if (param.indicator_type == "adjoint_residual")
    {
      std::cout<<"Using Adjoint Residual Error Estimator with Patch Recovery Weights"<<std::endl<<std::endl;

      AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;
      
      error_estimator.reset (adjoint_residual_estimator);
      
      adjoint_residual_estimator->adjoint_already_solved = true;

      adjoint_residual_estimator->error_plot_suffix = "error.gmv";
      
      PatchRecoveryErrorEstimator *p1 =
	new PatchRecoveryErrorEstimator;
      adjoint_residual_estimator->primal_error_estimator().reset(p1);
	   	   	   
      PatchRecoveryErrorEstimator *p2 =
	new PatchRecoveryErrorEstimator;
      adjoint_residual_estimator->dual_error_estimator().reset(p2);   

      adjoint_residual_estimator->primal_error_estimator()->error_norm.set_type(0, H1_SEMINORM);
	   
      adjoint_residual_estimator->dual_error_estimator()->error_norm.set_type(0, H1_SEMINORM);	   	   	   	   
    }
  else
    {
      std::cerr << "Unknown indicator_type" << std::endl;
      libmesh_error();
    }
  return error_estimator;
}

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip adaptive examples on a non-adaptive libMesh build
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_assert(false, "--enable-amr");
#else
  // Only our PETSc interface currently supports adjoint solves
  libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");


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
  GetPot infile_l_shaped("l-shaped.in");

  // Read in parameters from the input file
  FEMParameters param;
  param.read(infile);

  // Skip this default-2D example if libMesh was compiled as 1D-only.
  libmesh_example_assert(2 <= LIBMESH_DIM, "2D support");
  
  // Create a mesh.
  Mesh mesh;

  // And an object to refine it
  AutoPtr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  std::cout << "Reading in and building the mesh" << std::endl;

   // Read in the mesh
  mesh.read(param.domainfile.c_str());
  // Make all the elements of the mesh second order so we can compute
  // with a higher order basis
  mesh.all_second_order();

  // Create a mesh refinement object to do the initial uniform refinements
  // on the coarse grid read in from lshaped.xda
  MeshRefinement initial_uniform_refinements(mesh);
  initial_uniform_refinements.uniformly_refine(param.coarserefinements);
    
  std::cout << "Building system" << std::endl;

  // Build the FEMSystem
  LaplaceSystem &system = equation_systems.add_system<LaplaceSystem> ("LaplaceSystem");

  // Set its parameters
  set_system_parameters(system, param);
  
  std::cout << "Initializing systems" << std::endl;

  equation_systems.init ();

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  {	
    // Adaptively solve the timestep
    unsigned int a_step = 0;
    for (; a_step != param.max_adaptivesteps; ++a_step)
      {	    	    
	// We can't adapt to both a tolerance and a
	// target mesh size
	if (param.global_tolerance != 0.)
	  libmesh_assert (param.nelem_target == 0);
	// If we aren't adapting to a tolerance we need a
	// target mesh size
            else
              libmesh_assert (param.nelem_target > 0);
    
	// Solve the forward problem
	system.solve();

	// Write out the computed primal solution	
	write_output(equation_systems, a_step, "primal");
	
	// Get a pointer to the primal solution vector
	NumericVector<Number> &primal_solution = *system.solution;
	
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

	// A SensitivityData object to hold the qois and parameters
	SensitivityData sensitivities(qois, system, system.get_parameter_vector());

	// Make sure we get the contributions to the adjoint RHS from the sides
	system.assemble_qoi_sides = true;

	// Compute the sensitivities
	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);

	GetPot infile("l-shaped.in");

	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
	Number sensitivity_QoI_0_0_exact = infile("sensitivity_0_0", 0.0);
	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
	Number sensitivity_QoI_0_1_exact = infile("sensitivity_0_1", 0.0);
		
	std::cout << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
		  << " active dofs." << std::endl ;

	std::cout<<"Sensitivity of QoI one to Parameter one is "<<sensitivity_QoI_0_0_computed<<std::endl;
	std::cout<<"Sensitivity of QoI one to Parameter two is "<<sensitivity_QoI_0_1_computed<<std::endl;

	std::cout<< "The relative error in sensitivity QoI_0_0 is " 
		 << std::setprecision(17) 
                 << std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact) /
                    std::abs(sensitivity_QoI_0_0_exact) << std::endl;
	std::cout<< "The relative error in sensitivity QoI_0_1 is " 
                 << std::setprecision(17) 
                 << std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact) /
                    std::abs(sensitivity_QoI_0_1_exact) << std::endl << std::endl;						       			

	// Get a pointer to the solution vector of the adjoint problem for QoI 0
	NumericVector<Number> &dual_solution_0 = system.get_adjoint_solution(0);

	// Swap the primal and dual solutions so we can write out the adjoint solution
	primal_solution.swap(dual_solution_0);	    
	write_output(equation_systems, a_step, "adjoint_0");

	// Swap back
	primal_solution.swap(dual_solution_0);
					    										    	    
	// We have to refine either based on reaching an error tolerance or 
	// a number of elements target, which should be verified above
	// Otherwise we flag elements by error tolerance or nelem target	  

	// Uniform refinement
	if(param.refine_uniformly)
	  {	  			    
	    std::cout<<"Refining Uniformly"<<std::endl<<std::endl;

	    mesh_refinement->uniformly_refine(1);
	  }
	// Adaptively refine based on reaching an error tolerance
	else if(param.global_tolerance >= 0. && param.nelem_target == 0.)
	  {	      	    
	    // Now we construct the data structures for the mesh refinement process	
	    ErrorVector error;
	    
	    // Build an error estimator object
	    AutoPtr<ErrorEstimator> error_estimator =
	      build_error_estimator(param);

	    // Estimate the error in each element using the Adjoint Residual or Kelly error estimator
	    error_estimator->estimate_error(system, error);
	
	    mesh_refinement->flag_elements_by_error_tolerance (error);              
	    
	    mesh_refinement->refine_and_coarsen_elements();
	  }
	// Adaptively refine based on reaching a target number of elements
	else 
	  {
	    // Now we construct the data structures for the mesh refinement process	
	    ErrorVector error;
	    
	    // Build an error estimator object
	    AutoPtr<ErrorEstimator> error_estimator =
	      build_error_estimator(param);
	    
	    // Estimate the error in each element using the Adjoint Residual or Kelly error estimator
	    error_estimator->estimate_error(system, error);

	    if (mesh.n_active_elem() >= param.nelem_target)
	      {
		std::cout<<"We reached the target number of elements."<<std::endl <<std::endl;
		break;
	      }				
	    
	    mesh_refinement->flag_elements_by_nelem_target (error);
	    	    
	    mesh_refinement->refine_and_coarsen_elements();
	  }
            
	// Dont forget to reinit the system after each adaptive refinement !
	equation_systems.reinit();
      }

    // Do one last solve if necessary
    if (a_step == param.max_adaptivesteps)
      {	    
	system.solve();
	
	write_output(equation_systems, a_step, "primal");

	NumericVector<Number> &primal_solution = *system.solution;
				     	
	QoISet qois;
	
	std::vector<unsigned int> qoi_indices;
	qoi_indices.push_back(0);
	qoi_indices.push_back(1);
	qois.add_indices(qoi_indices);
	
	qois.set_weight(0, 0.5);
	qois.set_weight(1, 0.5);
	
	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
	
	system.assemble_qoi_sides = true;
	
	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
	
	GetPot infile("l-shaped.in");

	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
	Number sensitivity_QoI_0_0_exact = infile("sensitivity_0_0", 0.0);
	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
	Number sensitivity_QoI_0_1_exact = infile("sensitivity_0_1", 0.0);
		
	std::cout << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
		  << " active elements and "
		  << equation_systems.n_active_dofs()
		  << " active dofs." << std::endl ;	

	std::cout<<"Sensitivity of QoI one to Parameter one is "<<sensitivity_QoI_0_0_computed<<std::endl;
	std::cout<<"Sensitivity of QoI one to Parameter two is "<<sensitivity_QoI_0_1_computed<<std::endl;

	std::cout<< "The error in sensitivity QoI_0_0 is "
                 << std::setprecision(17)
                 << std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact)/sensitivity_QoI_0_0_exact << std::endl;
	std::cout<< "The error in sensitivity QoI_0_1 is "
                 << std::setprecision(17)
                 << std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact)/sensitivity_QoI_0_1_exact << std::endl << std::endl;
		
	NumericVector<Number> &dual_solution_0 = system.get_adjoint_solution(0);
	
	primal_solution.swap(dual_solution_0);	    
	write_output(equation_systems, a_step, "adjoint_0");

	primal_solution.swap(dual_solution_0);									
      }
  }

  std::cerr << '[' << libMesh::processor_id() 
            << "] Completing output." << std::endl; 
    
#endif // #ifndef LIBMESH_ENABLE_AMR
    
  // All done.  
  return 0;
}
