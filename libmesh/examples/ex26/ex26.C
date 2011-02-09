/* $Id$ */

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

 // <h1>Example 26 - Laplace Equation in the L-Shaped Domain with Adjoint based mesh refinement</h1>
 //
 // This example solves the Laplace equation on the classic "L-shaped"
 // domain with adaptive mesh refinement. The exact
 // solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). The kelly and 
 // adjoint residual error estimators are used to develop error indicators and 
 // guide mesh adaptation. Since we use the adjoint capabilities of libMesh in 
 // this example, we use the DiffSystem framework. This file (main.C) contains 
 // the declaration of mesh and equation system objects, L-shaped.C contains the 
 // assembly of the system, element_qoi_derivative.C and side_qoi_derivative.C
 // contain the RHS for the adjoint systems. Postprocessing to compute the 
 // the QoIs is done in element_postprocess.C and side_postprocess.C

 // The initial mesh contains three QUAD9 elements which represent the
 // standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
 // i.e.
 // Element 0: [-1,0]x[ 0,1]
 // Element 1: [ 0,1]x[ 0,1]
 // Element 2: [-1,0]x[-1,0]
 // The mesh is provided in the standard libMesh ASCII format file
 // named "lshaped.xda".  In addition, an input file named "ex26.in"
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
#include "mysystems.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Local function declarations

// Number output files, the files are give a prefix of primal or adjoint_i depending on
// whether the output is the primal solution or the dual solution for the ith QoI

// Write gmv output

void write_output(EquationSystems &es,		                   
		  unsigned int a_step,       // The adaptive step count
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

  system.time_solver =
      AutoPtr<TimeSolver>(new SteadySolver(system));

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
      std::cout<<"Using Adjoint Residual Error Estimator with Patch Recovery Weights"<<std::endl;

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
  FEMSystem &system = build_system(equation_systems, infile, param);

  // Set its parameters
  set_system_parameters(system, param);
  
  std::cout << "Initializing systems" << std::endl;

  equation_systems.init ();

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  // We catch any solver errors, so that we can write the proper
  // footers before closing
  try
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
	NumericVector<Number> &primal_solution =
          dynamic_cast<NumericVector<Number> &>(*system.solution);

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

	// Solve the adjoint system. This takes the transpose of the stiffness matrix and then
	// solves the resulting system
	system.adjoint_solve();

	// Get a pointer to the solution vector of the adjoint problem for QoI 0
	NumericVector<Number> &dual_solution_0 =
          dynamic_cast<NumericVector<Number> &>(system.get_adjoint_solution(0));

	// Swap the primal and dual solutions so we can write out the adjoint solution
	primal_solution.swap(dual_solution_0);	    
	write_output(equation_systems, a_step, "adjoint_0");

	// Swap back
	primal_solution.swap(dual_solution_0);
	
	// Get a pointer to the solution vector of the adjoint problem for QoI 0
	NumericVector<Number> &dual_solution_1 =
          dynamic_cast<NumericVector<Number> &>(system.get_adjoint_solution(1));
	
	// Swap again
	primal_solution.swap(dual_solution_1);	    
	write_output(equation_systems, a_step, "adjoint_1");

	// Swap back again
	primal_solution.swap(dual_solution_1);
			    
	std::cout << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
		  << " active dofs." << std::endl ;

	// Postprocess, compute the approximate QoIs and write them out to the console
	std::cout << "Postprocessing: " << std::endl; 
	system.postprocess_sides = true;
	system.postprocess();
	Number QoI_0_computed = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("computed", 0);
	Number QoI_0_exact = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("exact", 0);
	Number QoI_1_computed = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("computed", 1);
	Number QoI_1_exact = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("exact", 1);
		
	std::cout<< "The relative error in QoI 0 is " << std::setprecision(17)
                 << std::abs(QoI_0_computed - QoI_0_exact) /
                    std::abs(QoI_0_exact) << std::endl;
		
	std::cout<< "The relative error in QoI 1 is " << std::setprecision(17) 
                 << std::abs(QoI_1_computed - QoI_1_exact) /
                    std::abs(QoI_1_exact) << std::endl << std::endl;
	
	// Now we construct the data structures for the mesh refinement process	
	ErrorVector error;
	
	// Build an error estimator object
	AutoPtr<ErrorEstimator> error_estimator =
	  build_error_estimator(param);
	
	// Estimate the error in each element using the Adjoint Residual or Kelly error estimator
	error_estimator->estimate_error(system, error);
			    	    
	// We have to refine either based on reaching an error tolerance or 
	// a number of elements target, which should be verified above
	// Otherwise we flag elements by error tolerance or nelem target	  

	// Uniform refinement
	if(param.refine_uniformly)
	  {	  			    
	    mesh_refinement->uniformly_refine(1);
	  }
	// Adaptively refine based on reaching an error tolerance
	else if(param.global_tolerance >= 0. && param.nelem_target == 0.)
	  {	      	    
	    mesh_refinement->flag_elements_by_error_tolerance (error);              
	    
	    mesh_refinement->refine_and_coarsen_elements();
	  }
	// Adaptively refine based on reaching a target number of elements
	else 
	  {	      	    	    
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

	NumericVector<Number> &primal_solution =
          dynamic_cast<NumericVector<Number> &>(*system.solution);
	
	QoISet qois;
	std::vector<unsigned int> qoi_indices;

	qoi_indices.push_back(0);
	qoi_indices.push_back(1);
	qois.add_indices(qoi_indices);
	
	qois.set_weight(0, 0.5);
	qois.set_weight(1, 0.5);

	system.assemble_qoi_sides = true;
	system.adjoint_solve();
		
	NumericVector<Number> &dual_solution_0 =
          dynamic_cast<NumericVector<Number> &>(system.get_adjoint_solution(0));
	
	primal_solution.swap(dual_solution_0);	    
	write_output(equation_systems, a_step, "adjoint_0");

	primal_solution.swap(dual_solution_0);
	
	NumericVector<Number> &dual_solution_1 =
          dynamic_cast<NumericVector<Number> &>(system.get_adjoint_solution(1));
	
	primal_solution.swap(dual_solution_1);	    
	write_output(equation_systems, a_step, "adjoint_1");

	primal_solution.swap(dual_solution_1);
				
	std::cout << "Adaptive step " << a_step << ", we have " << mesh.n_active_elem()
		  << " active elements and "
		  << equation_systems.n_active_dofs()
		  << " active dofs." << std::endl ;	
	
	std::cout << "Postprocessing: " << std::endl; 
	system.postprocess_sides = true;
	system.postprocess();
	
	Number QoI_0_computed = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("computed", 0);
	Number QoI_0_exact = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("exact", 0);
	Number QoI_1_computed = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("computed", 1);
	Number QoI_1_exact = (dynamic_cast<LaplaceSystem&>(system)).get_QoI_value("exact", 1);
	
	std::cout<< "The relative error in QoI 0 is " << std::setprecision(17)
                 << std::abs(QoI_0_computed - QoI_0_exact) /
                    std::abs(QoI_0_exact) << std::endl;
	
	std::cout<< "The relative error in QoI 1 is " << std::setprecision(17)
                 << std::abs(QoI_1_computed - QoI_1_exact) /
                    std::abs(QoI_1_exact) << std::endl << std::endl;
	
	ErrorVector error;
	
	AutoPtr<ErrorEstimator> error_estimator =
	  build_error_estimator(param);
	
	error_estimator->estimate_error(system, error);
      }
  }
  catch (...)
  {
    std::cerr << '[' << libMesh::processor_id()
              << "] Caught exception; exiting early." << std::endl; 
  }

  std::cerr << '[' << libMesh::processor_id() 
            << "] Completing output." << std::endl; 
    
  // All done.  
  return 0;
}
