<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex27",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">

<br><br>This example solves the Laplace equation on the classic "L-shaped"
domain with adaptive mesh refinement. The exact solution is
u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). We scale this
exact solution by a combination of parameters, (\alpha_{1} + 2
*\alpha_{2}) to get u = (\alpha_{1} + 2 *\alpha_{2}) * r^{2/3} *
\sin ( (2/3) * \theta), which again satisfies the Laplace
Equation. We define the Quantity of Interest in element_qoi.C, and
compute the sensitivity of the QoI to \alpha_{1} and \alpha_{2}
using the adjoint sensitivity method. Since we use the adjoint
capabilities of libMesh in this example, we use the DiffSystem
framework. This file (main.C) contains the declaration of mesh and
equation system objects, L-shaped.C contains the assembly of the
system, element_qoi_derivative.C contains
the RHS for the adjoint system. Postprocessing to compute the the
QoIs is done in element_qoi.C


<br><br>The initial mesh contains three QUAD9 elements which represent the
standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
i.e.
Element 0: [-1,0]x[ 0,1]
Element 1: [ 0,1]x[ 0,1]
Element 2: [-1,0]x[-1,0]
The mesh is provided in the standard libMesh ASCII format file
named "lshaped.xda".  In addition, an input file named "general.in"
is provided which allows the user to set several parameters for
the solution so that the problem can be re-run without a
re-compile.  The solution technique employed is to have a
refinement loop with a linear (forward and adjoint) solve inside followed by a
refinement of the grid and projection of the solution to the new grid
In the final loop iteration, there is no additional
refinement after the solve.  In the input file "general.in", the variable
"max_adaptivesteps" controls the number of refinement steps, and
"refine_fraction" / "coarsen_fraction" determine the number of
elements which will be refined / coarsened at each step.


<br><br>C++ includes
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        
</pre>
</div>
<div class = "comment">
General libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "equation_systems.h"
        #include "error_vector.h"
        #include "mesh.h"
        #include "mesh_refinement.h"
        #include "newton_solver.h"
        #include "numeric_vector.h"
        #include "steady_solver.h"
        #include "system_norm.h"
        
</pre>
</div>
<div class = "comment">
Sensitivity Calculation related includes
</div>

<div class ="fragment">
<pre>
        #include "parameter_vector.h"
        #include "sensitivity_data.h"
        
</pre>
</div>
<div class = "comment">
Error Estimator includes
</div>

<div class ="fragment">
<pre>
        #include "kelly_error_estimator.h"
        #include "patch_recovery_error_estimator.h"
        
</pre>
</div>
<div class = "comment">
Adjoint Related includes
</div>

<div class ="fragment">
<pre>
        #include "adjoint_residual_error_estimator.h"
        #include "qoi_set.h"
        
</pre>
</div>
<div class = "comment">
libMesh I/O includes
</div>

<div class ="fragment">
<pre>
        #include "getpot.h"
        #include "gmv_io.h"
        
</pre>
</div>
<div class = "comment">
Some (older) compilers do not offer full stream
functionality, OStringStream works around this.
</div>

<div class ="fragment">
<pre>
        #include "o_string_stream.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
        #include "femparameters.h"
        #include "mysystems.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Local function declarations


<br><br>Number output files, the files are give a prefix of primal or adjoint_i depending on
whether the output is the primal solution or the dual solution for the ith QoI


<br><br>Write gmv output


<br><br></div>

<div class ="fragment">
<pre>
        void write_output(EquationSystems &es,                                   
        		  unsigned int a_step, // The adaptive step count
        		  std::string solution_type, // primal or adjoint solve
                          FEMParameters &param)
        {
          MeshBase &mesh = es.get_mesh();
        
        #ifdef LIBMESH_HAVE_GMV
          
          OStringStream file_name_gmv;
          file_name_gmv &lt;&lt; solution_type &lt;&lt; ".out.gmv.";      
          OSSRealzeroright(file_name_gmv,2,0,a_step);
          
          GMVIO(mesh).write_equation_systems
            (file_name_gmv.str(), es);
            
        #endif    
        }
        
</pre>
</div>
<div class = "comment">
Set the parameters for the nonlinear and linear solvers to be used during the simulation


<br><br></div>

<div class ="fragment">
<pre>
        void set_system_parameters(FEMSystem &system, FEMParameters &param)
        {
</pre>
</div>
<div class = "comment">
Verify analytic jacobians against numerical ones?
</div>

<div class ="fragment">
<pre>
          system.verify_analytic_jacobians = param.verify_analytic_jacobians;
        
</pre>
</div>
<div class = "comment">
More desperate debugging options
</div>

<div class ="fragment">
<pre>
          system.print_solution_norms = param.print_solution_norms;
          system.print_solutions      = param.print_solutions;
          system.print_residual_norms = param.print_residual_norms;
          system.print_residuals      = param.print_residuals;
          system.print_jacobian_norms = param.print_jacobian_norms;
          system.print_jacobians      = param.print_jacobians;
        
          system.time_solver =
              AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
          {
            NewtonSolver *solver = new NewtonSolver(system);
            system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
            
            solver-&gt;quiet                       = param.solver_quiet;
            solver-&gt;max_nonlinear_iterations    = param.max_nonlinear_iterations;
            solver-&gt;minsteplength               = param.min_step_length;
            solver-&gt;relative_step_tolerance     = param.relative_step_tolerance;
            solver-&gt;relative_residual_tolerance = param.relative_residual_tolerance;
            solver-&gt;require_residual_reduction  = param.require_residual_reduction;
            solver-&gt;linear_tolerance_multiplier = param.linear_tolerance_multiplier;
            if (system.time_solver-&gt;reduce_deltat_on_diffsolver_failure)
              {
        	solver-&gt;continue_after_max_iterations = true;
        	solver-&gt;continue_after_backtrack_failure = true;
              }
            
</pre>
</div>
<div class = "comment">
And the linear solver options
</div>

<div class ="fragment">
<pre>
            solver-&gt;max_linear_iterations       = param.max_linear_iterations;
            solver-&gt;initial_linear_tolerance    = param.initial_linear_tolerance;
            solver-&gt;minimum_linear_tolerance    = param.minimum_linear_tolerance;      
          }
        }
        
</pre>
</div>
<div class = "comment">
Build the mesh refinement object and set parameters for refining/coarsening etc


<br><br></div>

<div class ="fragment">
<pre>
        AutoPtr&lt;MeshRefinement&gt; build_mesh_refinement(MeshBase &mesh,
                                                      FEMParameters &param)
        {
          AutoPtr&lt;MeshRefinement&gt; mesh_refinement(new MeshRefinement(mesh));
          mesh_refinement-&gt;coarsen_by_parents() = true;
          mesh_refinement-&gt;absolute_global_tolerance() = param.global_tolerance;
          mesh_refinement-&gt;nelem_target()      = param.nelem_target;  
          mesh_refinement-&gt;refine_fraction()   = param.refine_fraction;
          mesh_refinement-&gt;coarsen_fraction()  = param.coarsen_fraction;  
          mesh_refinement-&gt;coarsen_threshold() = param.coarsen_threshold;
        
          return mesh_refinement;
        }
        
</pre>
</div>
<div class = "comment">
This is where we declare the error estimators to be built and used for
mesh refinement. The adjoint residual estimator needs two estimators.
One for the forward component of the estimate and one for the adjoint
weighting factor. Here we use the Patch Recovery indicator to estimate both the
forward and adjoint weights. The H1 seminorm component of the error is used
as dictated by the weak form the Laplace equation.


<br><br></div>

<div class ="fragment">
<pre>
        AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters &param)
        {
          AutoPtr&lt;ErrorEstimator&gt; error_estimator;
        
          if (param.indicator_type == "kelly")
            {
              std::cout&lt;&lt;"Using Kelly Error Estimator"&lt;&lt;std::endl;
        
              error_estimator.reset(new KellyErrorEstimator);
            }
          else if (param.indicator_type == "adjoint_residual")
            {
              std::cout&lt;&lt;"Using Adjoint Residual Error Estimator with Patch Recovery Weights"&lt;&lt;std::endl&lt;&lt;std::endl;
        
              AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;
              
              error_estimator.reset (adjoint_residual_estimator);
              
              adjoint_residual_estimator-&gt;adjoint_already_solved = true;
        
              adjoint_residual_estimator-&gt;error_plot_suffix = "error.gmv";
              
              PatchRecoveryErrorEstimator *p1 =
        	new PatchRecoveryErrorEstimator;
              adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
        	   	   	   
              PatchRecoveryErrorEstimator *p2 =
        	new PatchRecoveryErrorEstimator;
              adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);   
        
              adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);
        	   
              adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);	   	   	   	   
            }
          else
            {
              std::cerr &lt;&lt; "Unknown indicator_type" &lt;&lt; std::endl;
              libmesh_error();
            }
          return error_estimator;
        }
        
</pre>
</div>
<div class = "comment">
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
          std::cout &lt;&lt; "Started " &lt;&lt; argv[0] &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Make sure the general input file exists, and parse it
</div>

<div class ="fragment">
<pre>
          {
            std::ifstream i("general.in");
            if (!i)
              {
                std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id() 
                          &lt;&lt; "] Can't find general.in; exiting early." 
                          &lt;&lt; std::endl; 
                libmesh_error();
              }
          }
          GetPot infile("general.in");
          GetPot infile_l_shaped("l-shaped.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
          FEMParameters param;
          param.read(infile);
        
</pre>
</div>
<div class = "comment">
Skip this default-2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
          
</pre>
</div>
<div class = "comment">
Create a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
</pre>
</div>
<div class = "comment">
And an object to refine it
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;MeshRefinement&gt; mesh_refinement =
            build_mesh_refinement(mesh, param);
        
</pre>
</div>
<div class = "comment">
And an EquationSystems to run on it
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
          std::cout &lt;&lt; "Reading in and building the mesh" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Read in the mesh
</div>

<div class ="fragment">
<pre>
          mesh.read(param.domainfile.c_str());
</pre>
</div>
<div class = "comment">
Make all the elements of the mesh second order so we can compute
with a higher order basis
</div>

<div class ="fragment">
<pre>
          mesh.all_second_order();
        
</pre>
</div>
<div class = "comment">
Create a mesh refinement object to do the initial uniform refinements
on the coarse grid read in from lshaped.xda
</div>

<div class ="fragment">
<pre>
          MeshRefinement initial_uniform_refinements(mesh);
          initial_uniform_refinements.uniformly_refine(param.coarserefinements);
            
          std::cout &lt;&lt; "Building system" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Build the FEMSystem
</div>

<div class ="fragment">
<pre>
          FEMSystem &system = build_system(equation_systems, infile, param);
        
</pre>
</div>
<div class = "comment">
Set its parameters
</div>

<div class ="fragment">
<pre>
          set_system_parameters(system, param);
          
          std::cout &lt;&lt; "Initializing systems" &lt;&lt; std::endl;
        
          equation_systems.init ();
        
</pre>
</div>
<div class = "comment">
Print information about the mesh and system to the screen.
</div>

<div class ="fragment">
<pre>
          mesh.print_info();
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
We catch any solver errors, so that we can write the proper
footers before closing
</div>

<div class ="fragment">
<pre>
          try
          {	
</pre>
</div>
<div class = "comment">
Adaptively solve the timestep
</div>

<div class ="fragment">
<pre>
            unsigned int a_step = 0;
            for (; a_step != param.max_adaptivesteps; ++a_step)
              {	    	    
</pre>
</div>
<div class = "comment">
We can't adapt to both a tolerance and a
target mesh size
</div>

<div class ="fragment">
<pre>
                if (param.global_tolerance != 0.)
        	  libmesh_assert (param.nelem_target == 0);
</pre>
</div>
<div class = "comment">
If we aren't adapting to a tolerance we need a
target mesh size
</div>

<div class ="fragment">
<pre>
                    else
                      libmesh_assert (param.nelem_target &gt; 0);
            
</pre>
</div>
<div class = "comment">
Solve the forward problem
</div>

<div class ="fragment">
<pre>
                system.solve();
        
</pre>
</div>
<div class = "comment">
Write out the computed primal solution	
</div>

<div class ="fragment">
<pre>
                write_output(equation_systems, a_step, "primal", param);
        	
</pre>
</div>
<div class = "comment">
Get a pointer to the primal solution vector
</div>

<div class ="fragment">
<pre>
                NumericVector&lt;Number&gt; &primal_solution =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(*system.solution);
        	
</pre>
</div>
<div class = "comment">
Declare a QoISet object, we need this object to set weights for our QoI error contributions
</div>

<div class ="fragment">
<pre>
                QoISet qois;
        
</pre>
</div>
<div class = "comment">
Declare a qoi_indices vector, each index will correspond to a QoI
</div>

<div class ="fragment">
<pre>
                std::vector&lt;unsigned int&gt; qoi_indices;
        	qoi_indices.push_back(0);
        	qoi_indices.push_back(1);
        	qois.add_indices(qoi_indices);
        
</pre>
</div>
<div class = "comment">
Set weights for each index, these will weight the contribution of each QoI in the final error
estimate to be used for flagging elements for refinement
</div>

<div class ="fragment">
<pre>
                qois.set_weight(0, 0.5);
        	qois.set_weight(1, 0.5);
        
</pre>
</div>
<div class = "comment">
A SensitivityData object to hold the qois and parameters
</div>

<div class ="fragment">
<pre>
                SensitivityData sensitivities(qois, system, (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_parameter_vector());
        
</pre>
</div>
<div class = "comment">
Make sure we get the contributions to the adjoint RHS from the sides
</div>

<div class ="fragment">
<pre>
                system.assemble_qoi_sides = true;
        
</pre>
</div>
<div class = "comment">
Compute the sensitivities
</div>

<div class ="fragment">
<pre>
                system.adjoint_qoi_parameter_sensitivity(qois, (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_parameter_vector(), sensitivities);
        
        	GetPot infile("l-shaped.in");
        
        	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
        	Number sensitivity_QoI_0_0_exact = infile("sensitivity_0_0", 0.0);
        	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
        	Number sensitivity_QoI_0_1_exact = infile("sensitivity_0_1", 0.0);
        		
        	std::cout &lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; ", we have " &lt;&lt; mesh.n_active_elem()
                              &lt;&lt; " active elements and "
                              &lt;&lt; equation_systems.n_active_dofs()
        		  &lt;&lt; " active dofs." &lt;&lt; std::endl ;
        
        	std::cout&lt;&lt;"Sensitivity of QoI one to Parameter one is "&lt;&lt;sensitivity_QoI_0_0_computed&lt;&lt;std::endl;
        	std::cout&lt;&lt;"Sensitivity of QoI one to Parameter two is "&lt;&lt;sensitivity_QoI_0_1_computed&lt;&lt;std::endl;
        
        	std::cout&lt;&lt; "The relative error in sensitivity QoI_0_0 is " 
        		 &lt;&lt; std::setprecision(17) 
                         &lt;&lt; std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact) /
                            std::abs(sensitivity_QoI_0_0_exact) &lt;&lt; std::endl;
        	std::cout&lt;&lt; "The relative error in sensitivity QoI_0_1 is " 
                         &lt;&lt; std::setprecision(17) 
                         &lt;&lt; std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact) /
                            std::abs(sensitivity_QoI_0_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;						       			
        
</pre>
</div>
<div class = "comment">
Get a pointer to the solution vector of the adjoint problem for QoI 0
</div>

<div class ="fragment">
<pre>
                NumericVector&lt;Number&gt; &dual_solution_0 =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(0));
        
</pre>
</div>
<div class = "comment">
Swap the primal and dual solutions so we can write out the adjoint solution
</div>

<div class ="fragment">
<pre>
                primal_solution.swap(dual_solution_0);            
        	write_output(equation_systems, a_step, "adjoint_0", param);
        
</pre>
</div>
<div class = "comment">
Swap back
</div>

<div class ="fragment">
<pre>
                primal_solution.swap(dual_solution_0);
        					    										    	    
</pre>
</div>
<div class = "comment">
We have to refine either based on reaching an error tolerance or 
a number of elements target, which should be verified above
Otherwise we flag elements by error tolerance or nelem target	  


<br><br>Uniform refinement
</div>

<div class ="fragment">
<pre>
                if(param.refine_uniformly)
        	  {	  			    
        	    std::cout&lt;&lt;"Refining Uniformly"&lt;&lt;std::endl&lt;&lt;std::endl;
        
        	    mesh_refinement-&gt;uniformly_refine(1);
        	  }
</pre>
</div>
<div class = "comment">
Adaptively refine based on reaching an error tolerance
</div>

<div class ="fragment">
<pre>
                else if(param.global_tolerance &gt;= 0. && param.nelem_target == 0.)
        	  {	      	    
</pre>
</div>
<div class = "comment">
Now we construct the data structures for the mesh refinement process	
</div>

<div class ="fragment">
<pre>
                    ErrorVector error;
        	    
</pre>
</div>
<div class = "comment">
Build an error estimator object
</div>

<div class ="fragment">
<pre>
                    AutoPtr&lt;ErrorEstimator&gt; error_estimator =
        	      build_error_estimator(param);
        
</pre>
</div>
<div class = "comment">
Estimate the error in each element using the Adjoint Residual or Kelly error estimator
</div>

<div class ="fragment">
<pre>
                    error_estimator-&gt;estimate_error(system, error);
        	
        	    mesh_refinement-&gt;flag_elements_by_error_tolerance (error);              
        	    
        	    mesh_refinement-&gt;refine_and_coarsen_elements();
        	  }
</pre>
</div>
<div class = "comment">
Adaptively refine based on reaching a target number of elements
</div>

<div class ="fragment">
<pre>
                else 
        	  {
</pre>
</div>
<div class = "comment">
Now we construct the data structures for the mesh refinement process	
</div>

<div class ="fragment">
<pre>
                    ErrorVector error;
        	    
</pre>
</div>
<div class = "comment">
Build an error estimator object
</div>

<div class ="fragment">
<pre>
                    AutoPtr&lt;ErrorEstimator&gt; error_estimator =
        	      build_error_estimator(param);
        	    
</pre>
</div>
<div class = "comment">
Estimate the error in each element using the Adjoint Residual or Kelly error estimator
</div>

<div class ="fragment">
<pre>
                    error_estimator-&gt;estimate_error(system, error);
        
        	    if (mesh.n_active_elem() &gt;= param.nelem_target)
        	      {
        		std::cout&lt;&lt;"We reached the target number of elements."&lt;&lt;std::endl &lt;&lt;std::endl;
        		break;
        	      }				
        	    
        	    mesh_refinement-&gt;flag_elements_by_nelem_target (error);
        	    	    
        	    mesh_refinement-&gt;refine_and_coarsen_elements();
        	  }
                    
</pre>
</div>
<div class = "comment">
Dont forget to reinit the system after each adaptive refinement !
</div>

<div class ="fragment">
<pre>
                equation_systems.reinit();
              }
        
</pre>
</div>
<div class = "comment">
Do one last solve if necessary
</div>

<div class ="fragment">
<pre>
            if (a_step == param.max_adaptivesteps)
              {	    
        	system.solve();
        	
        	write_output(equation_systems, a_step, "primal", param);	    
        
        	NumericVector&lt;Number&gt; &primal_solution =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(*system.solution);
        				     	
        	QoISet qois;
        	
        	std::vector&lt;unsigned int&gt; qoi_indices;
        	qoi_indices.push_back(0);
        	qoi_indices.push_back(1);
        	qois.add_indices(qoi_indices);
        	
        	qois.set_weight(0, 0.5);
        	qois.set_weight(1, 0.5);
        	
        	SensitivityData sensitivities(qois, system, (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_parameter_vector());
        	
        	system.assemble_qoi_sides = true;
        	
        	system.adjoint_qoi_parameter_sensitivity(qois, (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_parameter_vector(), sensitivities);
        	
        	GetPot infile("l-shaped.in");
        
        	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
        	Number sensitivity_QoI_0_0_exact = infile("sensitivity_0_0", 0.0);
        	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
        	Number sensitivity_QoI_0_1_exact = infile("sensitivity_0_1", 0.0);
        		
        	std::cout &lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; ", we have " &lt;&lt; mesh.n_active_elem()
        		  &lt;&lt; " active elements and "
        		  &lt;&lt; equation_systems.n_active_dofs()
        		  &lt;&lt; " active dofs." &lt;&lt; std::endl ;	
        
        	std::cout&lt;&lt;"Sensitivity of QoI one to Parameter one is "&lt;&lt;sensitivity_QoI_0_0_computed&lt;&lt;std::endl;
        	std::cout&lt;&lt;"Sensitivity of QoI one to Parameter two is "&lt;&lt;sensitivity_QoI_0_1_computed&lt;&lt;std::endl;
        
        	std::cout&lt;&lt; "The error in sensitivity QoI_0_0 is "
                         &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact)/sensitivity_QoI_0_0_exact &lt;&lt; std::endl;
        	std::cout&lt;&lt; "The error in sensitivity QoI_0_1 is "
                         &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact)/sensitivity_QoI_0_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
        		
        	NumericVector&lt;Number&gt; &dual_solution_0 =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(0));
        	
        	primal_solution.swap(dual_solution_0);	    
        	write_output(equation_systems, a_step, "adjoint_0", param);
        
        	primal_solution.swap(dual_solution_0);									
              }
          }
          catch (...)
          {
            std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id()
                      &lt;&lt; "] Caught exception; exiting early." &lt;&lt; std::endl; 
          }
        
          std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id() 
                    &lt;&lt; "] Completing output." &lt;&lt; std::endl; 
            
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return 0;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;system_norm.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;parameter_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sensitivity_data.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;patch_recovery_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;adjoint_residual_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;qoi_set.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mysystems.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output(EquationSystems &amp;es,		                   
  		  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step, <I><FONT COLOR="#B22222">// The adaptive step count
</FONT></I>  		  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type, <I><FONT COLOR="#B22222">// primal or adjoint solve
</FONT></I>                    FEMParameters &amp;param)
  {
    MeshBase &amp;mesh = es.get_mesh();
  
  #ifdef LIBMESH_HAVE_GMV
    
    OStringStream file_name_gmv;
    file_name_gmv &lt;&lt; solution_type &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.out.gmv.&quot;</FONT></B>;      
    OSSRealzeroright(file_name_gmv,2,0,a_step);
    
    GMVIO(mesh).write_equation_systems
      (file_name_gmv.str(), es);
      
  #endif    
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> set_system_parameters(FEMSystem &amp;system, FEMParameters &amp;param)
  {
    system.verify_analytic_jacobians = param.verify_analytic_jacobians;
  
    system.print_solution_norms = param.print_solution_norms;
    system.print_solutions      = param.print_solutions;
    system.print_residual_norms = param.print_residual_norms;
    system.print_residuals      = param.print_residuals;
    system.print_jacobian_norms = param.print_jacobian_norms;
    system.print_jacobians      = param.print_jacobians;
  
    system.time_solver =
        AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
  
    {
      NewtonSolver *solver = <B><FONT COLOR="#A020F0">new</FONT></B> NewtonSolver(system);
      system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
      
      solver-&gt;quiet                       = param.solver_quiet;
      solver-&gt;max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver-&gt;minsteplength               = param.min_step_length;
      solver-&gt;relative_step_tolerance     = param.relative_step_tolerance;
      solver-&gt;relative_residual_tolerance = param.relative_residual_tolerance;
      solver-&gt;require_residual_reduction  = param.require_residual_reduction;
      solver-&gt;linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      <B><FONT COLOR="#A020F0">if</FONT></B> (system.time_solver-&gt;reduce_deltat_on_diffsolver_failure)
        {
  	solver-&gt;continue_after_max_iterations = true;
  	solver-&gt;continue_after_backtrack_failure = true;
        }
      
      solver-&gt;max_linear_iterations       = param.max_linear_iterations;
      solver-&gt;initial_linear_tolerance    = param.initial_linear_tolerance;
      solver-&gt;minimum_linear_tolerance    = param.minimum_linear_tolerance;      
    }
  }
  
  
  AutoPtr&lt;MeshRefinement&gt; build_mesh_refinement(MeshBase &amp;mesh,
                                                FEMParameters &amp;param)
  {
    AutoPtr&lt;MeshRefinement&gt; mesh_refinement(<B><FONT COLOR="#A020F0">new</FONT></B> MeshRefinement(mesh));
    mesh_refinement-&gt;coarsen_by_parents() = true;
    mesh_refinement-&gt;absolute_global_tolerance() = param.global_tolerance;
    mesh_refinement-&gt;nelem_target()      = param.nelem_target;  
    mesh_refinement-&gt;refine_fraction()   = param.refine_fraction;
    mesh_refinement-&gt;coarsen_fraction()  = param.coarsen_fraction;  
    mesh_refinement-&gt;coarsen_threshold() = param.coarsen_threshold;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> mesh_refinement;
  }
  
  
  AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters &amp;param)
  {
    AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.indicator_type == <B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Using Kelly Error Estimator&quot;</FONT></B>&lt;&lt;std::endl;
  
        error_estimator.reset(<B><FONT COLOR="#A020F0">new</FONT></B> KellyErrorEstimator);
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (param.indicator_type == <B><FONT COLOR="#BC8F8F">&quot;adjoint_residual&quot;</FONT></B>)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Using Adjoint Residual Error Estimator with Patch Recovery Weights&quot;</FONT></B>&lt;&lt;std::endl&lt;&lt;std::endl;
  
        AdjointResidualErrorEstimator *adjoint_residual_estimator = <B><FONT COLOR="#A020F0">new</FONT></B> AdjointResidualErrorEstimator;
        
        error_estimator.reset (adjoint_residual_estimator);
        
        adjoint_residual_estimator-&gt;adjoint_already_solved = true;
  
        adjoint_residual_estimator-&gt;error_plot_suffix = <B><FONT COLOR="#BC8F8F">&quot;error.gmv&quot;</FONT></B>;
        
        PatchRecoveryErrorEstimator *p1 =
  	<B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
        adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
  	   	   	   
        PatchRecoveryErrorEstimator *p2 =
  	<B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
        adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);   
  
        adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);
  	   
        adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);	   	   	   	   
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Unknown indicator_type&quot;</FONT></B> &lt;&lt; std::endl;
        libmesh_error();
      }
    <B><FONT COLOR="#A020F0">return</FONT></B> error_estimator;
  }
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Started &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; std::endl;
  
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::ifstream i(<B><FONT COLOR="#BC8F8F">&quot;general.in&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">if</FONT></B> (!i)
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id() 
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Can't find general.in; exiting early.&quot;</FONT></B> 
                    &lt;&lt; std::endl; 
          libmesh_error();
        }
    }
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;general.in&quot;</FONT></B>);
    GetPot infile_l_shaped(<B><FONT COLOR="#BC8F8F">&quot;l-shaped.in&quot;</FONT></B>);
  
    FEMParameters param;
    param.read(infile);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    Mesh mesh;
  
    AutoPtr&lt;MeshRefinement&gt; mesh_refinement =
      build_mesh_refinement(mesh, param);
  
    EquationSystems equation_systems (mesh);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Reading in and building the mesh&quot;</FONT></B> &lt;&lt; std::endl;
  
    mesh.read(param.domainfile.c_str());
    mesh.all_second_order();
  
    MeshRefinement initial_uniform_refinements(mesh);
    initial_uniform_refinements.uniformly_refine(param.coarserefinements);
      
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Building system&quot;</FONT></B> &lt;&lt; std::endl;
  
    FEMSystem &amp;system = build_system(equation_systems, infile, param);
  
    set_system_parameters(system, param);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initializing systems&quot;</FONT></B> &lt;&lt; std::endl;
  
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
  
    <B><FONT COLOR="#A020F0">try</FONT></B>
    {	
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step = 0;
      <B><FONT COLOR="#A020F0">for</FONT></B> (; a_step != param.max_adaptivesteps; ++a_step)
        {	    	    
  	<B><FONT COLOR="#A020F0">if</FONT></B> (param.global_tolerance != 0.)
  	  libmesh_assert (param.nelem_target == 0);
              <B><FONT COLOR="#A020F0">else</FONT></B>
                libmesh_assert (param.nelem_target &gt; 0);
      
  	system.solve();
  
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>, param);
  	
  	NumericVector&lt;Number&gt; &amp;primal_solution =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(*system.solution);
  	
  	QoISet qois;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	SensitivityData sensitivities(qois, system, (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_parameter_vector());
  
  	system.assemble_qoi_sides = true;
  
  	system.adjoint_qoi_parameter_sensitivity(qois, (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_parameter_vector(), sensitivities);
  
  	GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;l-shaped.in&quot;</FONT></B>);
  
  	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
  	Number sensitivity_QoI_0_0_exact = infile(<B><FONT COLOR="#BC8F8F">&quot;sensitivity_0_0&quot;</FONT></B>, 0.0);
  	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
  	Number sensitivity_QoI_0_1_exact = infile(<B><FONT COLOR="#BC8F8F">&quot;sensitivity_0_1&quot;</FONT></B>, 0.0);
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, we have &quot;</FONT></B> &lt;&lt; mesh.n_active_elem()
                        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                        &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl ;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Sensitivity of QoI one to Parameter one is &quot;</FONT></B>&lt;&lt;sensitivity_QoI_0_0_computed&lt;&lt;std::endl;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Sensitivity of QoI one to Parameter two is &quot;</FONT></B>&lt;&lt;sensitivity_QoI_0_1_computed&lt;&lt;std::endl;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in sensitivity QoI_0_0 is &quot;</FONT></B> 
  		 &lt;&lt; std::setprecision(17) 
                   &lt;&lt; std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(sensitivity_QoI_0_0_exact) &lt;&lt; std::endl;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in sensitivity QoI_0_1 is &quot;</FONT></B> 
                   &lt;&lt; std::setprecision(17) 
                   &lt;&lt; std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(sensitivity_QoI_0_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;						       			
  
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>, param);
  
  	primal_solution.swap(dual_solution_0);
  					    										    	    
  
  	<B><FONT COLOR="#A020F0">if</FONT></B>(param.refine_uniformly)
  	  {	  			    
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Refining Uniformly&quot;</FONT></B>&lt;&lt;std::endl&lt;&lt;std::endl;
  
  	    mesh_refinement-&gt;uniformly_refine(1);
  	  }
  	<B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B>(param.global_tolerance &gt;= 0. &amp;&amp; param.nelem_target == 0.)
  	  {	      	    
  	    ErrorVector error;
  	    
  	    AutoPtr&lt;ErrorEstimator&gt; error_estimator =
  	      build_error_estimator(param);
  
  	    error_estimator-&gt;estimate_error(system, error);
  	
  	    mesh_refinement-&gt;flag_elements_by_error_tolerance (error);              
  	    
  	    mesh_refinement-&gt;refine_and_coarsen_elements();
  	  }
  	<B><FONT COLOR="#A020F0">else</FONT></B> 
  	  {
  	    ErrorVector error;
  	    
  	    AutoPtr&lt;ErrorEstimator&gt; error_estimator =
  	      build_error_estimator(param);
  	    
  	    error_estimator-&gt;estimate_error(system, error);
  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (mesh.n_active_elem() &gt;= param.nelem_target)
  	      {
  		<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;We reached the target number of elements.&quot;</FONT></B>&lt;&lt;std::endl &lt;&lt;std::endl;
  		<B><FONT COLOR="#A020F0">break</FONT></B>;
  	      }				
  	    
  	    mesh_refinement-&gt;flag_elements_by_nelem_target (error);
  	    	    
  	    mesh_refinement-&gt;refine_and_coarsen_elements();
  	  }
              
  	equation_systems.reinit();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (a_step == param.max_adaptivesteps)
        {	    
  	system.solve();
  	
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>, param);	    
  
  	NumericVector&lt;Number&gt; &amp;primal_solution =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(*system.solution);
  				     	
  	QoISet qois;
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  	
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  	
  	SensitivityData sensitivities(qois, system, (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_parameter_vector());
  	
  	system.assemble_qoi_sides = true;
  	
  	system.adjoint_qoi_parameter_sensitivity(qois, (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_parameter_vector(), sensitivities);
  	
  	GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;l-shaped.in&quot;</FONT></B>);
  
  	Number sensitivity_QoI_0_0_computed = sensitivities[0][0];
  	Number sensitivity_QoI_0_0_exact = infile(<B><FONT COLOR="#BC8F8F">&quot;sensitivity_0_0&quot;</FONT></B>, 0.0);
  	Number sensitivity_QoI_0_1_computed = sensitivities[0][1];
  	Number sensitivity_QoI_0_1_exact = infile(<B><FONT COLOR="#BC8F8F">&quot;sensitivity_0_1&quot;</FONT></B>, 0.0);
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, we have &quot;</FONT></B> &lt;&lt; mesh.n_active_elem()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
  		  &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl ;	
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Sensitivity of QoI one to Parameter one is &quot;</FONT></B>&lt;&lt;sensitivity_QoI_0_0_computed&lt;&lt;std::endl;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Sensitivity of QoI one to Parameter two is &quot;</FONT></B>&lt;&lt;sensitivity_QoI_0_1_computed&lt;&lt;std::endl;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in sensitivity QoI_0_0 is &quot;</FONT></B>
                   &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(sensitivity_QoI_0_0_computed - sensitivity_QoI_0_0_exact)/sensitivity_QoI_0_0_exact &lt;&lt; std::endl;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in sensitivity QoI_0_1 is &quot;</FONT></B>
                   &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(sensitivity_QoI_0_1_computed - sensitivity_QoI_0_1_exact)/sensitivity_QoI_0_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
  		
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  	
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>, param);
  
  	primal_solution.swap(dual_solution_0);									
        }
    }
    <B><FONT COLOR="#A020F0">catch</FONT></B> (...)
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id()
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Caught exception; exiting early.&quot;</FONT></B> &lt;&lt; std::endl; 
    }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id() 
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Completing output.&quot;</FONT></B> &lt;&lt; std::endl; 
      
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) element_qoi.C...
Compiling C++ (in optimized mode) element_qoi_derivative.C...
Compiling C++ (in optimized mode) ex27.C...
Compiling C++ (in optimized mode) femparameters.C...
Compiling C++ (in optimized mode) L-shaped.C...
Linking ./ex27-opt...
***************************************************************
* Running Example  mpirun -np 2 ./ex27-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Started ./ex27-opt
Reading in and building the mesh
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=833
    n_local_nodes()=433
  n_elem()=255
    n_local_elem()=131
    n_active_elem()=192
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "LaplaceSystem"
    Type "Implicit"
    Variables="T" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="SECOND" 
    n_dofs()=833
    n_local_dofs()=433
    n_constrained_dofs()=0
    n_vectors()=1

Assembling the System
Nonlinear Residual: 1.50294e+10
Linear solve starting, tolerance 1e-12
Linear solve finished, step 14, residual 0.00965966
Trying full Newton step
  Current Residual: 0.00965963
  Nonlinear solver converged, step 0, residual reduction 6.42714e-13 < 1e-09
  Nonlinear solver relative step size inf > 1e-08
Adaptive step 0, we have 192 active elements and 833 active dofs.
Sensitivity of QoI one to Parameter one is 0.00543661
Sensitivity of QoI one to Parameter two is 0.0108732
The relative error in sensitivity QoI_0_0 is 9.9048031112916826e-05
The relative error in sensitivity QoI_0_1 is 9.9048172719716341e-05

Refining Uniformly

Assembling the System
Nonlinear Residual: 23210.667571448994
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 30, residual 1.1397939567548371e-08
Trying full Newton step
  Current Residual: 6.3851740204349951e-07
  Nonlinear solver converged, step 0, residual reduction 2.7509652623215706e-11 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.0006540817657206202 > 1e-08
Adaptive step 1, we have 768 active elements and 3201 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054368750920342828
Sensitivity of QoI one to Parameter two is 0.010873750183332958
The relative error in sensitivity QoI_0_0 is 5.0113869550874657e-05
The relative error in sensitivity QoI_0_1 is 5.0113937197032565e-05

Refining Uniformly

Assembling the System
Nonlinear Residual: 2112.7541760032022
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 59, residual 1.4150014043371439e-09
Trying full Newton step
  Current Residual: 4.4564694130789011e-07
  Nonlinear solver converged, step 0, residual reduction 2.1093175267126518e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.00021159976135482435 > 1e-08
Adaptive step 2, we have 3072 active elements and 12545 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054370392126232249
Sensitivity of QoI one to Parameter two is 0.010874078423072843
The relative error in sensitivity QoI_0_0 is 1.9928816290282997e-05
The relative error in sensitivity QoI_0_1 is 1.9929016174807365e-05

Refining Uniformly

Assembling the System
Nonlinear Residual: 188.85090644899773
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 99, residual 1.6750664498746172e-10
Trying full Newton step
  Current Residual: 2.8947774522029521e-07
  Nonlinear solver reached maximum step 0 with norm 402.47626074292833
  Continuing...
Adaptive step 3, we have 12288 active elements and 49665 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054371045782562856
Sensitivity of QoI one to Parameter two is 0.010874209154935239
The error in sensitivity QoI_0_0 is 7.9067712313990531e-06
The error in sensitivity QoI_0_1 is 7.9069162825058786e-06

[1] Completing output.
[0] Completing output.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex27-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Tue Feb 22 12:24:31 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.975e+00      1.03521   3.907e+00
Objects:              6.230e+02      1.00000   6.230e+02
Flops:                1.857e+09      1.01405   1.844e+09  3.688e+09
Flops/sec:            4.769e+08      1.02086   4.720e+08  9.440e+08
MPI Messages:         8.370e+02      1.00000   8.370e+02  1.674e+03
MPI Message Lengths:  5.854e+06      1.00848   6.964e+03  1.166e+07
MPI Reductions:       1.763e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.9072e+00 100.0%  3.6878e+09 100.0%  1.674e+03 100.0%  6.964e+03      100.0%  1.628e+03  92.3% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

KSPGMRESOrthog       362 1.0 1.2364e-01 1.1 2.82e+08 1.0 0.0e+00 0.0e+00 3.6e+02  3 15  0  0 21   3 15  0  0 22  4540
KSPSetup              16 1.0 4.5204e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               8 1.0 2.0021e+00 1.0 1.86e+09 1.0 7.5e+02 1.6e+03 7.7e+02 51100 45 11 44  51100 45 11 47  1841
PCSetUp               16 1.0 9.3896e-01 1.0 3.56e+08 1.0 0.0e+00 0.0e+00 2.4e+01 24 19  0  0  1  24 19  0  0  1   750
PCSetUpOnBlocks        8 1.0 9.3855e-01 1.0 3.56e+08 1.0 0.0e+00 0.0e+00 2.4e+01 24 19  0  0  1  24 19  0  0  1   750
PCApply              377 1.0 7.5246e-01 1.0 1.04e+09 1.0 0.0e+00 0.0e+00 0.0e+00 19 56  0  0  0  19 56  0  0  0  2737
VecDot                 8 1.0 1.3828e-04 1.0 1.33e+05 1.0 0.0e+00 0.0e+00 8.0e+00  0  0  0  0  0   0  0  0  0  0  1916
VecMDot              362 1.0 8.4022e-02 1.2 1.41e+08 1.0 0.0e+00 0.0e+00 3.6e+02  2  8  0  0 21   2  8  0  0 22  3341
VecNorm              401 1.0 4.0324e-02 1.0 1.05e+07 1.0 0.0e+00 0.0e+00 4.0e+02  1  1  0  0 23   1  1  0  0 25   520
VecScale             393 1.0 1.9016e-03 1.0 5.20e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  5442
VecCopy               59 1.0 8.3733e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               464 1.0 3.9530e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               42 1.0 5.8842e-04 1.1 1.02e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3435
VecMAXPY             377 1.0 4.2344e-02 1.0 1.51e+08 1.0 0.0e+00 0.0e+00 0.0e+00  1  8  0  0  0   1  8  0  0  0  7086
VecAssemblyBegin     151 1.0 2.0834e-02 1.8 0.00e+00 0.0 7.6e+01 8.6e+02 4.2e+02  0  0  5  1 24   0  0  5  1 26     0
VecAssemblyEnd       151 1.0 1.3804e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      523 1.0 3.7556e-03 1.3 0.00e+00 0.0 9.1e+02 1.9e+03 0.0e+00  0  0 54 15  0   0  0 54 15  0     0
VecScatterEnd        523 1.0 2.0218e-0213.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         377 1.0 3.8020e-02 1.0 1.52e+07 1.0 0.0e+00 0.0e+00 3.8e+02  1  1  0  0 21   1  1  0  0 23   796
MatMult              377 1.0 1.6577e-01 1.1 1.56e+08 1.0 7.5e+02 1.6e+03 0.0e+00  4  8 45 11  0   4  8 45 11  0  1867
MatSolve             377 1.0 7.4458e-01 1.0 1.04e+09 1.0 0.0e+00 0.0e+00 0.0e+00 19 56  0  0  0  19 56  0  0  0  2766
MatLUFactorNum         8 1.0 3.0688e-01 1.0 3.56e+08 1.0 0.0e+00 0.0e+00 0.0e+00  8 19  0  0  0   8 19  0  0  0  2295
MatILUFactorSym        8 1.0 6.3111e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 8.0e+00 16  0  0  0  0  16  0  0  0  0     0
MatAssemblyBegin      32 1.0 5.9409e-03 2.2 0.00e+00 0.0 6.0e+01 5.8e+03 6.4e+01  0  0  4  3  4   0  0  4  3  4     0
MatAssemblyEnd        32 1.0 2.5513e-02 1.0 0.00e+00 0.0 3.2e+01 3.0e+02 8.0e+01  1  0  2  0  5   1  0  2  0  5     0
MatGetRowIJ            8 1.0 3.8147e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         8 1.0 3.4499e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.6e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        20 1.0 3.2954e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatTranspose           4 1.0 4.0432e-02 1.0 0.00e+00 0.0 4.0e+01 2.4e+03 3.6e+01  1  0  2  1  2   1  0  2  1  2     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver    16             16       151040     0
      Preconditioner    16             16        11264     0
                 Vec   343            343     22369480     0
         Vec Scatter    88             88        76384     0
           Index Set   124            124       873200     0
   IS L to G Mapping     4              4       138008     0
              Matrix    32             32    102614112     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.7643e-06
Average time for zero size MPI_Send(): 6.91414e-06
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-pc_type bjacobi
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------
 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex27-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
</pre>
</div>
<?php make_footer() ?>
</body>
</html>
<?php if (0) { ?>
\#Local Variables:
\#mode: html
\#End:
<?php } ?>
