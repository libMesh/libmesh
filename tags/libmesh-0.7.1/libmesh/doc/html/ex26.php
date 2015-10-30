<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex26",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 26 - Laplace Equation in the L-Shaped Domain with Adjoint based mesh refinement</h1>

<br><br>This example solves the Laplace equation on the classic "L-shaped"
domain with adaptive mesh refinement. The exact
solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). The kelly and 
adjoint residual error estimators are used to develop error indicators and 
guide mesh adaptation. Since we use the adjoint capabilities of libMesh in 
this example, we use the DiffSystem framework. This file (main.C) contains 
the declaration of mesh and equation system objects, L-shaped.C contains the 
assembly of the system, element_qoi_derivative.C and side_qoi_derivative.C
contain the RHS for the adjoint systems. Postprocessing to compute the 
the QoIs is done in element_postprocess.C and side_postprocess.C


<br><br>The initial mesh contains three QUAD9 elements which represent the
standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
i.e.
Element 0: [-1,0]x[ 0,1]
Element 1: [ 0,1]x[ 0,1]
Element 2: [-1,0]x[-1,0]
The mesh is provided in the standard libMesh ASCII format file
named "lshaped.xda".  In addition, an input file named "ex26.in"
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
        		  unsigned int a_step,       // The adaptive step count
        		  std::string solution_type) // primal or adjoint solve
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
              std::cout&lt;&lt;"Using Adjoint Residual Error Estimator with Patch Recovery Weights"&lt;&lt;std::endl;
        
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
                write_output(equation_systems, a_step, "primal");
        	
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
Make sure we get the contributions to the adjoint RHS from the sides
</div>

<div class ="fragment">
<pre>
                system.assemble_qoi_sides = true;
        
</pre>
</div>
<div class = "comment">
Solve the adjoint system. This takes the transpose of the stiffness matrix and then
solves the resulting system
</div>

<div class ="fragment">
<pre>
                system.adjoint_solve();
        
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
        	write_output(equation_systems, a_step, "adjoint_0");
        
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
Get a pointer to the solution vector of the adjoint problem for QoI 0
</div>

<div class ="fragment">
<pre>
                NumericVector&lt;Number&gt; &dual_solution_1 =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(1));
        	
</pre>
</div>
<div class = "comment">
Swap again
</div>

<div class ="fragment">
<pre>
                primal_solution.swap(dual_solution_1);            
        	write_output(equation_systems, a_step, "adjoint_1");
        
</pre>
</div>
<div class = "comment">
Swap back again
</div>

<div class ="fragment">
<pre>
                primal_solution.swap(dual_solution_1);
        			    
        	std::cout &lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; ", we have " &lt;&lt; mesh.n_active_elem()
                              &lt;&lt; " active elements and "
                              &lt;&lt; equation_systems.n_active_dofs()
        		  &lt;&lt; " active dofs." &lt;&lt; std::endl ;
        
</pre>
</div>
<div class = "comment">
Postprocess, compute the approximate QoIs and write them out to the console
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; "Postprocessing: " &lt;&lt; std::endl; 
        	system.postprocess_sides = true;
        	system.postprocess();
        	Number QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 0);
        	Number QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 0);
        	Number QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 1);
        	Number QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 1);
        		
        	std::cout&lt;&lt; "The relative error in QoI 0 is " &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                            std::abs(QoI_0_exact) &lt;&lt; std::endl;
        		
        	std::cout&lt;&lt; "The relative error in QoI 1 is " &lt;&lt; std::setprecision(17) 
                         &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                            std::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
        	
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
        	
        	write_output(equation_systems, a_step, "primal");
        
        	NumericVector&lt;Number&gt; &primal_solution =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(*system.solution);
        	
        	QoISet qois;
        	std::vector&lt;unsigned int&gt; qoi_indices;
        
        	qoi_indices.push_back(0);
        	qoi_indices.push_back(1);
        	qois.add_indices(qoi_indices);
        	
        	qois.set_weight(0, 0.5);
        	qois.set_weight(1, 0.5);
        
        	system.assemble_qoi_sides = true;
        	system.adjoint_solve();
        		
        	NumericVector&lt;Number&gt; &dual_solution_0 =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(0));
        	
        	primal_solution.swap(dual_solution_0);	    
        	write_output(equation_systems, a_step, "adjoint_0");
        
        	primal_solution.swap(dual_solution_0);
        	
        	NumericVector&lt;Number&gt; &dual_solution_1 =
                  dynamic_cast&lt;NumericVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(1));
        	
        	primal_solution.swap(dual_solution_1);	    
        	write_output(equation_systems, a_step, "adjoint_1");
        
        	primal_solution.swap(dual_solution_1);
        				
        	std::cout &lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; ", we have " &lt;&lt; mesh.n_active_elem()
        		  &lt;&lt; " active elements and "
        		  &lt;&lt; equation_systems.n_active_dofs()
        		  &lt;&lt; " active dofs." &lt;&lt; std::endl ;	
        	
        	std::cout &lt;&lt; "Postprocessing: " &lt;&lt; std::endl; 
        	system.postprocess_sides = true;
        	system.postprocess();
        	
        	Number QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 0);
        	Number QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 0);
        	Number QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 1);
        	Number QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 1);
        	
        	std::cout&lt;&lt; "The relative error in QoI 0 is " &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                            std::abs(QoI_0_exact) &lt;&lt; std::endl;
        	
        	std::cout&lt;&lt; "The relative error in QoI 1 is " &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                            std::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
        	
        	ErrorVector error;
        	
        	AutoPtr&lt;ErrorEstimator&gt; error_estimator =
        	  build_error_estimator(param);
        	
        	error_estimator-&gt;estimate_error(system, error);
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
  		  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step,       <I><FONT COLOR="#B22222">// The adaptive step count
</FONT></I>  		  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type) <I><FONT COLOR="#B22222">// primal or adjoint solve
</FONT></I>  {
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
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Using Adjoint Residual Error Estimator with Patch Recovery Weights&quot;</FONT></B>&lt;&lt;std::endl;
  
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
  
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
  	
  	NumericVector&lt;Number&gt; &amp;primal_solution =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(*system.solution);
  
  	QoISet qois;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	system.assemble_qoi_sides = true;
  
  	system.adjoint_solve();
  
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	NumericVector&lt;Number&gt; &amp;dual_solution_1 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(1));
  	
  	primal_solution.swap(dual_solution_1);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_1&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_1);
  			    
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, we have &quot;</FONT></B> &lt;&lt; mesh.n_active_elem()
                        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                        &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl ;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Postprocessing: &quot;</FONT></B> &lt;&lt; std::endl; 
  	system.postprocess_sides = true;
  	system.postprocess();
  	Number QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Number QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Number QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Number QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_0_exact) &lt;&lt; std::endl;
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) 
                   &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
  	
  	ErrorVector error;
  	
  	AutoPtr&lt;ErrorEstimator&gt; error_estimator =
  	  build_error_estimator(param);
  	
  	error_estimator-&gt;estimate_error(system, error);
  			    	    
  
  	<B><FONT COLOR="#A020F0">if</FONT></B>(param.refine_uniformly)
  	  {	  			    
  	    mesh_refinement-&gt;uniformly_refine(1);
  	  }
  	<B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B>(param.global_tolerance &gt;= 0. &amp;&amp; param.nelem_target == 0.)
  	  {	      	    
  	    mesh_refinement-&gt;flag_elements_by_error_tolerance (error);              
  	    
  	    mesh_refinement-&gt;refine_and_coarsen_elements();
  	  }
  	<B><FONT COLOR="#A020F0">else</FONT></B> 
  	  {	      	    	    
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
  	
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
  
  	NumericVector&lt;Number&gt; &amp;primal_solution =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(*system.solution);
  	
  	QoISet qois;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  	
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	system.assemble_qoi_sides = true;
  	system.adjoint_solve();
  		
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  	
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	NumericVector&lt;Number&gt; &amp;dual_solution_1 =
            dynamic_cast&lt;NumericVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(1));
  	
  	primal_solution.swap(dual_solution_1);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_1&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_1);
  				
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, we have &quot;</FONT></B> &lt;&lt; mesh.n_active_elem()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
  		  &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl ;	
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Postprocessing: &quot;</FONT></B> &lt;&lt; std::endl; 
  	system.postprocess_sides = true;
  	system.postprocess();
  	
  	Number QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Number QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Number QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Number QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_0_exact) &lt;&lt; std::endl;
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
  	
  	ErrorVector error;
  	
  	AutoPtr&lt;ErrorEstimator&gt; error_estimator =
  	  build_error_estimator(param);
  	
  	error_estimator-&gt;estimate_error(system, error);
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
Compiling C++ (in optimized mode) element_postprocess.C...
Compiling C++ (in optimized mode) element_qoi_derivative.C...
Compiling C++ (in optimized mode) ex26.C...
Compiling C++ (in optimized mode) femparameters.C...
Compiling C++ (in optimized mode) L-shaped.C...
Compiling C++ (in optimized mode) side_postprocess.C...
Compiling C++ (in optimized mode) side_qoi_derivative.C...
Linking ./ex26-opt...
***************************************************************
* Running Example  mpirun -np 2 ./ex26-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Started ./ex26-opt
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
Nonlinear Residual: 5.00981e+09
Linear solve starting, tolerance 1e-12
Linear solve finished, step 13, residual 0.00220256
Trying full Newton step
  Current Residual: 0.00220255
  Nonlinear solver converged, step 0, residual reduction 4.39647e-13 < 1e-09
  Nonlinear solver relative step size inf > 1e-08
Adaptive step 0, we have 192 active elements and 833 active dofs.
Postprocessing: 
The relative error in QoI 0 is 0.00016365688099197935
The relative error in QoI 1 is 0.00032228775655961795

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 1893.5349258185458
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 20, residual 1.7729685303990687e-09
Trying full Newton step
  Current Residual: 3.2900050245111023e-07
  Nonlinear solver converged, step 0, residual reduction 1.7374937106527802e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.0011166144801248779 > 1e-08
Adaptive step 1, we have 222 active elements and 943 active dofs.
Postprocessing: 
The relative error in QoI 0 is 4.9719364089729268e-05
The relative error in QoI 1 is 7.8937527867902526e-05

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 1204.6689776172188
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 22, residual 9.1003101883228154e-10
Trying full Newton step
  Current Residual: 2.8710620380044727e-07
  Nonlinear solver converged, step 0, residual reduction 2.3832788021845676e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.00056434343309561452 > 1e-08
Adaptive step 2, we have 267 active elements and 1117 active dofs.
Postprocessing: 
The relative error in QoI 0 is 1.9351576491789132e-05
The relative error in QoI 1 is 2.6065414311629888e-07

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 957.76625412544399
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 22, residual 3.1763266653953041e-10
Trying full Newton step
  Current Residual: 3.0232149069470178e-07
  Nonlinear solver converged, step 0, residual reduction 3.1565268602072196e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.00034212214213161246 > 1e-08
Adaptive step 3, we have 312 active elements and 1287 active dofs.
Postprocessing: 
The relative error in QoI 0 is 7.1570030459360749e-06
The relative error in QoI 1 is 3.8377315458006266e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 1451.8907743524978
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 17, residual 1.0475874281079329e-09
Trying full Newton step
  Current Residual: 3.3350906241485033e-07
  Nonlinear solver converged, step 0, residual reduction 2.2970671644606732e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 0.00019305877998599563 > 1e-08
Adaptive step 4, we have 363 active elements and 1483 active dofs.
Postprocessing: 
The relative error in QoI 0 is 2.4148442282445063e-06
The relative error in QoI 1 is 9.8808967341529639e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 87.196932224487057
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 19, residual 7.0635037093672258e-11
Trying full Newton step
  Current Residual: 2.6136485000673542e-07
  Nonlinear solver reached maximum step 0 with norm 26.250204443587734
  Continuing...
Adaptive step 5, we have 420 active elements and 1701 active dofs.
Postprocessing: 
The relative error in QoI 0 is 5.2931276195391376e-07
The relative error in QoI 1 is 1.8488725104388126e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Assembling the System
Nonlinear Residual: 1515.0737761863265
Linear solve starting, tolerance 9.9999999999999998e-13
Linear solve finished, step 18, residual 8.4106935535299839e-10
Trying full Newton step
  Current Residual: 2.5087225679027682e-07
  Nonlinear solver converged, step 0, residual reduction 1.655841852281021e-10 < 1.0000000000000001e-09
  Nonlinear solver relative step size 6.6469881637661968e-05 > 1e-08
Adaptive step 6, we have 486 active elements and 1951 active dofs.
Postprocessing: 
The relative error in QoI 0 is 5.1247810322998242e-08
The relative error in QoI 1 is 1.5614371976136709e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
[1] Completing output.
[0] Completing output.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex26-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Tue Feb 22 12:24:15 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.601e+00      1.00264   1.599e+00
Objects:              1.281e+03      1.00000   1.281e+03
Flops:                2.533e+08      1.08533   2.434e+08  4.867e+08
Flops/sec:            1.582e+08      1.08248   1.522e+08  3.044e+08
MPI Messages:         1.208e+03      1.00000   1.208e+03  2.416e+03
MPI Message Lengths:  1.417e+06      1.01556   1.164e+03  2.812e+06
MPI Reductions:       3.057e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.5989e+00 100.0%  4.8672e+08 100.0%  2.416e+03 100.0%  1.164e+03      100.0%  2.673e+03  87.4% 

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

KSPGMRESOrthog       379 1.0 7.9799e-03 1.4 1.06e+07 1.0 0.0e+00 0.0e+00 3.8e+02  0  4  0  0 12   0  4  0  0 14  2602
KSPSetup              42 1.0 2.3699e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              21 1.0 4.6235e-01 1.0 2.53e+08 1.1 8.0e+02 4.2e+02 8.9e+02 29100 33 12 29  29100 33 12 33  1052
PCSetUp               42 1.0 3.7043e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 9.3e+01 22 61  0  0  3  22 61  0  0  3   797
PCSetUpOnBlocks       21 1.0 3.6981e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 9.3e+01 22 61  0  0  3  22 61  0  0  3   798
PCApply              400 1.0 3.9100e-02 1.1 7.86e+07 1.1 0.0e+00 0.0e+00 0.0e+00  2 31  0  0  0   2 31  0  0  0  3846
VecMDot              379 1.0 6.4294e-03 1.6 5.30e+06 1.0 0.0e+00 0.0e+00 3.8e+02  0  2  0  0 12   0  2  0  0 14  1614
VecNorm              449 1.0 6.5663e-03 1.0 6.60e+05 1.0 0.0e+00 0.0e+00 4.5e+02  0  0  0  0 15   0  0  0  0 17   197
VecScale             400 1.0 2.3532e-04 1.0 2.93e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2442
VecCopy               81 1.0 4.4584e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               590 1.0 2.3460e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               49 1.0 7.8201e-05 1.2 7.29e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1827
VecMAXPY             400 1.0 1.5182e-03 1.0 5.86e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0  7557
VecAssemblyBegin     329 1.0 5.8289e-03 1.6 0.00e+00 0.0 1.0e+02 2.1e+02 7.4e+02  0  0  4  1 24   0  0  4  1 27     0
VecAssemblyEnd       329 1.0 1.6284e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      644 1.0 1.1535e-03 1.0 0.00e+00 0.0 1.1e+03 1.0e+03 0.0e+00  0  0 46 40  0   0  0 46 40  0     0
VecScatterEnd        644 1.0 6.6599e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         400 1.0 3.3822e-03 1.0 8.80e+05 1.0 0.0e+00 0.0e+00 4.0e+02  0  0  0  0 13   0  0  0  0 15   510
MatMult              400 1.0 7.5426e-02 1.8 8.92e+06 1.1 8.0e+02 4.2e+02 0.0e+00  4  4 33 12  0   4  4 33 12  0   230
MatSolve             400 1.0 3.6803e-02 1.1 7.86e+07 1.1 0.0e+00 0.0e+00 0.0e+00  2 31  0  0  0   2 31  0  0  0  4086
MatLUFactorNum        21 1.0 1.2201e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 0.0e+00  7 61  0  0  0   7 61  0  0  0  2419
MatILUFactorSym       21 1.0 2.4646e-01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.1e+01 15  0  0  0  1  15  0  0  0  1     0
MatAssemblyBegin      63 1.0 3.0406e-03 1.7 0.00e+00 0.0 9.0e+01 1.6e+03 1.3e+02  0  0  4  5  4   0  0  4  5  5     0
MatAssemblyEnd        63 1.0 3.7935e-03 1.0 0.00e+00 0.0 5.6e+01 1.1e+02 1.5e+02  0  0  2  0  5   0  0  2  0  5     0
MatGetRowIJ           21 1.0 6.9141e-06 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        21 1.0 4.9257e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 7.2e+01  0  0  0  0  2   0  0  0  0  3     0
MatZeroEntries        28 1.0 1.7214e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatTranspose           7 1.0 5.2724e-03 1.0 0.00e+00 0.0 7.0e+01 7.8e+02 6.3e+01  0  0  3  2  2   0  0  3  2  2     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver    28             28       264320     0
      Preconditioner    28             28        19712     0
                 Vec   633            633      4477920     0
         Vec Scatter   190            190       164920     0
           Index Set   318            318       449484     0
   IS L to G Mapping    21             21        40904     0
              Matrix    63             63     28934468     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.62125e-06
Average time for zero size MPI_Send(): 5.48363e-06
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
* Done Running Example  mpirun -np 2 ./ex26-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
