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
        #include &lt;sys/time.h&gt;
        #include &lt;boost/bind.hpp&gt;
        #include &lt;boost/function.hpp&gt;
        
</pre>
</div>
<div class = "comment">
General libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "twostep_time_solver.h"
        #include "equation_systems.h"
        #include "error_vector.h"
        #include "euler_solver.h"
        #include "euler2_solver.h"
        #include "mesh.h"
        #include "serial_mesh.h"
        #include "mesh_tools.h"
        #include "mesh_base.h"
        #include "mesh_refinement.h"
        #include "newton_solver.h"
        #include "numeric_vector.h"
        #include "petsc_diff_solver.h"
        #include "steady_solver.h"
        #include "system_norm.h"
        #include "tecplot_io.h"
        #include "trilinos_epetra_matrix.h"
        #include "petsc_vector.h"
        
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


<br><br>Write tecplot, gmv and xda/xdr output


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
Create a mesh.
</div>

<div class ="fragment">
<pre>
          SerialMesh mesh (param.dimension);
        
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
          mesh.read("lshaped.xda");
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
                PetscVector&lt;Number&gt; &primal_solution = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(*system.solution);
        
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
                PetscVector&lt;Number&gt; &dual_solution_0 = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(0));
        
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
                PetscVector&lt;Number&gt; &dual_solution_1 = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(1));
        	
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
        	Real QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 0);
        	Real QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 0);
        	Real QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 1);
        	Real QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 1);
        		
        	std::cout&lt;&lt; "The error in QoI 0 is " &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_0_computed - QoI_0_exact)/QoI_0_exact &lt;&lt; std::endl;
        		
        	std::cout&lt;&lt; "The error in QoI 1 is " &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_1_computed - QoI_1_exact)/QoI_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
        	
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
        
        	PetscVector&lt;Number&gt; &primal_solution = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(*system.solution);
        	
        	QoISet qois;
        	std::vector&lt;unsigned int&gt; qoi_indices;
        
        	qoi_indices.push_back(0);
        	qoi_indices.push_back(1);
        	qois.add_indices(qoi_indices);
        	
        	qois.set_weight(0, 0.5);
        	qois.set_weight(1, 0.5);
        
        	system.assemble_qoi_sides = true;
        	system.adjoint_solve();
        		
        	PetscVector&lt;Number&gt; &dual_solution_0 = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(0));
        	
        	primal_solution.swap(dual_solution_0);	    
        	write_output(equation_systems, a_step, "adjoint_0");
        
        	primal_solution.swap(dual_solution_0);
        	
        	PetscVector&lt;Number&gt; &dual_solution_1 = dynamic_cast&lt;PetscVector&lt;Number&gt; &&gt;(system.get_adjoint_solution(1));
        	
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
        	
        	Real QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 0);
        	Real QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 0);
        	Real QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("computed", 1);
        	Real QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&&gt;(system)).get_QoI_value("exact", 1);
        	
        	std::cout&lt;&lt; "The error in QoI 0 is " &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_0_computed - QoI_0_exact)/QoI_0_exact &lt;&lt; std::endl;
        	
        	std::cout&lt;&lt; "The error in QoI 1 is " &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_1_computed - QoI_1_exact)/QoI_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
        	
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
  #include &lt;sys/time.h&gt;
  #include &lt;boost/bind.hpp&gt;
  #include &lt;boost/function.hpp&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;twostep_time_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;euler_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;euler2_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_tools.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;petsc_diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;system_norm.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;trilinos_epetra_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;petsc_vector.h&quot;</FONT></B>
  
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
  
    SerialMesh mesh (param.dimension);
  
    AutoPtr&lt;MeshRefinement&gt; mesh_refinement =
      build_mesh_refinement(mesh, param);
  
    EquationSystems equation_systems (mesh);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Reading in and building the mesh&quot;</FONT></B> &lt;&lt; std::endl;
  
    mesh.read(<B><FONT COLOR="#BC8F8F">&quot;lshaped.xda&quot;</FONT></B>);
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
  	
  	PetscVector&lt;Number&gt; &amp;primal_solution = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(*system.solution);
  
  	QoISet qois;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	system.assemble_qoi_sides = true;
  
  	system.adjoint_solve();
  
  	PetscVector&lt;Number&gt; &amp;dual_solution_0 = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	PetscVector&lt;Number&gt; &amp;dual_solution_1 = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(1));
  	
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
  	Real QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Real QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Real QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Real QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_0_computed - QoI_0_exact)/QoI_0_exact &lt;&lt; std::endl;
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_1_computed - QoI_1_exact)/QoI_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
  	
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
  
  	PetscVector&lt;Number&gt; &amp;primal_solution = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(*system.solution);
  	
  	QoISet qois;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  	
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	system.assemble_qoi_sides = true;
  	system.adjoint_solve();
  		
  	PetscVector&lt;Number&gt; &amp;dual_solution_0 = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(0));
  	
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	PetscVector&lt;Number&gt; &amp;dual_solution_1 = dynamic_cast&lt;PetscVector&lt;Number&gt; &amp;&gt;(system.get_adjoint_solution(1));
  	
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
  	
  	Real QoI_0_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Real QoI_0_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Real QoI_1_computed = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Real QoI_1_exact = (dynamic_cast&lt;LaplaceSystem&amp;&gt;(system)).get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_0_computed - QoI_0_exact)/QoI_0_exact &lt;&lt; std::endl;
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) &lt;&lt; fabs(QoI_1_computed - QoI_1_exact)/QoI_1_exact &lt;&lt; std::endl &lt;&lt; std::endl;
  	
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
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
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
The error in QoI 0 is 0.00016365688099197935
The error in QoI 1 is 0.00032228775655961795

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
The error in QoI 0 is 4.9719364089729268e-05
The error in QoI 1 is 7.8937527867902526e-05

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
The error in QoI 0 is 1.9351576491789132e-05
The error in QoI 1 is 2.6065414311629888e-07

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
The error in QoI 0 is 7.1570030459360749e-06
The error in QoI 1 is 3.8377315458006266e-06

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
The error in QoI 0 is 2.4148442282445063e-06
The error in QoI 1 is 9.8808967341529639e-06

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
The error in QoI 0 is 5.2931276195391376e-07
The error in QoI 1 is 1.8488725104388126e-06

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
The error in QoI 0 is 5.1247810322998242e-08
The error in QoI 1 is 1.5614371976136709e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
[1] Completing output.
[0] Completing output.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex26-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:10 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.843e+00      1.00291   1.841e+00
Objects:              1.281e+03      1.00000   1.281e+03
Flops:                2.533e+08      1.08533   2.434e+08  4.867e+08
Flops/sec:            1.374e+08      1.08218   1.322e+08  2.644e+08
MPI Messages:         1.208e+03      1.00000   1.208e+03  2.416e+03
MPI Message Lengths:  1.417e+06      1.01556   1.164e+03  2.812e+06
MPI Reductions:       3.057e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8406e+00 100.0%  4.8672e+08 100.0%  2.416e+03 100.0%  1.164e+03      100.0%  2.673e+03  87.4% 

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

KSPGMRESOrthog       379 1.0 7.8912e-03 1.4 1.06e+07 1.0 0.0e+00 0.0e+00 3.8e+02  0  4  0  0 12   0  4  0  0 14  2631
KSPSetup              42 1.0 2.3317e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              21 1.0 4.6056e-01 1.0 2.53e+08 1.1 8.0e+02 4.2e+02 8.9e+02 25100 33 12 29  25100 33 12 33  1057
PCSetUp               42 1.0 3.6896e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 9.3e+01 19 61  0  0  3  19 61  0  0  3   800
PCSetUpOnBlocks       21 1.0 3.6835e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 9.3e+01 19 61  0  0  3  19 61  0  0  3   801
PCApply              400 1.0 3.8940e-02 1.1 7.86e+07 1.1 0.0e+00 0.0e+00 0.0e+00  2 31  0  0  0   2 31  0  0  0  3862
VecMDot              379 1.0 6.3627e-03 1.6 5.30e+06 1.0 0.0e+00 0.0e+00 3.8e+02  0  2  0  0 12   0  2  0  0 14  1631
VecNorm              449 1.0 6.5522e-03 1.0 6.60e+05 1.0 0.0e+00 0.0e+00 4.5e+02  0  0  0  0 15   0  0  0  0 17   197
VecScale             400 1.0 2.2340e-04 1.0 2.93e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2572
VecCopy               81 1.0 4.7922e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               590 1.0 2.2936e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               49 1.0 8.2016e-05 1.1 7.29e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1742
VecMAXPY             400 1.0 1.4949e-03 1.0 5.86e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0  7675
VecAssemblyBegin     329 1.0 7.7894e-03 1.6 0.00e+00 0.0 1.0e+02 2.1e+02 7.4e+02  0  0  4  1 24   0  0  4  1 27     0
VecAssemblyEnd       329 1.0 1.5926e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      644 1.0 1.1880e-03 1.0 0.00e+00 0.0 1.1e+03 1.0e+03 0.0e+00  0  0 46 40  0   0  0 46 40  0     0
VecScatterEnd        644 1.0 6.4212e-02 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         400 1.0 3.3655e-03 1.1 8.80e+05 1.0 0.0e+00 0.0e+00 4.0e+02  0  0  0  0 13   0  0  0  0 15   512
MatMult              400 1.0 7.3036e-02 1.7 8.92e+06 1.1 8.0e+02 4.2e+02 0.0e+00  3  4 33 12  0   3  4 33 12  0   238
MatSolve             400 1.0 3.6645e-02 1.1 7.86e+07 1.1 0.0e+00 0.0e+00 0.0e+00  2 31  0  0  0   2 31  0  0  0  4104
MatLUFactorNum        21 1.0 1.2115e-01 1.1 1.54e+08 1.1 0.0e+00 0.0e+00 0.0e+00  6 61  0  0  0   6 61  0  0  0  2436
MatILUFactorSym       21 1.0 2.4588e-01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.1e+01 13  0  0  0  1  13  0  0  0  1     0
MatAssemblyBegin      63 1.0 3.6271e-03 1.4 0.00e+00 0.0 9.0e+01 1.6e+03 1.3e+02  0  0  4  5  4   0  0  4  5  5     0
MatAssemblyEnd        63 1.0 3.7966e-03 1.0 0.00e+00 0.0 5.6e+01 1.1e+02 1.5e+02  0  0  2  0  5   0  0  2  0  5     0
MatGetRowIJ           21 1.0 7.8678e-06 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        21 1.0 4.8995e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 7.2e+01  0  0  0  0  2   0  0  0  0  3     0
MatZeroEntries        28 1.0 1.5402e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatTranspose           7 1.0 5.2490e-03 1.0 0.00e+00 0.0 7.0e+01 7.8e+02 6.3e+01  0  0  3  2  2   0  0  3  2  2     0
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
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 2.38419e-06
Average time for zero size MPI_Send(): 6.07967e-06
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

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:10 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.85177, Active time=1.80323                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| AdjointResidualErrorEstimator                                                                               |
|   estimate_error()              7         0.0594      0.008489    1.0213      0.145901    3.30     56.64    |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  21        0.0020      0.000096    0.0022      0.000106    0.11     0.12     |
|   build_constraint_matrix()     5170      0.0070      0.000001    0.0070      0.000001    0.39     0.39     |
|   cnstrn_elem_mat_vec()         1034      0.0020      0.000002    0.0020      0.000002    0.11     0.11     |
|   compute_sparsity()            7         0.0063      0.000895    0.0073      0.001040    0.35     0.40     |
|   constrain_elem_matrix()       1034      0.0017      0.000002    0.0017      0.000002    0.10     0.10     |
|   constrain_elem_vector()       3102      0.0009      0.000000    0.0009      0.000000    0.05     0.05     |
|   create_dof_constraints()      21        0.0048      0.000226    0.0064      0.000303    0.26     0.35     |
|   distribute_dofs()             21        0.0071      0.000340    0.0183      0.000871    0.40     1.01     |
|   dof_indices()                 117334    0.0494      0.000000    0.0494      0.000000    2.74     2.74     |
|   enforce_constraints_exactly() 48        0.0041      0.000085    0.0041      0.000085    0.23     0.23     |
|   old_dof_indices()             6204      0.0025      0.000000    0.0025      0.000000    0.14     0.14     |
|   prepare_send_list()           21        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   reinit()                      21        0.0104      0.000495    0.0104      0.000495    0.58     0.58     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          95573     0.0958      0.000001    0.0958      0.000001    5.31     5.31     |
|   compute_face_map()            1810      0.0066      0.000004    0.0161      0.000009    0.37     0.89     |
|   compute_shape_functions()     95573     0.0552      0.000001    0.0552      0.000001    3.06     3.06     |
|   init_face_shape_functions()   440       0.0003      0.000001    0.0003      0.000001    0.02     0.02     |
|   init_shape_functions()        5235      0.0627      0.000012    0.0627      0.000012    3.48     3.48     |
|   inverse_map()                 16485     0.0225      0.000001    0.0225      0.000001    1.25     1.25     |
|                                                                                                             |
| FEMSystem                                                                                                   |
|   assemble_qoi_derivative()     7         0.0124      0.001775    0.0271      0.003874    0.69     1.50     |
|   assembly()                    7         0.0173      0.002474    0.0488      0.006966    0.96     2.70     |
|   assembly(get_jacobian)        7         0.0166      0.002375    0.0474      0.006765    0.92     2.63     |
|   assembly(get_residual)        7         0.0087      0.001249    0.0217      0.003100    0.48     1.20     |
|   numerical_elem_jacobian()     2260      0.0229      0.000010    0.0229      0.000010    1.27     1.27     |
|   numerical_side_jacobian()     724       0.0084      0.000012    0.0084      0.000012    0.47     0.47     |
|   postprocess()                 7         0.0057      0.000812    0.0164      0.002350    0.32     0.91     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            21        0.0768      0.003656    0.0768      0.003656    4.26     4.26     |
|                                                                                                             |
| ImplicitSystem                                                                                              |
|   adjoint_solve()               7         0.0073      0.001036    0.3900      0.055716    0.40     21.63    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        3220      0.0019      0.000001    0.0019      0.000001    0.11     0.11     |
|   init()                        15        0.0029      0.000193    0.0029      0.000193    0.16     0.16     |
|                                                                                                             |
| Mesh                                                                                                        |
|   all_first_order()             14        0.0018      0.000132    0.0018      0.000132    0.10     0.10     |
|   all_second_order()            1         0.0000      0.000036    0.0000      0.000036    0.00     0.00     |
|   contract()                    6         0.0001      0.000025    0.0004      0.000061    0.01     0.02     |
|   find_neighbors()              37        0.0236      0.000638    0.0242      0.000655    1.31     1.34     |
|   renumber_nodes_and_elem()     78        0.0024      0.000031    0.0024      0.000031    0.13     0.13     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   broadcast_bcs()               1         0.0000      0.000005    0.0000      0.000006    0.00     0.00     |
|   broadcast_mesh()              1         0.0001      0.000054    0.0001      0.000071    0.00     0.00     |
|   compute_hilbert_indices()     38        0.0179      0.000472    0.0179      0.000472    0.99     0.99     |
|   find_global_indices()         38        0.0033      0.000086    0.0237      0.000623    0.18     1.31     |
|   parallel_sort()               38        0.0015      0.000040    0.0018      0.000048    0.08     0.10     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           12        0.0002      0.000019    0.0003      0.000021    0.01     0.01     |
|   _refine_elements()            15        0.0036      0.000238    0.0084      0.000562    0.20     0.47     |
|   add_point()                   3220      0.0024      0.000001    0.0046      0.000001    0.13     0.25     |
|   make_coarsening_compatible()  15        0.0033      0.000218    0.0033      0.000218    0.18     0.18     |
|   make_refinement_compatible()  15        0.0003      0.000017    0.0003      0.000020    0.01     0.02     |
|                                                                                                             |
| MetisPartitioner                                                                                            |
|   partition()                   37        0.0161      0.000435    0.0398      0.001075    0.89     2.20     |
|                                                                                                             |
| NewtonSolver                                                                                                |
|   solve()                       7         0.0075      0.001065    0.2371      0.033878    0.41     13.15    |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   102       0.0003      0.000003    0.0003      0.000003    0.02     0.02     |
|   broadcast()                   9         0.0000      0.000005    0.0000      0.000003    0.00     0.00     |
|   gather()                      1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(bool)                     49        0.0001      0.000002    0.0001      0.000002    0.01     0.01     |
|   max(scalar)                   112       0.0010      0.000009    0.0010      0.000009    0.05     0.05     |
|   max(vector)                   38        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   min(bool)                     30        0.0001      0.000003    0.0001      0.000003    0.01     0.01     |
|   min(vector)                   44        0.0003      0.000007    0.0003      0.000007    0.02     0.02     |
|   probe()                       234       0.0003      0.000001    0.0003      0.000001    0.02     0.02     |
|   receive()                     234       0.0004      0.000002    0.0007      0.000003    0.02     0.04     |
|   send()                        234       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send_receive()                310       0.0005      0.000002    0.0016      0.000005    0.03     0.09     |
|   sum()                         218       0.0163      0.000075    0.0163      0.000075    0.90     0.90     |
|   wait()                        234       0.0001      0.000000    0.0001      0.000000    0.01     0.01     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      37        0.0072      0.000194    0.0077      0.000208    0.40     0.43     |
|   set_parent_processor_ids()    37        0.0012      0.000033    0.0012      0.000033    0.07     0.07     |
|                                                                                                             |
| PatchRecoveryErrorEstimator                                                                                 |
|   estimate_error()              21        0.6295      0.029978    0.8895      0.042357    34.91    49.33    |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       21        0.4648      0.022133    0.4648      0.022133    25.78    25.78    |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    18        0.0079      0.000440    0.0154      0.000854    0.44     0.85     |
|                                                                                                             |
| System                                                                                                      |
|   project_vector()              18        0.0068      0.000376    0.0251      0.001396    0.38     1.39     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         361017    1.8032                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
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
