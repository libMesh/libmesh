<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adjoints_ex2",$root)?>
 
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
        #include "L-shaped.h"
        
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
        void set_system_parameters(LaplaceSystem &system, FEMParameters &param)
        {
</pre>
</div>
<div class = "comment">
Use analytical jacobians?
</div>

<div class ="fragment">
<pre>
          system.analytic_jacobians() = param.analytic_jacobians;
        
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
Use the prescribed FE type
</div>

<div class ="fragment">
<pre>
          system.fe_family() = param.fe_family[0];
          system.fe_order() = param.fe_order[0];
        
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
        
</pre>
</div>
<div class = "comment">
No transient time solver
</div>

<div class ="fragment">
<pre>
          system.time_solver =
              AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
</pre>
</div>
<div class = "comment">
Nonlinear solver options
</div>

<div class ="fragment">
<pre>
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
        #ifdef LIBMESH_ENABLE_AMR
        
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
        
        #endif // LIBMESH_ENABLE_AMR
        
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
        
</pre>
</div>
<div class = "comment">
Skip adaptive examples on a non-adaptive libMesh build
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
</pre>
</div>
<div class = "comment">
Only our PETSc interface currently supports adjoint solves
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
        
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
          LaplaceSystem &system = equation_systems.add_system&lt;LaplaceSystem&gt; ("LaplaceSystem");
        
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
                NumericVector&lt;Number&gt; &primal_solution = *system.solution;
        	
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
                SensitivityData sensitivities(qois, system, system.get_parameter_vector());
        
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
                system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
        
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
                NumericVector&lt;Number&gt; &dual_solution_0 = system.get_adjoint_solution(0);
        
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
        
                std::cout &lt;&lt; "Refined mesh to "
                          &lt;&lt; mesh.n_active_elem()
                          &lt;&lt; " active elements and "
                          &lt;&lt; equation_systems.n_active_dofs()
                          &lt;&lt; " active dofs." &lt;&lt; std::endl;
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
        
        	NumericVector&lt;Number&gt; &primal_solution = *system.solution;
        				     	
        	QoISet qois;
        	
        	std::vector&lt;unsigned int&gt; qoi_indices;
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
        		
        	NumericVector&lt;Number&gt; &dual_solution_0 = system.get_adjoint_solution(0);
        	
        	primal_solution.swap(dual_solution_0);	    
        	write_output(equation_systems, a_step, "adjoint_0");
        
        	primal_solution.swap(dual_solution_0);									
              }
          }
        
          std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id() 
                    &lt;&lt; "] Completing output." &lt;&lt; std::endl; 
            
        #endif // #ifndef LIBMESH_ENABLE_AMR
            
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
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output(EquationSystems &amp;es,		                   
  		  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step, <I><FONT COLOR="#B22222">// The adaptive step count
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
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> set_system_parameters(LaplaceSystem &amp;system, FEMParameters &amp;param)
  {
    system.analytic_jacobians() = param.analytic_jacobians;
  
    system.verify_analytic_jacobians = param.verify_analytic_jacobians;
  
    system.fe_family() = param.fe_family[0];
    system.fe_order() = param.fe_order[0];
  
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
  
  
  #ifdef LIBMESH_ENABLE_AMR
  
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
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_ENABLE_AMR
</FONT></I>  
  
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
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
  
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
  
    LaplaceSystem &amp;system = equation_systems.add_system&lt;LaplaceSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;LaplaceSystem&quot;</FONT></B>);
  
    set_system_parameters(system, param);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initializing systems&quot;</FONT></B> &lt;&lt; std::endl;
  
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
  
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
  	
  	NumericVector&lt;Number&gt; &amp;primal_solution = *system.solution;
  	
  	QoISet qois;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  
  	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
  
  	system.assemble_qoi_sides = true;
  
  	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
  
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
  
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 = system.get_adjoint_solution(0);
  
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
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
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Refined mesh to &quot;</FONT></B>
                    &lt;&lt; mesh.n_active_elem()
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                    &lt;&lt; equation_systems.n_active_dofs()
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl;
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (a_step == param.max_adaptivesteps)
        {	    
  	system.solve();
  	
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
  
  	NumericVector&lt;Number&gt; &amp;primal_solution = *system.solution;
  				     	
  	QoISet qois;
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
  	qoi_indices.push_back(0);
  	qoi_indices.push_back(1);
  	qois.add_indices(qoi_indices);
  	
  	qois.set_weight(0, 0.5);
  	qois.set_weight(1, 0.5);
  	
  	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
  	
  	system.assemble_qoi_sides = true;
  	
  	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
  	
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
  		
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 = system.get_adjoint_solution(0);
  	
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);									
        }
    }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id() 
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Completing output.&quot;</FONT></B> &lt;&lt; std::endl; 
      
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>      
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) adjoints_ex2.C...
Compiling C++ (in optimized mode) element_qoi.C...
Compiling C++ (in optimized mode) element_qoi_derivative.C...
Compiling C++ (in optimized mode) femparameters.C...
Compiling C++ (in optimized mode) L-shaped.C...
Linking ./adjoints_ex2-opt...
***************************************************************
* Running Example  ./adjoints_ex2-opt
***************************************************************
 
Started ./adjoints_ex2-opt
Reading in and building the mesh
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=833
    n_local_nodes()=833
  n_elem()=255
    n_local_elem()=255
    n_active_elem()=192
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "LaplaceSystem"
    Type "Implicit"
    Variables="T" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=833
    n_local_dofs()=833
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 15.06
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 25
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  Nonlinear solver converged, step 0, residual reduction 1.18002e-13 < 1e-09
Adaptive step 0, we have 192 active elements and 833 active dofs.
Sensitivity of QoI one to Parameter one is 0.00543647
Sensitivity of QoI one to Parameter two is 0.0108729
The relative error in sensitivity QoI_0_0 is 0.00012463024667683093
The relative error in sensitivity QoI_0_1 is 0.00012463023657107274

Refining Uniformly

Refined mesh to 768 active elements and 3201 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.537150574096903e-11 < 1.0000000000000001e-09
Adaptive step 1, we have 768 active elements and 3201 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054368750523949884
Sensitivity of QoI one to Parameter two is 0.010873750103319355
The relative error in sensitivity QoI_0_0 is 5.0121160008591931e-05
The relative error in sensitivity QoI_0_1 is 5.0121295246669185e-05

Refining Uniformly

Refined mesh to 3072 active elements and 12545 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.1735620101414059e-10 < 1.0000000000000001e-09
Adaptive step 2, we have 3072 active elements and 12545 active dofs.
Sensitivity of QoI one to Parameter one is 0.005437039030023595
Sensitivity of QoI one to Parameter two is 0.010874078058810081
The relative error in sensitivity QoI_0_0 is 1.9962400008134908e-05
The relative error in sensitivity QoI_0_1 is 1.9962513772373264e-05

Refining Uniformly

Refined mesh to 12288 active elements and 49665 active dofs.
  Nonlinear solver reached maximum step 1, latest evaluated residual 2.9189471512256663e-05
  Continuing...
Adaptive step 3, we have 12288 active elements and 49665 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054371043992110322
Sensitivity of QoI one to Parameter two is 0.010874208797892675
The error in sensitivity QoI_0_0 is 7.9397012284596045e-06
The error in sensitivity QoI_0_1 is 7.9397499107471927e-06

[0] Completing output.

-------------------------------------------------------------------
| Time:           Sat Apr  7 15:57:41 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=31.0199, Active time=30.5695                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 4         0.0232      0.005793    0.0232      0.005793    0.08     0.08     |
|   compute_sparsity()           4         0.0789      0.019733    0.1073      0.026832    0.26     0.35     |
|   create_dof_constraints()     4         0.0299      0.007480    0.0299      0.007480    0.10     0.10     |
|   distribute_dofs()            4         0.0695      0.017371    0.1695      0.042376    0.23     0.55     |
|   dof_indices()                538176    0.9448      0.000002    0.9448      0.000002    3.09     3.09     |
|   old_dof_indices()            64512     0.0627      0.000001    0.0627      0.000001    0.21     0.21     |
|   prepare_send_list()          4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                     4         0.1000      0.025004    0.1000      0.025004    0.33     0.33     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      8         0.0749      0.009363    0.1265      0.015813    0.25     0.41     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         241920    0.3932      0.000002    0.3932      0.000002    1.29     1.29     |
|   compute_face_map()           13440     0.0833      0.000006    0.1897      0.000014    0.27     0.62     |
|   compute_shape_functions()    241920    0.2437      0.000001    0.2437      0.000001    0.80     0.80     |
|   init_face_shape_functions()  56        0.0002      0.000004    0.0002      0.000004    0.00     0.00     |
|   init_shape_functions()       13496     0.1276      0.000009    0.1276      0.000009    0.42     0.42     |
|   inverse_map()                171142    0.4912      0.000003    0.4912      0.000003    1.61     1.61     |
|                                                                                                            |
| FEMSystem                                                                                                  |
|   assemble_qoi()               20        0.3532      0.017658    1.1107      0.055535    1.16     3.63     |
|   assemble_qoi_derivative()    4         0.1243      0.031071    0.2659      0.066464    0.41     0.87     |
|   assembly()                   8         0.5871      0.073383    0.8331      0.104141    1.92     2.73     |
|   assembly(get_jacobian)       4         0.2520      0.062999    0.3695      0.092369    0.82     1.21     |
|   assembly(get_residual)       20        0.5869      0.029345    1.2550      0.062752    1.92     4.11     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           8         0.7992      0.099905    0.7992      0.099905    2.61     2.61     |
|                                                                                                            |
| ImplicitSystem                                                                                             |
|   adjoint_solve()              4         0.0003      0.000069    8.7111      2.177781    0.00     28.50    |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       81900     0.1718      0.000002    0.1718      0.000002    0.56     0.56     |
|   init()                       9         0.0482      0.005353    0.0482      0.005353    0.16     0.16     |
|                                                                                                            |
| Mesh                                                                                                       |
|   all_second_order()           1         0.0001      0.000118    0.0001      0.000118    0.00     0.00     |
|   contract()                   3         0.0020      0.000656    0.0112      0.003728    0.01     0.04     |
|   find_neighbors()             6         0.1082      0.018027    0.1082      0.018027    0.35     0.35     |
|   renumber_nodes_and_elem()    13        0.0285      0.002192    0.0285      0.002192    0.09     0.09     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     8         0.0002      0.000021    0.9259      0.115739    0.00     3.03     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _coarsen_elements()          3         0.0012      0.000416    0.0012      0.000416    0.00     0.00     |
|   _refine_elements()           9         0.1337      0.014851    0.4250      0.047222    0.44     1.39     |
|   add_point()                  81900     0.0929      0.000001    0.2782      0.000003    0.30     0.91     |
|   make_coarsening_compatible() 3         0.0380      0.012655    0.0380      0.012655    0.12     0.12     |
|   make_refinement_compatible() 3         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| NewtonSolver                                                                                               |
|   solve()                      4         0.0032      0.000806    16.7382     4.184555    0.01     54.75    |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  4         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           6         0.0056      0.000936    0.0056      0.000936    0.02     0.02     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      8         24.1441     3.018018    24.1441     3.018018    78.98    78.98    |
|                                                                                                            |
| ProjectVector                                                                                              |
|   operator()                   6         0.3163      0.052711    0.8268      0.137794    1.03     2.70     |
|                                                                                                            |
| System                                                                                                     |
|   project_vector()             6         0.0496      0.008264    0.9139      0.152323    0.16     2.99     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        1448654   30.5695                                         100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./adjoints_ex2-opt
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
