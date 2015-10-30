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
Linking ./adjoints_ex2-opt...
***************************************************************
* Running Example  mpirun -np 6 ./adjoints_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Started ./adjoints_ex2-opt
Reading in and building the mesh
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=833
    n_local_nodes()=157
  n_elem()=255
    n_local_elem()=47
    n_active_elem()=192
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=157
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 13.8662
      Average Off-Processor Bandwidth <= 1.47771
      Maximum  On-Processor Bandwidth <= 25
      Maximum Off-Processor Bandwidth <= 16
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  Nonlinear solver converged, step 0, residual reduction 9.24088e-13 < 1e-09
Adaptive step 0, we have 192 active elements and 833 active dofs.
Sensitivity of QoI one to Parameter one is 0.00560483
Sensitivity of QoI one to Parameter two is 0.0112097
The relative error in sensitivity QoI_0_0 is 0.030839617732778258
The relative error in sensitivity QoI_0_1 is 0.030839617735495619

Refining Uniformly

Refined mesh to 768 active elements and 3201 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.5906774900778154e-11 < 1.0000000000000001e-09
Adaptive step 1, we have 768 active elements and 3201 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054368759975284182
Sensitivity of QoI one to Parameter two is 0.01087375199426964
The relative error in sensitivity QoI_0_0 is 4.9947331099602652e-05
The relative error in sensitivity QoI_0_1 is 4.9947403489878758e-05

Refining Uniformly

Refined mesh to 3072 active elements and 12545 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.2275937434957728e-10 < 1.0000000000000001e-09
Adaptive step 2, we have 3072 active elements and 12545 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054370393125588775
Sensitivity of QoI one to Parameter two is 0.010874078624478088
The relative error in sensitivity QoI_0_0 is 1.9910436128373573e-05
The relative error in sensitivity QoI_0_1 is 1.9910494951833968e-05

Refining Uniformly

Refined mesh to 12288 active elements and 49665 active dofs.
  Nonlinear solver reached maximum step 1, latest evaluated residual 3.2753467194096035e-05
  Continuing...
Adaptive step 3, we have 12288 active elements and 49665 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054371045617331174
Sensitivity of QoI one to Parameter two is 0.010874209121387543
The error in sensitivity QoI_0_0 is 7.9098101719541892e-06
The error in sensitivity QoI_0_1 is 7.9100013280604072e-06

[[5] Completing output.
[2] Completing output.
[3] Completing output.
[1] Completing output.
4] Completing output.
[0] Completing output.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./adjoints_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:20:58 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           6.072e+00      1.02328   5.964e+00
Objects:              4.540e+02      1.00889   4.507e+02
Flops:                5.893e+08      1.05462   5.752e+08  3.451e+09
Flops/sec:            9.800e+07      1.04133   9.645e+07  5.787e+08
MPI Messages:         2.138e+03      2.01271   1.546e+03  9.276e+03
MPI Message Lengths:  1.565e+06      1.97702   8.099e+02  7.512e+06
MPI Reductions:       1.604e+03      1.00250

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.9641e+00 100.0%  3.4514e+09 100.0%  9.276e+03 100.0%  8.099e+02      100.0%  1.498e+03  93.4% 

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

KSPGMRESOrthog       409 1.0 5.9332e-01 2.0 9.23e+07 1.0 0.0e+00 0.0e+00 4.1e+02  7 16  0  0 25   7 16  0  0 27   915
KSPSetup              16 1.0 4.3488e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               4 1.0 9.3763e-01 1.0 3.00e+08 1.1 3.0e+03 7.8e+02 4.0e+02 16 51 33 32 25  16 51 33 32 27  1878
PCSetUp               16 1.0 5.2967e-01 1.7 1.09e+08 1.1 0.0e+00 0.0e+00 2.5e+01  6 18  0  0  2   6 18  0  0  2  1202
PCSetUpOnBlocks        4 1.0 2.6952e-01 1.8 5.47e+07 1.1 0.0e+00 0.0e+00 1.2e+01  3  9  0  0  1   3  9  0  0  1  1181
PCApply              427 1.0 7.3746e-01 1.4 3.81e+08 1.1 0.0e+00 0.0e+00 1.2e+01 11 65  0  0  1  11 65  0  0  1  3023
VecDot                 8 1.0 7.4172e-04 1.1 4.52e+04 1.1 0.0e+00 0.0e+00 8.0e+00  0  0  0  0  0   0  0  0  0  1   357
VecMDot              409 1.0 5.6058e-01 2.1 4.62e+07 1.0 0.0e+00 0.0e+00 4.1e+02  7  8  0  0 25   7  8  0  0 27   484
VecNorm              451 1.0 4.9220e-01 5.1 3.51e+06 1.0 0.0e+00 0.0e+00 4.5e+02  5  1  0  0 28   5  1  0  0 30    42
VecScale             443 1.0 2.3344e-03 2.5 1.73e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  4357
VecCopy              138 1.0 1.0333e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               981 1.0 5.6994e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               48 1.0 5.2309e-04 1.9 3.24e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3635
VecMAXPY             427 1.0 3.8892e-02 1.8 4.94e+07 1.0 0.0e+00 0.0e+00 0.0e+00  1  8  0  0  0   1  8  0  0  0  7471
VecAssemblyBegin     165 1.0 1.1018e+00 2.9 0.00e+00 0.0 6.1e+02 4.6e+02 4.6e+02 15  0  7  4 29  15  0  7  4 31     0
VecAssemblyEnd       165 1.0 3.2496e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      503 1.0 4.0421e-03 1.5 0.00e+00 0.0 7.9e+03 7.7e+02 0.0e+00  0  0 85 81  0   0  0 85 81  0     0
VecScatterEnd        503 1.0 6.0518e-01 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  7  0  0  0  0   7  0  0  0  0     0
VecNormalize         427 1.0 4.8029e-01 5.7 5.05e+06 1.0 0.0e+00 0.0e+00 4.3e+02  4  1  0  0 27   4  1  0  0 29    62
MatMult              195 1.0 2.7669e-01 5.3 2.66e+07 1.1 3.0e+03 7.8e+02 0.0e+00  2  5 33 32  0   2  5 33 32  0   562
MatMultTranspose     232 1.0 2.6678e-01 2.3 2.63e+07 1.1 3.7e+03 6.6e+02 0.0e+00  4  4 39 32  0   4  4 39 32  0   574
MatSolve             195 1.0 2.4672e-01 1.7 1.67e+08 1.1 0.0e+00 0.0e+00 0.0e+00  3 28  0  0  0   3 28  0  0  0  3977
MatSolveTranspos     232 1.0 3.1684e-01 1.6 1.59e+08 1.1 0.0e+00 0.0e+00 0.0e+00  4 27  0  0  0   4 27  0  0  0  2935
MatLUFactorNum         8 1.0 1.6264e-01 1.7 1.09e+08 1.1 0.0e+00 0.0e+00 0.0e+00  2 18  0  0  0   2 18  0  0  0  3915
MatILUFactorSym        8 1.0 3.6604e-01 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 8.0e+00  4  0  0  0  0   4  0  0  0  1     0
MatAssemblyBegin      20 1.0 7.5204e-02 4.5 0.00e+00 0.0 2.9e+02 3.3e+03 4.0e+01  1  0  3 12  2   1  0  3 12  3     0
MatAssemblyEnd        20 1.0 1.1942e-02 1.3 0.00e+00 0.0 1.3e+02 1.4e+02 4.4e+01  0  0  1  0  3   0  0  1  0  3     0
MatGetRowIJ            8 1.0 1.0014e-05 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         8 1.0 2.1410e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.7e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        20 1.0 2.0044e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver    16             16       151040     0
      Preconditioner    16             16        11264     0
                 Vec   334            334      8753856     0
         Vec Scatter    14             14        12152     0
           Index Set    50             50       307824     0
   IS L to G Mapping     4              4        51128     0
              Matrix    20             20     28643232     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.282e-05
Average time for zero size MPI_Send(): 4.65314e-05
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
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:20:58 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=6.17442, Active time=5.74662                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     4         0.0048      0.001201    0.0051      0.001276    0.08     0.09     |
|   build_sparsity()                 4         0.0362      0.009040    0.0552      0.013798    0.63     0.96     |
|   create_dof_constraints()         4         0.0323      0.008086    0.0323      0.008086    0.56     0.56     |
|   distribute_dofs()                4         0.0377      0.009417    0.1492      0.037310    0.66     2.60     |
|   dof_indices()                    103826    0.0641      0.000001    0.0641      0.000001    1.12     1.12     |
|   old_dof_indices()                10752     0.0052      0.000000    0.0052      0.000000    0.09     0.09     |
|   prepare_send_list()              4         0.0001      0.000025    0.0001      0.000025    0.00     0.00     |
|   reinit()                         4         0.0340      0.008493    0.0340      0.008493    0.59     0.59     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0110      0.001379    0.1484      0.018547    0.19     2.58     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        80276     0.7264      0.000009    0.7264      0.000009    12.64    12.64    |
|   init_shape_functions()           4228      0.0417      0.000010    0.0417      0.000010    0.72     0.72     |
|   inverse_map()                    34622     0.0711      0.000002    0.0711      0.000002    1.24     1.24     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assemble_qoi()                   20        0.1260      0.006301    0.6539      0.032697    2.19     11.38    |
|   assemble_qoi_derivative()        4         0.0744      0.018600    0.1358      0.033954    1.29     2.36     |
|   assembly()                       8         0.1258      0.015723    0.3001      0.037513    2.19     5.22     |
|   assembly(get_jacobian)           4         0.0661      0.016518    0.1568      0.039205    1.15     2.73     |
|   assembly(get_residual)           20        0.1823      0.009115    0.5505      0.027525    3.17     9.58     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             80276     0.1272      0.000002    0.1272      0.000002    2.21     2.21     |
|   compute_face_map()               4116      0.0311      0.000008    0.0612      0.000015    0.54     1.07     |
|   init_face_shape_functions()      112       0.0002      0.000002    0.0002      0.000002    0.00     0.00     |
|   init_reference_to_physical_map() 4228      0.0285      0.000007    0.0285      0.000007    0.50     0.50     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               8         0.5705      0.071310    0.5705      0.071310    9.93     9.93     |
|                                                                                                                |
| ImplicitSystem                                                                                                 |
|   adjoint_solve()                  4         0.0006      0.000149    1.3443      0.336076    0.01     23.39    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           81900     0.0441      0.000001    0.0441      0.000001    0.77     0.77     |
|   init()                           9         0.0401      0.004461    0.0401      0.004461    0.70     0.70     |
|                                                                                                                |
| Mesh                                                                                                           |
|   all_second_order()               1         0.0000      0.000034    0.0000      0.000034    0.00     0.00     |
|   contract()                       3         0.0008      0.000274    0.0056      0.001865    0.01     0.10     |
|   find_neighbors()                 6         0.0441      0.007348    0.0974      0.016237    0.77     1.70     |
|   renumber_nodes_and_elem()        13        0.0225      0.001734    0.0225      0.001734    0.39     0.39     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0001      0.000063    0.0006      0.000597    0.00     0.01     |
|   compute_hilbert_indices()        7         0.0257      0.003668    0.0257      0.003668    0.45     0.45     |
|   find_global_indices()            7         0.0036      0.000521    0.0662      0.009453    0.06     1.15     |
|   parallel_sort()                  7         0.0066      0.000939    0.0328      0.004688    0.11     0.57     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0001      0.000013    0.7190      0.089870    0.00     12.51    |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              3         0.0007      0.000245    0.0036      0.001196    0.01     0.06     |
|   _refine_elements()               9         0.1110      0.012331    0.3650      0.040556    1.93     6.35     |
|   add_point()                      81900     0.0871      0.000001    0.1368      0.000002    1.52     2.38     |
|   make_coarsening_compatible()     3         0.0257      0.008562    0.0257      0.008562    0.45     0.45     |
|   make_refinement_compatible()     3         0.0000      0.000003    0.0064      0.002140    0.00     0.11     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      6         0.0346      0.005759    0.1162      0.019367    0.60     2.02     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          4         0.0778      0.019454    1.3036      0.325895    1.35     22.68    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      27        0.0306      0.001131    0.0306      0.001131    0.53     0.53     |
|   broadcast()                      10        0.0001      0.000008    0.0000      0.000005    0.00     0.00     |
|   gather()                         1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(bool)                        19        0.1223      0.006436    0.1223      0.006436    2.13     2.13     |
|   max(scalar)                      11        0.0600      0.005455    0.0600      0.005455    1.04     1.04     |
|   max(vector)                      7         0.0004      0.000056    0.0004      0.000056    0.01     0.01     |
|   min(bool)                        6         0.0163      0.002718    0.0163      0.002718    0.28     0.28     |
|   min(vector)                      7         0.0168      0.002396    0.0168      0.002396    0.29     0.29     |
|   probe()                          210       0.0671      0.000320    0.0671      0.000320    1.17     1.17     |
|   receive()                        210       0.0189      0.000090    0.0861      0.000410    0.33     1.50     |
|   send()                           210       0.0005      0.000003    0.0005      0.000003    0.01     0.01     |
|   send_receive()                   224       0.0005      0.000002    0.0953      0.000425    0.01     1.66     |
|   sum()                            95        0.3440      0.003621    0.3440      0.003621    5.99     5.99     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           210       0.0082      0.000039    0.0082      0.000039    0.14     0.14     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         6         0.0165      0.002744    0.0603      0.010053    0.29     1.05     |
|   set_parent_processor_ids()       6         0.0034      0.000564    0.0034      0.000564    0.06     0.06     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          8         1.9949      0.249368    1.9949      0.249368    34.72    34.72    |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       6         0.0421      0.007011    0.0914      0.015229    0.73     1.59     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 6         0.0822      0.013697    0.1769      0.029478    1.43     3.08     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            487705    5.7466                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./adjoints_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
