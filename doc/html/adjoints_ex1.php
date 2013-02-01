<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adjoints_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file femparameters.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __fem_parameters_h__
        #define __fem_parameters_h__
        
        #include &lt;limits&gt;
        #include &lt;string&gt;
        
        #include "libmesh/libmesh_common.h"
        #include "libmesh/getpot.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        class FEMParameters
        {
        public:
            FEMParameters() :  
              domainfile("lshaped.xda"),
              coarserefinements(0),            
              solver_quiet(true),
        	reuse_preconditioner(false),
              require_residual_reduction(true),
              min_step_length(1e-5),
              max_linear_iterations(200000), max_nonlinear_iterations(20),
              relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
              initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
              linear_tolerance_multiplier(1.e-3),
              nelem_target(30000), global_tolerance(0.0),
              refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
              refine_uniformly(false),
              max_adaptivesteps(1),      
        	indicator_type("kelly"), patch_reuse(true),
              fe_family(1, "LAGRANGE"), fe_order(1, 1),	
              analytic_jacobians(true), verify_analytic_jacobians(0.0),
              print_solution_norms(false), print_solutions(false),
              print_residual_norms(false), print_residuals(false),
              print_jacobian_norms(false), print_jacobians(false) {}
        
            void read(GetPot &input);
           
            std::string domainfile;
            unsigned int coarserefinements;
                    
            bool solver_quiet, reuse_preconditioner, require_residual_reduction;
            Real min_step_length;
            unsigned int max_linear_iterations, max_nonlinear_iterations;
            Real relative_step_tolerance, relative_residual_tolerance,
        	 initial_linear_tolerance, minimum_linear_tolerance,
        	 linear_tolerance_multiplier;
        
            unsigned int nelem_target;
            Real global_tolerance;
            Real refine_fraction, coarsen_fraction, coarsen_threshold;
            bool refine_uniformly;
            unsigned int max_adaptivesteps;
            
            std::string indicator_type;    
            bool patch_reuse;
        
            std::vector&lt;std::string&gt; fe_family;
            std::vector&lt;unsigned int&gt; fe_order;
        
            bool analytic_jacobians;
            Real verify_analytic_jacobians;
        
            bool print_solution_norms, print_solutions,
                 print_residual_norms, print_residuals,
                 print_jacobian_norms, print_jacobians;
        };
        
        #endif // __fem_parameters_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file L-shaped.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/enum_fe_family.h"
        #include "libmesh/fem_system.h"
        #include "libmesh/qoi_set.h"
        #include "libmesh/system.h"
        
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
FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
but we must specify element residuals
</div>

<div class ="fragment">
<pre>
        class LaplaceSystem : public FEMSystem
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          LaplaceSystem(EquationSystems& es,
                       const std::string& name,
                       const unsigned int number)
          : FEMSystem(es, name, number),
            _fe_family("LAGRANGE"), _fe_order(1),
            _analytic_jacobians(true) { qoi.resize(2); }
        
          std::string & fe_family() { return _fe_family;  }
          unsigned int & fe_order() { return _fe_order;  }
          bool & analytic_jacobians() { return _analytic_jacobians; }
        
</pre>
</div>
<div class = "comment">
Postprocessing function which we are going to override for this application


<br><br></div>

<div class ="fragment">
<pre>
          virtual void postprocess(void);
        
          Number &get_QoI_value(std::string type, unsigned int QoI_index)
            {
              if(type == "exact")
        	{
        	  return exact_QoI[QoI_index];
        	}
              else
        	{
        	  return computed_QoI[QoI_index];
        	}
            }
        
          protected:
</pre>
</div>
<div class = "comment">
System initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_data ();
        
</pre>
</div>
<div class = "comment">
Context initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_context (DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Element residual and jacobian calculations
Time dependent parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_time_derivative (bool request_jacobian,
        					DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Constraint parts
</div>

<div class ="fragment">
<pre>
          virtual bool side_constraint (bool request_jacobian,
        				DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Overloading the postprocess function


<br><br></div>

<div class ="fragment">
<pre>
          virtual void element_postprocess(DiffContext &context);  
        
          virtual void side_postprocess(DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Overloading the qoi function on elements


<br><br></div>

<div class ="fragment">
<pre>
          virtual void element_qoi_derivative
            (DiffContext &context,
             const QoISet & qois);  
        
</pre>
</div>
<div class = "comment">
Overloading the qoi function on sides


<br><br></div>

<div class ="fragment">
<pre>
          virtual void side_qoi_derivative
            (DiffContext &context,
             const QoISet & qois);  
        
          Number exact_solution (const Point&);
        
</pre>
</div>
<div class = "comment">
Variables to hold the computed QoIs
  

<br><br></div>

<div class ="fragment">
<pre>
          Number computed_QoI[2];
        
</pre>
</div>
<div class = "comment">
Variables to read in the exact QoIs from l-shaped.in
 

<br><br></div>

<div class ="fragment">
<pre>
          Number exact_QoI[2];
          
</pre>
</div>
<div class = "comment">
The FE type to use
</div>

<div class ="fragment">
<pre>
          std::string _fe_family;
          unsigned int _fe_order;
        
</pre>
</div>
<div class = "comment">
Calculate Jacobians analytically or not?
</div>

<div class ="fragment">
<pre>
          bool _analytic_jacobians;
        };
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file adjoints_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Adjoints Example 1 - Laplace Equation in the L-Shaped Domain with Adjoint based mesh refinement</h1>

<br><br>This example solves the Laplace equation on the classic "L-shaped"
domain with adaptive mesh refinement. The exact
solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta). The kelly and 
adjoint residual error estimators are used to develop error indicators and 
guide mesh adaptation. Since we use the adjoint capabilities of libMesh in 
this example, we use the DiffSystem framework. This file (adjoints_ex1.C)
contains the declaration of mesh and equation system objects, L-shaped.C
contains the assembly of the system, element_qoi_derivative.C and 
side_qoi_derivative.C contain the RHS for the adjoint systems. 
Postprocessing to compute the QoIs is done in element_postprocess.C and 
side_postprocess.C.


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
        #include &lt;iomanip&gt;
        
</pre>
</div>
<div class = "comment">
General libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        #include "libmesh/linear_solver.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/newton_solver.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/steady_solver.h"
        #include "libmesh/system_norm.h"
        
</pre>
</div>
<div class = "comment">
Error Estimator includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/patch_recovery_error_estimator.h"
        
</pre>
</div>
<div class = "comment">
Adjoint Related includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/adjoint_residual_error_estimator.h"
        #include "libmesh/qoi_set.h"
        
</pre>
</div>
<div class = "comment">
libMesh I/O includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/getpot.h"
        #include "libmesh/gmv_io.h"
        
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
        		  unsigned int a_step,       // The adaptive step count
        		  std::string solution_type) // primal or adjoint solve
        {
          MeshBase &mesh = es.get_mesh();
        
        #ifdef LIBMESH_HAVE_GMV
          std::ostringstream file_name_gmv;
          file_name_gmv &lt;&lt; solution_type
                        &lt;&lt; ".out.gmv."
                        &lt;&lt; std::setw(2)
                        &lt;&lt; std::setfill('0')
                        &lt;&lt; std::right
                        &lt;&lt; a_step;
          
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
        AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters &param, QoISet &qois)
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
              
              adjoint_residual_estimator-&gt;qoi_set() = qois;            
        
              adjoint_residual_estimator-&gt;error_plot_suffix = "error.gmv";
              
              PatchRecoveryErrorEstimator *p1 =
        	new PatchRecoveryErrorEstimator;
              adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
        	   	   	   
              PatchRecoveryErrorEstimator *p2 =
        	new PatchRecoveryErrorEstimator;
              adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);   
        
              adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);
              p1-&gt;set_patch_reuse(param.patch_reuse);   
        
              adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);	 
              p2-&gt;set_patch_reuse(param.patch_reuse);
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
          LinearSolver&lt;Number&gt; *linear_solver = system.get_linear_solver();
        
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
        	  libmesh_assert_equal_to (param.nelem_target, 0);
</pre>
</div>
<div class = "comment">
If we aren't adapting to a tolerance we need a
target mesh size
</div>

<div class ="fragment">
<pre>
                    else
                      libmesh_assert_greater (param.nelem_target, 0);
            	
        	linear_solver-&gt;reuse_preconditioner(false);
        
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
Make sure we get the contributions to the adjoint RHS from the sides
</div>

<div class ="fragment">
<pre>
                system.assemble_qoi_sides = true;
        
</pre>
</div>
<div class = "comment">
We are about to solve the adjoint system, but before we do this we see the same preconditioner
flag to reuse the preconditioner from the forward solver		
</div>

<div class ="fragment">
<pre>
                linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
        	
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
Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(true);
        
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
Get a pointer to the solution vector of the adjoint problem for QoI 0
</div>

<div class ="fragment">
<pre>
                NumericVector&lt;Number&gt; &dual_solution_1 = system.get_adjoint_solution(1);
        	
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
        	Number QoI_0_computed = system.get_QoI_value("computed", 0);
        	Number QoI_0_exact = system.get_QoI_value("exact", 0);
        	Number QoI_1_computed = system.get_QoI_value("computed", 1);
        	Number QoI_1_exact = system.get_QoI_value("exact", 1);
        		
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
        	  build_error_estimator(param, qois);
        	
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
        	linear_solver-&gt;reuse_preconditioner(false);	
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
        
        	system.assemble_qoi_sides = true;	
        	linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
        	system.adjoint_solve();
        		
</pre>
</div>
<div class = "comment">
Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(true);
        
        	NumericVector&lt;Number&gt; &dual_solution_0 = system.get_adjoint_solution(0);
        	
        	primal_solution.swap(dual_solution_0);	    
        	write_output(equation_systems, a_step, "adjoint_0");
        
        	primal_solution.swap(dual_solution_0);
        	
        	NumericVector&lt;Number&gt; &dual_solution_1 = system.get_adjoint_solution(1);
        	
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
        	
        	Number QoI_0_computed = system.get_QoI_value("computed", 0);
        	Number QoI_0_exact = system.get_QoI_value("exact", 0);
        	Number QoI_1_computed = system.get_QoI_value("computed", 1);
        	Number QoI_1_exact = system.get_QoI_value("exact", 1);
        	
        	std::cout&lt;&lt; "The relative error in QoI 0 is " &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                            std::abs(QoI_0_exact) &lt;&lt; std::endl;
        	
        	std::cout&lt;&lt; "The relative error in QoI 1 is " &lt;&lt; std::setprecision(17)
                         &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                            std::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
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

<a name="comments"></a> 
<br><br><br> <h1> The source file element_postprocess.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        #include "libmesh/elem.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/point.h"
        #include "libmesh/quadrature.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
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
Last modified: Feb 13, 2009


<br><br>Define the postprocess function to compute QoI 0, the integral of the the solution 
over a subdomain


<br><br></div>

<div class ="fragment">
<pre>
        void LaplaceSystem::element_postprocess (DiffContext &context)
        
        {  
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
                    
</pre>
</div>
<div class = "comment">
Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.element_fe_var[0]-&gt;get_JxW();
        
          const std::vector&lt;Point&gt; &xyz = c.element_fe_var[0]-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable


<br><br></div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
</pre>
</div>
<div class = "comment">
The function R = int_{omega} T dR 
omega is a subset of Omega (the whole domain), omega = [0.75, 1.0] x [0.0, 0.25]


<br><br></div>

<div class ="fragment">
<pre>
          Number dQoI_0 = 0.;
          
</pre>
</div>
<div class = "comment">
Loop over quadrature points  
  

<br><br></div>

<div class ="fragment">
<pre>
          for (unsigned int qp = 0; qp != n_qpoints; qp++)    
            {
</pre>
</div>
<div class = "comment">
Get co-ordinate locations of the current quadrature point
</div>

<div class ="fragment">
<pre>
              const Real x = xyz[qp](0);
              const Real y = xyz[qp](1);
        
</pre>
</div>
<div class = "comment">
If in the sub-domain omega, add the contribution to the integral R
</div>

<div class ="fragment">
<pre>
              if(fabs(x - 0.875) &lt;= 0.125 && fabs(y - 0.125) &lt;= 0.125)
          	{      
</pre>
</div>
<div class = "comment">
Get the solution value at the quadrature point
</div>

<div class ="fragment">
<pre>
                    Number T = c.interior_value(0, qp);
        
</pre>
</div>
<div class = "comment">
Update the elemental increment dR for each qp
</div>

<div class ="fragment">
<pre>
                    dQoI_0 += JxW[qp] * T;
          	}
            }
        
</pre>
</div>
<div class = "comment">
Update the computed value of the global functional R, by adding the contribution from this element


<br><br></div>

<div class ="fragment">
<pre>
          computed_QoI[0] = computed_QoI[0] + dQoI_0;    
          
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file element_qoi_derivative.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        #include "libmesh/elem.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/point.h"
        #include "libmesh/quadrature.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
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
We only have one QoI, so we don't bother checking the qois argument
to see if it was requested from us
</div>

<div class ="fragment">
<pre>
        void LaplaceSystem::element_qoi_derivative (DiffContext &context,
                                                    const QoISet & /* qois */)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
          
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.element_fe_var[0]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The basis functions for the element
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;          &phi = c.element_fe_var[0]-&gt;get_phi();
          
</pre>
</div>
<div class = "comment">
The element quadrature points
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point &gt; &q_point = c.element_fe_var[0]-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_T_dofs = c.dof_indices_var[0].size();
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();  
          
</pre>
</div>
<div class = "comment">
Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
we fill in the [0][i] subderivatives, i corresponding to the variable index.
Our system has only one variable, so we only have to fill the [0][0] subderivative
</div>

<div class ="fragment">
<pre>
          DenseSubVector&lt;Number&gt; &Q = *c.elem_qoi_subderivatives[0][0];
              
</pre>
</div>
<div class = "comment">
Loop over the qps
</div>

<div class ="fragment">
<pre>
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              const Real x = q_point[qp](0);
              const Real y = q_point[qp](1);
                          
</pre>
</div>
<div class = "comment">
If in the sub-domain over which QoI 0 is supported, add contributions
to the adjoint rhs
</div>

<div class ="fragment">
<pre>
              if(fabs(x - 0.875) &lt;= 0.125 && fabs(y - 0.125) &lt;= 0.125)
              	{  	        	 	  
        	  for (unsigned int i=0; i != n_T_dofs; i++)
        	    Q(i) += JxW[qp] *phi[i][qp] ;
              	}
                    
            } // end of the quadrature point qp-loop
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file femparameters.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "femparameters.h"
        
        #define GETPOT_INPUT(A) { A = input(#A, A);\
          std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); };
        #define GETPOT_INT_INPUT(A) { A = input(#A, (int)A);\
          std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); };
        #define GETPOT_REGISTER(A) { \
          std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); };
        
        
        void FEMParameters::read(GetPot &input)
        {
            std::vector&lt;std::string&gt; variable_names;
           
            
            GETPOT_INT_INPUT(coarserefinements);  
            GETPOT_INPUT(domainfile);
        
            GETPOT_INPUT(solver_quiet);
            GETPOT_INPUT(reuse_preconditioner);
            GETPOT_INPUT(require_residual_reduction);
            GETPOT_INPUT(min_step_length);
            GETPOT_INT_INPUT(max_linear_iterations);
            GETPOT_INT_INPUT(max_nonlinear_iterations);
            GETPOT_INPUT(relative_step_tolerance);
            GETPOT_INPUT(relative_residual_tolerance);
            GETPOT_INPUT(initial_linear_tolerance);
            GETPOT_INPUT(minimum_linear_tolerance);
            GETPOT_INPUT(linear_tolerance_multiplier);
            GETPOT_INT_INPUT(nelem_target);
            GETPOT_INPUT(global_tolerance);
            GETPOT_INPUT(refine_fraction);
            GETPOT_INPUT(coarsen_fraction);
            GETPOT_INPUT(coarsen_threshold);
            GETPOT_INT_INPUT(max_adaptivesteps);   
            GETPOT_INPUT(refine_uniformly);
            GETPOT_INPUT(indicator_type);
            GETPOT_INPUT(patch_reuse);
        
            GETPOT_REGISTER(fe_family);
            const unsigned int n_fe_family =
              std::max(1u, input.vector_variable_size("fe_family"));
            fe_family.resize(n_fe_family, "LAGRANGE");
            for (unsigned int i=0; i != n_fe_family; ++i)
              fe_family[i]              = input("fe_family", fe_family[i].c_str(), i);
            GETPOT_REGISTER(fe_order);
            const unsigned int n_fe_order =
              input.vector_variable_size("fe_order");
            fe_order.resize(n_fe_order, 1);
            for (unsigned int i=0; i != n_fe_order; ++i)
              fe_order[i]               = input("fe_order", (int)fe_order[i], i);
        
            GETPOT_INPUT(analytic_jacobians);
            GETPOT_INPUT(verify_analytic_jacobians);
            GETPOT_INPUT(print_solution_norms);
            GETPOT_INPUT(print_solutions);
            GETPOT_INPUT(print_residual_norms);
            GETPOT_INPUT(print_residuals);
            GETPOT_INPUT(print_jacobian_norms);
            GETPOT_INPUT(print_jacobians);
        
          std::vector&lt;std::string&gt; bad_variables =
            input.unidentified_arguments(variable_names);
        
          if (libMesh::processor_id() == 0 && !bad_variables.empty())
            {
              std::cerr &lt;&lt; "ERROR: Unrecognized variables:" &lt;&lt; std::endl;
              for (unsigned int i = 0; i != bad_variables.size(); ++i)
                std::cerr &lt;&lt; bad_variables[i] &lt;&lt; std::endl;
              std::cerr &lt;&lt; "not found among recognized variables." &lt;&lt; std::endl;
              for (unsigned int i = 0; i != variable_names.size(); ++i)
                std::cerr &lt;&lt; variable_names[i] &lt;&lt; std::endl;
              libmesh_error();
            }
         }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file L-shaped.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/getpot.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/parallel.h"
        #include "libmesh/fem_context.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
        #include "L-shaped.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        void LaplaceSystem::init_data ()
        {
          this-&gt;add_variable ("T", static_cast&lt;Order&gt;(_fe_order),
                              Utility::string_to_enum&lt;FEFamily&gt;(_fe_family));
        
          GetPot infile("l-shaped.in");
          exact_QoI[0] = infile("QoI_0", 0.0);
          exact_QoI[1] = infile("QoI_1", 0.0);
          
</pre>
</div>
<div class = "comment">
Do the parent's initialization after variables are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        
          this-&gt;time_evolving(0);
        }
        
        void LaplaceSystem::init_context(DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Now make sure we have requested all the data
we need to build the linear system.
</div>

<div class ="fragment">
<pre>
          c.element_fe_var[0]-&gt;get_JxW();
          c.element_fe_var[0]-&gt;get_phi();
          c.element_fe_var[0]-&gt;get_dphi();
        
          c.side_fe_var[0]-&gt;get_JxW();
          c.side_fe_var[0]-&gt;get_phi();
          c.side_fe_var[0]-&gt;get_dphi();
        }
        
        #define optassert(X) {if (!(X)) libmesh_error();}
        
</pre>
</div>
<div class = "comment">
Assemble the element contributions to the stiffness matrix
</div>

<div class ="fragment">
<pre>
        bool LaplaceSystem::element_time_derivative (bool request_jacobian,
        					  DiffContext &context)
        {
</pre>
</div>
<div class = "comment">
Are the jacobians specified analytically ?
</div>

<div class ="fragment">
<pre>
          bool compute_jacobian = request_jacobian && _analytic_jacobians;
        
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.element_fe_var[0]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
Element basis functions
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &dphi = c.element_fe_var[0]-&gt;get_dphi();
          
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_T_dofs = c.dof_indices_var[0].size(); 
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &K = *c.elem_subjacobians[0][0];
          DenseSubVector&lt;Number&gt; &F = *c.elem_subresiduals[0];
        
</pre>
</div>
<div class = "comment">
Now we will build the element Jacobian and residual.
Constructing the residual requires the solution and its
gradient from the previous timestep.  This must be
calculated at each quadrature point by summing the
solution degree-of-freedom values by the appropriate
weight functions.
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {     	
</pre>
</div>
<div class = "comment">
Compute the solution gradient at the Newton iterate
</div>

<div class ="fragment">
<pre>
              Gradient grad_T = c.interior_gradient(0, qp);
        
</pre>
</div>
<div class = "comment">
The residual contribution from this element
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_T_dofs; i++)
                F(i) += JxW[qp] * ( grad_T * dphi[i][qp] ) ;
              if (compute_jacobian)
                for (unsigned int i=0; i != n_T_dofs; i++)
                  for (unsigned int j=0; j != n_T_dofs; ++j)
</pre>
</div>
<div class = "comment">
The analytic jacobian
</div>

<div class ="fragment">
<pre>
                    K(i,j) += JxW[qp] * ( dphi[i][qp] * dphi[j][qp] );
            } // end of the quadrature point qp-loop
          
          return compute_jacobian;
        }
        
</pre>
</div>
<div class = "comment">
Set Dirichlet bcs, side contributions to global stiffness matrix
</div>

<div class ="fragment">
<pre>
        bool LaplaceSystem::side_constraint (bool request_jacobian,
        				  DiffContext &context)
        {
</pre>
</div>
<div class = "comment">
Are the jacobians specified analytically ? 
</div>

<div class ="fragment">
<pre>
          bool compute_jacobian = request_jacobian && _analytic_jacobians;
        
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
          
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.side_fe_var[0]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
Side basis functions
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &phi = c.side_fe_var[0]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
Side Quadrature points
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point &gt; &qside_point = c.side_fe_var[0]-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_T_dofs = c.dof_indices_var[0].size(); 
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &K = *c.elem_subjacobians[0][0];
          DenseSubVector&lt;Number&gt; &F = *c.elem_subresiduals[0];
              
          unsigned int n_qpoints = (c.get_side_qrule())-&gt;n_points();
        
          const Real penalty = 1./(TOLERANCE*TOLERANCE);
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
</pre>
</div>
<div class = "comment">
Compute the solution at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Number T = c.side_value(0, qp);
                    
</pre>
</div>
<div class = "comment">
We get the Dirichlet bcs from the exact solution
</div>

<div class ="fragment">
<pre>
              Number u_dirichlet = exact_solution (qside_point[qp]);
        
</pre>
</div>
<div class = "comment">
The residual from the boundary terms, penalize non-zero temperature
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_T_dofs; i++)
        	F(i) += JxW[qp] * penalty * ( T - u_dirichlet) * phi[i][qp];
              if (compute_jacobian)
        	for (unsigned int i=0; i != n_T_dofs; i++)
        	  for (unsigned int j=0; j != n_T_dofs; ++j)
</pre>
</div>
<div class = "comment">
The analytic jacobian
</div>

<div class ="fragment">
<pre>
                    K(i,j) += JxW[qp] * penalty * phi[i][qp] * phi[j][qp];    
        	
            } // end of the quadrature point qp-loop
          
          return compute_jacobian;
        }
        
</pre>
</div>
<div class = "comment">
Override the default DiffSystem postprocess function to compute the 
approximations to the QoIs 
</div>

<div class ="fragment">
<pre>
        void LaplaceSystem::postprocess()
        {  
</pre>
</div>
<div class = "comment">
Reset the array holding the computed QoIs
</div>

<div class ="fragment">
<pre>
          computed_QoI[0] = 0.0;
          computed_QoI[1] = 0.0;
        
          FEMSystem::postprocess();
        
          Parallel::sum(computed_QoI[0]);
        
          Parallel::sum(computed_QoI[1]);
             
        }
        
</pre>
</div>
<div class = "comment">
The exact solution to the singular problem,
u_exact = r^(2/3)*sin(2*theta/3). We use this to set the Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
        Number LaplaceSystem::exact_solution(const Point& p)// xyz location   
        {
          const Real x1 = p(0);
          const Real x2 = p(1);
        
          Real theta = atan2(x2,x1);
        
</pre>
</div>
<div class = "comment">
Make sure 0 <= theta <= 2*pi
</div>

<div class ="fragment">
<pre>
          if (theta &lt; 0)
            theta += 2. * libMesh::pi;
            
          return pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
                   
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file side_postprocess.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        #include "libmesh/elem.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/point.h"
        #include "libmesh/quadrature.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
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
Define the postprocess function to compute QoI 1, the integral of the the normal
derivative of the solution over part of the boundary


<br><br></div>

<div class ="fragment">
<pre>
        void LaplaceSystem::side_postprocess(DiffContext &context)
        {
          
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
          
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.side_fe_var[0]-&gt;get_JxW();
        
          const std::vector&lt;Point &gt; &q_point = c.side_fe_var[0]-&gt;get_xyz();
        
          const std::vector&lt;Point&gt; &face_normals = c.side_fe_var[0]-&gt;get_normals();
        
          unsigned int n_qpoints = (c.get_side_qrule())-&gt;n_points();  
          
          Number dQoI_1 = 0. ;
        
</pre>
</div>
<div class = "comment">
Loop over qp's, compute the function at each qp and add
to get the QoI
</div>

<div class ="fragment">
<pre>
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              const Real x = q_point[qp](0);
              const Real y = q_point[qp](1);
        
              const Real TOL = 1.e-5;
              
</pre>
</div>
<div class = "comment">
If on the bottom horizontal bdry (y = -1)
</div>

<div class ="fragment">
<pre>
              if(fabs(y - 1.0) &lt;= TOL && x &gt; 0.0)
        	{
</pre>
</div>
<div class = "comment">
Get the value of the gradient at this point
</div>

<div class ="fragment">
<pre>
                  const Gradient grad_T = c.side_gradient(0,qp);          
        
</pre>
</div>
<div class = "comment">
Add the contribution of this qp to the integral QoI
</div>

<div class ="fragment">
<pre>
                  dQoI_1 += JxW[qp] * (grad_T * face_normals[qp]) ;          
        	}
              
            } // end of the quadrature point qp-loop
          
          computed_QoI[1] = computed_QoI[1] + dQoI_1; 
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file side_qoi_derivative.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        #include "libmesh/elem.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/point.h"
        #include "libmesh/quadrature.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
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
We only have one QoI, so we don't bother checking the qois argument
to see if it was requested from us
</div>

<div class ="fragment">
<pre>
        void LaplaceSystem::side_qoi_derivative (DiffContext &context,
                                                 const QoISet & /* qois */)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
          
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.side_fe_var[0]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The basis functions for the side
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &dphi = c.side_fe_var[0]-&gt;get_dphi();
          
</pre>
</div>
<div class = "comment">
The side quadrature points
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point &gt; &q_point = c.side_fe_var[0]-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
Get the normal to the side at each qp
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt; &face_normals = c.side_fe_var[0]-&gt;get_normals();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_T_dofs = c.dof_indices_var[0].size();
          unsigned int n_qpoints = (c.get_side_qrule())-&gt;n_points();  
        
</pre>
</div>
<div class = "comment">
Fill the QoI RHS corresponding to this QoI. Since this is QoI 1
we fill in the [1][i] subderivatives, i corresponding to the variable index.
Our system has only one variable, so we only have to fill the [1][0] subderivative
</div>

<div class ="fragment">
<pre>
          DenseSubVector&lt;Number&gt; &Q = *c.elem_qoi_subderivatives[1][0];
        
          const Real TOL = 1.e-5;
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              const Real x = q_point[qp](0);
              const Real y = q_point[qp](1);                  
              
</pre>
</div>
<div class = "comment">
If on the sides where the boundary QoI is supported, add contributions
to the adjoint rhs
</div>

<div class ="fragment">
<pre>
              if(fabs(y - 1.0) &lt;= TOL && x &gt; 0.0)
          	{  	  	  	  
          	  for (unsigned int i=0; i != n_T_dofs; i++)
          	    Q(i) += JxW[qp] * (dphi[i][qp] * face_normals[qp]);  	    	    	    	    	  
          	}      
             
            } // end of the quadrature point qp-loop
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.h without comments: </h1> 
<pre> 
  #ifndef __fem_parameters_h__
  #define __fem_parameters_h__
  
  #include &lt;limits&gt;
  #include &lt;string&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> FEMParameters
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
      FEMParameters() :  
        domainfile(<B><FONT COLOR="#BC8F8F">&quot;lshaped.xda&quot;</FONT></B>),
        coarserefinements(0),            
        solver_quiet(true),
  	reuse_preconditioner(false),
        require_residual_reduction(true),
        min_step_length(1e-5),
        max_linear_iterations(200000), max_nonlinear_iterations(20),
        relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
        initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
        linear_tolerance_multiplier(1.e-3),
        nelem_target(30000), global_tolerance(0.0),
        refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
        refine_uniformly(false),
        max_adaptivesteps(1),      
  	indicator_type(<B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>), patch_reuse(true),
        fe_family(1, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), fe_order(1, 1),	
        analytic_jacobians(true), verify_analytic_jacobians(0.0),
        print_solution_norms(false), print_solutions(false),
        print_residual_norms(false), print_residuals(false),
        print_jacobian_norms(false), print_jacobians(false) {}
  
      <B><FONT COLOR="#228B22">void</FONT></B> read(GetPot &amp;input);
     
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string domainfile;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarserefinements;
              
      <B><FONT COLOR="#228B22">bool</FONT></B> solver_quiet, reuse_preconditioner, require_residual_reduction;
      Real min_step_length;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_linear_iterations, max_nonlinear_iterations;
      Real relative_step_tolerance, relative_residual_tolerance,
  	 initial_linear_tolerance, minimum_linear_tolerance,
  	 linear_tolerance_multiplier;
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nelem_target;
      Real global_tolerance;
      Real refine_fraction, coarsen_fraction, coarsen_threshold;
      <B><FONT COLOR="#228B22">bool</FONT></B> refine_uniformly;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_adaptivesteps;
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string indicator_type;    
      <B><FONT COLOR="#228B22">bool</FONT></B> patch_reuse;
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; fe_family;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; fe_order;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> analytic_jacobians;
      Real verify_analytic_jacobians;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> print_solution_norms, print_solutions,
           print_residual_norms, print_residuals,
           print_jacobian_norms, print_jacobians;
  };
  
  #endif <I><FONT COLOR="#B22222">// __fem_parameters_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file L-shaped.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/enum_fe_family.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/qoi_set.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> LaplaceSystem : <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    LaplaceSystem(EquationSystems&amp; es,
                 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
                 <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : FEMSystem(es, name, number),
      _fe_family(<B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), _fe_order(1),
      _analytic_jacobians(true) { qoi.resize(2); }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string &amp; fe_family() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_family;  }
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; fe_order() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_order;  }
    <B><FONT COLOR="#228B22">bool</FONT></B> &amp; analytic_jacobians() { <B><FONT COLOR="#A020F0">return</FONT></B> _analytic_jacobians; }
  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> postprocess(<B><FONT COLOR="#228B22">void</FONT></B>);
  
    Number &amp;get_QoI_value(std::string type, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> QoI_index)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B>(type == <B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>)
  	{
  	  <B><FONT COLOR="#A020F0">return</FONT></B> exact_QoI[QoI_index];
  	}
        <B><FONT COLOR="#A020F0">else</FONT></B>
  	{
  	  <B><FONT COLOR="#A020F0">return</FONT></B> computed_QoI[QoI_index];
  	}
      }
  
    <B><FONT COLOR="#228B22">protected</FONT></B>:
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context (DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  					DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> side_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  				DiffContext &amp;context);
  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> element_postprocess(DiffContext &amp;context);  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> side_postprocess(DiffContext &amp;context);
  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> element_qoi_derivative
      (DiffContext &amp;context,
       <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois);  
  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> side_qoi_derivative
      (DiffContext &amp;context,
       <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois);  
  
    Number exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;);
  
    
    Number computed_QoI[2];
  
   
    Number exact_QoI[2];
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string _fe_family;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _fe_order;
  
    <B><FONT COLOR="#228B22">bool</FONT></B> _analytic_jacobians;
  };
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoints_ex1.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;iomanip&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system_norm.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/patch_recovery_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/adjoint_residual_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/qoi_set.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output(EquationSystems &amp;es,		                   
  		  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step,       <I><FONT COLOR="#B22222">// The adaptive step count
</FONT></I>  		  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type) <I><FONT COLOR="#B22222">// primal or adjoint solve
</FONT></I>  {
    MeshBase &amp;mesh = es.get_mesh();
  
  #ifdef LIBMESH_HAVE_GMV
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name_gmv;
    file_name_gmv &lt;&lt; solution_type
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.out.gmv.&quot;</FONT></B>
                  &lt;&lt; std::setw(2)
                  &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                  &lt;&lt; std::right
                  &lt;&lt; a_step;
    
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
  
  AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters &amp;param, QoISet &amp;qois)
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
        
        adjoint_residual_estimator-&gt;qoi_set() = qois;            
  
        adjoint_residual_estimator-&gt;error_plot_suffix = <B><FONT COLOR="#BC8F8F">&quot;error.gmv&quot;</FONT></B>;
        
        PatchRecoveryErrorEstimator *p1 =
  	<B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
        adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
  	   	   	   
        PatchRecoveryErrorEstimator *p2 =
  	<B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
        adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);   
  
        adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);
        p1-&gt;set_patch_reuse(param.patch_reuse);   
  
        adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(0, H1_SEMINORM);	 
        p2-&gt;set_patch_reuse(param.patch_reuse);
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
  
    LaplaceSystem &amp;system = equation_systems.add_system&lt;LaplaceSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;LaplaceSystem&quot;</FONT></B>);
  
    set_system_parameters(system, param);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initializing systems&quot;</FONT></B> &lt;&lt; std::endl;
  
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
    LinearSolver&lt;Number&gt; *linear_solver = system.get_linear_solver();
  
    {	
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step = 0;
      <B><FONT COLOR="#A020F0">for</FONT></B> (; a_step != param.max_adaptivesteps; ++a_step)
        {	    	    
  	<B><FONT COLOR="#A020F0">if</FONT></B> (param.global_tolerance != 0.)
  	  libmesh_assert_equal_to (param.nelem_target, 0);
              <B><FONT COLOR="#A020F0">else</FONT></B>
                libmesh_assert_greater (param.nelem_target, 0);
      	
  	linear_solver-&gt;reuse_preconditioner(false);
  
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
  
  	system.assemble_qoi_sides = true;
  
  	linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
  	
  	system.adjoint_solve();
  	
  	system.set_adjoint_already_solved(true);
  
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 = system.get_adjoint_solution(0);
  
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	NumericVector&lt;Number&gt; &amp;dual_solution_1 = system.get_adjoint_solution(1);
  	
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
  	Number QoI_0_computed = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Number QoI_0_exact = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Number QoI_1_computed = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Number QoI_1_exact = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_0_exact) &lt;&lt; std::endl;
  		
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17) 
                   &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
  	
  	ErrorVector error;
  	
  	AutoPtr&lt;ErrorEstimator&gt; error_estimator =
  	  build_error_estimator(param, qois);
  	
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
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Refined mesh to &quot;</FONT></B>
                    &lt;&lt; mesh.n_active_elem()
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                    &lt;&lt; equation_systems.n_active_dofs()
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl;
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (a_step == param.max_adaptivesteps)
        {	 
  	linear_solver-&gt;reuse_preconditioner(false);	
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
  
  	system.assemble_qoi_sides = true;	
  	linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
  	system.adjoint_solve();
  		
  	system.set_adjoint_already_solved(true);
  
  	NumericVector&lt;Number&gt; &amp;dual_solution_0 = system.get_adjoint_solution(0);
  	
  	primal_solution.swap(dual_solution_0);	    
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint_0&quot;</FONT></B>);
  
  	primal_solution.swap(dual_solution_0);
  	
  	NumericVector&lt;Number&gt; &amp;dual_solution_1 = system.get_adjoint_solution(1);
  	
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
  	
  	Number QoI_0_computed = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 0);
  	Number QoI_0_exact = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 0);
  	Number QoI_1_computed = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;computed&quot;</FONT></B>, 1);
  	Number QoI_1_exact = system.get_QoI_value(<B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>, 1);
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 0 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_0_computed - QoI_0_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_0_exact) &lt;&lt; std::endl;
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The relative error in QoI 1 is &quot;</FONT></B> &lt;&lt; std::setprecision(17)
                   &lt;&lt; std::abs(QoI_1_computed - QoI_1_exact) /
                      <B><FONT COLOR="#5F9EA0">std</FONT></B>::abs(QoI_1_exact) &lt;&lt; std::endl &lt;&lt; std::endl;
        }
    }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id() 
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Completing output.&quot;</FONT></B> &lt;&lt; std::endl; 
  
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>      
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file element_postprocess.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::element_postprocess (DiffContext &amp;context)
  
  {  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
              
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;xyz = c.element_fe_var[0]-&gt;get_xyz();
  
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
  
    Number dQoI_0 = 0.;
    
    
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp = 0; qp != n_qpoints; qp++)    
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real x = xyz[qp](0);
        <B><FONT COLOR="#228B22">const</FONT></B> Real y = xyz[qp](1);
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(fabs(x - 0.875) &lt;= 0.125 &amp;&amp; fabs(y - 0.125) &lt;= 0.125)
    	{      
    	  Number T = c.interior_value(0, qp);
  
    	  dQoI_0 += JxW[qp] * T;
    	}
      }
  
  
    computed_QoI[0] = computed_QoI[0] + dQoI_0;    
    
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file element_qoi_derivative.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::element_qoi_derivative (DiffContext &amp;context,
                                              <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
    
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;          &amp;phi = c.element_fe_var[0]-&gt;get_phi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;q_point = c.element_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_T_dofs = c.dof_indices_var[0].size();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();  
    
    DenseSubVector&lt;Number&gt; &amp;Q = *c.elem_qoi_subderivatives[0][0];
        
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
        <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);
                    
        <B><FONT COLOR="#A020F0">if</FONT></B>(fabs(x - 0.875) &lt;= 0.125 &amp;&amp; fabs(y - 0.125) &lt;= 0.125)
        	{  	        	 	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
  	    Q(i) += JxW[qp] *phi[i][qp] ;
        	}
              
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  
  #define GETPOT_INPUT(A) { A = input(#A, A);\
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); };
  #define GETPOT_INT_INPUT(A) { A = input(#A, (<B><FONT COLOR="#228B22">int</FONT></B>)A);\
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); };
  #define GETPOT_REGISTER(A) { \
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); };
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> FEMParameters::read(GetPot &amp;input)
  {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; variable_names;
     
      
      GETPOT_INT_INPUT(coarserefinements);  
      GETPOT_INPUT(domainfile);
  
      GETPOT_INPUT(solver_quiet);
      GETPOT_INPUT(reuse_preconditioner);
      GETPOT_INPUT(require_residual_reduction);
      GETPOT_INPUT(min_step_length);
      GETPOT_INT_INPUT(max_linear_iterations);
      GETPOT_INT_INPUT(max_nonlinear_iterations);
      GETPOT_INPUT(relative_step_tolerance);
      GETPOT_INPUT(relative_residual_tolerance);
      GETPOT_INPUT(initial_linear_tolerance);
      GETPOT_INPUT(minimum_linear_tolerance);
      GETPOT_INPUT(linear_tolerance_multiplier);
      GETPOT_INT_INPUT(nelem_target);
      GETPOT_INPUT(global_tolerance);
      GETPOT_INPUT(refine_fraction);
      GETPOT_INPUT(coarsen_fraction);
      GETPOT_INPUT(coarsen_threshold);
      GETPOT_INT_INPUT(max_adaptivesteps);   
      GETPOT_INPUT(refine_uniformly);
      GETPOT_INPUT(indicator_type);
      GETPOT_INPUT(patch_reuse);
  
      GETPOT_REGISTER(fe_family);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_fe_family =
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::max(1u, input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;fe_family&quot;</FONT></B>));
      fe_family.resize(n_fe_family, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_fe_family; ++i)
        fe_family[i]              = input(<B><FONT COLOR="#BC8F8F">&quot;fe_family&quot;</FONT></B>, fe_family[i].c_str(), i);
      GETPOT_REGISTER(fe_order);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_fe_order =
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;fe_order&quot;</FONT></B>);
      fe_order.resize(n_fe_order, 1);
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_fe_order; ++i)
        fe_order[i]               = input(<B><FONT COLOR="#BC8F8F">&quot;fe_order&quot;</FONT></B>, (<B><FONT COLOR="#228B22">int</FONT></B>)fe_order[i], i);
  
      GETPOT_INPUT(analytic_jacobians);
      GETPOT_INPUT(verify_analytic_jacobians);
      GETPOT_INPUT(print_solution_norms);
      GETPOT_INPUT(print_solutions);
      GETPOT_INPUT(print_residual_norms);
      GETPOT_INPUT(print_residuals);
      GETPOT_INPUT(print_jacobian_norms);
      GETPOT_INPUT(print_jacobians);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; bad_variables =
      input.unidentified_arguments(variable_names);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0 &amp;&amp; !bad_variables.empty())
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Unrecognized variables:&quot;</FONT></B> &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != bad_variables.size(); ++i)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; bad_variables[i] &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;not found among recognized variables.&quot;</FONT></B> &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != variable_names.size(); ++i)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; variable_names[i] &lt;&lt; std::endl;
        libmesh_error();
      }
   }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file L-shaped.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parallel.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::init_data ()
  {
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;T&quot;</FONT></B>, static_cast&lt;Order&gt;(_fe_order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(_fe_family));
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;l-shaped.in&quot;</FONT></B>);
    exact_QoI[0] = infile(<B><FONT COLOR="#BC8F8F">&quot;QoI_0&quot;</FONT></B>, 0.0);
    exact_QoI[1] = infile(<B><FONT COLOR="#BC8F8F">&quot;QoI_1&quot;</FONT></B>, 0.0);
    
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(0);
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::init_context(DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    c.element_fe_var[0]-&gt;get_JxW();
    c.element_fe_var[0]-&gt;get_phi();
    c.element_fe_var[0]-&gt;get_dphi();
  
    c.side_fe_var[0]-&gt;get_JxW();
    c.side_fe_var[0]-&gt;get_phi();
    c.side_fe_var[0]-&gt;get_dphi();
  }
  
  #define optassert(X) {<B><FONT COLOR="#A020F0">if</FONT></B> (!(X)) libmesh_error();}
  
  <B><FONT COLOR="#228B22">bool</FONT></B> LaplaceSystem::element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  					  DiffContext &amp;context)
  {
    <B><FONT COLOR="#228B22">bool</FONT></B> compute_jacobian = request_jacobian &amp;&amp; _analytic_jacobians;
  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &amp;dphi = c.element_fe_var[0]-&gt;get_dphi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_T_dofs = c.dof_indices_var[0].size(); 
  
    DenseSubMatrix&lt;Number&gt; &amp;K = *c.elem_subjacobians[0][0];
    DenseSubVector&lt;Number&gt; &amp;F = *c.elem_subresiduals[0];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {     	
        Gradient grad_T = c.interior_gradient(0, qp);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
          F(i) += JxW[qp] * ( grad_T * dphi[i][qp] ) ;
        <B><FONT COLOR="#A020F0">if</FONT></B> (compute_jacobian)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_T_dofs; ++j)
              K(i,j) += JxW[qp] * ( dphi[i][qp] * dphi[j][qp] );
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> compute_jacobian;
  }
  
  <B><FONT COLOR="#228B22">bool</FONT></B> LaplaceSystem::side_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  				  DiffContext &amp;context)
  {
    <B><FONT COLOR="#228B22">bool</FONT></B> compute_jacobian = request_jacobian &amp;&amp; _analytic_jacobians;
  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
    
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.side_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;phi = c.side_fe_var[0]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;qside_point = c.side_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_T_dofs = c.dof_indices_var[0].size(); 
  
    DenseSubMatrix&lt;Number&gt; &amp;K = *c.elem_subjacobians[0][0];
    DenseSubVector&lt;Number&gt; &amp;F = *c.elem_subresiduals[0];
        
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_side_qrule())-&gt;n_points();
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1./(TOLERANCE*TOLERANCE);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Number T = c.side_value(0, qp);
              
        Number u_dirichlet = exact_solution (qside_point[qp]);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
  	F(i) += JxW[qp] * penalty * ( T - u_dirichlet) * phi[i][qp];
        <B><FONT COLOR="#A020F0">if</FONT></B> (compute_jacobian)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_T_dofs; ++j)
  	    K(i,j) += JxW[qp] * penalty * phi[i][qp] * phi[j][qp];    
  	
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> compute_jacobian;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::postprocess()
  {  
    computed_QoI[0] = 0.0;
    computed_QoI[1] = 0.0;
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::postprocess();
  
    <B><FONT COLOR="#5F9EA0">Parallel</FONT></B>::sum(computed_QoI[0]);
  
    <B><FONT COLOR="#5F9EA0">Parallel</FONT></B>::sum(computed_QoI[1]);
       
  }
  
  Number LaplaceSystem::exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)<I><FONT COLOR="#B22222">// xyz location   
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x1 = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real x2 = p(1);
  
    Real theta = atan2(x2,x1);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
      theta += 2. * libMesh::pi;
      
    <B><FONT COLOR="#A020F0">return</FONT></B> pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
             
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file side_postprocess.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::side_postprocess(DiffContext &amp;context)
  {
    
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
    
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.side_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;q_point = c.side_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;face_normals = c.side_fe_var[0]-&gt;get_normals();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_side_qrule())-&gt;n_points();  
    
    Number dQoI_1 = 0. ;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
        <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real TOL = 1.e-5;
        
        <B><FONT COLOR="#A020F0">if</FONT></B>(fabs(y - 1.0) &lt;= TOL &amp;&amp; x &gt; 0.0)
  	{
  	  <B><FONT COLOR="#228B22">const</FONT></B> Gradient grad_T = c.side_gradient(0,qp);	  
  
  	  dQoI_1 += JxW[qp] * (grad_T * face_normals[qp]) ;	  
  	}
        
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>    
    computed_QoI[1] = computed_QoI[1] + dQoI_1; 
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file side_qoi_derivative.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::side_qoi_derivative (DiffContext &amp;context,
                                           <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
    
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.side_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &amp;dphi = c.side_fe_var[0]-&gt;get_dphi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;q_point = c.side_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;face_normals = c.side_fe_var[0]-&gt;get_normals();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_T_dofs = c.dof_indices_var[0].size();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_side_qrule())-&gt;n_points();  
  
    DenseSubVector&lt;Number&gt; &amp;Q = *c.elem_qoi_subderivatives[1][0];
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real TOL = 1.e-5;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
        <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);                  
        
        <B><FONT COLOR="#A020F0">if</FONT></B>(fabs(y - 1.0) &lt;= TOL &amp;&amp; x &gt; 0.0)
    	{  	  	  	  
    	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
    	    Q(i) += JxW[qp] * (dphi[i][qp] * face_normals[qp]);  	    	    	    	    	  
    	}      
       
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example adjoints_ex1:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Started /workspace/libmesh/examples/adjoints/adjoints_ex1/.libs/lt-example-devel
Reading in and building the mesh
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=833
    n_local_nodes()=83
  n_elem()=255
    n_local_elem()=24
    n_active_elem()=192
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=83
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 12.898
      Average Off-Processor Bandwidth <= 2.94358
      Maximum  On-Processor Bandwidth <= 26
      Maximum Off-Processor Bandwidth <= 20
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  Nonlinear solver converged, step 0, residual reduction 9.3416e-11 < 1e-09
Adaptive step 0, we have 192 active elements and 833 active dofs.
Postprocessing: 
The relative error in QoI 0 is 0.00012462746176825957
The relative error in QoI 1 is 0.00014255932627301295

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 222 active elements and 941 active dofs.
  Nonlinear solver converged, step 0, residual reduction 1.7777534988354176e-10 < 1.0000000000000001e-09
Adaptive step 1, we have 222 active elements and 941 active dofs.
Postprocessing: 
The relative error in QoI 0 is 4.9921625044337115e-05
The relative error in QoI 1 is 3.8093557244440054e-05

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 258 active elements and 1077 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.0095509241827086e-10 < 1.0000000000000001e-09
Adaptive step 2, we have 258 active elements and 1077 active dofs.
Postprocessing: 
The relative error in QoI 0 is 2.000891334291647e-05
The relative error in QoI 1 is 4.5494817875652178e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 303 active elements and 1243 active dofs.
  Nonlinear solver converged, step 0, residual reduction 1.7734566692867963e-10 < 1.0000000000000001e-09
Adaptive step 3, we have 303 active elements and 1243 active dofs.
Postprocessing: 
The relative error in QoI 0 is 7.9821047379054138e-06
The relative error in QoI 1 is 1.3954930350315467e-05

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 351 active elements and 1425 active dofs.
  Nonlinear solver converged, step 0, residual reduction 1.7742142618023534e-10 < 1.0000000000000001e-09
Adaptive step 4, we have 351 active elements and 1425 active dofs.
Postprocessing: 
The relative error in QoI 0 is 6.1717814821293693e-06
The relative error in QoI 1 is 9.3546529622144245e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 414 active elements and 1657 active dofs.
  Nonlinear solver reached maximum step 1, latest evaluated residual 2.7625286229497285e-05
  Continuing...
Adaptive step 5, we have 414 active elements and 1657 active dofs.
Postprocessing: 
The relative error in QoI 0 is 2.7129633668951609e-06
The relative error in QoI 1 is 6.5267095938953275e-06

Using Adjoint Residual Error Estimator with Patch Recovery Weights
Refined mesh to 495 active elements and 1969 active dofs.
  Nonlinear solver reached maximum step 1, latest evaluated residual 2.5272056423555221e-05
  Continuing...
Adaptive step 6, we have 495 active elements and 1969 active dofs.
Postprocessing: 
The relative error in QoI 0 is 1.2918613719670313e-06
The relative error in QoI 1 is 1.5052253967931045e-05

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adjoints/adjoints_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:01:48 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           5.051e+00      1.00000   5.051e+00
Objects:              1.125e+03      1.00357   1.125e+03
Flops:                3.473e+07      1.47859   2.933e+07  3.519e+08
Flops/sec:            6.876e+06      1.47859   5.806e+06  6.968e+07
MPI Messages:         7.034e+03      1.44832   6.167e+03  7.401e+04
MPI Message Lengths:  1.483e+06      1.17258   2.293e+02  1.697e+07
MPI Reductions:       3.684e+03      1.00109

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.0506e+00 100.0%  3.5191e+08 100.0%  7.401e+04 100.0%  2.293e+02      100.0%  3.683e+03 100.0% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

KSPGMRESOrthog       947 1.0 2.0976e-02 1.5 6.96e+06 1.3 0.0e+00 0.0e+00 9.5e+02  0 21  0  0 26   0 21  0  0 26  3554
KSPSetUp              35 1.0 4.5276e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               7 1.0 4.8899e-02 1.0 1.77e+07 1.5 2.1e+04 1.2e+02 1.0e+03  1 51 29 16 27   1 51 29 16 27  3665
PCSetUp               28 1.0 2.9595e-02 1.6 4.78e+06 1.9 0.0e+00 0.0e+00 1.0e+02  0 13  0  0  3   0 13  0  0  3  1488
PCSetUpOnBlocks        7 1.0 1.3483e-02 1.7 2.39e+06 1.9 0.0e+00 0.0e+00 3.7e+01  0  6  0  0  1   0  6  0  0  1  1632
PCApply             1009 1.0 3.5747e-02 1.5 2.05e+07 1.6 0.0e+00 0.0e+00 3.7e+01  1 58  0  0  1   1 58  0  0  1  5730
VecMDot              947 1.0 1.8362e-02 1.7 3.47e+06 1.3 0.0e+00 0.0e+00 9.5e+02  0 11  0  0 26   0 11  0  0 26  2025
VecNorm             1037 1.0 4.3935e-02 1.3 2.75e+05 1.3 0.0e+00 0.0e+00 1.0e+03  1  1  0  0 28   1  1  0  0 28    67
VecScale             988 1.0 6.0654e-04 1.5 1.31e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2310
VecCopy              144 1.0 1.4329e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              2322 1.0 9.8157e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               89 1.0 1.4162e-04 1.5 2.33e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1764
VecMAXPY             988 1.0 2.4393e-03 1.3 3.74e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0 11  0  0  0   0 11  0  0  0 16413
VecAssemblyBegin     343 1.0 2.3293e-01 2.0 0.00e+00 0.0 2.1e+03 7.7e+01 7.8e+02  4  0  3  1 21   4  0  3  1 21     0
VecAssemblyEnd       343 1.0 3.5429e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1139 1.0 4.4320e-03 1.3 0.00e+00 0.0 5.5e+04 2.3e+02 0.0e+00  0  0 74 75  0   0  0 74 75  0     0
VecScatterEnd       1139 1.0 1.0535e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         988 1.0 2.2455e-02 1.8 3.93e+05 1.3 0.0e+00 0.0e+00 9.9e+02  0  1  0  0 27   0  1  0  0 27   187
MatMult              479 1.0 1.1908e-02 1.7 2.08e+06 1.4 2.1e+04 1.2e+02 0.0e+00  0  6 29 16  0   0  6 29 16  0  1791
MatMultTranspose     509 1.0 6.5222e-03 1.2 2.12e+06 1.4 2.2e+04 1.2e+02 0.0e+00  0  6 30 16  0   0  6 30 16  0  3363
MatSolve             486 1.0 6.9845e-03 1.5 9.30e+06 1.6 0.0e+00 0.0e+00 0.0e+00  0 26  0  0  0   0 26  0  0  0 13315
MatSolveTranspos     523 1.0 9.3865e-03 1.4 8.86e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0 26  0  0  0   0 26  0  0  0  9568
MatLUFactorNum        14 1.0 6.5181e-03 1.9 4.78e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0  6754
MatILUFactorSym       14 1.0 1.8974e-02 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 4.2e+01  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin      42 1.0 2.7864e-02 2.1 0.00e+00 0.0 9.7e+02 6.0e+02 8.4e+01  0  0  1  3  2   0  0  1  3  2     0
MatAssemblyEnd        42 1.0 2.1966e-03 1.2 0.00e+00 0.0 6.2e+02 3.3e+01 5.6e+01  0  0  1  0  2   0  0  1  0  2     0
MatGetRowIJ           14 1.0 1.9073e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        14 1.0 5.6624e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 3.2e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        28 1.0 8.7261e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver    28             28       271040     0
      Preconditioner    28             28        24976     0
              Vector   758            758      2760880     0
      Vector Scatter    92             92        95312     0
           Index Set   164            164       136952     0
   IS L to G Mapping    19             19        10716     0
              Matrix    35             35      1878696     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.3869e-06
Average time for zero size MPI_Send(): 1.43449e-05
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:01:49 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ---------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=5.18094, Active time=4.93375                                                        |
 ---------------------------------------------------------------------------------------------------------------------
| Event                                   nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                   w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------------|
|                                                                                                                     |
|                                                                                                                     |
| AdjointResidualErrorEstimator                                                                                       |
|   estimate_error()                      6         0.2008      0.033466    1.9832      0.330542    4.07     40.20    |
|                                                                                                                     |
| DofMap                                                                                                              |
|   add_neighbors_to_send_list()          19        0.0344      0.001811    0.0672      0.003539    0.70     1.36     |
|   build_constraint_matrix()             850       0.0157      0.000018    0.0157      0.000018    0.32     0.32     |
|   build_sparsity()                      7         0.0168      0.002404    0.0496      0.007083    0.34     1.00     |
|   cnstrn_elem_mat_vec()                 170       0.0023      0.000014    0.0023      0.000014    0.05     0.05     |
|   constrain_elem_matrix()               170       0.0019      0.000011    0.0019      0.000011    0.04     0.04     |
|   constrain_elem_vector()               510       0.0026      0.000005    0.0026      0.000005    0.05     0.05     |
|   create_dof_constraints()              19        0.1174      0.006180    0.2304      0.012125    2.38     4.67     |
|   distribute_dofs()                     19        0.1935      0.010186    0.5857      0.030827    3.92     11.87    |
|   dof_indices()                         10382     0.6737      0.000065    0.6737      0.000065    13.65    13.65    |
|   enforce_constraints_exactly()         48        0.0171      0.000356    0.0171      0.000356    0.35     0.35     |
|   old_dof_indices()                     1020      0.1186      0.000116    0.1186      0.000116    2.40     2.40     |
|   prepare_send_list()                   19        0.0007      0.000035    0.0007      0.000035    0.01     0.01     |
|   reinit()                              19        0.3302      0.017378    0.3302      0.017378    6.69     6.69     |
|                                                                                                                     |
| EquationSystems                                                                                                     |
|   build_discontinuous_solution_vector() 12        0.0210      0.001749    0.0877      0.007305    0.43     1.78     |
|   build_solution_vector()               21        0.0133      0.000634    0.0903      0.004298    0.27     1.83     |
|                                                                                                                     |
| FE                                                                                                                  |
|   compute_shape_functions()             3461      0.2825      0.000082    0.2825      0.000082    5.72     5.72     |
|   init_shape_functions()                630       0.0424      0.000067    0.0424      0.000067    0.86     0.86     |
|   inverse_map()                         8407      0.0675      0.000008    0.0675      0.000008    1.37     1.37     |
|                                                                                                                     |
| FEMSystem                                                                                                           |
|   assemble_qoi_derivative()             7         0.0468      0.006688    0.1740      0.024852    0.95     3.53     |
|   assembly()                            7         0.0251      0.003583    0.2426      0.034664    0.51     4.92     |
|   assembly(get_jacobian)                7         0.0227      0.003247    0.2390      0.034138    0.46     4.84     |
|   assembly(get_residual)                7         0.0204      0.002914    0.1447      0.020673    0.41     2.93     |
|   numerical_elem_jacobian()             372       0.1547      0.000416    0.1547      0.000416    3.14     3.14     |
|   numerical_side_jacobian()             106       0.0135      0.000127    0.0135      0.000127    0.27     0.27     |
|   postprocess()                         7         0.0134      0.001908    0.1335      0.019072    0.27     2.71     |
|                                                                                                                     |
| FEMap                                                                                                               |
|   compute_affine_map()                  3461      0.0620      0.000018    0.0620      0.000018    1.26     1.26     |
|   compute_face_map()                    530       0.0180      0.000034    0.0447      0.000084    0.36     0.91     |
|   init_face_shape_functions()           70        0.0011      0.000015    0.0011      0.000015    0.02     0.02     |
|   init_reference_to_physical_map()      630       0.0398      0.000063    0.0398      0.000063    0.81     0.81     |
|                                                                                                                     |
| GMVIO                                                                                                               |
|   write_nodal_data()                    21        0.2066      0.009838    0.2092      0.009961    4.19     4.24     |
|                                                                                                                     |
| ImplicitSystem                                                                                                      |
|   adjoint_solve()                       7         0.0025      0.000357    0.4954      0.070772    0.05     10.04    |
|                                                                                                                     |
| LocationMap                                                                                                         |
|   find()                                3280      0.0208      0.000006    0.0208      0.000006    0.42     0.42     |
|   init()                                15        0.0305      0.002035    0.0305      0.002035    0.62     0.62     |
|                                                                                                                     |
| Mesh                                                                                                                |
|   all_first_order()                     12        0.0642      0.005353    0.0642      0.005353    1.30     1.30     |
|   all_second_order()                    1         0.0004      0.000426    0.0004      0.000426    0.01     0.01     |
|   contract()                            6         0.0020      0.000332    0.0136      0.002260    0.04     0.27     |
|   find_neighbors()                      33        0.2768      0.008387    0.3281      0.009944    5.61     6.65     |
|   renumber_nodes_and_elem()             6         0.0116      0.001928    0.0116      0.001928    0.23     0.23     |
|                                                                                                                     |
| MeshCommunication                                                                                                   |
|   broadcast()                           1         0.0011      0.001062    0.0016      0.001581    0.02     0.03     |
|   compute_hilbert_indices()             22        0.0464      0.002111    0.0464      0.002111    0.94     0.94     |
|   find_global_indices()                 22        0.0217      0.000986    0.1100      0.005000    0.44     2.23     |
|   parallel_sort()                       22        0.0169      0.000768    0.0267      0.001212    0.34     0.54     |
|                                                                                                                     |
| MeshOutput                                                                                                          |
|   write_equation_systems()              21        0.0011      0.000052    0.3030      0.014430    0.02     6.14     |
|                                                                                                                     |
| MeshRefinement                                                                                                      |
|   _coarsen_elements()                   12        0.0016      0.000132    0.0018      0.000152    0.03     0.04     |
|   _refine_elements()                    15        0.0265      0.001766    0.0716      0.004772    0.54     1.45     |
|   add_point()                           3280      0.0191      0.000006    0.0405      0.000012    0.39     0.82     |
|   make_coarsening_compatible()          27        0.0496      0.001837    0.0496      0.001837    1.01     1.01     |
|   make_refinement_compatible()          27        0.0031      0.000116    0.0035      0.000128    0.06     0.07     |
|                                                                                                                     |
| MetisPartitioner                                                                                                    |
|   partition()                           21        0.6181      0.029432    0.7328      0.034894    12.53    14.85    |
|                                                                                                                     |
| NewtonSolver                                                                                                        |
|   solve()                               7         0.0800      0.011426    0.5321      0.076016    1.62     10.79    |
|                                                                                                                     |
| Parallel                                                                                                            |
|   allgather()                           115       0.0288      0.000250    0.0294      0.000256    0.58     0.60     |
|   broadcast()                           9         0.0004      0.000045    0.0003      0.000036    0.01     0.01     |
|   max(bool)                             68        0.0048      0.000071    0.0048      0.000071    0.10     0.10     |
|   max(scalar)                           2940      0.0282      0.000010    0.0282      0.000010    0.57     0.57     |
|   max(vector)                           691       0.0104      0.000015    0.0272      0.000039    0.21     0.55     |
|   max(vector<bool>)                     12        0.0010      0.000085    0.0014      0.000120    0.02     0.03     |
|   min(bool)                             3594      0.0313      0.000009    0.0313      0.000009    0.63     0.63     |
|   min(scalar)                           2871      0.2694      0.000094    0.2694      0.000094    5.46     5.46     |
|   min(vector)                           697       0.0108      0.000016    0.0296      0.000043    0.22     0.60     |
|   probe()                               2200      0.0275      0.000013    0.0275      0.000013    0.56     0.56     |
|   receive()                             2200      0.0126      0.000006    0.0406      0.000018    0.26     0.82     |
|   send()                                2200      0.0065      0.000003    0.0065      0.000003    0.13     0.13     |
|   send_receive()                        2244      0.0150      0.000007    0.0672      0.000030    0.30     1.36     |
|   sum()                                 244       0.0345      0.000141    0.0421      0.000173    0.70     0.85     |
|                                                                                                                     |
| Parallel::Request                                                                                                   |
|   wait()                                2200      0.0038      0.000002    0.0038      0.000002    0.08     0.08     |
|                                                                                                                     |
| Partitioner                                                                                                         |
|   set_node_processor_ids()              33        0.0474      0.001435    0.1325      0.004016    0.96     2.69     |
|   set_parent_processor_ids()            21        0.0166      0.000789    0.0166      0.000789    0.34     0.34     |
|                                                                                                                     |
| PatchRecoveryErrorEstimator                                                                                         |
|   estimate_error()                      18        0.1308      0.007269    0.4322      0.024010    2.65     8.76     |
|                                                                                                                     |
| PetscLinearSolver                                                                                                   |
|   solve()                               21        0.1251      0.005958    0.1251      0.005958    2.54     2.54     |
|                                                                                                                     |
| ProjectVector                                                                                                       |
|   operator()                            18        0.0103      0.000575    0.1375      0.007640    0.21     2.79     |
|                                                                                                                     |
| System                                                                                                              |
|   project_vector()                      18        0.0479      0.002664    0.2556      0.014197    0.97     5.18     |
 ---------------------------------------------------------------------------------------------------------------------
| Totals:                                 60265     4.9338                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adjoints_ex1:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
