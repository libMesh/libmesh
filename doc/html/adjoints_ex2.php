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
              indicator_type("kelly"),
              fe_family(1, "LAGRANGE"), fe_order(1, 1),	
              analytic_jacobians(true), verify_analytic_jacobians(0.0),
              print_solution_norms(false), print_solutions(false),
              print_residual_norms(false), print_residuals(false),
              print_jacobian_norms(false), print_jacobians(false) {}
        
            void read(GetPot &input);
           
            unsigned int dimension;  
            Real elementorder;    
            std::string domainfile;
            unsigned int coarserefinements;
                    
            bool solver_quiet, require_residual_reduction;
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
<br><br><br> <h1> The source file L-qoi.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef L_QOI_H
        #define L_QOI_H
        
        #include "libmesh/libmesh_common.h"
        #include "libmesh/elem.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/point.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/diff_qoi.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        class LaplaceQoI : public DifferentiableQoI
        {
        public:
          LaplaceQoI(){}
          virtual ~LaplaceQoI(){} 
        
          virtual void init_qoi( std::vector&lt;Number&gt;& sys_qoi );
        
          virtual void postprocess( ){} 
          
          virtual void element_qoi_derivative(DiffContext &context, const QoISet & qois);  
        
          virtual void element_qoi (DiffContext &context, const QoISet & qois); 
        
          virtual AutoPtr&lt;DifferentiableQoI&gt; clone( ) {
            return AutoPtr&lt;DifferentiableQoI&gt; ( new LaplaceQoI(*this) );
          }
        
        };
        #endif // L_QOI_H
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
        #include "libmesh/parameter_vector.h"
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
            _analytic_jacobians(true) { }
        
          std::string & fe_family() { return _fe_family;  }
          unsigned int & fe_order() { return _fe_order;  }
          bool & analytic_jacobians() { return _analytic_jacobians; }
        
          Number &get_parameter_value(unsigned int parameter_index)
            {
              return parameters[parameter_index];
            }
        
          ParameterVector &get_parameter_vector()
            {
              parameter_vector.resize(parameters.size());
              for(unsigned int i = 0; i != parameters.size(); ++i)
        	{
        	  parameter_vector[i] = &parameters[i];
        	}
              
              return parameter_vector;
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
         
          Number exact_solution (const Point&);
        
</pre>
</div>
<div class = "comment">
Parameters associated with the system
</div>

<div class ="fragment">
<pre>
          std::vector&lt;Number&gt; parameters;
          
</pre>
</div>
<div class = "comment">
The ParameterVector object that will contain pointers to
the system parameters
</div>

<div class ="fragment">
<pre>
          ParameterVector parameter_vector;
          
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
<br><br><br> <h1> The source file adjoints_ex2.C with comments: </h1> 
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
        #include &lt;iomanip&gt;
        
</pre>
</div>
<div class = "comment">
General libMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
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
Sensitivity Calculation related includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/parameter_vector.h"
        #include "libmesh/sensitivity_data.h"
        
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
        #include "L-qoi.h"
        
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
        
          QoISet qois;
        	
          std::vector&lt;unsigned int&gt; qoi_indices;
          qoi_indices.push_back(0);
          qois.add_indices(qoi_indices);
          
          qois.set_weight(0, 0.5);
        
</pre>
</div>
<div class = "comment">
Put some scope here to test that the cloning is working right
</div>

<div class ="fragment">
<pre>
          {
            LaplaceQoI qoi;
            system.attach_qoi( &qoi );
          }
        
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
Here we solve the adjoint problem inside the adjoint_qoi_parameter_sensitivity
function, so we have to set the adjoint_already_solved boolean to false
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(false);
        
</pre>
</div>
<div class = "comment">
Compute the sensitivities
</div>

<div class ="fragment">
<pre>
                system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
        
</pre>
</div>
<div class = "comment">
Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(true);
        
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
        	
        	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
        	
        	system.assemble_qoi_sides = true;
        	
</pre>
</div>
<div class = "comment">
Here we solve the adjoint problem inside the adjoint_qoi_parameter_sensitivity
function, so we have to set the adjoint_already_solved boolean to false
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(false);
        
        	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
        	
</pre>
</div>
<div class = "comment">
Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
</div>

<div class ="fragment">
<pre>
                system.set_adjoint_already_solved(true);
        
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
<br><br><br> <h1> The source file L-qoi.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "L-qoi.h"
        
        using namespace libMesh;
        
        void LaplaceQoI::init_qoi( std::vector&lt;Number&gt;& sys_qoi )
        {
</pre>
</div>
<div class = "comment">
Only 1 qoi to worry about
</div>

<div class ="fragment">
<pre>
          sys_qoi.resize(1);
          return;
        }
        
</pre>
</div>
<div class = "comment">
We only have one QoI, so we don't bother checking the qois argument
to see if it was requested from us
</div>

<div class ="fragment">
<pre>
        void LaplaceQoI::element_qoi (DiffContext &context,
        			      const QoISet & /* qois */ )
        
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
        
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
          
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
              const Real xf = xyz[qp](0);
              const Real yf = xyz[qp](1);
        
</pre>
</div>
<div class = "comment">
If in the sub-domain omega, add the contribution to the integral R
</div>

<div class ="fragment">
<pre>
              if(fabs(xf - 0.875) &lt;= 0.125 && fabs(yf - 0.125) &lt;= 0.125)
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
          c.elem_qoi[0] = c.elem_qoi[0] + dQoI_0;    
          
        }
        
</pre>
</div>
<div class = "comment">
We only have one QoI, so we don't bother checking the qois argument
to see if it was requested from us
</div>

<div class ="fragment">
<pre>
        void LaplaceQoI::element_qoi_derivative (DiffContext &context,
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
          parameters.push_back(infile("alpha_1", 1.0));
          parameters.push_back(infile("alpha_2", 1.0));
          
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
Element basis function gradients
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
                F(i) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( grad_T * dphi[i][qp] ) ;
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
                    K(i,j) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( dphi[i][qp] * dphi[j][qp] );
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
            
          return (parameters[0] + (2.*parameters[1])) * pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
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
        indicator_type(<B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>),
        fe_family(1, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), fe_order(1, 1),	
        analytic_jacobians(true), verify_analytic_jacobians(0.0),
        print_solution_norms(false), print_solutions(false),
        print_residual_norms(false), print_residuals(false),
        print_jacobian_norms(false), print_jacobians(false) {}
  
      <B><FONT COLOR="#228B22">void</FONT></B> read(GetPot &amp;input);
     
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dimension;  
      Real elementorder;    
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string domainfile;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarserefinements;
              
      <B><FONT COLOR="#228B22">bool</FONT></B> solver_quiet, require_residual_reduction;
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
<br><br><br> <h1> The source file L-qoi.h without comments: </h1> 
<pre> 
  #ifndef L_QOI_H
  #define L_QOI_H
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/diff_qoi.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> LaplaceQoI : <B><FONT COLOR="#228B22">public</FONT></B> DifferentiableQoI
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    LaplaceQoI(){}
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~LaplaceQoI(){} 
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_qoi( std::vector&lt;Number&gt;&amp; sys_qoi );
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> postprocess( ){} 
    
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> element_qoi_derivative(DiffContext &amp;context, <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois);  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> element_qoi (DiffContext &amp;context, <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois); 
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;DifferentiableQoI&gt; clone( ) {
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;DifferentiableQoI&gt; ( <B><FONT COLOR="#A020F0">new</FONT></B> LaplaceQoI(*<B><FONT COLOR="#A020F0">this</FONT></B>) );
    }
  
  };
  #endif <I><FONT COLOR="#B22222">// L_QOI_H
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file L-shaped.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/enum_fe_family.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameter_vector.h&quot;</FONT></B>
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
      _analytic_jacobians(true) { }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string &amp; fe_family() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_family;  }
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; fe_order() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_order;  }
    <B><FONT COLOR="#228B22">bool</FONT></B> &amp; analytic_jacobians() { <B><FONT COLOR="#A020F0">return</FONT></B> _analytic_jacobians; }
  
    Number &amp;get_parameter_value(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> parameter_index)
      {
        <B><FONT COLOR="#A020F0">return</FONT></B> parameters[parameter_index];
      }
  
    ParameterVector &amp;get_parameter_vector()
      {
        parameter_vector.resize(parameters.size());
        <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != parameters.size(); ++i)
  	{
  	  parameter_vector[i] = &amp;parameters[i];
  	}
        
        <B><FONT COLOR="#A020F0">return</FONT></B> parameter_vector;
      }
   
    <B><FONT COLOR="#228B22">protected</FONT></B>:
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context (DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  					DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> side_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
  				DiffContext &amp;context);
   
    Number exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; parameters;
    
    ParameterVector parameter_vector;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string _fe_family;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _fe_order;
  
    <B><FONT COLOR="#228B22">bool</FONT></B> _analytic_jacobians;
  };
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoints_ex2.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;iomanip&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system_norm.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameter_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sensitivity_data.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/patch_recovery_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/adjoint_residual_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/qoi_set.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;L-shaped.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;L-qoi.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output(EquationSystems &amp;es,		                   
  		  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step, <I><FONT COLOR="#B22222">// The adaptive step count
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
  
    QoISet qois;
  	
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; qoi_indices;
    qoi_indices.push_back(0);
    qois.add_indices(qoi_indices);
    
    qois.set_weight(0, 0.5);
  
    {
      LaplaceQoI qoi;
      system.attach_qoi( &amp;qoi );
    }
  
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
  	  libmesh_assert_equal_to (param.nelem_target, 0);
              <B><FONT COLOR="#A020F0">else</FONT></B>
                libmesh_assert_greater (param.nelem_target, 0);
      
  	system.solve();
  
  	write_output(equation_systems, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
  	
  	NumericVector&lt;Number&gt; &amp;primal_solution = *system.solution;
  
  	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
  
  	system.assemble_qoi_sides = true;
  
  	system.set_adjoint_already_solved(false);
  
  	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
  
  	system.set_adjoint_already_solved(true);
  
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
  	
  	SensitivityData sensitivities(qois, system, system.get_parameter_vector());
  	
  	system.assemble_qoi_sides = true;
  	
  	system.set_adjoint_already_solved(false);
  
  	system.adjoint_qoi_parameter_sensitivity(qois, system.get_parameter_vector(), sensitivities);
  	
  	system.set_adjoint_already_solved(true);
  
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
<br><br><br> <h1> The source file L-qoi.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;L-qoi.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceQoI::init_qoi( std::vector&lt;Number&gt;&amp; sys_qoi )
  {
    sys_qoi.resize(1);
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceQoI::element_qoi (DiffContext &amp;context,
  			      <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I> )
  
  {  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
              
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;xyz = c.element_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
    
    Number dQoI_0 = 0.;
    
    
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp = 0; qp != n_qpoints; qp++)    
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real xf = xyz[qp](0);
        <B><FONT COLOR="#228B22">const</FONT></B> Real yf = xyz[qp](1);
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(fabs(xf - 0.875) &lt;= 0.125 &amp;&amp; fabs(yf - 0.125) &lt;= 0.125)
        	{      
        	  Number T = c.interior_value(0, qp);
  
        	  dQoI_0 += JxW[qp] * T;
        	}
  
      }
  
  
    c.elem_qoi[0] = c.elem_qoi[0] + dQoI_0;    
    
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceQoI::element_qoi_derivative (DiffContext &amp;context,
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
    parameters.push_back(infile(<B><FONT COLOR="#BC8F8F">&quot;alpha_1&quot;</FONT></B>, 1.0));
    parameters.push_back(infile(<B><FONT COLOR="#BC8F8F">&quot;alpha_2&quot;</FONT></B>, 1.0));
    
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
          F(i) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( grad_T * dphi[i][qp] ) ;
        <B><FONT COLOR="#A020F0">if</FONT></B> (compute_jacobian)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_T_dofs; ++j)
              K(i,j) += JxW[qp] * (parameters[0] + (2.*parameters[1])) * ( dphi[i][qp] * dphi[j][qp] );
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
  
  Number LaplaceSystem::exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)<I><FONT COLOR="#B22222">// xyz location   
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x1 = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real x2 = p(1);
  
    Real theta = atan2(x2,x1);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
      theta += 2. * libMesh::pi;
      
    <B><FONT COLOR="#A020F0">return</FONT></B> (parameters[0] + (2.*parameters[1])) * pow(x1*x1 + x2*x2, 1./3.)*sin(2./3.*theta);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example adjoints_ex2:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Started /workspace/libmesh/examples/adjoints/adjoints_ex2/.libs/lt-example-devel
Reading in and building the mesh
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=225
    n_local_nodes()=25
  n_elem()=63
    n_local_elem()=6
    n_active_elem()=48
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
    n_dofs()=225
    n_local_dofs()=25
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.8444
      Average Off-Processor Bandwidth <= 4.74667
      Maximum  On-Processor Bandwidth <= 26
      Maximum Off-Processor Bandwidth <= 20
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  Nonlinear solver converged, step 0, residual reduction 8.53523e-14 < 1e-09
Adaptive step 0, we have 48 active elements and 225 active dofs.
Sensitivity of QoI one to Parameter one is 0.00543563
Sensitivity of QoI one to Parameter two is 0.0108713
The relative error in sensitivity QoI_0_0 is 0.00027992885259054269
The relative error in sensitivity QoI_0_1 is 0.00027992872361382737

Refining Uniformly

Refined mesh to 192 active elements and 833 active dofs.
  Nonlinear solver converged, step 0, residual reduction 3.9494998974879854e-12 < 1.0000000000000001e-09
Adaptive step 1, we have 192 active elements and 833 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054364699551922796
Sensitivity of QoI one to Parameter two is 0.010872939909142702
The relative error in sensitivity QoI_0_0 is 0.0001246266240597838
The relative error in sensitivity QoI_0_1 is 0.00012462673826064247

Refining Uniformly

Refined mesh to 768 active elements and 3201 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.2957525992439153e-11 < 1.0000000000000001e-09
Adaptive step 2, we have 768 active elements and 3201 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054368750660834285
Sensitivity of QoI one to Parameter two is 0.010873750133095982
The relative error in sensitivity QoI_0_0 is 5.0118642431136876e-05
The relative error in sensitivity QoI_0_1 is 5.0118556988513685e-05

Refining Uniformly

Refined mesh to 3072 active elements and 12545 active dofs.
  Nonlinear solver converged, step 0, residual reduction 2.1791378693485278e-10 < 1.0000000000000001e-09
Adaptive step 3, we have 3072 active elements and 12545 active dofs.
Sensitivity of QoI one to Parameter one is 0.0054370389039995774
Sensitivity of QoI one to Parameter two is 0.010874077807717411
The error in sensitivity QoI_0_0 is 1.9985578341263046e-05
The error in sensitivity QoI_0_1 is 1.9985604250060418e-05

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adjoints/adjoints_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:02:38 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           6.166e+00      1.00038   6.166e+00
Objects:              4.860e+02      1.01674   4.820e+02
Flops:                8.585e+07      1.19725   7.845e+07  9.413e+08
Flops/sec:            1.392e+07      1.19726   1.272e+07  1.527e+08
MPI Messages:         3.241e+03      2.11761   2.490e+03  2.988e+04
MPI Message Lengths:  7.581e+05      1.88936   2.251e+02  6.727e+06
MPI Reductions:       1.781e+03      1.00451

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 6.1658e+00 100.0%  9.4134e+08 100.0%  2.988e+04 100.0%  2.251e+02      100.0%  1.776e+03  99.9% 

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

KSPGMRESOrthog       507 1.0 2.6292e-02 1.6 1.66e+07 1.2 0.0e+00 0.0e+00 5.1e+02  0 20  0  0 29   0 20  0  0 29  7006
KSPSetUp              16 1.0 2.6131e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               4 1.0 1.0999e-01 1.0 6.08e+07 1.2 1.6e+04 2.2e+02 7.6e+02  2 71 53 51 42   2 71 53 51 43  6067
PCSetUp               16 1.0 4.4907e-02 1.2 1.06e+07 1.2 0.0e+00 0.0e+00 6.0e+01  1 12  0  0  3   1 12  0  0  3  2577
PCSetUpOnBlocks        4 1.0 2.2063e-02 1.2 5.30e+06 1.2 0.0e+00 0.0e+00 2.2e+01  0  6  0  0  1   0  6  0  0  1  2623
PCApply              537 1.0 9.8479e-02 1.1 5.34e+07 1.2 0.0e+00 0.0e+00 2.2e+01  1 62  0  0  1   1 62  0  0  1  5923
VecDot                 8 1.0 1.8096e-04 1.6 6.11e+03 1.2 0.0e+00 0.0e+00 8.0e+00  0  0  0  0  0   0  0  0  0  0   371
VecMDot              507 1.0 2.0796e-02 1.9 8.31e+06 1.2 0.0e+00 0.0e+00 5.1e+02  0 10  0  0 29   0 10  0  0 29  4427
VecNorm              553 1.0 1.9994e-02 1.3 6.02e+05 1.2 0.0e+00 0.0e+00 5.5e+02  0  1  0  0 31   0  1  0  0 31   334
VecScale             545 1.0 6.2132e-04 1.2 2.98e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  5313
VecCopy              106 1.0 1.6904e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               969 1.0 9.6512e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               56 1.0 2.1434e-04 1.6 5.42e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2792
VecMAXPY             529 1.0 5.2593e-03 1.2 8.88e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0 10  0  0  0   0 10  0  0  0 18703
VecAssemblyBegin     165 1.0 4.2658e-01 3.6 0.00e+00 0.0 1.6e+03 1.3e+02 4.6e+02  5  0  5  3 26   5  0  5  3 26     0
VecAssemblyEnd       165 1.0 3.1042e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      605 1.0 3.1927e-03 1.3 0.00e+00 0.0 2.6e+04 2.2e+02 0.0e+00  0  0 88 85  0   0  0 88 85  0     0
VecScatterEnd        605 1.0 6.6600e-03 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         529 1.0 1.0015e-02 1.8 8.76e+05 1.2 0.0e+00 0.0e+00 5.3e+02  0  1  0  0 30   0  1  0  0 30   969
MatMult              368 1.0 1.7569e-02 1.3 6.61e+06 1.2 1.6e+04 2.2e+02 0.0e+00  0  8 53 51  0   0  8 53 51  0  4111
MatMultTranspose     161 1.0 5.4486e-03 1.1 2.47e+06 1.2 7.0e+03 1.9e+02 0.0e+00  0  3 23 20  0   0  3 23 20  0  4937
MatSolve             372 1.0 4.9938e-02 1.1 3.53e+07 1.2 0.0e+00 0.0e+00 0.0e+00  1 41  0  0  0   1 41  0  0  0  7740
MatSolveTranspos     165 1.0 2.1765e-02 1.3 1.27e+07 1.2 0.0e+00 0.0e+00 0.0e+00  0 15  0  0  0   0 15  0  0  0  6382
MatLUFactorNum         8 1.0 1.3646e-02 1.2 1.06e+07 1.2 0.0e+00 0.0e+00 0.0e+00  0 12  0  0  0   0 12  0  0  0  8481
MatILUFactorSym        8 1.0 2.8668e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 2.4e+01  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin      20 1.0 1.7528e-0220.4 0.00e+00 0.0 7.4e+02 9.7e+02 4.0e+01  0  0  2 11  2   0  0  2 11  2     0
MatAssemblyEnd        20 1.0 1.8444e-03 1.4 0.00e+00 0.0 3.5e+02 4.4e+01 3.2e+01  0  0  1  0  2   0  0  1  0  2     0
MatGetRowIJ            8 1.0 2.4796e-0513.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         8 1.0 4.3297e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        20 1.0 2.7728e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver    16             16       154880     0
      Preconditioner    16             16        14272     0
              Vector   361            361      1846624     0
      Vector Scatter    14             14        14504     0
           Index Set    50             50        52696     0
   IS L to G Mapping     4              4         2256     0
              Matrix    20             20      3285596     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.00679e-06
Average time for zero size MPI_Send(): 1.37488e-05
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
| Time:           Thu Jan 31 22:02:38 2013                                                                             |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=6.19578, Active time=5.99688                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     4         0.0429      0.010731    0.0636      0.015904    0.72     1.06     |
|   build_sparsity()                 4         0.0251      0.006282    0.0737      0.018417    0.42     1.23     |
|   create_dof_constraints()         4         0.0199      0.004971    0.0199      0.004971    0.33     0.33     |
|   distribute_dofs()                4         0.2044      0.051109    0.5671      0.141771    3.41     9.46     |
|   dof_indices()                    11651     1.3774      0.000118    1.3774      0.000118    22.97    22.97    |
|   old_dof_indices()                1376      0.1605      0.000117    0.1605      0.000117    2.68     2.68     |
|   prepare_send_list()              4         0.0004      0.000105    0.0004      0.000105    0.01     0.01     |
|   reinit()                         4         0.3374      0.084355    0.3374      0.084355    5.63     5.63     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0146      0.001821    0.1033      0.012917    0.24     1.72     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        10432     1.4099      0.000135    1.4099      0.000135    23.51    23.51    |
|   init_shape_functions()           800       0.0535      0.000067    0.0535      0.000067    0.89     0.89     |
|   inverse_map()                    5062      0.0803      0.000016    0.0803      0.000016    1.34     1.34     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assemble_qoi()                   20        0.0646      0.003231    1.1959      0.059793    1.08     19.94    |
|   assemble_qoi_derivative()        4         0.0340      0.008505    0.2591      0.064775    0.57     4.32     |
|   assembly()                       8         0.1257      0.015716    0.5163      0.064531    2.10     8.61     |
|   assembly(get_jacobian)           4         0.0611      0.015263    0.2549      0.063728    1.02     4.25     |
|   assembly(get_residual)           20        0.1321      0.006604    1.1109      0.055546    2.20     18.52    |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             10432     0.1836      0.000018    0.1836      0.000018    3.06     3.06     |
|   compute_face_map()               688       0.0228      0.000033    0.0568      0.000083    0.38     0.95     |
|   init_face_shape_functions()      64        0.0009      0.000014    0.0009      0.000014    0.01     0.01     |
|   init_reference_to_physical_map() 800       0.0446      0.000056    0.0446      0.000056    0.74     0.74     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               8         0.2316      0.028947    0.2325      0.029064    3.86     3.88     |
|                                                                                                                |
| ImplicitSystem                                                                                                 |
|   adjoint_solve()                  4         0.0006      0.000139    0.5761      0.144034    0.01     9.61     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           20460     0.1295      0.000006    0.1295      0.000006    2.16     2.16     |
|   init()                           8         0.0367      0.004589    0.0367      0.004589    0.61     0.61     |
|                                                                                                                |
| Mesh                                                                                                           |
|   all_second_order()               1         0.0004      0.000424    0.0004      0.000424    0.01     0.01     |
|   contract()                       3         0.0053      0.001757    0.0286      0.009523    0.09     0.48     |
|   find_neighbors()                 6         0.1347      0.022445    0.1367      0.022791    2.25     2.28     |
|   renumber_nodes_and_elem()        3         0.0233      0.007766    0.0233      0.007766    0.39     0.39     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0011      0.001073    0.0016      0.001603    0.02     0.03     |
|   compute_hilbert_indices()        7         0.0335      0.004784    0.0335      0.004784    0.56     0.56     |
|   find_global_indices()            7         0.0148      0.002120    0.0603      0.008615    0.25     1.01     |
|   parallel_sort()                  7         0.0055      0.000790    0.0075      0.001068    0.09     0.12     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0005      0.000059    0.3373      0.042161    0.01     5.62     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              3         0.0017      0.000560    0.0017      0.000582    0.03     0.03     |
|   _refine_elements()               8         0.1577      0.019711    0.4237      0.052965    2.63     7.07     |
|   add_point()                      20460     0.1185      0.000006    0.2524      0.000012    1.98     4.21     |
|   make_coarsening_compatible()     6         0.0401      0.006678    0.0401      0.006678    0.67     0.67     |
|   make_refinement_compatible()     6         0.0000      0.000004    0.0001      0.000012    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      6         0.2496      0.041597    0.3094      0.051562    4.16     5.16     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          4         0.0448      0.011202    0.6417      0.160432    0.75     10.70    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      28        0.0105      0.000375    0.0106      0.000380    0.17     0.18     |
|   broadcast()                      9         0.0004      0.000045    0.0003      0.000036    0.01     0.01     |
|   max(bool)                        19        0.0105      0.000552    0.0105      0.000552    0.17     0.17     |
|   max(scalar)                      865       0.0061      0.000007    0.0061      0.000007    0.10     0.10     |
|   max(vector)                      201       0.0031      0.000015    0.0072      0.000036    0.05     0.12     |
|   min(bool)                        1056      0.0070      0.000007    0.0070      0.000007    0.12     0.12     |
|   min(scalar)                      850       0.0438      0.000052    0.0438      0.000052    0.73     0.73     |
|   min(vector)                      201       0.0028      0.000014    0.0073      0.000036    0.05     0.12     |
|   probe()                          550       0.0121      0.000022    0.0121      0.000022    0.20     0.20     |
|   receive()                        550       0.0035      0.000006    0.0157      0.000028    0.06     0.26     |
|   send()                           550       0.0019      0.000003    0.0019      0.000003    0.03     0.03     |
|   send_receive()                   564       0.0043      0.000008    0.0233      0.000041    0.07     0.39     |
|   sum()                            105       0.0044      0.000042    0.0162      0.000154    0.07     0.27     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           550       0.0011      0.000002    0.0011      0.000002    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         6         0.0282      0.004700    0.0373      0.006212    0.47     0.62     |
|   set_parent_processor_ids()       6         0.0129      0.002146    0.0129      0.002146    0.21     0.21     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          8         0.1730      0.021631    0.1730      0.021631    2.89     2.89     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       6         0.0279      0.004642    0.2370      0.039492    0.46     3.95     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 6         0.0178      0.002967    0.3366      0.056106    0.30     5.61     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            88539     5.9969                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adjoints_ex2:
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
