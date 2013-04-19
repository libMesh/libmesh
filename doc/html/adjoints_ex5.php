<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adjoints_ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file adjoint_initial.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/parameters.h"
        #include "libmesh/point.h"
        #include "libmesh/vector_value.h"
        
        using namespace libMesh;
        
        void adjoint_read_initial_parameters();
        void adjoint_finish_initialization();
        
        Number adjoint_initial_value(const Point& p,
                             const Parameters&,
                             const std::string&,
                             const std::string&);
        
        Gradient adjoint_initial_grad(const Point& p,
                              const Parameters&,
                              const std::string&,
                              const std::string&);
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file femparameters.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __fem_parameters_h__
        #define __fem_parameters_h__
        
        #include "libmesh/libmesh_common.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/enum_norm_type.h"
        #include "libmesh/function_base.h"
        #include "libmesh/getpot.h"
        #include "libmesh/id_types.h"
        #include "libmesh/periodic_boundaries.h"
        
        #include &lt;limits&gt;
        #include &lt;map&gt;
        #include &lt;string&gt;
        #include &lt;vector&gt;
        
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
            FEMParameters();
        
            ~FEMParameters();
        
            void read(GetPot &input);
        
</pre>
</div>
<div class = "comment">
Parameters applicable to entire EquationSystems:


<br><br></div>

<div class ="fragment">
<pre>
            unsigned int initial_timestep, n_timesteps;
            bool transient;
            unsigned int deltat_reductions;
            std::string timesolver_core;
            Real end_time, deltat, timesolver_theta,
        	 timesolver_maxgrowth, timesolver_tolerance,
        	 timesolver_upper_tolerance, steadystate_tolerance;
            std::vector&lt;FEMNormType&gt; timesolver_norm;
        
</pre>
</div>
<div class = "comment">
Mesh generation


<br><br></div>

<div class ="fragment">
<pre>
            unsigned int dimension;
            std::string domaintype, domainfile, elementtype;
            Real elementorder;
            Real domain_xmin, domain_ymin, domain_zmin;
            Real domain_edge_width, domain_edge_length, domain_edge_height;
            unsigned int coarsegridx, coarsegridy, coarsegridz;
            unsigned int coarserefinements, extrarefinements;
        
</pre>
</div>
<div class = "comment">
Mesh refinement


<br><br></div>

<div class ="fragment">
<pre>
            unsigned int nelem_target;
            Real global_tolerance;
            Real refine_fraction, coarsen_fraction, coarsen_threshold;
            unsigned int max_adaptivesteps;
            unsigned int initial_adaptivesteps;
        
</pre>
</div>
<div class = "comment">
Output


<br><br></div>

<div class ="fragment">
<pre>
            unsigned int write_interval;
            bool write_gmv_error, write_tecplot_error, output_xda, output_xdr,
        	 output_bz2, output_gz, output_gmv, output_tecplot;
        
</pre>
</div>
<div class = "comment">
Types of Systems to create


<br><br></div>

<div class ="fragment">
<pre>
            std::vector&lt;std::string&gt; system_types;
        
</pre>
</div>
<div class = "comment">
Parameters applicable to each system:


<br><br>Boundary and initial conditions


<br><br>std::vector<PeriodicBoundary> periodic_boundaries;


<br><br></div>

<div class ="fragment">
<pre>
            std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt; initial_conditions;
            std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;  dirichlet_conditions,
                                                                neumann_conditions;
            std::map&lt;boundary_id_type, std::vector&lt;unsigned int&gt; &gt; dirichlet_condition_variables,
                                                                   neumann_condition_variables;
            std::map&lt;int, std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt; &gt;
              other_interior_functions;
            std::map&lt;int, std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt; &gt;
              other_boundary_functions;
        
</pre>
</div>
<div class = "comment">
Execution type


<br><br></div>

<div class ="fragment">
<pre>
            bool run_simulation, run_postprocess;
        
</pre>
</div>
<div class = "comment">
Approximation type


<br><br></div>

<div class ="fragment">
<pre>
            std::vector&lt;std::string&gt; fe_family;
            std::vector&lt;unsigned int&gt; fe_order;
            int extra_quadrature_order;
        
</pre>
</div>
<div class = "comment">
Assembly options


<br><br></div>

<div class ="fragment">
<pre>
            bool analytic_jacobians;
            Real verify_analytic_jacobians;
            Real numerical_jacobian_h;
        
            bool print_solution_norms, print_solutions,
                 print_residual_norms, print_residuals,
                 print_jacobian_norms, print_jacobians,
                 print_element_jacobians;
        
</pre>
</div>
<div class = "comment">
Solver options


<br><br></div>

<div class ="fragment">
<pre>
            bool use_petsc_snes;
            bool time_solver_quiet, solver_quiet, solver_verbose, require_residual_reduction;
            Real min_step_length;
            unsigned int max_linear_iterations, max_nonlinear_iterations;
            Real relative_step_tolerance, relative_residual_tolerance,
        	 initial_linear_tolerance, minimum_linear_tolerance,
        	 linear_tolerance_multiplier;
        
</pre>
</div>
<div class = "comment">
Initialization


<br><br></div>

<div class ="fragment">
<pre>
            unsigned int initial_sobolev_order;
            unsigned int initial_extra_quadrature;
        
</pre>
</div>
<div class = "comment">
Error indicators


<br><br></div>

<div class ="fragment">
<pre>
            std::string indicator_type;
            unsigned int sobolev_order;
        
</pre>
</div>
<div class = "comment">
System-specific parameters file


<br><br></div>

<div class ="fragment">
<pre>
            std::string system_config_file;
        };
        
        #endif // __fem_parameters_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file heatsystem.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/enum_fe_family.h"
        #include "libmesh/fem_system.h"
        #include "libmesh/parameter_vector.h"
        
</pre>
</div>
<div class = "comment">
FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
but we must specify element residuals
</div>

<div class ="fragment">
<pre>
        class HeatSystem : public FEMSystem
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          HeatSystem(EquationSystems& es,
                       const std::string& name_in,
                       const unsigned int number_in)
          : FEMSystem(es, name_in, number_in),
            _k(1.0),
            _fe_family("LAGRANGE"), _fe_order(1),
            _analytic_jacobians(true), R_plus_dp(0.0), R_minus_dp(0.0), dp(1.e-6) { qoi.resize(1); }
        
          std::string & fe_family() { return _fe_family;  }
          unsigned int & fe_order() { return _fe_order;  }
          Real & k() { return _k;  }
          bool & analytic_jacobians() { return _analytic_jacobians; }
        
</pre>
</div>
<div class = "comment">
A function to compute and accumulate residuals
</div>

<div class ="fragment">
<pre>
          void perturb_accumulate_residuals(const ParameterVector& parameters);
        
</pre>
</div>
<div class = "comment">
Sensitivity Calculation
</div>

<div class ="fragment">
<pre>
          Number & compute_final_sensitivity()
            {
              final_sensitivity = -(R_plus_dp - R_minus_dp)/(2*dp);
        
              return final_sensitivity;
            }
        
            void set_tf(Real val)
          {
            tf = val;
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
        
          Number &get_QoI_value(unsigned int QoI_index)
            {
              return computed_QoI[QoI_index];
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
virtual bool side_constraint (bool request_jacobian,
DiffContext &context);


<br><br>RHS for adjoint problem
</div>

<div class ="fragment">
<pre>
          virtual void element_qoi_derivative (DiffContext &context,
        				       const QoISet & /* qois */);
        
</pre>
</div>
<div class = "comment">
virtual void element_qoi (DiffContext &context, const QoISet & qois);


<br><br>Parameters associated with the system
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
The parameters to solve for
</div>

<div class ="fragment">
<pre>
          Real _k;
        
</pre>
</div>
<div class = "comment">
The final time parameter
</div>

<div class ="fragment">
<pre>
          Real tf;
        
</pre>
</div>
<div class = "comment">
Variables to hold the computed QoIs
</div>

<div class ="fragment">
<pre>
          Number computed_QoI[1];
        
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
Index for T variable
</div>

<div class ="fragment">
<pre>
          unsigned int T_var;
        
</pre>
</div>
<div class = "comment">
Calculate Jacobians analytically or not?
</div>

<div class ="fragment">
<pre>
          bool _analytic_jacobians;
        
</pre>
</div>
<div class = "comment">
Variables to hold the perturbed residuals
</div>

<div class ="fragment">
<pre>
          Number R_plus_dp;
          Number R_minus_dp;
        
</pre>
</div>
<div class = "comment">
Perturbation parameter
</div>

<div class ="fragment">
<pre>
          Real dp;
        
</pre>
</div>
<div class = "comment">
The final computed sensitivity
</div>

<div class ="fragment">
<pre>
          Number final_sensitivity;
        
        };
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file initial.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/parameters.h"
        #include "libmesh/point.h"
        #include "libmesh/vector_value.h"
        
        using namespace libMesh;
        
        void read_initial_parameters();
        void finish_initialization();
        
        Number initial_value(const Point& p,
                             const Parameters&,
                             const std::string&,
                             const std::string&);
        
        Gradient initial_grad(const Point& p,
                              const Parameters&,
                              const std::string&,
                              const std::string&);
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file adjoint_initial.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "adjoint_initial.h"
        
        using namespace libMesh;
        
        void adjoint_read_initial_parameters()
        {
        }
        
        void adjoint_finish_initialization()
        {
        }
        
        
        
</pre>
</div>
<div class = "comment">
Initial conditions
</div>

<div class ="fragment">
<pre>
        Number adjoint_initial_value(const Point& p,
                             const Parameters&,
                             const std::string&,
                             const std::string&)
        {
          Real x = p(0), y = p(1);
        
          Number val = 0.;
        
          val = sin(M_PI * x) * sin(M_PI * y);
        
          return val;
        }
        
        
        
        Gradient adjoint_initial_grad(const Point& p,
                              const Parameters&,
                              const std::string&,
                              const std::string&)
        {
          Real x = p(0), y = p(1);
        
          return Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
          M_PI*sin(M_PI * x) * cos(M_PI * y));
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file adjoints_ex5.C with comments: </h1> 
<div class = "comment">
partial(T)/partial(t) - K Laplacian(T) = 0
with initial condition:
T(x,y;0) = sin(pi*x) sin(pi*y) and boundary conditions:
T(boundary;t) = 0


<br><br>For these initial and boundary conditions, the exact solution
u = exp(-K pi^2 t) * sin(pi*x) * sin(pi*y)


<br><br>We specify our Quantity of Interest (QoI) as
Q(u) = int_{domain} u(x,y;1) sin(pi*x) sin(pi*y) dx dy, and
are interested in computing the sensitivity dQ/dK


<br><br>For this QoI, the continuous adjoint problem reads,
-partial(z)/partial(t) - K Laplacian(z) = 0
with initial condition:
T(x,y;1) = sin(pi*x) sin(pi*y)
and boundary condition:
T(boundary;t) = 0


<br><br>which has the exact solution,
z = exp(-K pi^2 (1 - t)) * sin(pi*x) * sin(pi*y)
which is the mirror image in time of the forward solution


<br><br>For an adjoint consistent space-time formulation, the discrete
adjoint can be obtained by marching backwards from the adjoint
initial condition and solving the transpose of the discrete primal
problem at the last nonlinear solve of the corresponding primal
timestep. This necessitates the storage of the primal solution at
all timesteps, which is accomplished here using a
MemorySolutionHistory object. As the name suggests, this object
simply stores the primal solution (and other vectors we may choose
to save), so that we can retrieve them later, whenever necessary.


<br><br>The discrete adjoint system for implicit time steppers requires the
localization of vectors other than system.solution, which is
accomplished using the localize_vectors method. In this particular
example, we use the localized adjoint solution to assemble the
residual contribution for the current adjoint timestep from the last
computed adjoint timestep.


<br><br>Finally, The adjoint_advance_timestep method, the backwards time
analog of advance_timestep prepares the time solver for solving the
adjoint system, while the retrieve_timestep method retrieves the
saved solutions at the current system.time, so that the adjoint
sensitivity contribution for the current time can be computed.


<br><br>Local includes
</div>

<div class ="fragment">
<pre>
        #include "initial.h"
        #include "adjoint_initial.h"
        #include "femparameters.h"
        #include "heatsystem.h"
        
</pre>
</div>
<div class = "comment">
Libmesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/system_norm.h"
        #include "libmesh/numeric_vector.h"
        
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/mesh_refinement.h"
        
        #include "libmesh/petsc_diff_solver.h"
        #include "libmesh/steady_solver.h"
        #include "libmesh/euler_solver.h"
        #include "libmesh/euler2_solver.h"
        #include "libmesh/newton_solver.h"
        #include "libmesh/twostep_time_solver.h"
        
        #include "libmesh/getpot.h"
        #include "libmesh/tecplot_io.h"
        #include "libmesh/gmv_io.h"
        
</pre>
</div>
<div class = "comment">
SolutionHistory Includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/solution_history.h"
        #include "libmesh/memory_solution_history.h"
        
</pre>
</div>
<div class = "comment">
C++ includes
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;sys/time.h&gt;
        #include &lt;iomanip&gt;
        
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
        
        void set_system_parameters(HeatSystem &system, FEMParameters &param)
        {
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
          system.numerical_jacobian_h = param.numerical_jacobian_h;
        
</pre>
</div>
<div class = "comment">
More desperate debugging options
</div>

<div class ="fragment">
<pre>
          system.print_solution_norms    = param.print_solution_norms;
          system.print_solutions         = param.print_solutions;
          system.print_residual_norms    = param.print_residual_norms;
          system.print_residuals         = param.print_residuals;
          system.print_jacobian_norms    = param.print_jacobian_norms;
          system.print_jacobians         = param.print_jacobians;
          system.print_element_jacobians = param.print_element_jacobians;
        
</pre>
</div>
<div class = "comment">
Solve this as a time-dependent or steady system
</div>

<div class ="fragment">
<pre>
          if (param.transient)
            {
              UnsteadySolver *innersolver;
               if (param.timesolver_core == "euler")
                {
                  EulerSolver *eulersolver =
                    new EulerSolver(system);
        
                  eulersolver-&gt;theta = param.timesolver_theta;
                  innersolver = eulersolver;
                }
              else
                {
                  std::cerr &lt;&lt; "This example (and unsteady adjoints in libMesh) only support Backward Euler and explicit methods."
                            &lt;&lt;  std::endl;
                  libmesh_error();
                }
               system.time_solver =
                  AutoPtr&lt;TimeSolver&gt;(innersolver);
            }
          else
            system.time_solver =
              AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
</pre>
</div>
<div class = "comment">
The Memory Solution History object we will set the system SolutionHistory object to
</div>

<div class ="fragment">
<pre>
          MemorySolutionHistory heatsystem_solution_history(system);
          system.time_solver-&gt;set_solution_history(heatsystem_solution_history);
        
          system.time_solver-&gt;reduce_deltat_on_diffsolver_failure =
                                                param.deltat_reductions;
          system.time_solver-&gt;quiet           = param.time_solver_quiet;
        
</pre>
</div>
<div class = "comment">
Create any Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
          typedef
            std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::
            const_iterator Iter;
        
          for (Iter i = param.dirichlet_conditions.begin();
               i != param.dirichlet_conditions.end(); ++i)
            {
              boundary_id_type b = i-&gt;first;
              FunctionBase&lt;Number&gt; *f = i-&gt;second;
              std::set&lt;boundary_id_type&gt; bdys; bdys.insert(b);
              system.get_dof_map().add_dirichlet_boundary
                (DirichletBoundary
                   (bdys, param.dirichlet_condition_variables[b], f));
              std::cout &lt;&lt; "Added Dirichlet boundary " &lt;&lt; b &lt;&lt; " for variables ";
              for (unsigned int vi=0; vi !=
                   param.dirichlet_condition_variables[b].size(); ++vi)
                std::cout &lt;&lt; param.dirichlet_condition_variables[b][vi];
              std::cout &lt;&lt; std::endl;
            }
        
</pre>
</div>
<div class = "comment">
Set the time stepping options
</div>

<div class ="fragment">
<pre>
          system.deltat = param.deltat;
        
</pre>
</div>
<div class = "comment">
And the integration options
</div>

<div class ="fragment">
<pre>
          system.extra_quadrature_order = param.extra_quadrature_order;
        
</pre>
</div>
<div class = "comment">
And the nonlinear solver options
</div>

<div class ="fragment">
<pre>
          if (param.use_petsc_snes)
            {
        #ifdef LIBMESH_HAVE_PETSC
              PetscDiffSolver *solver = new PetscDiffSolver(system);
              system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
        #else
              libmesh_error();
        #endif
            }
          else
            {
              NewtonSolver *solver = new NewtonSolver(system);
              system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
        
              solver-&gt;quiet                       = param.solver_quiet;
              solver-&gt;verbose                     = param.solver_verbose;
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
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
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
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
        
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
Create a mesh with the given dimension, distributed
across the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm(), param.dimension);
        
</pre>
</div>
<div class = "comment">
And an object to refine it
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;MeshRefinement&gt; mesh_refinement(new MeshRefinement(mesh));
        
</pre>
</div>
<div class = "comment">
And an EquationSystems to run on it
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
          std::cout &lt;&lt; "Building mesh" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Build a unit square
</div>

<div class ="fragment">
<pre>
          ElemType elemtype;
        
          if (param.elementtype == "tri" ||
              param.elementtype == "unstructured")
            elemtype = TRI3;
          else
            elemtype = QUAD4;
        
          MeshTools::Generation::build_square
            (mesh, param.coarsegridx, param.coarsegridy,
             param.domain_xmin, param.domain_xmin + param.domain_edge_width,
             param.domain_ymin, param.domain_ymin + param.domain_edge_length,
             elemtype);
        
          std::cout &lt;&lt; "Building system" &lt;&lt; std::endl;
        
          HeatSystem &system = equation_systems.add_system&lt;HeatSystem&gt; ("HeatSystem");
        
          set_system_parameters(system, param);
        
          std::cout &lt;&lt; "Initializing systems" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Initialize the system
</div>

<div class ="fragment">
<pre>
          equation_systems.init ();
        
</pre>
</div>
<div class = "comment">
Refine the grid again if requested
</div>

<div class ="fragment">
<pre>
          for (unsigned int i=0; i != param.extrarefinements; ++i)
            {
              mesh_refinement-&gt;uniformly_refine(1);
              equation_systems.reinit();
            }
        
          std::cout&lt;&lt;"Setting primal initial conditions"&lt;&lt;std::endl;
        
          read_initial_parameters();
        
          system.project_solution(initial_value, initial_grad,
                                          equation_systems.parameters);
        
</pre>
</div>
<div class = "comment">
Output the H1 norm of the initial conditions
</div>

<div class ="fragment">
<pre>
          libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
Add an adjoint vector, this will be computed after the forward
time stepping is complete

<br><br>Tell the library not to save adjoint solutions during the forward
solve

<br><br>Tell the library not to project this vector, and hence, memory
solution history to not save it.

<br><br>Make this vector ghosted so we can localize it to each element
later.
</div>

<div class ="fragment">
<pre>
          const std::string & adjoint_solution_name = "adjoint_solution0";
          system.add_vector("adjoint_solution0", false, GHOSTED);
        
</pre>
</div>
<div class = "comment">
Close up any resources initial.C needed
</div>

<div class ="fragment">
<pre>
          finish_initialization();
        
</pre>
</div>
<div class = "comment">
Plot the initial conditions
</div>

<div class ="fragment">
<pre>
          write_output(equation_systems, 0, "primal");
        
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
In optimized mode we catch any solver errors, so that we can
write the proper footers before closing.  In debug mode we just
let the exception throw so that gdb can grab it.
</div>

<div class ="fragment">
<pre>
        #ifdef NDEBUG
          try
          {
        #endif
</pre>
</div>
<div class = "comment">
Now we begin the timestep loop to compute the time-accurate
solution of the equations.
</div>

<div class ="fragment">
<pre>
            for (unsigned int t_step=param.initial_timestep;
                 t_step != param.initial_timestep + param.n_timesteps; ++t_step)
              {
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; " Solving time step " &lt;&lt; t_step &lt;&lt; ", time = "
                          &lt;&lt; system.time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Solve the forward problem at time t, to obtain the solution at time t + dt
</div>

<div class ="fragment">
<pre>
                system.solve();
        
</pre>
</div>
<div class = "comment">
Output the H1 norm of the computed solution
</div>

<div class ="fragment">
<pre>
                libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time + system.deltat&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Advance to the next timestep in a transient problem
</div>

<div class ="fragment">
<pre>
                std::cout&lt;&lt;"Advancing timestep"&lt;&lt;std::endl&lt;&lt;std::endl;
                system.time_solver-&gt;advance_timestep();
        
</pre>
</div>
<div class = "comment">
Write out this timestep
</div>

<div class ="fragment">
<pre>
                  write_output(equation_systems, t_step+1, "primal");
              }
</pre>
</div>
<div class = "comment">
End timestep loop


<br><br>/////////////// Now for the Adjoint Solution //////////////////////////////////////


<br><br>Now we will solve the backwards in time adjoint problem
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; std::endl &lt;&lt; "Solving the adjoint problem" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
We need to tell the library that it needs to project the adjoint, so
MemorySolutionHistory knows it has to save it


<br><br>Tell the library to project the adjoint vector, and hence, memory solution history to
save it
</div>

<div class ="fragment">
<pre>
          system.set_vector_preservation(adjoint_solution_name, true);
        
          std::cout &lt;&lt; "Setting adjoint initial conditions Z("&lt;&lt;system.time&lt;&lt;")"&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
Need to call adjoint_advance_timestep once for the initial condition setup
</div>

<div class ="fragment">
<pre>
          std::cout&lt;&lt;"Retrieving solutions at time t="&lt;&lt;system.time&lt;&lt;std::endl;
          system.time_solver-&gt;adjoint_advance_timestep();
        
</pre>
</div>
<div class = "comment">
Output the H1 norm of the retrieved solutions (u^i and u^i+1)
</div>

<div class ="fragment">
<pre>
          libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time + system.deltat&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
        
          libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
The first thing we have to do is to apply the adjoint initial
condition. The user should supply these. Here they are specified
in the functions adjoint_initial_value and adjoint_initial_gradient
</div>

<div class ="fragment">
<pre>
          system.project_vector(adjoint_initial_value, adjoint_initial_grad, equation_systems.parameters, system.get_adjoint_solution(0));
        
</pre>
</div>
<div class = "comment">
Since we have specified an adjoint solution for the current time (T), set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the adjoint sensitivity method
</div>

<div class ="fragment">
<pre>
          system.set_adjoint_already_solved(true);
        
          libMesh::out &lt;&lt; "|Z(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(system.get_adjoint_solution(), 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
        
          write_output(equation_systems, param.n_timesteps, "dual");
        
</pre>
</div>
<div class = "comment">
Now that the adjoint initial condition is set, we will start the
backwards in time adjoint integration


<br><br>For loop stepping backwards in time
</div>

<div class ="fragment">
<pre>
          for (unsigned int t_step=param.initial_timestep;
               t_step != param.initial_timestep + param.n_timesteps; ++t_step)
            {
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; " Solving adjoint time step " &lt;&lt; t_step &lt;&lt; ", time = "
        		&lt;&lt; system.time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
The adjoint_advance_timestep
function calls the retrieve function of the memory_solution_history
class via the memory_solution_history object we declared earlier.
The retrieve function sets the system primal vectors to their values
at the current timestep
</div>

<div class ="fragment">
<pre>
              std::cout&lt;&lt;"Retrieving solutions at time t="&lt;&lt;system.time&lt;&lt;std::endl;
              system.time_solver-&gt;adjoint_advance_timestep();
        
</pre>
</div>
<div class = "comment">
Output the H1 norm of the retrieved solution
</div>

<div class ="fragment">
<pre>
              libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time + system.deltat &lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
        
              libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) &lt;&lt; std::endl;
        
              system.set_adjoint_already_solved(false);
        
              system.adjoint_solve();
        
</pre>
</div>
<div class = "comment">
Now that we have solved the adjoint, set the adjoint_already_solved boolean to true, so we dont solve unneccesarily in the error estimator
</div>

<div class ="fragment">
<pre>
              system.set_adjoint_already_solved(true);
        
              libMesh::out &lt;&lt; "|Z(" &lt;&lt;system.time&lt;&lt; ")|= "&lt;&lt; system.calculate_norm(system.get_adjoint_solution(), 0, H1) &lt;&lt; std::endl &lt;&lt; std::endl;
        
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
        
              write_output(equation_systems, param.n_timesteps - (t_step + 1), "dual");
        
</pre>
</div>
<div class = "comment">
Swap back
</div>

<div class ="fragment">
<pre>
              primal_solution.swap(dual_solution_0);
            }
</pre>
</div>
<div class = "comment">
End adjoint timestep loop


<br><br>Now that we have computed both the primal and adjoint solutions, we compute the sensitivties to the parameter p
dQ/dp = partialQ/partialp - partialR/partialp
partialQ/partialp = (Q(p+dp) - Q(p-dp))/(2*dp), this is not supported by the library yet
partialR/partialp = (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp), where
R(u,z;p+dp) = int_{0}^{T} f(z;p+dp) - <partialu/partialt, z>(p+dp) - <g(u),z>(p+dp)
To do this we need to step forward in time, and compute the perturbed R at each time step and accumulate it
Then once all time steps are over, we can compute (R(u,z;p+dp) - R(u,z;p-dp))/(2*dp)


<br><br>Now we begin the timestep loop to compute the time-accurate
adjoint sensitivities
</div>

<div class ="fragment">
<pre>
            for (unsigned int t_step=param.initial_timestep;
                 t_step != param.initial_timestep + param.n_timesteps; ++t_step)
              {
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; "Retrieving " &lt;&lt; t_step &lt;&lt; ", time = "
                          &lt;&lt; system.time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Retrieve the primal and adjoint solutions at the current timestep
</div>

<div class ="fragment">
<pre>
                system.time_solver-&gt;retrieve_timestep();
        
        	libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time + system.deltat &lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
        
              libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) &lt;&lt; std::endl;
        
              libMesh::out &lt;&lt; "|Z(" &lt;&lt;system.time&lt;&lt; ")|= "&lt;&lt; system.calculate_norm(system.get_adjoint_solution(0), 0, H1) &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Call the postprocess function which we have overloaded to compute
accumulate the perturbed residuals
</div>

<div class ="fragment">
<pre>
                (dynamic_cast&lt;HeatSystem&&gt;(system)).perturb_accumulate_residuals(dynamic_cast&lt;HeatSystem&&gt;(system).get_parameter_vector());
        
</pre>
</div>
<div class = "comment">
Move the system time forward (retrieve_timestep does not do this)
</div>

<div class ="fragment">
<pre>
                system.time += system.deltat;
              }
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; "Retrieving " &lt;&lt; " final time = "
        	      &lt;&lt; system.time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Retrieve the primal and adjoint solutions at the current timestep
</div>

<div class ="fragment">
<pre>
            system.time_solver-&gt;retrieve_timestep();
        
            libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time + system.deltat &lt;&lt; ")|= " &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
        
              libMesh::out &lt;&lt; "|U(" &lt;&lt;system.time&lt;&lt; ")|= " &lt;&lt; system.calculate_norm(system.get_vector("_old_nonlinear_solution"), 0, H1) &lt;&lt; std::endl;
        
              libMesh::out &lt;&lt; "|Z(" &lt;&lt;system.time&lt;&lt; ")|= "&lt;&lt; system.calculate_norm(system.get_adjoint_solution(0), 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
Call the postprocess function which we have overloaded to compute
accumulate the perturbed residuals
</div>

<div class ="fragment">
<pre>
            (dynamic_cast&lt;HeatSystem&&gt;(system)).perturb_accumulate_residuals(dynamic_cast&lt;HeatSystem&&gt;(system).get_parameter_vector());
        
</pre>
</div>
<div class = "comment">
Now that we computed the accumulated, perturbed residuals, we can compute the
approximate sensitivity
</div>

<div class ="fragment">
<pre>
            Number sensitivity_0_0 = (dynamic_cast&lt;HeatSystem&&gt;(system)).compute_final_sensitivity();
        
</pre>
</div>
<div class = "comment">
Print it out
</div>

<div class ="fragment">
<pre>
            std::cout&lt;&lt;"Sensitivity of QoI 0 w.r.t parameter 0 is: " &lt;&lt; sensitivity_0_0 &lt;&lt; std::endl;
        
        #ifdef NDEBUG
        }
          catch (...)
          {
            std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id()
                      &lt;&lt; "] Caught exception; exiting early." &lt;&lt; std::endl;
          }
        #endif
        
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
        
        #endif // LIBMESH_ENABLE_AMR
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
        #include "heatsystem.h"
        
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
        void HeatSystem::element_qoi_derivative (DiffContext &context,
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
A reference to the system context is built with
</div>

<div class ="fragment">
<pre>
          const System & sys = c.get_system();
        
</pre>
</div>
<div class = "comment">
Get a pointer to the adjoint solution vector
</div>

<div class ="fragment">
<pre>
          NumericVector&lt;Number&gt; &adjoint_solution = const_cast&lt;System &&gt;(sys).get_adjoint_solution(0);
        
</pre>
</div>
<div class = "comment">
Get the previous adjoint solution values at all the qps


<br><br></div>

<div class ="fragment">
<pre>
          std::vector&lt;Number&gt; old_adjoint (n_qpoints, 0);
        
          c.interior_values&lt;Number&gt;(0, adjoint_solution, old_adjoint);
        
</pre>
</div>
<div class = "comment">
Our QoI depends solely on the final time, so there are no QoI contributions.
However, there is a contribution from the adjoint solution timestep, for the
time part of the residual of the adjoint problem
Loop over the qps
</div>

<div class ="fragment">
<pre>
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              for (unsigned int i=0; i != n_T_dofs; i++)
        	{
        	  Q(i) += -JxW[qp] * old_adjoint[qp] * phi[i][qp] ;
              	}
        
            } // end of the quadrature point qp-loop
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file factoryfunction.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dense_vector.h"
        #include "libmesh/factory.h"
        #include "libmesh/function_base.h"
        
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
An example of a hand-coded function


<br><br></div>

<div class ="fragment">
<pre>
        class ExampleOneFunction : public FunctionBase&lt;Number&gt;
        {
          virtual Number operator() (const Point&  /*p*/,
                                     const Real /*time*/)
            {
              return 1;
            }
        
          virtual void operator() (const Point&  /*p*/,
                                   const Real /*time*/,
                                   DenseVector&lt;Number&gt;& output)
            {
              for (unsigned int i=0; i != output.size(); ++i)
                output(i) = 1;
            }
        
          virtual void init() {}
          virtual void clear() {}
          virtual AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() const {
            return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
              (new ExampleOneFunction());
          }
        };
        
        #ifdef LIBMESH_USE_SEPARATE_NAMESPACE
        namespace libMesh {
        #endif
        
</pre>
</div>
<div class = "comment">
-------------------------------------------------
Full specialization for the Factory<FunctionBase<Number> >
So we can look up hand-coded functions by name string
</div>

<div class ="fragment">
<pre>
        template&lt;&gt;
        std::map&lt;std::string, Factory&lt;FunctionBase&lt;Number&gt; &gt;*&gt;&
        Factory&lt;FunctionBase&lt;Number&gt; &gt;::factory_map()
        {
          static std::map&lt;std::string, Factory&lt;FunctionBase&lt;Number&gt; &gt;*&gt; _map;
          return _map;
        }
        
        FactoryImp&lt;ExampleOneFunction, FunctionBase&lt;Number&gt; &gt; example_one_factory ("example_one");
        
        #ifdef LIBMESH_USE_SEPARATE_NAMESPACE
        } // namespace libMesh
        #endif
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file femparameters.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "femparameters.h"
        
        #include "libmesh/factory.h"
        #include "libmesh/parsed_function.h"
        #include "libmesh/zero_function.h"
        #include "libmesh/periodic_boundaries.h"
        #include "libmesh/vector_value.h" // RealVectorValue
        
        #define GETPOT_INPUT(A) { A = input(#A, A);\
          const std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); }
        #define GETPOT_INT_INPUT(A) { A = input(#A, (int)A);\
          const std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); }
        
        #define GETPOT_FUNCTION_INPUT(A) {\
          const std::string type  = input(#A "_type", "zero");\
          const std::string value = input(#A "_value", "");\
          A = new_function_base(type, value); }
        #define GETPOT_REGISTER(A) { \
          std::string stringval = input(#A, std::string());\
          variable_names.push_back(std::string(#A "=") + stringval); }
        
        FEMParameters::FEMParameters() :
              initial_timestep(0), n_timesteps(100),
              transient(true),
              deltat_reductions(0),
              timesolver_core("euler"),
              end_time(std::numeric_limits&lt;Real&gt;::max()),
              deltat(0.0001), timesolver_theta(0.5),
              timesolver_maxgrowth(0.), timesolver_tolerance(0.),
              timesolver_upper_tolerance(0.),
              steadystate_tolerance(0.),
              timesolver_norm(0,L2),
        
              dimension(2),
              domaintype("square"), domainfile("mesh.xda"), elementtype("quad"),
              elementorder(2),
              domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
              domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
              coarsegridx(1), coarsegridy(1), coarsegridz(1),
              coarserefinements(4), extrarefinements(0),
        
              nelem_target(8000), global_tolerance(0.0),
              refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
              max_adaptivesteps(1),
              initial_adaptivesteps(0),
        
              write_interval(10),
              write_gmv_error(false), write_tecplot_error(false),
              output_xda(false), output_xdr(false),
              output_bz2(true), output_gz(true),
              output_gmv(false), output_tecplot(false),
        
              system_types(0),
        
</pre>
</div>
<div class = "comment">
periodic_boundaries(0),


<br><br></div>

<div class ="fragment">
<pre>
              run_simulation(true), run_postprocess(false),
        
              fe_family(1, "LAGRANGE"), fe_order(1, 1),
              extra_quadrature_order(0),
        
              analytic_jacobians(true), verify_analytic_jacobians(0.0),
              numerical_jacobian_h(TOLERANCE),
              print_solution_norms(false), print_solutions(false),
              print_residual_norms(false), print_residuals(false),
              print_jacobian_norms(false), print_jacobians(false),
              print_element_jacobians(false),
        
              use_petsc_snes(false),
              time_solver_quiet(true), solver_quiet(true), solver_verbose(false),
              require_residual_reduction(true),
              min_step_length(1e-5),
              max_linear_iterations(200000), max_nonlinear_iterations(20),
              relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
              initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
              linear_tolerance_multiplier(1.e-3),
        
              initial_sobolev_order(1),
              initial_extra_quadrature(0),
              indicator_type("kelly"), sobolev_order(1)
        {
        }
        
        FEMParameters::~FEMParameters()
        {
          for (std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
               i = initial_conditions.begin(); i != initial_conditions.end();
               ++i)
            delete i-&gt;second;
        
          for (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
               i = dirichlet_conditions.begin(); i != dirichlet_conditions.end();
               ++i)
            delete i-&gt;second;
        
          for (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
               i = neumann_conditions.begin(); i != neumann_conditions.end();
               ++i)
            delete i-&gt;second;
        
          for (std::map&lt;int,
        	 std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;
               &gt;::iterator
               i = other_boundary_functions.begin(); i != other_boundary_functions.end();
               ++i)
            for (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
                 j = i-&gt;second.begin(); j != i-&gt;second.end();
                 ++j)
              delete j-&gt;second;
        
          for (std::map&lt;int,
        	 std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;
               &gt;::iterator
               i = other_interior_functions.begin(); i != other_interior_functions.end();
               ++i)
            for (std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
                 j = i-&gt;second.begin(); j != i-&gt;second.end();
                 ++j)
              delete j-&gt;second;
        }
        
        
        AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; new_function_base(const std::string& func_type,
                                                const std::string& func_value)
        {
          if (func_type == "parsed")
            return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
              (new ParsedFunction&lt;Number&gt;(func_value));
          else if (func_type == "factory")
            return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
              (Factory&lt;FunctionBase&lt;Number&gt; &gt;::build(func_value).release());
          else if (func_type == "zero")
            return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
              (new ZeroFunction&lt;Number&gt;);
          else
            libmesh_not_implemented();
        
          return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;();
        }
        
        
        void FEMParameters::read(GetPot &input)
        {
            std::vector&lt;std::string&gt; variable_names;
        
            GETPOT_INT_INPUT(initial_timestep);
            GETPOT_INT_INPUT(n_timesteps);
            GETPOT_INPUT(transient);
            GETPOT_INT_INPUT(deltat_reductions);
            GETPOT_INPUT(timesolver_core);
            GETPOT_INPUT(end_time);
            GETPOT_INPUT(deltat);
            GETPOT_INPUT(timesolver_theta);
            GETPOT_INPUT(timesolver_maxgrowth);
            GETPOT_INPUT(timesolver_tolerance);
            GETPOT_INPUT(timesolver_upper_tolerance);
            GETPOT_INPUT(steadystate_tolerance);
        
            GETPOT_REGISTER(timesolver_norm);
            const unsigned int n_timesolver_norm = input.vector_variable_size("timesolver_norm");
            timesolver_norm.resize(n_timesolver_norm, L2);
            for (unsigned int i=0; i != n_timesolver_norm; ++i)
              {
                int current_norm = 0; // L2
                if (timesolver_norm[i] == H1)
                  current_norm = 1;
                if (timesolver_norm[i] == H2)
                  current_norm = 2;
                current_norm = input("timesolver_norm", current_norm, i);
                if (current_norm == 0)
                  timesolver_norm[i] = L2;
                else if (current_norm == 1)
                  timesolver_norm[i] = H1;
                else if (current_norm == 2)
                  timesolver_norm[i] = H2;
                else
                  timesolver_norm[i] = DISCRETE_L2;
              }
        
        
            GETPOT_INT_INPUT(dimension);
            GETPOT_INPUT(domaintype);
            GETPOT_INPUT(domainfile);
            GETPOT_INPUT(elementtype);
            GETPOT_INPUT(elementorder);
            GETPOT_INPUT(domain_xmin);
            GETPOT_INPUT(domain_ymin);
            GETPOT_INPUT(domain_zmin);
            GETPOT_INPUT(domain_edge_width);
            GETPOT_INPUT(domain_edge_length);
            GETPOT_INPUT(domain_edge_height);
            GETPOT_INT_INPUT(coarsegridx);
            GETPOT_INT_INPUT(coarsegridy);
            GETPOT_INT_INPUT(coarsegridz);
            GETPOT_INT_INPUT(coarserefinements);
            GETPOT_INT_INPUT(extrarefinements);
        
        
            GETPOT_INT_INPUT(nelem_target);
            GETPOT_INPUT(global_tolerance);
            GETPOT_INPUT(refine_fraction);
            GETPOT_INPUT(coarsen_fraction);
            GETPOT_INPUT(coarsen_threshold);
            GETPOT_INT_INPUT(max_adaptivesteps);
            GETPOT_INT_INPUT(initial_adaptivesteps);
        
        
            GETPOT_INT_INPUT(write_interval);
            GETPOT_INPUT(output_xda);
            GETPOT_INPUT(output_xdr);
            GETPOT_INPUT(output_gz);
        #ifndef LIBMESH_HAVE_GZSTREAM
            output_gz                   = false;
        #endif
            GETPOT_INPUT(output_bz2);
        #ifndef LIBMESH_HAVE_BZ2
            output_bz2                  = false;
        #endif
            GETPOT_INPUT(output_gmv);
            GETPOT_INPUT(write_gmv_error);
        #ifndef LIBMESH_HAVE_GMV
            output_gmv                  = false;
            write_gmv_error             = false;
        #endif
            GETPOT_INPUT(output_tecplot);
            GETPOT_INPUT(write_tecplot_error);
        #ifndef LIBMESH_HAVE_TECPLOT_API
            output_tecplot              = false;
            write_tecplot_error         = false;
        #endif
        
        
            GETPOT_REGISTER(system_types);
            const unsigned int n_system_types =
              input.vector_variable_size("system_types");
            if (n_system_types)
              {
                system_types.resize(n_system_types, "");
                for (unsigned int i=0; i != n_system_types; ++i)
                  {
                    system_types[i] = input("system_types", system_types[i], i);
                  }
              }
        
        
            GETPOT_REGISTER(periodic_boundaries);
            const unsigned int n_periodic_bcs =
              input.vector_variable_size("periodic_boundaries");
        
            if (n_periodic_bcs)
              {
                if (domaintype != "square" &&
                    domaintype != "cylinder" &&
                    domaintype != "file" &&
                    domaintype != "od2")
                  {
                    libMesh::out &lt;&lt; "Periodic boundaries need rectilinear domains" &lt;&lt; std::endl;;
                    libmesh_error();
                  }
                for (unsigned int i=0; i != n_periodic_bcs; ++i)
                  {
                    unsigned int myboundary =
                      input("periodic_boundaries", -1, i);
</pre>
</div>
<div class = "comment">
unsigned int pairedboundary = 0;
</div>

<div class ="fragment">
<pre>
                    RealVectorValue translation_vector;
                    if (dimension == 2)
                      switch (myboundary)
                      {
                      case 0:
</pre>
</div>
<div class = "comment">
pairedboundary = 2;
</div>

<div class ="fragment">
<pre>
                        translation_vector = RealVectorValue(0., domain_edge_length);
                        break;
                      case 1:
</pre>
</div>
<div class = "comment">
pairedboundary = 3;
</div>

<div class ="fragment">
<pre>
                        translation_vector = RealVectorValue(-domain_edge_width, 0);
                        break;
                      default:
                        libMesh::out &lt;&lt; "Unrecognized periodic boundary id " &lt;&lt;
                                        myboundary &lt;&lt; std::endl;;
                        libmesh_error();
                      }
                    else if (dimension == 3)
                      switch (myboundary)
                      {
                      case 0:
</pre>
</div>
<div class = "comment">
pairedboundary = 5;
</div>

<div class ="fragment">
<pre>
                        translation_vector = RealVectorValue((Real) 0., (Real) 0., domain_edge_height);
                        break;
                      case 1:
</pre>
</div>
<div class = "comment">
pairedboundary = 3;
</div>

<div class ="fragment">
<pre>
                        translation_vector = RealVectorValue((Real) 0., domain_edge_length, (Real) 0.);
                        break;
                      case 2:
</pre>
</div>
<div class = "comment">
pairedboundary = 4;
</div>

<div class ="fragment">
<pre>
                        translation_vector = RealVectorValue(-domain_edge_width, (Real) 0., (Real) 0.);
                        break;
                      default:
                        libMesh::out &lt;&lt; "Unrecognized periodic boundary id " &lt;&lt;
                                        myboundary &lt;&lt; std::endl;;
                        libmesh_error();
                      }
</pre>
</div>
<div class = "comment">
periodic_boundaries.push_back(PeriodicBoundary(translation_vector));
periodic_boundaries[i].myboundary = myboundary;
periodic_boundaries[i].pairedboundary = pairedboundary;
</div>

<div class ="fragment">
<pre>
                  }
              }
        
</pre>
</div>
<div class = "comment">
Use std::string inputs so GetPot doesn't have to make a bunch
of internal C string copies
</div>

<div class ="fragment">
<pre>
            std::string zero_string = "zero";
            std::string empty_string = "";
        
            GETPOT_REGISTER(dirichlet_condition_types);
            GETPOT_REGISTER(dirichlet_condition_values);
            GETPOT_REGISTER(dirichlet_condition_boundaries);
            GETPOT_REGISTER(dirichlet_condition_variables);
        
            const unsigned int n_dirichlet_conditions=
              input.vector_variable_size("dirichlet_condition_types");
        
            if (n_dirichlet_conditions !=
                input.vector_variable_size("dirichlet_condition_values"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_dirichlet_conditions
                             &lt;&lt; " Dirichlet condition types does not match "
                             &lt;&lt; input.vector_variable_size("dirichlet_condition_values")
                             &lt;&lt; " Dirichlet condition values." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_dirichlet_conditions !=
                input.vector_variable_size("dirichlet_condition_boundaries"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_dirichlet_conditions
                             &lt;&lt; " Dirichlet condition types does not match "
                             &lt;&lt; input.vector_variable_size("dirichlet_condition_boundaries")
                             &lt;&lt; " Dirichlet condition boundaries." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_dirichlet_conditions !=
                input.vector_variable_size("dirichlet_condition_variables"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_dirichlet_conditions
                             &lt;&lt; " Dirichlet condition types does not match "
                             &lt;&lt; input.vector_variable_size("dirichlet_condition_variables")
                             &lt;&lt; " Dirichlet condition variables sets." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            for (unsigned int i=0; i != n_dirichlet_conditions; ++i)
              {
                const std::string func_type =
                  input("dirichlet_condition_types", zero_string, i);
        
                const std::string func_value =
                  input("dirichlet_condition_values", empty_string, i);
        
                const boundary_id_type func_boundary =
                  input("dirichlet_condition_boundaries", boundary_id_type(0), i);
        
                dirichlet_conditions[func_boundary] =
                  (new_function_base(func_type, func_value).release());
        
                const std::string variable_set =
                  input("dirichlet_condition_variables", empty_string, i);
        
                for (unsigned int j=0; j != variable_set.size(); ++j)
                  {
                    if (variable_set[j] == '1')
                      dirichlet_condition_variables[func_boundary].push_back(j);
                    else if (variable_set[j] != '0')
                      {
                        libMesh::out &lt;&lt; "Unable to understand Dirichlet variable set"
                                     &lt;&lt; variable_set &lt;&lt; std::endl;
                        libmesh_error();
                      }
                  }
              }
        
            GETPOT_REGISTER(neumann_condition_types);
            GETPOT_REGISTER(neumann_condition_values);
            GETPOT_REGISTER(neumann_condition_boundaries);
            GETPOT_REGISTER(neumann_condition_variables);
        
            const unsigned int n_neumann_conditions=
              input.vector_variable_size("neumann_condition_types");
        
            if (n_neumann_conditions !=
                input.vector_variable_size("neumann_condition_values"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_neumann_conditions
                             &lt;&lt; " Neumann condition types does not match "
                             &lt;&lt; input.vector_variable_size("neumann_condition_values")
                             &lt;&lt; " Neumann condition values." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_neumann_conditions !=
                input.vector_variable_size("neumann_condition_boundaries"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_neumann_conditions
                             &lt;&lt; " Neumann condition types does not match "
                             &lt;&lt; input.vector_variable_size("neumann_condition_boundaries")
                             &lt;&lt; " Neumann condition boundaries." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_neumann_conditions !=
                input.vector_variable_size("neumann_condition_variables"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_neumann_conditions
                             &lt;&lt; " Neumann condition types does not match "
                             &lt;&lt; input.vector_variable_size("neumann_condition_variables")
                             &lt;&lt; " Neumann condition variables sets." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            for (unsigned int i=0; i != n_neumann_conditions; ++i)
              {
                const std::string func_type =
                  input("neumann_condition_types", zero_string, i);
        
                const std::string func_value =
                  input("neumann_condition_values", empty_string, i);
        
                const boundary_id_type func_boundary =
                  input("neumann_condition_boundaries", boundary_id_type(0), i);
        
                neumann_conditions[func_boundary] =
                  (new_function_base(func_type, func_value).release());
        
                const std::string variable_set =
                  input("neumann_condition_variables", empty_string, i);
        
                for (unsigned int j=0; j != variable_set.size(); ++j)
                  {
                    if (variable_set[j] == '1')
                      neumann_condition_variables[func_boundary].push_back(j);
                    else if (variable_set[j] != '0')
                      {
                        libMesh::out &lt;&lt; "Unable to understand Neumann variable set"
                                     &lt;&lt; variable_set &lt;&lt; std::endl;
                        libmesh_error();
                      }
                  }
              }
        
            GETPOT_REGISTER(initial_condition_types);
            GETPOT_REGISTER(initial_condition_values);
            GETPOT_REGISTER(initial_condition_subdomains);
        
            const unsigned int n_initial_conditions=
              input.vector_variable_size("initial_condition_types");
        
            if (n_initial_conditions !=
                input.vector_variable_size("initial_condition_values"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_initial_conditions
                             &lt;&lt; " initial condition types does not match "
                             &lt;&lt; input.vector_variable_size("initial_condition_values")
                             &lt;&lt; " initial condition values." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_initial_conditions !=
                input.vector_variable_size("initial_condition_subdomains"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_initial_conditions
                             &lt;&lt; " initial condition types does not match "
                             &lt;&lt; input.vector_variable_size("initial_condition_subdomains")
                             &lt;&lt; " initial condition subdomains." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            for (unsigned int i=0; i != n_initial_conditions; ++i)
              {
                const std::string func_type =
                  input("initial_condition_types", zero_string, i);
        
                const std::string func_value =
                  input("initial_condition_values", empty_string, i);
        
                const subdomain_id_type func_subdomain =
                  input("initial_condition_subdomains", subdomain_id_type(0), i);
        
                initial_conditions[func_subdomain] =
                  (new_function_base(func_type, func_value).release());
              }
        
            GETPOT_REGISTER(other_interior_function_types);
            GETPOT_REGISTER(other_interior_function_values);
            GETPOT_REGISTER(other_interior_function_subdomains);
            GETPOT_REGISTER(other_interior_function_ids);
        
            const unsigned int n_other_interior_functions =
              input.vector_variable_size("other_interior_function_types");
        
            if (n_other_interior_functions !=
                input.vector_variable_size("other_interior_function_values"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_interior_functions
                             &lt;&lt; " other interior function types does not match "
                             &lt;&lt; input.vector_variable_size("other_interior_function_values")
                             &lt;&lt; " other interior function values." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_other_interior_functions !=
                input.vector_variable_size("other_interior_function_subdomains"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_interior_functions
                             &lt;&lt; " other interior function types does not match "
                             &lt;&lt; input.vector_variable_size("other_interior_function_subdomains")
                             &lt;&lt; " other interior function subdomains." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_other_interior_functions !=
                input.vector_variable_size("other_interior_function_ids"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_interior_functions
                             &lt;&lt; " other interior function types does not match "
                             &lt;&lt; input.vector_variable_size("other_interior_function_ids")
                             &lt;&lt; " other interior function ids." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            for (unsigned int i=0; i != n_other_interior_functions; ++i)
              {
                const std::string func_type =
                  input("other_interior_function_types", zero_string, i);
        
                const std::string func_value =
                  input("other_interior_function_values", empty_string, i);
        
                const subdomain_id_type func_subdomain =
                  input("other_interior_condition_subdomains", subdomain_id_type(0), i);
        
                const int func_id =
                  input("other_interior_condition_ids", int(0), i);
        
                other_interior_functions[func_id][func_subdomain] =
                  (new_function_base(func_type, func_value).release());
              }
        
            GETPOT_REGISTER(other_boundary_function_types);
            GETPOT_REGISTER(other_boundary_function_values);
            GETPOT_REGISTER(other_boundary_function_boundaries);
            GETPOT_REGISTER(other_boundary_function_ids);
        
            const unsigned int n_other_boundary_functions =
              input.vector_variable_size("other_boundary_function_types");
        
            if (n_other_boundary_functions !=
                input.vector_variable_size("other_boundary_function_values"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_boundary_functions
                             &lt;&lt; " other boundary function types does not match "
                             &lt;&lt; input.vector_variable_size("other_boundary_function_values")
                             &lt;&lt; " other boundary function values." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_other_boundary_functions !=
                input.vector_variable_size("other_boundary_function_boundaries"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_boundary_functions
                             &lt;&lt; " other boundary function types does not match "
                             &lt;&lt; input.vector_variable_size("other_boundary_function_boundaries")
                             &lt;&lt; " other boundary function boundaries." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            if (n_other_boundary_functions !=
                input.vector_variable_size("other_boundary_function_ids"))
              {
                libMesh::out &lt;&lt; "Error: " &lt;&lt; n_other_boundary_functions
                             &lt;&lt; " other boundary function types does not match "
                             &lt;&lt; input.vector_variable_size("other_boundary_function_ids")
                             &lt;&lt; " other boundary function ids." &lt;&lt; std::endl;
        
                libmesh_error();
              }
        
            for (unsigned int i=0; i != n_other_boundary_functions; ++i)
              {
                const std::string func_type =
                  input("other_boundary_function_types", zero_string, i);
        
                const std::string func_value =
                  input("other_boundary_function_values", empty_string, i);
        
                const boundary_id_type func_boundary =
                  input("other_boundary_function_boundaries", boundary_id_type(0), i);
        
                const int func_id =
                  input("other_boundary_function_ids", int(0), i);
        
                other_boundary_functions[func_id][func_boundary] =
                  (new_function_base(func_type, func_value).release());
              }
        
            GETPOT_INPUT(run_simulation);
            GETPOT_INPUT(run_postprocess);
        
        
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
            GETPOT_INPUT(extra_quadrature_order);
        
        
            GETPOT_INPUT(analytic_jacobians);
            GETPOT_INPUT(verify_analytic_jacobians);
            GETPOT_INPUT(numerical_jacobian_h);
            GETPOT_INPUT(print_solution_norms);
            GETPOT_INPUT(print_solutions);
            GETPOT_INPUT(print_residual_norms);
            GETPOT_INPUT(print_residuals);
            GETPOT_INPUT(print_jacobian_norms);
            GETPOT_INPUT(print_jacobians);
            GETPOT_INPUT(print_element_jacobians);
        
        
            GETPOT_INPUT(use_petsc_snes);
            GETPOT_INPUT(time_solver_quiet);
            GETPOT_INPUT(solver_quiet);
            GETPOT_INPUT(solver_verbose);
            GETPOT_INPUT(require_residual_reduction);
            GETPOT_INPUT(min_step_length);
            GETPOT_INT_INPUT(max_linear_iterations);
            GETPOT_INT_INPUT(max_nonlinear_iterations);
            GETPOT_INPUT(relative_step_tolerance);
            GETPOT_INPUT(relative_residual_tolerance);
            GETPOT_INPUT(initial_linear_tolerance);
            GETPOT_INPUT(minimum_linear_tolerance);
            GETPOT_INPUT(linear_tolerance_multiplier);
        
        
            GETPOT_INT_INPUT(initial_sobolev_order);
            GETPOT_INT_INPUT(initial_extra_quadrature);
            GETPOT_INPUT(indicator_type);
            GETPOT_INT_INPUT(sobolev_order);
        
            GETPOT_INPUT(system_config_file);
        
          std::vector&lt;std::string&gt; bad_variables =
            input.unidentified_arguments(variable_names);
        
          if (libMesh::processor_id() == 0 && !bad_variables.empty())
            {
              std::cerr &lt;&lt; "ERROR: Unrecognized variables:" &lt;&lt; std::endl;
              for (unsigned int i = 0; i != bad_variables.size(); ++i)
                std::cerr &lt;&lt; bad_variables[i] &lt;&lt; std::endl;
              std::cerr &lt;&lt; "Not found among recognized variables:" &lt;&lt; std::endl;
              for (unsigned int i = 0; i != variable_names.size(); ++i)
                std::cerr &lt;&lt; variable_names[i] &lt;&lt; std::endl;
              libmesh_error();
            }
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file heatsystem.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/getpot.h"
        
        #include "heatsystem.h"
        
        #include "libmesh/fe_base.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/mesh.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/system.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/zero_function.h"
        #include "libmesh/boundary_info.h"
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/dof_map.h"
        
        void HeatSystem::init_data ()
        {
          T_var = this-&gt;add_variable ("T", static_cast&lt;Order&gt;(_fe_order),
                              Utility::string_to_enum&lt;FEFamily&gt;(_fe_family));
        
</pre>
</div>
<div class = "comment">
Make sure the input file heat.in exists, and parse it.
</div>

<div class ="fragment">
<pre>
          {
            std::ifstream i("heat.in");
            if (!i)
              {
                std::cerr &lt;&lt; '[' &lt;&lt; libMesh::processor_id()
                          &lt;&lt; "] Can't find heat.in; exiting early."
                          &lt;&lt; std::endl;
                libmesh_error();
              }
          }
          GetPot infile("heat.in");
          _k = infile("k", 1.0);
          _analytic_jacobians = infile("analytic_jacobians", true);
        
          parameters.push_back(_k);
        
</pre>
</div>
<div class = "comment">
Set equation system parameters _k and theta, so they can be read by the exact solution
</div>

<div class ="fragment">
<pre>
          this-&gt;get_equation_systems().parameters.set&lt;Real&gt;("_k") = _k;
        
          this-&gt;time_evolving(T_var);
        
          const boundary_id_type all_ids[4] = {0, 1, 2, 3};
          std::set&lt;boundary_id_type&gt; all_bdys(all_ids, all_ids+(2*2));
        
          std::vector&lt;unsigned int&gt; T_only(1, T_var);
        
          ZeroFunction&lt;Number&gt; zero;
        
          this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (all_bdys, T_only, &zero));
        
          FEMSystem::init_data();
        }
        
        
        
        void HeatSystem::init_context(DiffContext &context)
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
          c.element_fe_var[0]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
We'll have a more automatic solution to preparing adjoint
solutions for time integration, eventually...
</div>

<div class ="fragment">
<pre>
          if (c.is_adjoint())
            {
</pre>
</div>
<div class = "comment">
A reference to the system context is built with
</div>

<div class ="fragment">
<pre>
              const System & sys = c.get_system();
        
</pre>
</div>
<div class = "comment">
Get a pointer to the adjoint solution vector
</div>

<div class ="fragment">
<pre>
              NumericVector&lt;Number&gt; &adjoint_solution =
                const_cast&lt;System &&gt;(sys).get_adjoint_solution(0);
        
</pre>
</div>
<div class = "comment">
Add this adjoint solution to the vectors that diff context should localize
</div>

<div class ="fragment">
<pre>
              c.add_localized_vector(adjoint_solution, sys);
            }
        
          FEMSystem::init_context(context);
        }
        
        #define optassert(X) {if (!(X)) libmesh_error();}
</pre>
</div>
<div class = "comment">
#define optassert(X) libmesh_assert(X);


<br><br></div>

<div class ="fragment">
<pre>
        bool HeatSystem::element_time_derivative (bool request_jacobian,
                                                  DiffContext &context)
        {
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
const std::vector<std::vector<Real> >          &phi = c.element_fe_var[0]->get_phi();
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &dphi = c.element_fe_var[0]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
Workaround for weird FC6 bug
</div>

<div class ="fragment">
<pre>
          optassert(c.dof_indices_var.size() &gt; 0);
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[0].size();
        
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
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
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
        
              for (unsigned int i=0; i != n_u_dofs; i++)
                F(i) += JxW[qp] * -parameters[0] * (grad_T * dphi[i][qp]);
              if (compute_jacobian)
                for (unsigned int i=0; i != n_u_dofs; i++)
                  for (unsigned int j=0; j != n_u_dofs; ++j)
                    K(i,j) += JxW[qp] * -parameters[0] * (dphi[i][qp] * dphi[j][qp]);
            } // end of the quadrature point qp-loop
        
          return compute_jacobian;
        }
        
</pre>
</div>
<div class = "comment">
Perturb and accumulate dual weighted residuals
</div>

<div class ="fragment">
<pre>
        void HeatSystem::perturb_accumulate_residuals(const ParameterVector& parameters_in)
        {
          const unsigned int Np = parameters_in.size();
        
          for (unsigned int j=0; j != Np; ++j)
            {
              Number old_parameter = *parameters_in[j];
        
              *parameters_in[j] = old_parameter - dp;
        
              this-&gt;assembly(true, false);
        
              this-&gt;rhs-&gt;close();
        
              AutoPtr&lt;NumericVector&lt;Number&gt; &gt; R_minus = this-&gt;rhs-&gt;clone();
        
</pre>
</div>
<div class = "comment">
The contribution at a single time step would be [f(z;p+dp) - <partialu/partialt, z>(p+dp) - <g(u),z>(p+dp)] * dt
But since we compute the residual already scaled by dt, there is no need for the * dt
</div>

<div class ="fragment">
<pre>
              R_minus_dp += -R_minus-&gt;dot(this-&gt;get_adjoint_solution(0));
        
              *parameters_in[j] = old_parameter + dp;
        
              this-&gt;assembly(true, false);
        
              this-&gt;rhs-&gt;close();
        
              AutoPtr&lt;NumericVector&lt;Number&gt; &gt; R_plus = this-&gt;rhs-&gt;clone();
        
              R_plus_dp += -R_plus-&gt;dot(this-&gt;get_adjoint_solution(0));
        
              *parameters_in[j] = old_parameter;
        
            }
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file initial.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "initial.h"
        
        using namespace libMesh;
        
        void read_initial_parameters()
        {
        }
        
        void finish_initialization()
        {
        }
        
        
        
</pre>
</div>
<div class = "comment">
Initial conditions
</div>

<div class ="fragment">
<pre>
        Number initial_value(const Point& p,
                             const Parameters&,
                             const std::string&,
                             const std::string&)
        {
          Real x = p(0), y = p(1);
        
          return sin(M_PI * x) * sin(M_PI * y);
        }
        
        
        
        Gradient initial_grad(const Point& p,
                              const Parameters&,
                              const std::string&,
                              const std::string&)
        {
          Real x = p(0), y = p(1);
        
          return Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
          M_PI*sin(M_PI * x) * cos(M_PI * y));
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoint_initial.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> adjoint_read_initial_parameters();
  <B><FONT COLOR="#228B22">void</FONT></B> adjoint_finish_initialization();
  
  Number adjoint_initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Gradient adjoint_initial_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.h without comments: </h1> 
<pre> 
  #ifndef __fem_parameters_h__
  #define __fem_parameters_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/enum_norm_type.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/function_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/id_types.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundaries.h&quot;</FONT></B>
  
  #include &lt;limits&gt;
  #include &lt;map&gt;
  #include &lt;string&gt;
  #include &lt;vector&gt;
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> FEMParameters
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
      FEMParameters();
  
      ~FEMParameters();
  
      <B><FONT COLOR="#228B22">void</FONT></B> read(GetPot &amp;input);
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_timestep, n_timesteps;
      <B><FONT COLOR="#228B22">bool</FONT></B> transient;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> deltat_reductions;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string timesolver_core;
      Real end_time, deltat, timesolver_theta,
  	 timesolver_maxgrowth, timesolver_tolerance,
  	 timesolver_upper_tolerance, steadystate_tolerance;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FEMNormType&gt; timesolver_norm;
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dimension;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string domaintype, domainfile, elementtype;
      Real elementorder;
      Real domain_xmin, domain_ymin, domain_zmin;
      Real domain_edge_width, domain_edge_length, domain_edge_height;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarsegridx, coarsegridy, coarsegridz;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarserefinements, extrarefinements;
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nelem_target;
      Real global_tolerance;
      Real refine_fraction, coarsen_fraction, coarsen_threshold;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_adaptivesteps;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_adaptivesteps;
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> write_interval;
      <B><FONT COLOR="#228B22">bool</FONT></B> write_gmv_error, write_tecplot_error, output_xda, output_xdr,
  	 output_bz2, output_gz, output_gmv, output_tecplot;
  
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; system_types;
  
  
  
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt; initial_conditions;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;  dirichlet_conditions,
                                                          neumann_conditions;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;boundary_id_type, std::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; &gt; dirichlet_condition_variables,
                                                             neumann_condition_variables;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">int</FONT></B>, std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt; &gt;
        other_interior_functions;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;<B><FONT COLOR="#228B22">int</FONT></B>, std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt; &gt;
        other_boundary_functions;
  
  
      <B><FONT COLOR="#228B22">bool</FONT></B> run_simulation, run_postprocess;
  
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; fe_family;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; fe_order;
      <B><FONT COLOR="#228B22">int</FONT></B> extra_quadrature_order;
  
  
      <B><FONT COLOR="#228B22">bool</FONT></B> analytic_jacobians;
      Real verify_analytic_jacobians;
      Real numerical_jacobian_h;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> print_solution_norms, print_solutions,
           print_residual_norms, print_residuals,
           print_jacobian_norms, print_jacobians,
           print_element_jacobians;
  
  
      <B><FONT COLOR="#228B22">bool</FONT></B> use_petsc_snes;
      <B><FONT COLOR="#228B22">bool</FONT></B> time_solver_quiet, solver_quiet, solver_verbose, require_residual_reduction;
      Real min_step_length;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_linear_iterations, max_nonlinear_iterations;
      Real relative_step_tolerance, relative_residual_tolerance,
  	 initial_linear_tolerance, minimum_linear_tolerance,
  	 linear_tolerance_multiplier;
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_sobolev_order;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_extra_quadrature;
  
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string indicator_type;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> sobolev_order;
  
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string system_config_file;
  };
  
  #endif <I><FONT COLOR="#B22222">// __fem_parameters_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file heatsystem.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/enum_fe_family.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameter_vector.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">class</FONT></B> HeatSystem : <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    HeatSystem(EquationSystems&amp; es,
                 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                 <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : FEMSystem(es, name_in, number_in),
      _k(1.0),
      _fe_family(<B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), _fe_order(1),
      _analytic_jacobians(true), R_plus_dp(0.0), R_minus_dp(0.0), dp(1.e-6) { qoi.resize(1); }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string &amp; fe_family() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_family;  }
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; fe_order() { <B><FONT COLOR="#A020F0">return</FONT></B> _fe_order;  }
    Real &amp; k() { <B><FONT COLOR="#A020F0">return</FONT></B> _k;  }
    <B><FONT COLOR="#228B22">bool</FONT></B> &amp; analytic_jacobians() { <B><FONT COLOR="#A020F0">return</FONT></B> _analytic_jacobians; }
  
    <B><FONT COLOR="#228B22">void</FONT></B> perturb_accumulate_residuals(<B><FONT COLOR="#228B22">const</FONT></B> ParameterVector&amp; parameters);
  
    Number &amp; compute_final_sensitivity()
      {
        final_sensitivity = -(R_plus_dp - R_minus_dp)/(2*dp);
  
        <B><FONT COLOR="#A020F0">return</FONT></B> final_sensitivity;
      }
  
      <B><FONT COLOR="#228B22">void</FONT></B> set_tf(Real val)
    {
      tf = val;
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
  
    Number &amp;get_QoI_value(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> QoI_index)
      {
        <B><FONT COLOR="#A020F0">return</FONT></B> computed_QoI[QoI_index];
      }
  
  <B><FONT COLOR="#228B22">protected</FONT></B>:
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context (DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                          DiffContext &amp;context);
  
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> element_qoi_derivative (DiffContext &amp;context,
  				       <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>);
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; parameters;
  
    ParameterVector parameter_vector;
  
    Real _k;
  
    Real tf;
  
    Number computed_QoI[1];
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string _fe_family;
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _fe_order;
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> T_var;
  
    <B><FONT COLOR="#228B22">bool</FONT></B> _analytic_jacobians;
  
    Number R_plus_dp;
    Number R_minus_dp;
  
    Real dp;
  
    Number final_sensitivity;
  
  };
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file initial.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> read_initial_parameters();
  <B><FONT COLOR="#228B22">void</FONT></B> finish_initialization();
  
  Number initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Gradient initial_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoint_initial.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;adjoint_initial.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> adjoint_read_initial_parameters()
  {
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> adjoint_finish_initialization()
  {
  }
  
  
  
  Number adjoint_initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Real x = p(0), y = p(1);
  
    Number val = 0.;
  
    val = sin(M_PI * x) * sin(M_PI * y);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> val;
  }
  
  
  
  Gradient adjoint_initial_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Real x = p(0), y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
    M_PI*sin(M_PI * x) * cos(M_PI * y));
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoints_ex5.C without comments: </h1> 
<pre> 
  
  
  
  
  
  
  
  
  #include <B><FONT COLOR="#BC8F8F">&quot;initial.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;adjoint_initial.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;heatsystem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system_norm.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/petsc_diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/euler_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/euler2_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/twostep_time_solver.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/solution_history.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/memory_solution_history.h&quot;</FONT></B>
  
  #include &lt;iostream&gt;
  #include &lt;sys/time.h&gt;
  #include &lt;iomanip&gt;
  
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
  
  <B><FONT COLOR="#228B22">void</FONT></B> set_system_parameters(HeatSystem &amp;system, FEMParameters &amp;param)
  {
    system.fe_family() = param.fe_family[0];
    system.fe_order() = param.fe_order[0];
  
    system.analytic_jacobians() = param.analytic_jacobians;
  
    system.verify_analytic_jacobians = param.verify_analytic_jacobians;
    system.numerical_jacobian_h = param.numerical_jacobian_h;
  
    system.print_solution_norms    = param.print_solution_norms;
    system.print_solutions         = param.print_solutions;
    system.print_residual_norms    = param.print_residual_norms;
    system.print_residuals         = param.print_residuals;
    system.print_jacobian_norms    = param.print_jacobian_norms;
    system.print_jacobians         = param.print_jacobians;
    system.print_element_jacobians = param.print_element_jacobians;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.transient)
      {
        UnsteadySolver *innersolver;
         <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_core == <B><FONT COLOR="#BC8F8F">&quot;euler&quot;</FONT></B>)
          {
            EulerSolver *eulersolver =
              <B><FONT COLOR="#A020F0">new</FONT></B> EulerSolver(system);
  
            eulersolver-&gt;theta = param.timesolver_theta;
            innersolver = eulersolver;
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;This example (and unsteady adjoints in libMesh) only support Backward Euler and explicit methods.&quot;</FONT></B>
                      &lt;&lt;  std::endl;
            libmesh_error();
          }
         system.time_solver =
            AutoPtr&lt;TimeSolver&gt;(innersolver);
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      system.time_solver =
        AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
  
    MemorySolutionHistory heatsystem_solution_history(system);
    system.time_solver-&gt;set_solution_history(heatsystem_solution_history);
  
    system.time_solver-&gt;reduce_deltat_on_diffsolver_failure =
                                          param.deltat_reductions;
    system.time_solver-&gt;quiet           = param.time_solver_quiet;
  
    <B><FONT COLOR="#228B22">typedef</FONT></B>
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::
      const_iterator Iter;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (Iter i = param.dirichlet_conditions.begin();
         i != param.dirichlet_conditions.end(); ++i)
      {
        boundary_id_type b = i-&gt;first;
        FunctionBase&lt;Number&gt; *f = i-&gt;second;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; bdys; bdys.insert(b);
        system.get_dof_map().add_dirichlet_boundary
          (DirichletBoundary
             (bdys, param.dirichlet_condition_variables[b], f));
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Added Dirichlet boundary &quot;</FONT></B> &lt;&lt; b &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; for variables &quot;</FONT></B>;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> vi=0; vi !=
             param.dirichlet_condition_variables[b].size(); ++vi)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; param.dirichlet_condition_variables[b][vi];
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl;
      }
  
    system.deltat = param.deltat;
  
    system.extra_quadrature_order = param.extra_quadrature_order;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.use_petsc_snes)
      {
  #ifdef LIBMESH_HAVE_PETSC
        PetscDiffSolver *solver = <B><FONT COLOR="#A020F0">new</FONT></B> PetscDiffSolver(system);
        system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
        libmesh_error();
  #endif
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        NewtonSolver *solver = <B><FONT COLOR="#A020F0">new</FONT></B> NewtonSolver(system);
        system.time_solver-&gt;diff_solver() = AutoPtr&lt;DiffSolver&gt;(solver);
  
        solver-&gt;quiet                       = param.solver_quiet;
        solver-&gt;verbose                     = param.solver_verbose;
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
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
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
  
    Mesh mesh(init.comm(), param.dimension);
  
    AutoPtr&lt;MeshRefinement&gt; mesh_refinement(<B><FONT COLOR="#A020F0">new</FONT></B> MeshRefinement(mesh));
  
    EquationSystems equation_systems (mesh);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Building mesh&quot;</FONT></B> &lt;&lt; std::endl;
  
    ElemType elemtype;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.elementtype == <B><FONT COLOR="#BC8F8F">&quot;tri&quot;</FONT></B> ||
        param.elementtype == <B><FONT COLOR="#BC8F8F">&quot;unstructured&quot;</FONT></B>)
      elemtype = TRI3;
    <B><FONT COLOR="#A020F0">else</FONT></B>
      elemtype = QUAD4;
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square
      (mesh, param.coarsegridx, param.coarsegridy,
       param.domain_xmin, param.domain_xmin + param.domain_edge_width,
       param.domain_ymin, param.domain_ymin + param.domain_edge_length,
       elemtype);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Building system&quot;</FONT></B> &lt;&lt; std::endl;
  
    HeatSystem &amp;system = equation_systems.add_system&lt;HeatSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;HeatSystem&quot;</FONT></B>);
  
    set_system_parameters(system, param);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initializing systems&quot;</FONT></B> &lt;&lt; std::endl;
  
    equation_systems.init ();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != param.extrarefinements; ++i)
      {
        mesh_refinement-&gt;uniformly_refine(1);
        equation_systems.reinit();
      }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Setting primal initial conditions&quot;</FONT></B>&lt;&lt;std::endl;
  
    read_initial_parameters();
  
    system.project_solution(initial_value, initial_grad,
                                    equation_systems.parameters);
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::string &amp; adjoint_solution_name = <B><FONT COLOR="#BC8F8F">&quot;adjoint_solution0&quot;</FONT></B>;
    system.add_vector(<B><FONT COLOR="#BC8F8F">&quot;adjoint_solution0&quot;</FONT></B>, false, GHOSTED);
  
    finish_initialization();
  
    write_output(equation_systems, 0, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
  
    mesh.print_info();
    equation_systems.print_info();
  
  #ifdef NDEBUG
    <B><FONT COLOR="#A020F0">try</FONT></B>
    {
  #endif
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B>
                    &lt;&lt; system.time &lt;&lt; std::endl;
  
  	system.solve();
  
  	<B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time + system.deltat&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Advancing timestep&quot;</FONT></B>&lt;&lt;std::endl&lt;&lt;std::endl;
          system.time_solver-&gt;advance_timestep();
  
            write_output(equation_systems, t_step+1, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>);
        }
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solving the adjoint problem&quot;</FONT></B> &lt;&lt; std::endl;
  
  
    system.set_vector_preservation(adjoint_solution_name, true);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Setting adjoint initial conditions Z(&quot;</FONT></B>&lt;&lt;system.time&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;)&quot;</FONT></B>&lt;&lt;std::endl;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Retrieving solutions at time t=&quot;</FONT></B>&lt;&lt;system.time&lt;&lt;std::endl;
    system.time_solver-&gt;adjoint_advance_timestep();
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time + system.deltat&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;_old_nonlinear_solution&quot;</FONT></B>), 0, H1) &lt;&lt; std::endl;
  
    system.project_vector(adjoint_initial_value, adjoint_initial_grad, equation_systems.parameters, system.get_adjoint_solution(0));
  
    system.set_adjoint_already_solved(true);
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|Z(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(system.get_adjoint_solution(), 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
  
    write_output(equation_systems, param.n_timesteps, <B><FONT COLOR="#BC8F8F">&quot;dual&quot;</FONT></B>);
  
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=param.initial_timestep;
         t_step != param.initial_timestep + param.n_timesteps; ++t_step)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving adjoint time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B>
  		&lt;&lt; system.time &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Retrieving solutions at time t=&quot;</FONT></B>&lt;&lt;system.time&lt;&lt;std::endl;
        system.time_solver-&gt;adjoint_advance_timestep();
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time + system.deltat &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;_old_nonlinear_solution&quot;</FONT></B>), 0, H1) &lt;&lt; std::endl;
  
        system.set_adjoint_already_solved(false);
  
        system.adjoint_solve();
  
        system.set_adjoint_already_solved(true);
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|Z(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B>&lt;&lt; system.calculate_norm(system.get_adjoint_solution(), 0, H1) &lt;&lt; std::endl &lt;&lt; std::endl;
  
        NumericVector&lt;Number&gt; &amp;primal_solution = *system.solution;
  
        NumericVector&lt;Number&gt; &amp;dual_solution_0 = system.get_adjoint_solution(0);
  
        primal_solution.swap(dual_solution_0);
  
        write_output(equation_systems, param.n_timesteps - (t_step + 1), <B><FONT COLOR="#BC8F8F">&quot;dual&quot;</FONT></B>);
  
        primal_solution.swap(dual_solution_0);
      }
  
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=param.initial_timestep;
           t_step != param.initial_timestep + param.n_timesteps; ++t_step)
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Retrieving &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B>
                    &lt;&lt; system.time &lt;&lt; std::endl;
  
  	system.time_solver-&gt;retrieve_timestep();
  
  	<B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time + system.deltat &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;_old_nonlinear_solution&quot;</FONT></B>), 0, H1) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|Z(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B>&lt;&lt; system.calculate_norm(system.get_adjoint_solution(0), 0, H1) &lt;&lt; std::endl &lt;&lt; std::endl;
  
  	(dynamic_cast&lt;HeatSystem&amp;&gt;(system)).perturb_accumulate_residuals(dynamic_cast&lt;HeatSystem&amp;&gt;(system).get_parameter_vector());
  
  	system.time += system.deltat;
        }
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Retrieving &quot;</FONT></B> &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; final time = &quot;</FONT></B>
  	      &lt;&lt; system.time &lt;&lt; std::endl;
  
      system.time_solver-&gt;retrieve_timestep();
  
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time + system.deltat &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(*system.solution, 0, H1) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|U(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B> &lt;&lt; system.calculate_norm(system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;_old_nonlinear_solution&quot;</FONT></B>), 0, H1) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;|Z(&quot;</FONT></B> &lt;&lt;system.time&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)|= &quot;</FONT></B>&lt;&lt; system.calculate_norm(system.get_adjoint_solution(0), 0, H1) &lt;&lt; std::endl&lt;&lt;std::endl;
  
      (dynamic_cast&lt;HeatSystem&amp;&gt;(system)).perturb_accumulate_residuals(dynamic_cast&lt;HeatSystem&amp;&gt;(system).get_parameter_vector());
  
      Number sensitivity_0_0 = (dynamic_cast&lt;HeatSystem&amp;&gt;(system)).compute_final_sensitivity();
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Sensitivity of QoI 0 w.r.t parameter 0 is: &quot;</FONT></B> &lt;&lt; sensitivity_0_0 &lt;&lt; std::endl;
  
  #ifdef NDEBUG
  }
    <B><FONT COLOR="#A020F0">catch</FONT></B> (...)
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id()
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Caught exception; exiting early.&quot;</FONT></B> &lt;&lt; std::endl;
    }
  #endif
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id()
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Completing output.&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_ENABLE_AMR
</FONT></I>  }
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;heatsystem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> HeatSystem::element_qoi_derivative (DiffContext &amp;context,
                                              <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;          &amp;phi = c.element_fe_var[0]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_T_dofs = c.dof_indices_var[0].size();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    DenseSubVector&lt;Number&gt; &amp;Q = *c.elem_qoi_subderivatives[0][0];
  
    <B><FONT COLOR="#228B22">const</FONT></B> System &amp; sys = c.get_system();
  
    NumericVector&lt;Number&gt; &amp;adjoint_solution = const_cast&lt;System &amp;&gt;(sys).get_adjoint_solution(0);
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; old_adjoint (n_qpoints, 0);
  
    c.interior_values&lt;Number&gt;(0, adjoint_solution, old_adjoint);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_T_dofs; i++)
  	{
  	  Q(i) += -JxW[qp] * old_adjoint[qp] * phi[i][qp] ;
        	}
  
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file factoryfunction.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/factory.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/function_base.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> ExampleOneFunction : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Number&gt;
  {
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;  <I><FONT COLOR="#B22222">/*p*/</FONT></I>,
                               <B><FONT COLOR="#228B22">const</FONT></B> Real <I><FONT COLOR="#B22222">/*time*/</FONT></I>)
      {
        <B><FONT COLOR="#A020F0">return</FONT></B> 1;
      }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;  <I><FONT COLOR="#B22222">/*p*/</FONT></I>,
                             <B><FONT COLOR="#228B22">const</FONT></B> Real <I><FONT COLOR="#B22222">/*time*/</FONT></I>,
                             DenseVector&lt;Number&gt;&amp; output)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != output.size(); ++i)
          output(i) = 1;
      }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init() {}
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> clear() {}
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() <B><FONT COLOR="#228B22">const</FONT></B> {
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
        (<B><FONT COLOR="#A020F0">new</FONT></B> ExampleOneFunction());
    }
  };
  
  #ifdef LIBMESH_USE_SEPARATE_NAMESPACE
  namespace libMesh {
  #endif
  
  <B><FONT COLOR="#228B22">template</FONT></B>&lt;&gt;
  <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;std::string, Factory&lt;FunctionBase&lt;Number&gt; &gt;*&gt;&amp;
  Factory&lt;FunctionBase&lt;Number&gt; &gt;::factory_map()
  {
    <B><FONT COLOR="#228B22">static</FONT></B> std::map&lt;std::string, Factory&lt;FunctionBase&lt;Number&gt; &gt;*&gt; _map;
    <B><FONT COLOR="#A020F0">return</FONT></B> _map;
  }
  
  FactoryImp&lt;ExampleOneFunction, FunctionBase&lt;Number&gt; &gt; example_one_factory (<B><FONT COLOR="#BC8F8F">&quot;example_one&quot;</FONT></B>);
  
  #ifdef LIBMESH_USE_SEPARATE_NAMESPACE
  } <I><FONT COLOR="#B22222">// namespace libMesh
</FONT></I>  #endif
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/factory.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parsed_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/zero_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B> <I><FONT COLOR="#B22222">// RealVectorValue
</FONT></I>  
  #define GETPOT_INPUT(A) { A = input(#A, A);\
    <B><FONT COLOR="#228B22">const</FONT></B> std::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); }
  #define GETPOT_INT_INPUT(A) { A = input(#A, (<B><FONT COLOR="#228B22">int</FONT></B>)A);\
    <B><FONT COLOR="#228B22">const</FONT></B> std::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); }
  
  #define GETPOT_FUNCTION_INPUT(A) {\
    <B><FONT COLOR="#228B22">const</FONT></B> std::string type  = input(#A <B><FONT COLOR="#BC8F8F">&quot;_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;zero&quot;</FONT></B>);\
    <B><FONT COLOR="#228B22">const</FONT></B> std::string value = input(#A <B><FONT COLOR="#BC8F8F">&quot;_value&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>);\
    A = new_function_base(type, value); }
  #define GETPOT_REGISTER(A) { \
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string stringval = input(#A, std::string());\
    variable_names.push_back(std::string(#A <B><FONT COLOR="#BC8F8F">&quot;=&quot;</FONT></B>) + stringval); }
  
  <B><FONT COLOR="#5F9EA0">FEMParameters</FONT></B>::FEMParameters() :
        initial_timestep(0), n_timesteps(100),
        transient(true),
        deltat_reductions(0),
        timesolver_core(<B><FONT COLOR="#BC8F8F">&quot;euler&quot;</FONT></B>),
        end_time(std::numeric_limits&lt;Real&gt;::max()),
        deltat(0.0001), timesolver_theta(0.5),
        timesolver_maxgrowth(0.), timesolver_tolerance(0.),
        timesolver_upper_tolerance(0.),
        steadystate_tolerance(0.),
        timesolver_norm(0,L2),
  
        dimension(2),
        domaintype(<B><FONT COLOR="#BC8F8F">&quot;square&quot;</FONT></B>), domainfile(<B><FONT COLOR="#BC8F8F">&quot;mesh.xda&quot;</FONT></B>), elementtype(<B><FONT COLOR="#BC8F8F">&quot;quad&quot;</FONT></B>),
        elementorder(2),
        domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
        domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
        coarsegridx(1), coarsegridy(1), coarsegridz(1),
        coarserefinements(4), extrarefinements(0),
  
        nelem_target(8000), global_tolerance(0.0),
        refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
        max_adaptivesteps(1),
        initial_adaptivesteps(0),
  
        write_interval(10),
        write_gmv_error(false), write_tecplot_error(false),
        output_xda(false), output_xdr(false),
        output_bz2(true), output_gz(true),
        output_gmv(false), output_tecplot(false),
  
        system_types(0),
  
  
        run_simulation(true), run_postprocess(false),
  
        fe_family(1, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), fe_order(1, 1),
        extra_quadrature_order(0),
  
        analytic_jacobians(true), verify_analytic_jacobians(0.0),
        numerical_jacobian_h(TOLERANCE),
        print_solution_norms(false), print_solutions(false),
        print_residual_norms(false), print_residuals(false),
        print_jacobian_norms(false), print_jacobians(false),
        print_element_jacobians(false),
  
        use_petsc_snes(false),
        time_solver_quiet(true), solver_quiet(true), solver_verbose(false),
        require_residual_reduction(true),
        min_step_length(1e-5),
        max_linear_iterations(200000), max_nonlinear_iterations(20),
        relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
        initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
        linear_tolerance_multiplier(1.e-3),
  
        initial_sobolev_order(1),
        initial_extra_quadrature(0),
        indicator_type(<B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>), sobolev_order(1)
  {
  }
  
  <B><FONT COLOR="#5F9EA0">FEMParameters</FONT></B>::~FEMParameters()
  {
    <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
         i = initial_conditions.begin(); i != initial_conditions.end();
         ++i)
      <B><FONT COLOR="#A020F0">delete</FONT></B> i-&gt;second;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
         i = dirichlet_conditions.begin(); i != dirichlet_conditions.end();
         ++i)
      <B><FONT COLOR="#A020F0">delete</FONT></B> i-&gt;second;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
         i = neumann_conditions.begin(); i != neumann_conditions.end();
         ++i)
      <B><FONT COLOR="#A020F0">delete</FONT></B> i-&gt;second;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;<B><FONT COLOR="#228B22">int</FONT></B>,
  	 <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;
         &gt;::iterator
         i = other_boundary_functions.begin(); i != other_boundary_functions.end();
         ++i)
      <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;boundary_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
           j = i-&gt;second.begin(); j != i-&gt;second.end();
           ++j)
        <B><FONT COLOR="#A020F0">delete</FONT></B> j-&gt;second;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;<B><FONT COLOR="#228B22">int</FONT></B>,
  	 <B><FONT COLOR="#5F9EA0">std</FONT></B>::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;
         &gt;::iterator
         i = other_interior_functions.begin(); i != other_interior_functions.end();
         ++i)
      <B><FONT COLOR="#A020F0">for</FONT></B> (std::map&lt;subdomain_id_type, FunctionBase&lt;Number&gt; *&gt;::iterator
           j = i-&gt;second.begin(); j != i-&gt;second.end();
           ++j)
        <B><FONT COLOR="#A020F0">delete</FONT></B> j-&gt;second;
  }
  
  
  AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; new_function_base(<B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; func_type,
                                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; func_value)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (func_type == <B><FONT COLOR="#BC8F8F">&quot;parsed&quot;</FONT></B>)
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
        (<B><FONT COLOR="#A020F0">new</FONT></B> ParsedFunction&lt;Number&gt;(func_value));
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (func_type == <B><FONT COLOR="#BC8F8F">&quot;factory&quot;</FONT></B>)
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
        (Factory&lt;FunctionBase&lt;Number&gt; &gt;::build(func_value).release());
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (func_type == <B><FONT COLOR="#BC8F8F">&quot;zero&quot;</FONT></B>)
      <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
        (<B><FONT COLOR="#A020F0">new</FONT></B> ZeroFunction&lt;Number&gt;);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      libmesh_not_implemented();
  
    <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;();
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> FEMParameters::read(GetPot &amp;input)
  {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; variable_names;
  
      GETPOT_INT_INPUT(initial_timestep);
      GETPOT_INT_INPUT(n_timesteps);
      GETPOT_INPUT(transient);
      GETPOT_INT_INPUT(deltat_reductions);
      GETPOT_INPUT(timesolver_core);
      GETPOT_INPUT(end_time);
      GETPOT_INPUT(deltat);
      GETPOT_INPUT(timesolver_theta);
      GETPOT_INPUT(timesolver_maxgrowth);
      GETPOT_INPUT(timesolver_tolerance);
      GETPOT_INPUT(timesolver_upper_tolerance);
      GETPOT_INPUT(steadystate_tolerance);
  
      GETPOT_REGISTER(timesolver_norm);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesolver_norm = input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;timesolver_norm&quot;</FONT></B>);
      timesolver_norm.resize(n_timesolver_norm, L2);
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_timesolver_norm; ++i)
        {
          <B><FONT COLOR="#228B22">int</FONT></B> current_norm = 0; <I><FONT COLOR="#B22222">// L2
</FONT></I>          <B><FONT COLOR="#A020F0">if</FONT></B> (timesolver_norm[i] == H1)
            current_norm = 1;
          <B><FONT COLOR="#A020F0">if</FONT></B> (timesolver_norm[i] == H2)
            current_norm = 2;
          current_norm = input(<B><FONT COLOR="#BC8F8F">&quot;timesolver_norm&quot;</FONT></B>, current_norm, i);
          <B><FONT COLOR="#A020F0">if</FONT></B> (current_norm == 0)
            timesolver_norm[i] = L2;
          <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (current_norm == 1)
            timesolver_norm[i] = H1;
          <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (current_norm == 2)
            timesolver_norm[i] = H2;
          <B><FONT COLOR="#A020F0">else</FONT></B>
            timesolver_norm[i] = DISCRETE_L2;
        }
  
  
      GETPOT_INT_INPUT(dimension);
      GETPOT_INPUT(domaintype);
      GETPOT_INPUT(domainfile);
      GETPOT_INPUT(elementtype);
      GETPOT_INPUT(elementorder);
      GETPOT_INPUT(domain_xmin);
      GETPOT_INPUT(domain_ymin);
      GETPOT_INPUT(domain_zmin);
      GETPOT_INPUT(domain_edge_width);
      GETPOT_INPUT(domain_edge_length);
      GETPOT_INPUT(domain_edge_height);
      GETPOT_INT_INPUT(coarsegridx);
      GETPOT_INT_INPUT(coarsegridy);
      GETPOT_INT_INPUT(coarsegridz);
      GETPOT_INT_INPUT(coarserefinements);
      GETPOT_INT_INPUT(extrarefinements);
  
  
      GETPOT_INT_INPUT(nelem_target);
      GETPOT_INPUT(global_tolerance);
      GETPOT_INPUT(refine_fraction);
      GETPOT_INPUT(coarsen_fraction);
      GETPOT_INPUT(coarsen_threshold);
      GETPOT_INT_INPUT(max_adaptivesteps);
      GETPOT_INT_INPUT(initial_adaptivesteps);
  
  
      GETPOT_INT_INPUT(write_interval);
      GETPOT_INPUT(output_xda);
      GETPOT_INPUT(output_xdr);
      GETPOT_INPUT(output_gz);
  #ifndef LIBMESH_HAVE_GZSTREAM
      output_gz                   = false;
  #endif
      GETPOT_INPUT(output_bz2);
  #ifndef LIBMESH_HAVE_BZ2
      output_bz2                  = false;
  #endif
      GETPOT_INPUT(output_gmv);
      GETPOT_INPUT(write_gmv_error);
  #ifndef LIBMESH_HAVE_GMV
      output_gmv                  = false;
      write_gmv_error             = false;
  #endif
      GETPOT_INPUT(output_tecplot);
      GETPOT_INPUT(write_tecplot_error);
  #ifndef LIBMESH_HAVE_TECPLOT_API
      output_tecplot              = false;
      write_tecplot_error         = false;
  #endif
  
  
      GETPOT_REGISTER(system_types);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_system_types =
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;system_types&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_system_types)
        {
          system_types.resize(n_system_types, <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>);
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_system_types; ++i)
            {
              system_types[i] = input(<B><FONT COLOR="#BC8F8F">&quot;system_types&quot;</FONT></B>, system_types[i], i);
            }
        }
  
  
      GETPOT_REGISTER(periodic_boundaries);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_periodic_bcs =
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;periodic_boundaries&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_periodic_bcs)
        {
          <B><FONT COLOR="#A020F0">if</FONT></B> (domaintype != <B><FONT COLOR="#BC8F8F">&quot;square&quot;</FONT></B> &amp;&amp;
              domaintype != <B><FONT COLOR="#BC8F8F">&quot;cylinder&quot;</FONT></B> &amp;&amp;
              domaintype != <B><FONT COLOR="#BC8F8F">&quot;file&quot;</FONT></B> &amp;&amp;
              domaintype != <B><FONT COLOR="#BC8F8F">&quot;od2&quot;</FONT></B>)
            {
              <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Periodic boundaries need rectilinear domains&quot;</FONT></B> &lt;&lt; std::endl;;
              libmesh_error();
            }
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_periodic_bcs; ++i)
            {
              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> myboundary =
                input(<B><FONT COLOR="#BC8F8F">&quot;periodic_boundaries&quot;</FONT></B>, -1, i);
              RealVectorValue translation_vector;
              <B><FONT COLOR="#A020F0">if</FONT></B> (dimension == 2)
                <B><FONT COLOR="#A020F0">switch</FONT></B> (myboundary)
                {
                <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
                  translation_vector = RealVectorValue(0., domain_edge_length);
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
                  translation_vector = RealVectorValue(-domain_edge_width, 0);
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                <B><FONT COLOR="#5F9EA0">default</FONT></B>:
                  <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Unrecognized periodic boundary id &quot;</FONT></B> &lt;&lt;
                                  myboundary &lt;&lt; std::endl;;
                  libmesh_error();
                }
              <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dimension == 3)
                <B><FONT COLOR="#A020F0">switch</FONT></B> (myboundary)
                {
                <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
                  translation_vector = RealVectorValue((Real) 0., (Real) 0., domain_edge_height);
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
                  translation_vector = RealVectorValue((Real) 0., domain_edge_length, (Real) 0.);
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
                  translation_vector = RealVectorValue(-domain_edge_width, (Real) 0., (Real) 0.);
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                <B><FONT COLOR="#5F9EA0">default</FONT></B>:
                  <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Unrecognized periodic boundary id &quot;</FONT></B> &lt;&lt;
                                  myboundary &lt;&lt; std::endl;;
                  libmesh_error();
                }
            }
        }
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string zero_string = <B><FONT COLOR="#BC8F8F">&quot;zero&quot;</FONT></B>;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string empty_string = <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>;
  
      GETPOT_REGISTER(dirichlet_condition_types);
      GETPOT_REGISTER(dirichlet_condition_values);
      GETPOT_REGISTER(dirichlet_condition_boundaries);
      GETPOT_REGISTER(dirichlet_condition_variables);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dirichlet_conditions=
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_types&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_dirichlet_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_values&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_dirichlet_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_values&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition values.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_dirichlet_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_boundaries&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_dirichlet_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_boundaries&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition boundaries.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_dirichlet_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_variables&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_dirichlet_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_variables&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Dirichlet condition variables sets.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_dirichlet_conditions; ++i)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_type =
            input(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_types&quot;</FONT></B>, zero_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_value =
            input(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_values&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type func_boundary =
            input(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_boundaries&quot;</FONT></B>, boundary_id_type(0), i);
  
          dirichlet_conditions[func_boundary] =
            (new_function_base(func_type, func_value).release());
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string variable_set =
            input(<B><FONT COLOR="#BC8F8F">&quot;dirichlet_condition_variables&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != variable_set.size(); ++j)
            {
              <B><FONT COLOR="#A020F0">if</FONT></B> (variable_set[j] == <B><FONT COLOR="#BC8F8F">'1'</FONT></B>)
                dirichlet_condition_variables[func_boundary].push_back(j);
              <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (variable_set[j] != <B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                {
                  <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Unable to understand Dirichlet variable set&quot;</FONT></B>
                               &lt;&lt; variable_set &lt;&lt; std::endl;
                  libmesh_error();
                }
            }
        }
  
      GETPOT_REGISTER(neumann_condition_types);
      GETPOT_REGISTER(neumann_condition_values);
      GETPOT_REGISTER(neumann_condition_boundaries);
      GETPOT_REGISTER(neumann_condition_variables);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_neumann_conditions=
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_types&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_neumann_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_values&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_neumann_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_values&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition values.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_neumann_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_boundaries&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_neumann_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_boundaries&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition boundaries.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_neumann_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_variables&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_neumann_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_variables&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Neumann condition variables sets.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_neumann_conditions; ++i)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_type =
            input(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_types&quot;</FONT></B>, zero_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_value =
            input(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_values&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type func_boundary =
            input(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_boundaries&quot;</FONT></B>, boundary_id_type(0), i);
  
          neumann_conditions[func_boundary] =
            (new_function_base(func_type, func_value).release());
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string variable_set =
            input(<B><FONT COLOR="#BC8F8F">&quot;neumann_condition_variables&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != variable_set.size(); ++j)
            {
              <B><FONT COLOR="#A020F0">if</FONT></B> (variable_set[j] == <B><FONT COLOR="#BC8F8F">'1'</FONT></B>)
                neumann_condition_variables[func_boundary].push_back(j);
              <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (variable_set[j] != <B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                {
                  <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Unable to understand Neumann variable set&quot;</FONT></B>
                               &lt;&lt; variable_set &lt;&lt; std::endl;
                  libmesh_error();
                }
            }
        }
  
      GETPOT_REGISTER(initial_condition_types);
      GETPOT_REGISTER(initial_condition_values);
      GETPOT_REGISTER(initial_condition_subdomains);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_initial_conditions=
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_types&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_initial_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_values&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_initial_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; initial condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_values&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; initial condition values.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_initial_conditions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_subdomains&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_initial_conditions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; initial condition types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_subdomains&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; initial condition subdomains.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_initial_conditions; ++i)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_type =
            input(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_types&quot;</FONT></B>, zero_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_value =
            input(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_values&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> subdomain_id_type func_subdomain =
            input(<B><FONT COLOR="#BC8F8F">&quot;initial_condition_subdomains&quot;</FONT></B>, subdomain_id_type(0), i);
  
          initial_conditions[func_subdomain] =
            (new_function_base(func_type, func_value).release());
        }
  
      GETPOT_REGISTER(other_interior_function_types);
      GETPOT_REGISTER(other_interior_function_values);
      GETPOT_REGISTER(other_interior_function_subdomains);
      GETPOT_REGISTER(other_interior_function_ids);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_other_interior_functions =
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_types&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_interior_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_values&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_interior_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_values&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function values.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_interior_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_subdomains&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_interior_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_subdomains&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function subdomains.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_interior_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_ids&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_interior_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_ids&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other interior function ids.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_other_interior_functions; ++i)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_type =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_types&quot;</FONT></B>, zero_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_value =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_interior_function_values&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> subdomain_id_type func_subdomain =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_interior_condition_subdomains&quot;</FONT></B>, subdomain_id_type(0), i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> func_id =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_interior_condition_ids&quot;</FONT></B>, <B><FONT COLOR="#228B22">int</FONT></B>(0), i);
  
          other_interior_functions[func_id][func_subdomain] =
            (new_function_base(func_type, func_value).release());
        }
  
      GETPOT_REGISTER(other_boundary_function_types);
      GETPOT_REGISTER(other_boundary_function_values);
      GETPOT_REGISTER(other_boundary_function_boundaries);
      GETPOT_REGISTER(other_boundary_function_ids);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_other_boundary_functions =
        input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_types&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_boundary_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_values&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_boundary_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_values&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function values.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_boundary_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_boundaries&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_boundary_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_boundaries&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function boundaries.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (n_other_boundary_functions !=
          input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_ids&quot;</FONT></B>))
        {
          <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; n_other_boundary_functions
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function types does not match &quot;</FONT></B>
                       &lt;&lt; input.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_ids&quot;</FONT></B>)
                       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; other boundary function ids.&quot;</FONT></B> &lt;&lt; std::endl;
  
          libmesh_error();
        }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_other_boundary_functions; ++i)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_type =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_types&quot;</FONT></B>, zero_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> std::string func_value =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_values&quot;</FONT></B>, empty_string, i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type func_boundary =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_boundaries&quot;</FONT></B>, boundary_id_type(0), i);
  
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> func_id =
            input(<B><FONT COLOR="#BC8F8F">&quot;other_boundary_function_ids&quot;</FONT></B>, <B><FONT COLOR="#228B22">int</FONT></B>(0), i);
  
          other_boundary_functions[func_id][func_boundary] =
            (new_function_base(func_type, func_value).release());
        }
  
      GETPOT_INPUT(run_simulation);
      GETPOT_INPUT(run_postprocess);
  
  
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
      GETPOT_INPUT(extra_quadrature_order);
  
  
      GETPOT_INPUT(analytic_jacobians);
      GETPOT_INPUT(verify_analytic_jacobians);
      GETPOT_INPUT(numerical_jacobian_h);
      GETPOT_INPUT(print_solution_norms);
      GETPOT_INPUT(print_solutions);
      GETPOT_INPUT(print_residual_norms);
      GETPOT_INPUT(print_residuals);
      GETPOT_INPUT(print_jacobian_norms);
      GETPOT_INPUT(print_jacobians);
      GETPOT_INPUT(print_element_jacobians);
  
  
      GETPOT_INPUT(use_petsc_snes);
      GETPOT_INPUT(time_solver_quiet);
      GETPOT_INPUT(solver_quiet);
      GETPOT_INPUT(solver_verbose);
      GETPOT_INPUT(require_residual_reduction);
      GETPOT_INPUT(min_step_length);
      GETPOT_INT_INPUT(max_linear_iterations);
      GETPOT_INT_INPUT(max_nonlinear_iterations);
      GETPOT_INPUT(relative_step_tolerance);
      GETPOT_INPUT(relative_residual_tolerance);
      GETPOT_INPUT(initial_linear_tolerance);
      GETPOT_INPUT(minimum_linear_tolerance);
      GETPOT_INPUT(linear_tolerance_multiplier);
  
  
      GETPOT_INT_INPUT(initial_sobolev_order);
      GETPOT_INT_INPUT(initial_extra_quadrature);
      GETPOT_INPUT(indicator_type);
      GETPOT_INT_INPUT(sobolev_order);
  
      GETPOT_INPUT(system_config_file);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; bad_variables =
      input.unidentified_arguments(variable_names);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0 &amp;&amp; !bad_variables.empty())
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Unrecognized variables:&quot;</FONT></B> &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != bad_variables.size(); ++i)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; bad_variables[i] &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Not found among recognized variables:&quot;</FONT></B> &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != variable_names.size(); ++i)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; variable_names[i] &lt;&lt; std::endl;
        libmesh_error();
      }
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file heatsystem.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;heatsystem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/zero_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">void</FONT></B> HeatSystem::init_data ()
  {
    T_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;T&quot;</FONT></B>, static_cast&lt;Order&gt;(_fe_order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(_fe_family));
  
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::ifstream i(<B><FONT COLOR="#BC8F8F">&quot;heat.in&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">if</FONT></B> (!i)
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">'['</FONT></B> &lt;&lt; libMesh::processor_id()
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;] Can't find heat.in; exiting early.&quot;</FONT></B>
                    &lt;&lt; std::endl;
          libmesh_error();
        }
    }
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;heat.in&quot;</FONT></B>);
    _k = infile(<B><FONT COLOR="#BC8F8F">&quot;k&quot;</FONT></B>, 1.0);
    _analytic_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;analytic_jacobians&quot;</FONT></B>, true);
  
    parameters.push_back(_k);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;_k&quot;</FONT></B>) = _k;
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(T_var);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type all_ids[4] = {0, 1, 2, 3};
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; all_bdys(all_ids, all_ids+(2*2));
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; T_only(1, T_var);
  
    ZeroFunction&lt;Number&gt; zero;
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (all_bdys, T_only, &amp;zero));
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> HeatSystem::init_context(DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    c.element_fe_var[0]-&gt;get_JxW();
    c.element_fe_var[0]-&gt;get_dphi();
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (c.is_adjoint())
      {
        <B><FONT COLOR="#228B22">const</FONT></B> System &amp; sys = c.get_system();
  
        NumericVector&lt;Number&gt; &amp;adjoint_solution =
          const_cast&lt;System &amp;&gt;(sys).get_adjoint_solution(0);
  
        c.add_localized_vector(adjoint_solution, sys);
      }
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_context(context);
  }
  
  #define optassert(X) {<B><FONT COLOR="#A020F0">if</FONT></B> (!(X)) libmesh_error();}
  
  <B><FONT COLOR="#228B22">bool</FONT></B> HeatSystem::element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                            DiffContext &amp;context)
  {
    <B><FONT COLOR="#228B22">bool</FONT></B> compute_jacobian = request_jacobian &amp;&amp; _analytic_jacobians;
  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &amp;dphi = c.element_fe_var[0]-&gt;get_dphi();
  
    optassert(c.dof_indices_var.size() &gt; 0);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[0].size();
  
    DenseSubMatrix&lt;Number&gt; &amp;K = *c.elem_subjacobians[0][0];
    DenseSubVector&lt;Number&gt; &amp;F = *c.elem_subresiduals[0];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Gradient grad_T = c.interior_gradient(0, qp);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          F(i) += JxW[qp] * -parameters[0] * (grad_T * dphi[i][qp]);
        <B><FONT COLOR="#A020F0">if</FONT></B> (compute_jacobian)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; ++j)
              K(i,j) += JxW[qp] * -parameters[0] * (dphi[i][qp] * dphi[j][qp]);
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> compute_jacobian;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> HeatSystem::perturb_accumulate_residuals(<B><FONT COLOR="#228B22">const</FONT></B> ParameterVector&amp; parameters_in)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Np = parameters_in.size();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != Np; ++j)
      {
        Number old_parameter = *parameters_in[j];
  
        *parameters_in[j] = old_parameter - dp;
  
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;assembly(true, false);
  
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;rhs-&gt;close();
  
        AutoPtr&lt;NumericVector&lt;Number&gt; &gt; R_minus = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;rhs-&gt;clone();
  
        R_minus_dp += -R_minus-&gt;dot(<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_adjoint_solution(0));
  
        *parameters_in[j] = old_parameter + dp;
  
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;assembly(true, false);
  
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;rhs-&gt;close();
  
        AutoPtr&lt;NumericVector&lt;Number&gt; &gt; R_plus = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;rhs-&gt;clone();
  
        R_plus_dp += -R_plus-&gt;dot(<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_adjoint_solution(0));
  
        *parameters_in[j] = old_parameter;
  
      }
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file initial.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;initial.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> read_initial_parameters()
  {
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> finish_initialization()
  {
  }
  
  
  
  Number initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Real x = p(0), y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> sin(M_PI * x) * sin(M_PI * y);
  }
  
  
  
  Gradient initial_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Real x = p(0), y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> Gradient(M_PI*cos(M_PI * x) * sin(M_PI * y),
    M_PI*sin(M_PI * x) * cos(M_PI * y));
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex5'
***************************************************************
* Running Example adjoints_ex5:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Started /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex5/.libs/lt-example-devel
Building mesh
Building system
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/memory_solution_history.h, line 41, compiled Apr 19 2013 at 11:48:31 ***
Initializing systems
Setting primal initial conditions
|U(0)|= 2.27632

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=1099
  n_elem()=5456
    n_local_elem()=1395
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "HeatSystem"
    Type "Implicit"
    Variables="T" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=1099
    n_constrained_dofs()=256
    n_local_constrained_dofs()=63
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.71976
      Average Off-Processor Bandwidth <= 0.205917
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 256
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Solving time step 0, time = 0
  Nonlinear solver converged, step 0, residual reduction 9.43741e-13 < 1e-09
|U(0.1)|= 2.27183
Advancing timestep

 Solving time step 1, time = 0.1
  Nonlinear solver converged, step 0, residual reduction 9.41233e-13 < 1e-09
|U(0.2)|= 2.26736
Advancing timestep

 Solving time step 2, time = 0.2
  Nonlinear solver converged, step 0, residual reduction 9.36696e-13 < 1e-09
|U(0.3)|= 2.26289
Advancing timestep

 Solving time step 3, time = 0.3
  Nonlinear solver converged, step 0, residual reduction 9.36474e-13 < 1e-09
|U(0.4)|= 2.25843
Advancing timestep

 Solving time step 4, time = 0.4
  Nonlinear solver converged, step 0, residual reduction 9.35902e-13 < 1e-09
|U(0.5)|= 2.25398
Advancing timestep

 Solving time step 5, time = 0.5
  Nonlinear solver converged, step 0, residual reduction 9.33977e-13 < 1e-09
|U(0.6)|= 2.24954
Advancing timestep

 Solving time step 6, time = 0.6
  Nonlinear solver converged, step 0, residual reduction 9.31019e-13 < 1e-09
|U(0.7)|= 2.24511
Advancing timestep

 Solving time step 7, time = 0.7
  Nonlinear solver converged, step 0, residual reduction 9.3018e-13 < 1e-09
|U(0.8)|= 2.24068
Advancing timestep

 Solving time step 8, time = 0.8
  Nonlinear solver converged, step 0, residual reduction 9.2767e-13 < 1e-09
|U(0.9)|= 2.23627
Advancing timestep

 Solving time step 9, time = 0.9
  Nonlinear solver converged, step 0, residual reduction 9.2706e-13 < 1e-09
|U(1)|= 2.23186
Advancing timestep


Solving the adjoint problem
Setting adjoint initial conditions Z(1)
Retrieving solutions at time t=1
|U(1.1)|= 2.23186
|U(1)|= 2.23186
|Z(1)|= 2.27632

 Solving adjoint time step 0, time = 1
Retrieving solutions at time t=1
|U(1)|= 2.23186
|U(0.9)|= 2.23627
|Z(0.9)|= 2.27183

 Solving adjoint time step 1, time = 0.9
Retrieving solutions at time t=0.9
|U(0.9)|= 2.23627
|U(0.8)|= 2.24068
|Z(0.8)|= 2.26736

 Solving adjoint time step 2, time = 0.8
Retrieving solutions at time t=0.8
|U(0.8)|= 2.24068
|U(0.7)|= 2.24511
|Z(0.7)|= 2.26289

 Solving adjoint time step 3, time = 0.7
Retrieving solutions at time t=0.7
|U(0.7)|= 2.24511
|U(0.6)|= 2.24954
|Z(0.6)|= 2.25843

 Solving adjoint time step 4, time = 0.6
Retrieving solutions at time t=0.6
|U(0.6)|= 2.24954
|U(0.5)|= 2.25398
|Z(0.5)|= 2.25398

 Solving adjoint time step 5, time = 0.5
Retrieving solutions at time t=0.5
|U(0.5)|= 2.25398
|U(0.4)|= 2.25843
|Z(0.4)|= 2.24954

 Solving adjoint time step 6, time = 0.4
Retrieving solutions at time t=0.4
|U(0.4)|= 2.25843
|U(0.3)|= 2.26289
|Z(0.3)|= 2.24511

 Solving adjoint time step 7, time = 0.3
Retrieving solutions at time t=0.3
|U(0.3)|= 2.26289
|U(0.2)|= 2.26736
|Z(0.2)|= 2.24068

 Solving adjoint time step 8, time = 0.2
Retrieving solutions at time t=0.2
|U(0.2)|= 2.26736
|U(0.1)|= 2.27183
|Z(0.1)|= 2.23627

 Solving adjoint time step 9, time = 0.1
Retrieving solutions at time t=0.1
|U(0.1)|= 2.27183
|U(2.77556e-17)|= 2.27632
|Z(2.77556e-17)|= 2.23186

Retrieving 0, time = 2.77556e-17
|U(0.1)|= 2.27183
|U(2.77556e-17)|= 2.27632
|Z(2.77556e-17)|= 2.23186

Retrieving 1, time = 0.1
|U(0.2)|= 2.26736
|U(0.1)|= 2.27183
|Z(0.1)|= 2.23627

Retrieving 2, time = 0.2
|U(0.3)|= 2.26289
|U(0.2)|= 2.26736
|Z(0.2)|= 2.24068

Retrieving 3, time = 0.3
|U(0.4)|= 2.25843
|U(0.3)|= 2.26289
|Z(0.3)|= 2.24511

Retrieving 4, time = 0.4
|U(0.5)|= 2.25398
|U(0.4)|= 2.25843
|Z(0.4)|= 2.24954

Retrieving 5, time = 0.5
|U(0.6)|= 2.24954
|U(0.5)|= 2.25398
|Z(0.5)|= 2.25398

Retrieving 6, time = 0.6
|U(0.7)|= 2.24511
|U(0.6)|= 2.24954
|Z(0.6)|= 2.25843

Retrieving 7, time = 0.7
|U(0.8)|= 2.24068
|U(0.7)|= 2.24511
|Z(0.7)|= 2.26289

Retrieving 8, time = 0.8
|U(0.9)|= 2.23627
|U(0.8)|= 2.24068
|Z(0.8)|= 2.26736

Retrieving 9, time = 0.9
|U(1)|= 2.23186
|U(0.9)|= 2.23627
|Z(0.9)|= 2.27183

Retrieving  final time = 1
|U(1.1)|= 2.23186
|U(1)|= 2.23186
|Z(1)|= 2.27632

Sensitivity of QoI 0 w.r.t parameter 0 is: -5.30953

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:49:06 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=9.21305, Active time=7.81359                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     5         0.0050      0.001002    0.0063      0.001253    0.06     0.08     |
|   build_constraint_matrix()        63488     0.0654      0.000001    0.0654      0.000001    0.84     0.84     |
|   build_sparsity()                 5         0.0039      0.000784    0.0141      0.002822    0.05     0.18     |
|   cnstrn_elem_mat_vec()            10240     0.0155      0.000002    0.0155      0.000002    0.20     0.20     |
|   constrain_elem_matrix()          10240     0.0108      0.000001    0.0108      0.000001    0.14     0.14     |
|   constrain_elem_vector()          43008     0.0398      0.000001    0.0398      0.000001    0.51     0.51     |
|   create_dof_constraints()         5         0.0192      0.003849    0.0454      0.009088    0.25     0.58     |
|   distribute_dofs()                5         0.0108      0.002164    0.0412      0.008248    0.14     0.53     |
|   dof_indices()                    285053    0.8995      0.000003    0.8995      0.000003    11.51    11.51    |
|   enforce_constraints_exactly()    50        0.0119      0.000239    0.0119      0.000239    0.15     0.15     |
|   old_dof_indices()                5440      0.0190      0.000003    0.0190      0.000003    0.24     0.24     |
|   prepare_send_list()              5         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   reinit()                         5         0.0164      0.003283    0.0164      0.003283    0.21     0.21     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          22        0.0438      0.001993    0.1992      0.009056    0.56     2.55     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        425264    1.1779      0.000003    1.1779      0.000003    15.07    15.07    |
|   init_shape_functions()           6649      0.0265      0.000004    0.0265      0.000004    0.34     0.34     |
|   inverse_map()                    15890     0.0402      0.000003    0.0402      0.000003    0.51     0.51     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assemble_qoi_derivative()        10        0.6396      0.063960    0.8595      0.085948    8.19     11.00    |
|   assembly()                       10        0.2145      0.021448    0.6289      0.062888    2.74     8.05     |
|   assembly(get_jacobian)           10        0.2191      0.021910    0.6592      0.065919    2.80     8.44     |
|   assembly(get_residual)           32        0.5898      0.018432    1.9361      0.060502    7.55     24.78    |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             425264    0.6627      0.000002    0.6627      0.000002    8.48     8.48     |
|   compute_face_map()               6448      0.0264      0.000004    0.0601      0.000009    0.34     0.77     |
|   init_face_shape_functions()      104       0.0003      0.000003    0.0003      0.000003    0.00     0.00     |
|   init_reference_to_physical_map() 6649      0.0206      0.000003    0.0206      0.000003    0.26     0.26     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               22        0.3610      0.016410    0.3627      0.016485    4.62     4.64     |
|                                                                                                                |
| ImplicitSystem                                                                                                 |
|   adjoint_solve()                  10        0.0006      0.000057    1.8579      0.185793    0.01     23.78    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           16320     0.0202      0.000001    0.0202      0.000001    0.26     0.26     |
|   init()                           8         0.0025      0.000312    0.0025      0.000312    0.03     0.03     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       4         0.0008      0.000209    0.0016      0.000400    0.01     0.02     |
|   find_neighbors()                 5         0.0281      0.005616    0.0285      0.005697    0.36     0.36     |
|   renumber_nodes_and_elem()        14        0.0027      0.000194    0.0027      0.000194    0.03     0.03     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        6         0.0228      0.003807    0.0228      0.003807    0.29     0.29     |
|   find_global_indices()            6         0.0034      0.000563    0.0284      0.004729    0.04     0.36     |
|   parallel_sort()                  6         0.0009      0.000146    0.0014      0.000225    0.01     0.02     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         22        0.0006      0.000029    0.5641      0.025641    0.01     7.22     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              4         0.0007      0.000182    0.0011      0.000285    0.01     0.01     |
|   _refine_elements()               8         0.0298      0.003724    0.0904      0.011300    0.38     1.16     |
|   add_point()                      16320     0.0305      0.000002    0.0526      0.000003    0.39     0.67     |
|   make_coarsening_compatible()     8         0.0098      0.001230    0.0098      0.001230    0.13     0.13     |
|   make_flags_parallel_consistent() 8         0.0070      0.000880    0.0132      0.001645    0.09     0.17     |
|   make_refinement_compatible()     8         0.0000      0.000002    0.0000      0.000006    0.00     0.00     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000091    0.0001      0.000091    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      5         0.0669      0.013372    0.0951      0.019025    0.86     1.22     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          10        0.9284      0.092839    2.2089      0.220890    11.88    28.27    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      29        0.0098      0.000338    0.0103      0.000356    0.13     0.13     |
|   max(bool)                        21        0.0077      0.000367    0.0077      0.000367    0.10     0.10     |
|   max(scalar)                      1736      0.0075      0.000004    0.0075      0.000004    0.10     0.10     |
|   max(vector)                      414       0.0027      0.000006    0.0078      0.000019    0.03     0.10     |
|   min(bool)                        2154      0.0150      0.000007    0.0150      0.000007    0.19     0.19     |
|   min(scalar)                      1722      0.1287      0.000075    0.1287      0.000075    1.65     1.65     |
|   min(vector)                      414       0.0031      0.000008    0.0083      0.000020    0.04     0.11     |
|   probe()                          252       0.0073      0.000029    0.0073      0.000029    0.09     0.09     |
|   receive()                        252       0.0008      0.000003    0.0081      0.000032    0.01     0.10     |
|   send()                           252       0.0005      0.000002    0.0005      0.000002    0.01     0.01     |
|   send_receive()                   264       0.0011      0.000004    0.0101      0.000038    0.01     0.13     |
|   sum()                            222       0.5955      0.002683    0.6735      0.003034    7.62     8.62     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           252       0.0003      0.000001    0.0003      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         5         0.0046      0.000926    0.0151      0.003029    0.06     0.19     |
|   set_parent_processor_ids()       5         0.0022      0.000431    0.0022      0.000431    0.03     0.03     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          20        0.3934      0.019670    0.3934      0.019670    5.03     5.03     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       8         0.0097      0.001212    0.0375      0.004688    0.12     0.48     |
|                                                                                                                |
| System                                                                                                         |
|   calculate_norm()                 77        0.2823      0.003666    1.4223      0.018472    3.61     18.20    |
|   project_vector()                 10        0.0457      0.004569    0.1026      0.010262    0.58     1.31     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1344805   7.8136                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adjoints_ex5:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex5'
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
