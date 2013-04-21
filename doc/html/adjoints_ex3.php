<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adjoints_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file coupled_system.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __coupled_system_h__
        #define __coupled_system_h__
        
</pre>
</div>
<div class = "comment">
DiffSystem framework files
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fem_function_base.h"
        #include "libmesh/fem_system.h"
        #include "libmesh/libmesh_common.h"
        #include "libmesh/parameter_vector.h"
        
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
The Navier-Stokes system class.
FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
but we must specify element residuals
</div>

<div class ="fragment">
<pre>
        class CoupledSystem : public FEMSystem
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          CoupledSystem(EquationSystems& es,
                       const std::string& name_in,
                       const unsigned int number_in)
            : FEMSystem(es, name_in, number_in), Peclet(1.) {qoi.resize(1);}
        
</pre>
</div>
<div class = "comment">
Function to get computed QoI values


<br><br></div>

<div class ="fragment">
<pre>
          Number &get_QoI_value()
            {
              return computed_QoI;
            }
        
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
        
          Real &get_Pe()
            {
              return Peclet;
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
          virtual void init_context(DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Element residual and jacobian calculations
Time dependent parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_time_derivative (bool request_jacobian,
                                                DiffContext& context);
        
</pre>
</div>
<div class = "comment">
Constraint parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_constraint (bool request_jacobian,
                                           DiffContext& context);
        
</pre>
</div>
<div class = "comment">
Postprocessed output
</div>

<div class ="fragment">
<pre>
          virtual void postprocess ();
        
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
Indices for each variable;
</div>

<div class ="fragment">
<pre>
          unsigned int p_var, u_var, v_var, C_var;
        
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
The Peclet number for the species transport
</div>

<div class ="fragment">
<pre>
          Real Peclet;
        
</pre>
</div>
<div class = "comment">
The functionals to be computed as QoIs
</div>

<div class ="fragment">
<pre>
          Number computed_QoI;
        
        };
        
        
        class CoupledFEMFunctionsx : public FEMFunctionBase&lt;Number&gt;
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          CoupledFEMFunctionsx(System& /* sys */, unsigned int var_number) {var = var_number;}
        
</pre>
</div>
<div class = "comment">
Destructor
</div>

<div class ="fragment">
<pre>
          virtual ~CoupledFEMFunctionsx () {}
        
          virtual AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt; clone () const
          {return AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt;( new CoupledFEMFunctionsx(*this) ); }
        
          virtual void operator() (const FEMContext&, const Point&,
        			   const Real, DenseVector&lt;Number&gt;&)
          {libmesh_error();}
        
          virtual Number operator() (const FEMContext&, const Point& p,
        			     const Real time = 0.);
        
         private:
        
          unsigned int var;
        
        };
        
        
        class CoupledFEMFunctionsy : public FEMFunctionBase&lt;Number&gt;
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          CoupledFEMFunctionsy(System& /* sys */, unsigned int var_number) {var = var_number;}
        
</pre>
</div>
<div class = "comment">
Destructor
</div>

<div class ="fragment">
<pre>
          virtual ~CoupledFEMFunctionsy () {}
        
          virtual AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt; clone () const
          {return AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt;( new CoupledFEMFunctionsy(*this) ); }
        
          virtual void operator() (const FEMContext&, const Point&,
        			   const Real,
        			   DenseVector&lt;Number&gt;&)
          {libmesh_error();}
        
          virtual Number operator() (const FEMContext&, const Point& p,
        			     const Real time = 0.);
        
         private:
        
          unsigned int var;
        
        };
        
        #endif //__coupled_system_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file domain.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "femparameters.h"
        #include "libmesh/mesh.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        class FEMParameters;
        
        void build_domain (Mesh &mesh, FEMParameters &param);
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
        
        #include &lt;limits&gt;
        
        #include "libmesh/libmesh_common.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/enum_norm_type.h"
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
        	fine_mesh_file_primal("fine_mesh.xda"), fine_mesh_soln_primal("fine_mesh_soln.xda"),
        	fine_mesh_file_adjoint("fine_mesh.xda"), fine_mesh_soln_adjoint("fine_mesh_soln.xda"),
              elementorder(2),
              domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
              domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
              coarsegridx(1), coarsegridy(1), coarsegridz(1),
        	coarserefinements(0), coarsecoarsenings(0), extrarefinements(0),
              use_petsc_snes(false),
              time_solver_quiet(true), solver_quiet(true),
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
              initial_adaptivesteps(0), initial_sobolev_order(1),
        	initial_extra_quadrature(0),
        	indicator_type("kelly"), patch_reuse(true), sobolev_order(1), adjoint_residual_type("patchpatch"), alternate_with_uniform_steps("false"), alternate_step_number(2),
        	component_wise_error("false"), compare_to_fine_solution_primal(false), compare_to_fine_solution_adjoint(false), compute_sensitivities(false), do_forward_sensitivity(true), do_adjoint_sensitivity(false),
        	postprocess_adjoint(false),
              write_interval(10),
              write_gmv_error(false), write_tecplot_error(false),
              output_xda(false), output_xdr(false),
              output_bz2(true), output_gz(true),
              output_gmv(false), output_tecplot(false),
              run_simulation(true), run_postprocess(false),
              fe_family(1, "LAGRANGE"), fe_order(1, 1),
              extra_quadrature_order(0),
              analytic_jacobians(true), verify_analytic_jacobians(0.0),
              print_solution_norms(false), print_solutions(false),
              print_residual_norms(false), print_residuals(false),
              print_jacobian_norms(false), print_jacobians(false) {}
        
            void read(GetPot &input);
        
        
        
            unsigned int initial_timestep, n_timesteps;
            bool transient;
            unsigned int deltat_reductions;
            std::string timesolver_core;
            Real end_time, deltat, timesolver_theta,
        	 timesolver_maxgrowth, timesolver_tolerance,
        	 timesolver_upper_tolerance, steadystate_tolerance;
            std::vector&lt;FEMNormType&gt; timesolver_norm;
        
            unsigned int dimension;
            std::string domaintype, domainfile, elementtype, fine_mesh_file_primal, fine_mesh_soln_primal, fine_mesh_file_adjoint, fine_mesh_soln_adjoint;
            Real elementorder;
            Real domain_xmin, domain_ymin, domain_zmin;
            Real domain_edge_width, domain_edge_length, domain_edge_height;
            unsigned int coarsegridx, coarsegridy, coarsegridz;
            unsigned int coarserefinements, coarsecoarsenings, extrarefinements;
        
            bool use_petsc_snes;
            bool time_solver_quiet, solver_quiet, reuse_preconditioner, require_residual_reduction;
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
            unsigned int initial_adaptivesteps;
            unsigned int initial_sobolev_order;
            unsigned int initial_extra_quadrature;
        
            std::string indicator_type;
            bool patch_reuse;
            unsigned int sobolev_order;
            std::string adjoint_residual_type;
            bool alternate_with_uniform_steps;
            unsigned int alternate_step_number;
            bool component_wise_error;
            bool compare_to_fine_solution_primal, compare_to_fine_solution_adjoint;
        
            bool compute_sensitivities;
            bool do_forward_sensitivity;
            bool do_adjoint_sensitivity;
        
            bool postprocess_adjoint;
        
            unsigned int write_interval;
            bool write_gmv_error, write_tecplot_error, output_xda, output_xdr,
        	 output_bz2, output_gz, output_gmv, output_tecplot;
            bool run_simulation, run_postprocess;
        
            std::vector&lt;std::string&gt; fe_family;
            std::vector&lt;unsigned int&gt; fe_order;
            int extra_quadrature_order;
        
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
<br><br><br> <h1> The source file H-qoi.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef H_QOI_H
        #define H_QOI_H
        
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
        
        class CoupledSystemQoI : public DifferentiableQoI
        {
        public:
          CoupledSystemQoI(){}
          virtual ~CoupledSystemQoI(){}
        
          virtual void init_qoi( std::vector&lt;Number&gt;& sys_qoi);
          virtual void postprocess( ){}
        
          virtual void side_qoi_derivative(DiffContext &context, const QoISet & qois);
        
          virtual void side_qoi(DiffContext &context, const QoISet & qois);
        
          virtual AutoPtr&lt;DifferentiableQoI&gt; clone( )
          { AutoPtr&lt;DifferentiableQoI&gt; my_clone( new CoupledSystemQoI );
            *my_clone = *this;
            return my_clone;
          }
        
        };
        #endif // H_QOI_H
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
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
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
<br><br><br> <h1> The source file mysystems.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "coupled_system.h"
        #include "femparameters.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        FEMSystem &build_system(EquationSystems &es, GetPot &, FEMParameters& /* param */)
        {
          CoupledSystem &system = es.add_system&lt;CoupledSystem&gt; ("CoupledSystem");
        
          return system;
        }
        
        Number exact_value(const Point& /* p */,       // xyz location
                           const Parameters& /* param */,  // EquationSystem parameters
                           const std::string&, // sys_name
                           const std::string&) // unknown_name
        {
          std::cout&lt;&lt;"Warning ! No exact value function specified ! Returning 0." &lt;&lt; std::endl;
        
          return 0;
        }
        
        Gradient exact_grad(const Point& /* p */,       // xyz location
                            const Parameters& /* param */,  // EquationSystems parameters
                            const std::string&, // sys_name
                            const std::string&) // unknown_name
        {
          std::cout&lt;&lt;"Warning ! No exact value function specified ! Returning 0." &lt;&lt; std::endl;
        
          return 0;
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file output.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include &lt;string&gt;
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        void start_output(unsigned int timesteps,
                          std::string filename,
                          std::string varname);
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file adjoints_ex3.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;sys/time.h&gt;
        #include &lt;iomanip&gt;
        
</pre>
</div>
<div class = "comment">
Libmesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        
        #include "libmesh/twostep_time_solver.h"
        #include "libmesh/euler_solver.h"
        #include "libmesh/euler2_solver.h"
        #include "libmesh/steady_solver.h"
        
        #include "libmesh/newton_solver.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/petsc_diff_solver.h"
        
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_tools.h"
        #include "libmesh/mesh_base.h"
        #include "libmesh/mesh_refinement.h"
        
        #include "libmesh/system.h"
        #include "libmesh/system_norm.h"
        
        #include "libmesh/adjoint_residual_error_estimator.h"
        #include "libmesh/const_fem_function.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/fem_function_base.h"
        #include "libmesh/getpot.h"
        #include "libmesh/gmv_io.h"
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/parameter_vector.h"
        #include "libmesh/patch_recovery_error_estimator.h"
        #include "libmesh/petsc_vector.h"
        #include "libmesh/sensitivity_data.h"
        #include "libmesh/tecplot_io.h"
        #include "libmesh/uniform_refinement_estimator.h"
        #include "libmesh/qoi_set.h"
        #include "libmesh/weighted_patch_recovery_error_estimator.h"
        
</pre>
</div>
<div class = "comment">
Local includes
</div>

<div class ="fragment">
<pre>
        #include "coupled_system.h"
        #include "domain.h"
        #include "initial.h"
        #include "femparameters.h"
        #include "mysystems.h"
        #include "output.h"
        #include "H-qoi.h"
        
</pre>
</div>
<div class = "comment">
We solve a coupled Stokes + Convection Diffusion system in an H channel geometry
with 2 inlets and 2 outlets. The QoI is the species
flux from left outlet


<br><br>Channel Geometry:
Wall
----------------------------------------------------------------
0                                                                1
----------------------------     -------------------------------
|     |
Wall |     | Wall
|     |
-----------------------------     -------------------------------
2                                                               2
-----------------------------------------------------------------
Wall


<br><br>The equations governing this flow are:
Stokes: -VectorLaplacian(velocity) + grad(pressure) = vector(0)
Convection-Diffusion: - dot(velocity, grad(concentration) ) + Laplacian(concentration) = 0


<br><br>The boundary conditions are:
u_1(0) = -(y-2)*(y-3), u_2(0) = 0 ; u_1(1) = (y-2)*(y-3), u_2(1) = 0 ;
u_1(walls) = 0, u_2(walls) = 0;
C(0) = 1 ; C(1) = 0;
grad(C) dot n (walls) = 0 ;
grad(C) dot n (2) = 0 ; grad(C) dot n (3) = 0


<br><br>The QoI is:
Q((u,v), C) = integral_{left_outlet}  - u * C ds


<br><br>The complete equal order adjoint QoI error estimate is: (Notation for derivatives: grad(C) = C,1 + C,2)
Q(u) - Q(u_h) \leq
|e(u_1)|_{H1} |e(u_1^*)|_{H1} +  |e(u_2)|_{H1} |e(u_2^*)|_{H1} +  (1/Pe) * |e(C)|_{H1} |e(C^*)|_{H1}  +
||e(u_1,1)||_{L2} ||e(p^*)||_{L2} + ||e(u_2,2)||_{L2} ||e(p^*)||_{L2} + ||e(u_1,1^*)||_{L2} ||e(p)||_{L2} + ||e(u_2,2^*)||_{L2} ||e(p)||_{L2}  +
||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} + ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2}
||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} + ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2}
= error_non_pressure + error_with_pressure + error_convection_diffusion_x + error_convection_diffusion_y


<br><br>Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Number output files


<br><br></div>

<div class ="fragment">
<pre>
        std::string numbered_filename(unsigned int t_step, // The timestep count
                                      unsigned int a_step, // The adaptive step count
                                      std::string solution_type, // primal or adjoint solve
                                      std::string type,
                                      std::string extension,
                                      FEMParameters &param)
        {
          std::ostringstream file_name;
          file_name &lt;&lt; solution_type
                    &lt;&lt; ".out."
                    &lt;&lt; type
                    &lt;&lt; '.'
                    &lt;&lt; std::setw(3)
                    &lt;&lt; std::setfill('0')
                    &lt;&lt; std::right
                    &lt;&lt; t_step
                    &lt;&lt; '.'
                    &lt;&lt; std::setw(2)
                    &lt;&lt; std::setfill('0')
                    &lt;&lt; std::right
                    &lt;&lt; a_step
                    &lt;&lt; '.'
                    &lt;&lt; extension;
        
          if (param.output_bz2)
            file_name &lt;&lt; ".bz2";
          else if (param.output_gz)
            file_name &lt;&lt; ".gz";
          return file_name.str();
        }
        
</pre>
</div>
<div class = "comment">
Write tecplot, gmv and xda/xdr output


<br><br></div>

<div class ="fragment">
<pre>
        void write_output(EquationSystems &es,
                          unsigned int t_step, // The timestep count
                          unsigned int a_step, // The adaptive step count
                          std::string solution_type, // primal or adjoint solve
                          FEMParameters &param)
        {
          MeshBase &mesh = es.get_mesh();
        
        #ifdef LIBMESH_HAVE_GMV
          if (param.output_gmv)
            {
              std::ostringstream file_name_gmv;
              file_name_gmv &lt;&lt; solution_type
                            &lt;&lt; ".out.gmv."
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_step
                            &lt;&lt; '.'
                            &lt;&lt; std::setw(2)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; a_step;
        
              GMVIO(mesh).write_equation_systems
                (file_name_gmv.str(), es);
            }
        #endif
        
        #ifdef LIBMESH_HAVE_TECPLOT_API
          if (param.output_tecplot)
            {
              std::ostringstream file_name_tecplot;
              file_name_tecplot &lt;&lt; solution_type
                                &lt;&lt; ".out."
                                &lt;&lt; std::setw(3)
                                &lt;&lt; std::setfill('0')
                                &lt;&lt; std::right
                                &lt;&lt; t_step
                                &lt;&lt; '.'
                                &lt;&lt; std::setw(2)
                                &lt;&lt; std::setfill('0')
                                &lt;&lt; std::right
                                &lt;&lt; a_step
                                &lt;&lt; ".plt";
        
              TecplotIO(mesh).write_equation_systems
                (file_name_tecplot.str(), es);
            }
        #endif
        
          if (param.output_xda || param.output_xdr)
            mesh.renumber_nodes_and_elements();
          if (param.output_xda)
            {
              mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param));
              es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xda", param),
                       libMeshEnums::WRITE, EquationSystems::WRITE_DATA |
                       EquationSystems::WRITE_ADDITIONAL_DATA);
            }
          if (param.output_xdr)
            {
              mesh.write(numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param));
              es.write(numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param),
                       libMeshEnums::WRITE, EquationSystems::WRITE_DATA |
                       EquationSystems::WRITE_ADDITIONAL_DATA);
            }
        }
        
        void write_output_headers(FEMParameters &param)
        {
</pre>
</div>
<div class = "comment">
Only one processor needs to take care of headers.
</div>

<div class ="fragment">
<pre>
          if (libMesh::processor_id() != 0)
            return;
        
          start_output(param.initial_timestep, "out_clocktime.m", "vector_clocktime");
        
          if (param.run_simulation)
            {
              start_output(param.initial_timestep, "out_activemesh.m", "table_activemesh");
              start_output(param.initial_timestep, "out_solvesteps.m", "table_solvesteps");
        
              if (param.timesolver_tolerance)
                {
                  start_output(param.initial_timestep, "out_time.m", "vector_time");
                  start_output(param.initial_timestep, "out_timesteps.m", "vector_timesteps");
                }
              if (param.steadystate_tolerance)
                start_output(param.initial_timestep, "out_changerate.m", "vector_changerate");
            }
        }
        
        void write_output_solvedata(EquationSystems &es,
                                    unsigned int a_step,
                                    unsigned int newton_steps,
                                    unsigned int krylov_steps,
                                    unsigned int tv_sec,
                                    unsigned int tv_usec)
        {
          MeshBase &mesh = es.get_mesh();
          unsigned int n_active_elem = mesh.n_active_elem();
          unsigned int n_active_dofs = es.n_active_dofs();
        
          if (libMesh::processor_id() == 0)
            {
</pre>
</div>
<div class = "comment">
Write out the number of elements/dofs used
</div>

<div class ="fragment">
<pre>
              std::ofstream activemesh ("out_activemesh.m",
                std::ios_base::app | std::ios_base::out);
              activemesh.precision(17);
              activemesh &lt;&lt; (a_step + 1) &lt;&lt; ' '
                         &lt;&lt; n_active_elem &lt;&lt; ' '
                         &lt;&lt; n_active_dofs &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Write out the number of solver steps used
</div>

<div class ="fragment">
<pre>
              std::ofstream solvesteps ("out_solvesteps.m",
                std::ios_base::app | std::ios_base::out);
              solvesteps.precision(17);
              solvesteps &lt;&lt; newton_steps &lt;&lt; ' '
                         &lt;&lt; krylov_steps &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Write out the clock time
</div>

<div class ="fragment">
<pre>
              std::ofstream clocktime ("out_clocktime.m",
                std::ios_base::app | std::ios_base::out);
              clocktime.precision(17);
              clocktime &lt;&lt; tv_sec &lt;&lt; '.' &lt;&lt; tv_usec &lt;&lt; std::endl;
            }
        }
        
        void write_output_footers(FEMParameters &param)
        {
</pre>
</div>
<div class = "comment">
Write footers on output .m files
</div>

<div class ="fragment">
<pre>
          if (libMesh::processor_id() == 0)
            {
              std::ofstream clocktime ("out_clocktime.m",
                std::ios_base::app | std::ios_base::out);
              clocktime &lt;&lt; "];" &lt;&lt; std::endl;
        
              if (param.run_simulation)
                {
                  std::ofstream activemesh ("out_activemesh.m",
                    std::ios_base::app | std::ios_base::out);
                  activemesh &lt;&lt; "];" &lt;&lt; std::endl;
        
                  std::ofstream solvesteps ("out_solvesteps.m",
                    std::ios_base::app | std::ios_base::out);
                  solvesteps &lt;&lt; "];" &lt;&lt; std::endl;
        
                  if (param.timesolver_tolerance)
                    {
                      std::ofstream times ("out_time.m",
                        std::ios_base::app | std::ios_base::out);
                      times &lt;&lt; "];" &lt;&lt; std::endl;
                      std::ofstream timesteps ("out_timesteps.m",
                        std::ios_base::app | std::ios_base::out);
                      timesteps &lt;&lt; "];" &lt;&lt; std::endl;
                    }
                  if (param.steadystate_tolerance)
                    {
                      std::ofstream changerate ("out_changerate.m",
                        std::ios_base::app | std::ios_base::out);
                      changerate &lt;&lt; "];" &lt;&lt; std::endl;
                    }
                }
        
</pre>
</div>
<div class = "comment">
We're done, so write out a file (for e.g. Make to check)
</div>

<div class ="fragment">
<pre>
              std::ofstream complete ("complete");
              complete &lt;&lt; "complete" &lt;&lt; std::endl;
            }
        }
        
        #if defined(LIBMESH_HAVE_GMV) || defined(LIBMESH_HAVE_TECPLOT_API)
        void write_error(EquationSystems &es,
                         ErrorVector &error,
                         unsigned int t_number,
                         unsigned int a_number,
                         FEMParameters &param,
                         std::string error_type)
        #else
        void write_error(EquationSystems &,
                         ErrorVector &,
                         unsigned int,
                         unsigned int,
                         FEMParameters &,
                         std::string)
        #endif
        {
        #ifdef LIBMESH_HAVE_GMV
          if (param.write_gmv_error)
            {
              std::ostringstream error_gmv;
              error_gmv &lt;&lt; "error.gmv."
                        &lt;&lt; std::setw(3)
                        &lt;&lt; std::setfill('0')
                        &lt;&lt; std::right
                        &lt;&lt; a_number
                        &lt;&lt; "."
                        &lt;&lt; std::setw(2)
                        &lt;&lt; std::setfill('0')
                        &lt;&lt; std::right
                        &lt;&lt; t_number
                        &lt;&lt; error_type;
        
              error.plot_error(error_gmv.str(), es.get_mesh());
            }
        #endif
        
        #ifdef LIBMESH_HAVE_TECPLOT_API
          if (param.write_tecplot_error)
            {
              std::ostringstream error_tecplot;
              error_tecplot &lt;&lt; "error.plt."
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; a_number
                            &lt;&lt; "."
                            &lt;&lt; std::setw(2)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_number
                            &lt;&lt; error_type;
        
              error.plot_error(error_tecplot.str(), es.get_mesh());
            }
        #endif
        }
        
        void read_output(EquationSystems &es,
                         unsigned int t_step,
                         unsigned int a_step,
                         std::string solution_type,
                         FEMParameters &param)
        {
          MeshBase &mesh = es.get_mesh();
        
          std::string file_name_mesh, file_name_soln;
</pre>
</div>
<div class = "comment">
Look for ASCII files first
</div>

<div class ="fragment">
<pre>
          if (param.output_xda)
            {
              file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xda", param);
              file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xda", param);
            }
          else if (param.output_xdr)
            {
              file_name_mesh = numbered_filename(t_step, a_step, solution_type, "mesh", "xdr", param);
              file_name_soln = numbered_filename(t_step, a_step, solution_type, "soln", "xdr", param);
            }
        
</pre>
</div>
<div class = "comment">
Read in the mesh
</div>

<div class ="fragment">
<pre>
          mesh.read(file_name_mesh);
        
</pre>
</div>
<div class = "comment">
And the stored solution
</div>

<div class ="fragment">
<pre>
          es.read(file_name_soln, libMeshEnums::READ,
                  EquationSystems::READ_HEADER |
                  EquationSystems::READ_DATA |
                  EquationSystems::READ_ADDITIONAL_DATA);
        
</pre>
</div>
<div class = "comment">
Put systems in a consistent state
</div>

<div class ="fragment">
<pre>
          for (unsigned int i = 0; i != es.n_systems(); ++i)
            es.get_system&lt;FEMSystem&gt;(i).update();
        
</pre>
</div>
<div class = "comment">
Figure out the current time
</div>

<div class ="fragment">
<pre>
          Real current_time = 0., current_timestep = 0.;
        
          if (param.timesolver_tolerance)
            {
              std::ifstream times ("out_time.m");
              std::ifstream timesteps ("out_timesteps.m");
              if (times.is_open() && timesteps.is_open())
                {
</pre>
</div>
<div class = "comment">
Read headers
</div>

<div class ="fragment">
<pre>
                  const unsigned int headersize = 25;
                  char header[headersize];
                  timesteps.getline (header, headersize);
                  if (strcmp(header, "vector_timesteps = [") != 0)
                    {
                      std::cout &lt;&lt; "Bad header in out_timesteps.m:" &lt;&lt; std::endl
                                &lt;&lt; header
                                &lt;&lt; std::endl;
                      libmesh_error();
                    }
        
                  times.getline (header, headersize);
                  if (strcmp(header, "vector_time = [") != 0)
                    {
                      std::cout &lt;&lt; "Bad header in out_time.m:" &lt;&lt; std::endl
                                &lt;&lt; header
                                &lt;&lt; std::endl;
                      libmesh_error();
                    }
        
</pre>
</div>
<div class = "comment">
Read each timestep
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i = 0; i != t_step; ++i)
                    {
                      if (!times.good())
                        libmesh_error();
                      times &gt;&gt; current_time;
                      timesteps &gt;&gt; current_timestep;
                    }
</pre>
</div>
<div class = "comment">
Remember to increment the last timestep; out_times.m
lists each *start* time
</div>

<div class ="fragment">
<pre>
                  current_time += current_timestep;
                }
              else
                libmesh_error();
            }
          else
            current_time = t_step * param.deltat;
        
          for (unsigned int i = 0; i != es.n_systems(); ++i)
            es.get_system&lt;FEMSystem&gt;(i).time = current_time;
        }
        
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
              if (param.timesolver_core == "euler2")
                {
                  Euler2Solver *euler2solver =
                    new Euler2Solver(system);
        
                  euler2solver-&gt;theta = param.timesolver_theta;
                  innersolver = euler2solver;
                }
              else if (param.timesolver_core == "euler")
                {
                  EulerSolver *eulersolver =
                    new EulerSolver(system);
        
                  eulersolver-&gt;theta = param.timesolver_theta;
                  innersolver = eulersolver;
                }
              else
                {
                  std::cerr &lt;&lt; "Don't recognize core TimeSolver type: "
                            &lt;&lt; param.timesolver_core &lt;&lt; std::endl;
                  libmesh_error();
                }
        
              if (param.timesolver_tolerance)
                {
                  TwostepTimeSolver *timesolver =
                    new TwostepTimeSolver(system);
        
                  timesolver-&gt;max_growth       = param.timesolver_maxgrowth;
                  timesolver-&gt;target_tolerance = param.timesolver_tolerance;
                  timesolver-&gt;upper_tolerance  = param.timesolver_upper_tolerance;
                  timesolver-&gt;component_norm   = SystemNorm(param.timesolver_norm);
        
                  timesolver-&gt;core_time_solver =
                    AutoPtr&lt;UnsteadySolver&gt;(innersolver);
                  system.time_solver =
                    AutoPtr&lt;UnsteadySolver&gt;(timesolver);
                }
              else
                system.time_solver =
                  AutoPtr&lt;TimeSolver&gt;(innersolver);
            }
          else
            system.time_solver =
              AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
          system.time_solver-&gt;reduce_deltat_on_diffsolver_failure =
                                                param.deltat_reductions;
          system.time_solver-&gt;quiet           = param.time_solver_quiet;
        
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
This function builds the Kelly error indicator. This indicator can be used
for comparisons of adjoint and non-adjoint based error indicators
</div>

<div class ="fragment">
<pre>
        AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters& /* param */)
        {
          AutoPtr&lt;ErrorEstimator&gt; error_estimator;
        
          error_estimator.reset(new KellyErrorEstimator);
        
          return error_estimator;
        }
        
</pre>
</div>
<div class = "comment">
Functions to build the adjoint based error indicators
The error_non_pressure and error_pressure constributions are estimated using
the build_error_estimator_component_wise function below
</div>

<div class ="fragment">
<pre>
        AutoPtr&lt;ErrorEstimator&gt;
        build_error_estimator_component_wise
          (FEMParameters &param,
           std::vector&lt;std::vector&lt;Real&gt; &gt; &term_weights,
           std::vector&lt;libMeshEnums::FEMNormType&gt; &primal_error_norm_type,
           std::vector&lt;libMeshEnums::FEMNormType&gt; &dual_error_norm_type)
        {
          AutoPtr&lt;ErrorEstimator&gt; error_estimator;
        
          AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;
        
          error_estimator.reset (adjoint_residual_estimator);
        
</pre>
</div>
<div class = "comment">
Both the primal and dual weights are going to be estimated using the patch recovery error estimator
</div>

<div class ="fragment">
<pre>
          PatchRecoveryErrorEstimator *p1 =
            new PatchRecoveryErrorEstimator;
          adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
        
          PatchRecoveryErrorEstimator *p2 =
            new PatchRecoveryErrorEstimator;
          adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);
        
</pre>
</div>
<div class = "comment">
Set the boolean for specifying whether we are reusing patches while building the patch recovery estimates
</div>

<div class ="fragment">
<pre>
          p1-&gt;set_patch_reuse(param.patch_reuse);
          p2-&gt;set_patch_reuse(param.patch_reuse);
        
</pre>
</div>
<div class = "comment">
Using the user filled error norm type vector, we pass the type of norm to be used for
the error in each variable, we can have different types of norms for the primal and
dual variables
</div>

<div class ="fragment">
<pre>
          unsigned int size = primal_error_norm_type.size();
        
          libmesh_assert_equal_to (size, dual_error_norm_type.size());
          for (unsigned int i = 0; i != size; ++i)
            {
              adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(i, primal_error_norm_type[i]);
              adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(i, dual_error_norm_type[i]);
            }
        
</pre>
</div>
<div class = "comment">
Now we set the right weights for each term in the error estimate, using the user provided
term_weights matrix
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (size, term_weights.size());
          for (unsigned int i = 0; i != term_weights.size(); ++i)
            {
              libmesh_assert_equal_to (size, term_weights[i].size());
              adjoint_residual_estimator-&gt;error_norm.set_weight(i, term_weights[i][i]);
              for (unsigned int j = 0; j != size; ++j)
                if (i != j)
                  adjoint_residual_estimator-&gt;error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
            }
        
          return error_estimator;
        }
        
</pre>
</div>
<div class = "comment">
The error_convection_diffusion_x and error_convection_diffusion_y are the nonlinear contributions which
are computed using the build_weighted_error_estimator_component_wise below
</div>

<div class ="fragment">
<pre>
        AutoPtr&lt;ErrorEstimator&gt;
        build_weighted_error_estimator_component_wise
          (FEMParameters &param,
           std::vector&lt;std::vector&lt;Real&gt; &gt; &term_weights,
           std::vector&lt;libMeshEnums::FEMNormType&gt; &primal_error_norm_type,
           std::vector&lt;libMeshEnums::FEMNormType&gt; &dual_error_norm_type,
           std::vector&lt;FEMFunctionBase&lt;Number&gt;*&gt; coupled_system_weight_functions)
        {
          AutoPtr&lt;ErrorEstimator&gt; error_estimator;
        
          AdjointResidualErrorEstimator *adjoint_residual_estimator = new AdjointResidualErrorEstimator;
        
          error_estimator.reset (adjoint_residual_estimator);
        
</pre>
</div>
<div class = "comment">
Using the user filled error norm type vector, we pass the type of norm to be used for
the error in each variable, we can have different types of norms for the primal and
dual variables


<br><br></div>

<div class ="fragment">
<pre>
          WeightedPatchRecoveryErrorEstimator *p1 =
            new WeightedPatchRecoveryErrorEstimator;
          adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
        
          PatchRecoveryErrorEstimator *p2 =
            new PatchRecoveryErrorEstimator;
          adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);
        
          p1-&gt;set_patch_reuse(param.patch_reuse);
          p2-&gt;set_patch_reuse(param.patch_reuse);
        
</pre>
</div>
<div class = "comment">
This is the critical difference with the build_error_estimate_component_wise
</div>

<div class ="fragment">
<pre>
          p1-&gt;weight_functions.clear();
        
</pre>
</div>
<div class = "comment">
We pass the pointers to the user specified weight functions to the patch recovery
error estimator objects declared above
</div>

<div class ="fragment">
<pre>
          unsigned int size = primal_error_norm_type.size();
        
          libmesh_assert(coupled_system_weight_functions.size() == size);
          libmesh_assert(dual_error_norm_type.size() == size);
        
          for (unsigned int i = 0; i != size; ++i)
            {
              p1-&gt;weight_functions.push_back(coupled_system_weight_functions[i]);
              adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(i, primal_error_norm_type[i]);
              adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(i, dual_error_norm_type[i]);
            }
        
</pre>
</div>
<div class = "comment">
Now we set the right weights for each term in the error estimate, using the user provided
term_weights matrix
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (size, term_weights.size());
          for (unsigned int i = 0; i != term_weights.size(); ++i)
            {
              libmesh_assert_equal_to (size, term_weights[i].size());
              adjoint_residual_estimator-&gt;error_norm.set_weight(i, term_weights[i][i]);
              for (unsigned int j = 0; j != size; ++j)
                if (i != j)
                  adjoint_residual_estimator-&gt;error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
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
This doesn't converge with Eigen BICGSTAB for some reason...
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc");
        
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
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
        
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
        
          std::cout &lt;&lt; "Building mesh" &lt;&lt; std::endl;
        
          if (!param.initial_timestep && param.run_simulation)
            build_domain(mesh, param);
        
          std::cout &lt;&lt; "Building system" &lt;&lt; std::endl;
        
          FEMSystem &system = build_system(equation_systems, infile, param);
        
          set_system_parameters(system, param);
        
          std::cout &lt;&lt; "Initializing systems" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Create a QoI object, we will need this to compute the QoI
</div>

<div class ="fragment">
<pre>
          {
            CoupledSystemQoI qoi;
</pre>
</div>
<div class = "comment">
Our QoI is computed on the side
</div>

<div class ="fragment">
<pre>
            qoi.assemble_qoi_sides = true;
            system.attach_qoi(&qoi);
          }
        
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
        
          std::cout&lt;&lt;"Starting adaptive loop"&lt;&lt;std::endl&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
Count the number of steps used
</div>

<div class ="fragment">
<pre>
          unsigned int a_step = 0;
        
</pre>
</div>
<div class = "comment">
Get the linear solver object to set the preconditioner reuse flag
</div>

<div class ="fragment">
<pre>
          LinearSolver&lt;Number&gt; *linear_solver = system.get_linear_solver();
        
</pre>
</div>
<div class = "comment">
Adaptively solve the timestep
</div>

<div class ="fragment">
<pre>
          for (; a_step &lt;= param.max_adaptivesteps; ++a_step)
            {
</pre>
</div>
<div class = "comment">
We have to ensure that we are refining to either an error
tolerance or to a target number of elements, and, not to a
non-zero tolerance and a target number of elements
simultaneously
</div>

<div class ="fragment">
<pre>
              if(param.global_tolerance &gt; 0 && param.nelem_target &gt; 0)
                {
                  std::cout&lt;&lt;"We cant refine to both a non-zero tolerance and a target number of elements, EXITING adaptive loop. "&lt;&lt;std::endl&lt;&lt;std::endl;
                  break;
                }
        
</pre>
</div>
<div class = "comment">
We dont want to reuse the preconditioner for the forward solve
</div>

<div class ="fragment">
<pre>
              linear_solver-&gt;reuse_preconditioner(false);
        
              std::cout&lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; std::endl;
        
              std::cout&lt;&lt; "Solving the forward problem" &lt;&lt;std::endl;
        
              std::cout &lt;&lt; "We have " &lt;&lt; mesh.n_active_elem()
                            &lt;&lt; " active elements and " &lt;&lt; equation_systems.n_active_dofs()
                        &lt;&lt; " active dofs." &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Solve the forward system
</div>

<div class ="fragment">
<pre>
              system.solve();
        
</pre>
</div>
<div class = "comment">
Write the output
</div>

<div class ="fragment">
<pre>
              write_output(equation_systems, 0, a_step, "primal", param);
        
</pre>
</div>
<div class = "comment">
Compute the QoI
</div>

<div class ="fragment">
<pre>
              system.assemble_qoi();
        
</pre>
</div>
<div class = "comment">
We just call a postprocess here to set the variable computed_QoI to the value computed by assemble_qoi
</div>

<div class ="fragment">
<pre>
              system.postprocess();
        
</pre>
</div>
<div class = "comment">
Get the value of the computed_QoI variable of the CoupledSystem class
</div>

<div class ="fragment">
<pre>
              Number QoI_0_computed = (dynamic_cast&lt;CoupledSystem&&gt;(system)).get_QoI_value();
        
              std::cout&lt;&lt; "The boundary QoI is " &lt;&lt; std::setprecision(17) &lt;&lt; QoI_0_computed &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
You dont need to compute error estimates and refine at the last
adaptive step, only before that
</div>

<div class ="fragment">
<pre>
              if(a_step!=param.max_adaptivesteps)
                {
                  NumericVector&lt;Number&gt; &primal_solution = (*system.solution);
        
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
We are about to solve the adjoint system, but before we
do this we see the same preconditioner flag to reuse the
preconditioner from the forward solver
</div>

<div class ="fragment">
<pre>
                  linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
        
</pre>
</div>
<div class = "comment">
Solve the adjoint system
</div>

<div class ="fragment">
<pre>
                  std::cout&lt;&lt; "Solving the adjoint problem" &lt;&lt;std::endl;
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
To plot the adjoint solution, we swap it with the primal solution
and use the write_output function
</div>

<div class ="fragment">
<pre>
                  NumericVector&lt;Number&gt; &dual_solution = system.get_adjoint_solution();
                  primal_solution.swap(dual_solution);
        
                  write_output(equation_systems, 0, a_step, "adjoint", param);
        
</pre>
</div>
<div class = "comment">
Swap back
</div>

<div class ="fragment">
<pre>
                  primal_solution.swap(dual_solution);
        
</pre>
</div>
<div class = "comment">
We need the values of the parameters Pe from the
system for the adjoint error estimate
</div>

<div class ="fragment">
<pre>
                  Real Pe = (dynamic_cast&lt;CoupledSystem&&gt;(system)).get_Pe();
        
</pre>
</div>
<div class = "comment">
The total error is the sum: error = error_non_pressure +
error_with_pressure + ...
error_estimator_convection_diffusion_x +
error_estimator_convection_diffusion_y
</div>

<div class ="fragment">
<pre>
                  ErrorVector error;
        
</pre>
</div>
<div class = "comment">
We first construct the non-pressure contributions
</div>

<div class ="fragment">
<pre>
                  ErrorVector error_non_pressure;
        
</pre>
</div>
<div class = "comment">
First we build the norm_type_vector_non_pressure vectors and
weights_matrix_non_pressure matrix for the non-pressure term
error contributions
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    primal_norm_type_vector_non_pressure;
                  primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
                  primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
                  primal_norm_type_vector_non_pressure.push_back(L2);
                  primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
        
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    dual_norm_type_vector_non_pressure;
                  dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
                  dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
                  dual_norm_type_vector_non_pressure.push_back(L2);
                  dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
        
                  std::vector&lt;std::vector&lt;Real&gt; &gt;
                    weights_matrix_non_pressure(system.n_vars(),
                      std::vector&lt;Real&gt;(system.n_vars(), 0.0));
                  weights_matrix_non_pressure[0][0] = 1.;
                  weights_matrix_non_pressure[1][1] = 1.;
                  weights_matrix_non_pressure[3][3] = 1./Pe;
        
</pre>
</div>
<div class = "comment">
We build the error estimator to estimate the contributions
to the QoI error from the non pressure term
</div>

<div class ="fragment">
<pre>
                  AutoPtr&lt;ErrorEstimator&gt; error_estimator_non_pressure =
                    build_error_estimator_component_wise
                      (param, weights_matrix_non_pressure,
                       primal_norm_type_vector_non_pressure,
                       dual_norm_type_vector_non_pressure);
        
</pre>
</div>
<div class = "comment">
Estimate the contributions to the QoI error from the non
pressure terms
</div>

<div class ="fragment">
<pre>
                  error_estimator_non_pressure-&gt;estimate_error(system, error_non_pressure);
        
</pre>
</div>
<div class = "comment">
Plot the estimated error from the non_pressure terms
</div>

<div class ="fragment">
<pre>
                  write_error(equation_systems, error_non_pressure, 0, a_step, param, "_non_pressure");
        
</pre>
</div>
<div class = "comment">
Now for the pressure contributions
</div>

<div class ="fragment">
<pre>
                  ErrorVector error_with_pressure;
        
</pre>
</div>
<div class = "comment">
Next we build the norm_type_vector_with_pressure vectors and
weights_matrix_with_pressure matrix for the pressure term
error contributions
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    primal_norm_type_vector_with_pressure;
                  primal_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
                  primal_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
                  primal_norm_type_vector_with_pressure.push_back(L2);
                  primal_norm_type_vector_with_pressure.push_back(L2);
        
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    dual_norm_type_vector_with_pressure;
                  dual_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
                  dual_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
                  dual_norm_type_vector_with_pressure.push_back(L2);
                  dual_norm_type_vector_with_pressure.push_back(L2);
        
                  std::vector&lt;std::vector&lt;Real&gt; &gt;
                    weights_matrix_with_pressure
                      (system.n_vars(),
                       std::vector&lt;Real&gt;(system.n_vars(), 0.0));
                  weights_matrix_with_pressure[0][2] = 1.;
        
                  weights_matrix_with_pressure[1][2] = 1.;
        
                  weights_matrix_with_pressure[2][0] = 1.;
                  weights_matrix_with_pressure[2][1] = 1.;
        
</pre>
</div>
<div class = "comment">
We build the error estimator to estimate the contributions
to the QoI error from the pressure term
</div>

<div class ="fragment">
<pre>
                  AutoPtr&lt;ErrorEstimator&gt; error_estimator_with_pressure =
                    build_error_estimator_component_wise
                      (param, weights_matrix_with_pressure,
                       primal_norm_type_vector_with_pressure,
                       dual_norm_type_vector_with_pressure);
        
</pre>
</div>
<div class = "comment">
Estimate the contributions to the QoI error from the pressure terms
</div>

<div class ="fragment">
<pre>
                  error_estimator_with_pressure-&gt;estimate_error(system, error_with_pressure);
        
</pre>
</div>
<div class = "comment">
Plot the error due to the pressure terms
</div>

<div class ="fragment">
<pre>
                  write_error(equation_systems, error_with_pressure, 0, a_step, param, "_with_pressure");
        
</pre>
</div>
<div class = "comment">
Now for the convection diffusion term errors (in the x and y directions)


<br><br></div>

<div class ="fragment">
<pre>
                  ErrorVector error_convection_diffusion_x;
        
</pre>
</div>
<div class = "comment">
The norm type vectors and weights matrix for the convection_diffusion_x errors
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    primal_norm_type_vector_convection_diffusion_x;
                  primal_norm_type_vector_convection_diffusion_x.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_x.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_x.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_x.push_back(H1_X_SEMINORM);
        
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    dual_norm_type_vector_convection_diffusion_x;
                  dual_norm_type_vector_convection_diffusion_x.push_back(L2);
                  dual_norm_type_vector_convection_diffusion_x.push_back(L2);
                  dual_norm_type_vector_convection_diffusion_x.push_back(L2);
</pre>
</div>
<div class = "comment">
Note that we need the error of the dual concentration in L2
</div>

<div class ="fragment">
<pre>
                  dual_norm_type_vector_convection_diffusion_x.push_back(L2);
        
                  std::vector&lt;std::vector&lt;Real&gt; &gt;
                    weights_matrix_convection_diffusion_x
                      (system.n_vars(),
                       std::vector&lt;Real&gt;(system.n_vars(), 0.0));
        	  weights_matrix_convection_diffusion_x[0][3] = 1.;
                  weights_matrix_convection_diffusion_x[3][3] = 1.;
        
</pre>
</div>
<div class = "comment">
We will also have to build and pass the weight functions to the weighted patch recovery estimators


<br><br>We pass the identity function as weights to error entries that the above matrix will scale to 0.
</div>

<div class ="fragment">
<pre>
                  ConstFEMFunction&lt;Number&gt; identity(1);
        
</pre>
</div>
<div class = "comment">
Declare object of class CoupledFEMFunctionsx, the definition of the function contains the weight
to be applied to the relevant terms
For ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2} term, returns C,1_h
</div>

<div class ="fragment">
<pre>
                  CoupledFEMFunctionsx convdiffx0(system, 0);
</pre>
</div>
<div class = "comment">
For ||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} term, returns (u_1)_h
</div>

<div class ="fragment">
<pre>
                  CoupledFEMFunctionsx convdiffx3(system, 3);
        
</pre>
</div>
<div class = "comment">
Make a vector of pointers to these objects
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;FEMFunctionBase&lt;Number&gt; *&gt; coupled_system_weight_functions_x;
                  coupled_system_weight_functions_x.push_back(&convdiffx0);
                  coupled_system_weight_functions_x.push_back(&identity);
                  coupled_system_weight_functions_x.push_back(&identity);
                  coupled_system_weight_functions_x.push_back(&convdiffx3);
        
</pre>
</div>
<div class = "comment">
Build the error estimator to estimate the contributions
to the QoI error from the convection diffusion x term
</div>

<div class ="fragment">
<pre>
                  AutoPtr&lt;ErrorEstimator&gt; error_estimator_convection_diffusion_x =
                  build_weighted_error_estimator_component_wise
                    (param, weights_matrix_convection_diffusion_x,
                     primal_norm_type_vector_convection_diffusion_x,
                     dual_norm_type_vector_convection_diffusion_x,
                     coupled_system_weight_functions_x);
        
</pre>
</div>
<div class = "comment">
Estimate the contributions to the QoI error from the
convection diffusion x term
</div>

<div class ="fragment">
<pre>
                  error_estimator_convection_diffusion_x-&gt;estimate_error
                    (system, error_convection_diffusion_x);
        
</pre>
</div>
<div class = "comment">
Plot this error
</div>

<div class ="fragment">
<pre>
                  write_error(equation_systems, error_convection_diffusion_x,
                              0, a_step, param, "_convection_diffusion_x");
        
</pre>
</div>
<div class = "comment">
Now for the y direction terms
</div>

<div class ="fragment">
<pre>
                  ErrorVector error_convection_diffusion_y;
        
</pre>
</div>
<div class = "comment">
The norm type vectors and weights matrix for the convection_diffusion_x errors
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    primal_norm_type_vector_convection_diffusion_y;
                  primal_norm_type_vector_convection_diffusion_y.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_y.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_y.push_back(L2);
                  primal_norm_type_vector_convection_diffusion_y.push_back(H1_Y_SEMINORM);
        
                  std::vector&lt;libMeshEnums::FEMNormType&gt;
                    dual_norm_type_vector_convection_diffusion_y;
                  dual_norm_type_vector_convection_diffusion_y.push_back(L2);
                  dual_norm_type_vector_convection_diffusion_y.push_back(L2);
                  dual_norm_type_vector_convection_diffusion_y.push_back(L2);
</pre>
</div>
<div class = "comment">
Note that we need the error of the dual concentration in L2
</div>

<div class ="fragment">
<pre>
                  dual_norm_type_vector_convection_diffusion_y.push_back(L2);
        
                  std::vector&lt;std::vector&lt;Real&gt; &gt;
                    weights_matrix_convection_diffusion_y
                      (system.n_vars(), std::vector&lt;Real&gt;(system.n_vars(), 0.0));
        	  weights_matrix_convection_diffusion_y[1][3] = 1.;
                  weights_matrix_convection_diffusion_y[3][3] = 1.;
        
</pre>
</div>
<div class = "comment">
For ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2} term, returns C,2_h
</div>

<div class ="fragment">
<pre>
                  CoupledFEMFunctionsy convdiffy1(system, 1);
</pre>
</div>
<div class = "comment">
For ||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} term, returns (u_2)_h
</div>

<div class ="fragment">
<pre>
                  CoupledFEMFunctionsy convdiffy3(system, 3);
        
</pre>
</div>
<div class = "comment">
Make a vector of pointers to these objects
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;FEMFunctionBase&lt;Number&gt; *&gt; coupled_system_weight_functions_y;
                  coupled_system_weight_functions_y.push_back(&identity);
                  coupled_system_weight_functions_y.push_back(&convdiffy1);
                  coupled_system_weight_functions_y.push_back(&identity);
                  coupled_system_weight_functions_y.push_back(&convdiffy3);
        
</pre>
</div>
<div class = "comment">
Build the error estimator to estimate the contributions
to the QoI error from the convection diffsion y term
</div>

<div class ="fragment">
<pre>
                  AutoPtr&lt;ErrorEstimator&gt; error_estimator_convection_diffusion_y =
                    build_weighted_error_estimator_component_wise
                      (param, weights_matrix_convection_diffusion_y,
                       primal_norm_type_vector_convection_diffusion_y,
                       dual_norm_type_vector_convection_diffusion_y,
                       coupled_system_weight_functions_y);
        
</pre>
</div>
<div class = "comment">
Estimate the contributions to the QoI error from the
convection diffusion y terms
</div>

<div class ="fragment">
<pre>
                  error_estimator_convection_diffusion_y-&gt;estimate_error(system, error_convection_diffusion_y);
        
</pre>
</div>
<div class = "comment">
Plot this error
</div>

<div class ="fragment">
<pre>
                  write_error(equation_systems, error_convection_diffusion_y, 0, a_step, param, "_convection_diffusion_y");
        
                  if(param.indicator_type == "adjoint_residual")
                    {
                      error.resize(error_non_pressure.size());
        
</pre>
</div>
<div class = "comment">
Now combine the contribs from the pressure and non
pressure terms to get the complete estimate
</div>

<div class ="fragment">
<pre>
                      for(unsigned int i = 0; i &lt; error.size(); i++)
                        {
                          error[i] = error_non_pressure[i] + error_with_pressure[i] + error_convection_diffusion_x[i] + error_convection_diffusion_y[i];
                        }
                    }
                  else
                    {
                      std::cout&lt;&lt;"Using Kelly Estimator"&lt;&lt;std::endl;
</pre>
</div>
<div class = "comment">
Build the Kelly error estimator
</div>

<div class ="fragment">
<pre>
                      AutoPtr&lt;ErrorEstimator&gt; error_estimator =  build_error_estimator(param);
        
</pre>
</div>
<div class = "comment">
Estimate the error
</div>

<div class ="fragment">
<pre>
                      error_estimator-&gt;estimate_error(system, error);
                    }
        
                  write_error(equation_systems, error, 0, a_step, param, "_total");
        
</pre>
</div>
<div class = "comment">
We have to refine either based on reaching an error
tolerance or a number of elements target, which should be
verified above


<br><br>Otherwise we flag elements by error tolerance or nelem
target


<br><br></div>

<div class ="fragment">
<pre>
                  if(param.refine_uniformly)
                    {
                      mesh_refinement-&gt;uniformly_refine(1);
                    }
                  else if(param.global_tolerance &gt;= 0. && param.nelem_target == 0.) // Refine based on reaching an error tolerance
                    {
                      mesh_refinement-&gt;flag_elements_by_error_tolerance (error);
        
</pre>
</div>
<div class = "comment">
Carry out the adaptive mesh refinement/coarsening
</div>

<div class ="fragment">
<pre>
                      mesh_refinement-&gt;refine_and_coarsen_elements();
                    }
                  else // Refine based on reaching a target number of elements
                    {
</pre>
</div>
<div class = "comment">
If we have reached the desired number of elements, we
dont do any more adaptive steps
</div>

<div class ="fragment">
<pre>
                      if (mesh.n_active_elem() &gt;= param.nelem_target)
                        {
                          std::cout&lt;&lt;"We reached the target number of elements."&lt;&lt;std::endl &lt;&lt;std::endl;
                          break;
                        }
        
                      mesh_refinement-&gt;flag_elements_by_nelem_target (error);
        
</pre>
</div>
<div class = "comment">
Carry out the adaptive mesh refinement/coarsening
</div>

<div class ="fragment">
<pre>
                      mesh_refinement-&gt;refine_and_coarsen_elements();
                    }
        
                  equation_systems.reinit();
        
                }
        
            }
        
        
          write_output_footers(param);
        
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
<br><br><br> <h1> The source file coupled_system.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/boundary_info.h"
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/getpot.h"
        #include "libmesh/mesh.h"
        #include "libmesh/parallel.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/zero_function.h"
        
        #include "coupled_system.h"
        
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
Function to set the Dirichlet boundary function for the inlet boundary velocities
</div>

<div class ="fragment">
<pre>
        class BdyFunction : public FunctionBase&lt;Number&gt;
        {
        public:
          BdyFunction (unsigned int u_var, unsigned int v_var, int sign)
            : _u_var(u_var), _v_var(v_var), _sign(sign)
            { this-&gt;_initialized = true; }
        
          virtual Number operator() (const Point&, const Real = 0)
            { libmesh_not_implemented(); }
        
          virtual void operator() (const Point& p,
                                   const Real,
                                   DenseVector&lt;Number&gt;& output)
            {
              output.resize(2);
              output.zero();
              const Real y=p(1);
</pre>
</div>
<div class = "comment">
Set the parabolic inflow boundary conditions at stations 0 & 1
</div>

<div class ="fragment">
<pre>
              output(_u_var) = (_sign)*((y-2) * (y-3));
              output(_v_var) = 0;
            }
        
          virtual AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() const
            { return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (new BdyFunction(_u_var, _v_var, _sign)); }
        
        private:
          const unsigned int _u_var, _v_var;
          const Real _sign;
        };
        
        
        void CoupledSystem::init_data ()
        {
</pre>
</div>
<div class = "comment">
Check the input file for Reynolds number, application type,
approximation type
</div>

<div class ="fragment">
<pre>
          GetPot infile("coupled_system.in");
          Peclet = infile("Peclet", 1.);
          unsigned int pressure_p = infile("pressure_p", 1);
          std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
        
</pre>
</div>
<div class = "comment">
LBB needs better-than-quadratic velocities for better-than-linear
pressures, and libMesh needs non-Lagrange elements for
better-than-quadratic velocities.
</div>

<div class ="fragment">
<pre>
          libmesh_assert((pressure_p == 1) || (fe_family != "LAGRANGE"));
        
          FEFamily fefamily = Utility::string_to_enum&lt;FEFamily&gt;(fe_family);
        
</pre>
</div>
<div class = "comment">
Add the velocity components "u" & "v".  They
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          u_var = this-&gt;add_variable ("u", static_cast&lt;Order&gt;(pressure_p+1),
        			      fefamily);
          v_var = this-&gt;add_variable ("v", static_cast&lt;Order&gt;(pressure_p+1),
        			      fefamily);
        
</pre>
</div>
<div class = "comment">
Add the pressure variable "p". This will
be approximated with a first-order basis,
providing an LBB-stable pressure-velocity pair.
</div>

<div class ="fragment">
<pre>
          p_var = this-&gt;add_variable ("p", static_cast&lt;Order&gt;(pressure_p),
        			      fefamily);
        
</pre>
</div>
<div class = "comment">
Add the Concentration variable "C". They will
be approximated using second-order approximation, the same as the velocity components
</div>

<div class ="fragment">
<pre>
          C_var = this-&gt;add_variable ("C", static_cast&lt;Order&gt;(pressure_p+1),
        			      fefamily);
        
</pre>
</div>
<div class = "comment">
Tell the system to march velocity forward in time, but
leave p as a constraint only
</div>

<div class ="fragment">
<pre>
          this-&gt;time_evolving(u_var);
          this-&gt;time_evolving(v_var);
          this-&gt;time_evolving(C_var);
        
</pre>
</div>
<div class = "comment">
Useful debugging options
</div>

<div class ="fragment">
<pre>
          this-&gt;verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
          this-&gt;print_jacobians = infile("print_jacobians", false);
          this-&gt;print_element_jacobians = infile("print_element_jacobians", false);
        
</pre>
</div>
<div class = "comment">
Set Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
          const boundary_id_type left_inlet_id = 0;
          std::set&lt;boundary_id_type&gt; left_inlet_bdy;
          left_inlet_bdy.insert(left_inlet_id);
        
          const boundary_id_type right_inlet_id = 1;
          std::set&lt;boundary_id_type&gt; right_inlet_bdy;
          right_inlet_bdy.insert(right_inlet_id);
        
          const boundary_id_type outlets_id = 2;
          std::set&lt;boundary_id_type&gt; outlets_bdy;
          outlets_bdy.insert(outlets_id);
        
          const boundary_id_type wall_id = 3;
          std::set&lt;boundary_id_type&gt; wall_bdy;
          wall_bdy.insert(wall_id);
        
</pre>
</div>
<div class = "comment">
The uv identifier for the setting the inlet and wall velocity boundary conditions
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; uv(1, u_var);
          uv.push_back(v_var);
</pre>
</div>
<div class = "comment">
The C_only identifier for setting the concentrations at the inlets
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; C_only(1, C_var);
        
</pre>
</div>
<div class = "comment">
The zero and constant functions
</div>

<div class ="fragment">
<pre>
          ZeroFunction&lt;Number&gt; zero;
          ConstFunction&lt;Number&gt; one(1);
        
</pre>
</div>
<div class = "comment">
We need two boundary functions for the inlets, because the signs on the velocities
will be different
</div>

<div class ="fragment">
<pre>
          int velocity_sign = 1;
          BdyFunction inflow_left(u_var, v_var, -velocity_sign);
          BdyFunction inflow_right(u_var, v_var, velocity_sign);
        
</pre>
</div>
<div class = "comment">
On the walls we will apply the no slip and no penetration boundary condition, u=0, v=0
</div>

<div class ="fragment">
<pre>
          this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (wall_bdy, uv, &zero));
        
</pre>
</div>
<div class = "comment">
On the inlet (left), we apply parabolic inflow boundary conditions for the velocity, u = - (y-2)*(y-3), v=0
and set C = 1
</div>

<div class ="fragment">
<pre>
          this-&gt;get_dof_map().add_dirichlet_boundary
            (DirichletBoundary (left_inlet_bdy, uv, &inflow_left));
          this-&gt;get_dof_map().add_dirichlet_boundary
            (DirichletBoundary (left_inlet_bdy, C_only, &one));
        
</pre>
</div>
<div class = "comment">
On the inlet (right), we apply parabolic inflow boundary conditions for the velocity, u = (y-2)*(y-3), v=0
and set C = 0
</div>

<div class ="fragment">
<pre>
          this-&gt;get_dof_map().add_dirichlet_boundary
            (DirichletBoundary (right_inlet_bdy, uv, &inflow_right));
          this-&gt;get_dof_map().add_dirichlet_boundary
            (DirichletBoundary (right_inlet_bdy, C_only, &zero));
        
</pre>
</div>
<div class = "comment">
Note that the remaining boundary conditions are the natural boundary conditions for the concentration
on the wall (grad(c) dot n = 0) and natural boundary conditions for the velocity and the concentration
on the outlets ((grad(velocity) dot n - p n) dot t = 0, grad(C) dot n = 0)


<br><br>Do the parent's initialization after variables and boundary constraints are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        }
        
        void CoupledSystem::init_context(DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
We should prerequest all the data
we will need to build the linear system.
Note that the concentration and velocity components
use the same basis.
</div>

<div class ="fragment">
<pre>
          c.element_fe_var[u_var]-&gt;get_JxW();
          c.element_fe_var[u_var]-&gt;get_phi();
          c.element_fe_var[u_var]-&gt;get_dphi();
          c.element_fe_var[u_var]-&gt;get_xyz();
        
          c.element_fe_var[p_var]-&gt;get_phi();
        
          c.side_fe_var[u_var]-&gt;get_JxW();
          c.side_fe_var[u_var]-&gt;get_phi();
          c.side_fe_var[u_var]-&gt;get_xyz();
        }
        
        
        bool CoupledSystem::element_time_derivative (bool request_jacobian,
                                                    DiffContext &context)
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
          const std::vector&lt;Real&gt; &JxW =
            c.element_fe_var[u_var]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi =
            c.element_fe_var[u_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
            c.element_fe_var[u_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi =
            c.element_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
          libmesh_assert_equal_to (n_u_dofs, c.dof_indices_var[v_var].size());
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &Kuu = *c.elem_subjacobians[u_var][u_var];
          DenseSubMatrix&lt;Number&gt; &Kup = *c.elem_subjacobians[u_var][p_var];
          DenseSubVector&lt;Number&gt; &Fu = *c.elem_subresiduals[u_var];
        
          DenseSubMatrix&lt;Number&gt; &Kvv = *c.elem_subjacobians[v_var][v_var];
          DenseSubMatrix&lt;Number&gt; &Kvp = *c.elem_subjacobians[v_var][p_var];
          DenseSubVector&lt;Number&gt; &Fv = *c.elem_subresiduals[v_var];
        
          DenseSubMatrix&lt;Number&gt; &KCu = *c.elem_subjacobians[C_var][u_var];
          DenseSubMatrix&lt;Number&gt; &KCv = *c.elem_subjacobians[C_var][v_var];
          DenseSubMatrix&lt;Number&gt; &KCC = *c.elem_subjacobians[C_var][C_var];
          DenseSubVector&lt;Number&gt; &FC = *c.elem_subresiduals[C_var];
        
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
Compute the solution & its gradient at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Number p = c.interior_value(p_var, qp),
        	u = c.interior_value(u_var, qp),
        	v = c.interior_value(v_var, qp);
              Gradient grad_u = c.interior_gradient(u_var, qp),
        	grad_v = c.interior_gradient(v_var, qp),
        	grad_C = c.interior_gradient(C_var, qp);
        
</pre>
</div>
<div class = "comment">
Definitions for convenience.  It is sometimes simpler to do a
dot product if you have the full vector at your disposal.
</div>

<div class ="fragment">
<pre>
              NumberVectorValue U     (u,     v);
              const Number C_x = grad_C(0);
              const Number C_y = grad_C(1);
        
</pre>
</div>
<div class = "comment">
First, an i-loop over the velocity degrees of freedom.
We know that n_u_dofs == n_v_dofs so we can compute contributions
for both at the same time.
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_u_dofs; i++)
                {
</pre>
</div>
<div class = "comment">
Stokes equations residuals
</div>

<div class ="fragment">
<pre>
                  Fu(i) += JxW[qp] *
                           (p*dphi[i][qp](0) -                // pressure term
        		    (grad_u*dphi[i][qp]));            // diffusion term
        
                  Fv(i) += JxW[qp] *
                           (p*dphi[i][qp](1) -                // pressure term
        		    (grad_v*dphi[i][qp]));            // diffusion term
        
</pre>
</div>
<div class = "comment">
Concentration Equation Residual
</div>

<div class ="fragment">
<pre>
                  FC(i) += JxW[qp] *
        	           ( (U*grad_C)*phi[i][qp] +                // convection term
        	            (1./Peclet)*(grad_C*dphi[i][qp]) );     // diffusion term
        
</pre>
</div>
<div class = "comment">
Note that the Fp block is identically zero unless we are using
some kind of artificial compressibility scheme...


<br><br></div>

<div class ="fragment">
<pre>
                  if (request_jacobian && c.elem_solution_derivative)
                    {
                      libmesh_assert (c.elem_solution_derivative == 1.0);
        
</pre>
</div>
<div class = "comment">
Matrix contributions for the uu and vv couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        {
                          Kuu(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); /* diffusion term  */
        
                          Kvv(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); /* diffusion term  */
        
        		  KCu(i,j) += JxW[qp]* ( (phi[j][qp]*C_x)*phi[i][qp] ); /* convection term */
        
        		  KCv(i,j) += JxW[qp]*( (phi[j][qp]*C_y)*phi[i][qp] );  /* convection term */
        
        		  KCC(i,j) += JxW[qp]*
        		              ( (U*dphi[j][qp])*phi[i][qp] +      /* nonlinear term (convection) */
        		              (1./Peclet)*(dphi[j][qp]*dphi[i][qp]) ); /* diffusion term */
        		}
        
</pre>
</div>
<div class = "comment">
Matrix contributions for the up and vp couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j != n_p_dofs; j++)
        		{
        		  Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
        		  Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
        		}
        	    }
        	}
        
            } // end of the quadrature point qp-loop
        
              return request_jacobian;
        }
        
        
        bool CoupledSystem::element_constraint (bool request_jacobian,
                                               DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weight for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.element_fe_var[u_var]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
            c.element_fe_var[u_var]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi =
            c.element_fe_var[p_var]-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
          const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &Kpu = *c.elem_subjacobians[p_var][u_var];
          DenseSubMatrix&lt;Number&gt; &Kpv = *c.elem_subjacobians[p_var][v_var];
          DenseSubVector&lt;Number&gt; &Fp = *c.elem_subresiduals[p_var];
        
</pre>
</div>
<div class = "comment">
Add the constraint given by the continuity equation
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
</pre>
</div>
<div class = "comment">
Compute the velocity gradient at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Gradient grad_u = c.interior_gradient(u_var, qp),
        	grad_v = c.interior_gradient(v_var, qp);
        
</pre>
</div>
<div class = "comment">
Now a loop over the pressure degrees of freedom.  This
computes the contributions of the continuity equation.
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_p_dofs; i++)
                {
                  Fp(i) += JxW[qp] * psi[i][qp] *
                           (grad_u(0) + grad_v(1));
        
                  if (request_jacobian && c.elem_solution_derivative)
                    {
                      libmesh_assert (c.elem_solution_derivative == 1.0);
        
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        {
                          Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                          Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                        }
                    }
                }
            } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
        void CoupledSystem::postprocess()
        {
</pre>
</div>
<div class = "comment">
We need to overload the postprocess function to set the computed_QoI variable of the CoupledSystem class
to the qoi value stored in System::qoi[0]


<br><br></div>

<div class ="fragment">
<pre>
          computed_QoI = 0.0;
        
          computed_QoI = System::qoi[0];
        }
        
</pre>
</div>
<div class = "comment">
These functions supply the nonlinear weighting for the adjoint residual error estimate which
arise due to the convection term in the convection-diffusion equation:
||e((u_1)_h C,1)||_{L2} ||e(C^*)||_{L2} + ||e(u1 C,1_h)||_{L2} ||e(C^*)||_{L2}
||e((u_2)_h C,2)||_{L2} ||e(C^*)||_{L2} + ||e(u2 C,2_h)||_{L2} ||e(C^*)||_{L2}
These functions compute (u_1)_h or C,1_h , and (u_2)_h or C,2_h , and supply it to the weighted patch recovery error estimator
In CoupledFEMFunctionsx, the object built with var = 0, returns the (u_1)_h weight, while
the object built with var = 1, returns the C,1_h weight. The switch statment
distinguishes the behavior of the two objects
Same thing for CoupledFEMFunctionsy
</div>

<div class ="fragment">
<pre>
        Number CoupledFEMFunctionsx::operator()(const FEMContext& c, const Point& p,
        					const Real /* time */)
        {
          Number weight = 0.0;
        
          switch(var)
            {
            case 0:
              {
        	Gradient grad_C = c.point_gradient(3, p);
        
        	weight = grad_C(0);
              }
              break;
        
            case 3:
              {
        	Number u = c.point_value(0, p);
        
        	weight = u;
              }
              break;
        
            default:
              {
        	std::cout&lt;&lt;"Wrong variable number"&lt;&lt;var&lt;&lt;" passed to CoupledFEMFunctionsx object ! Quitting !"&lt;&lt;std::endl;
        	libmesh_error();
              }
        
            }
        
          return weight;
        }
        
        Number CoupledFEMFunctionsy::operator()(const FEMContext& c, const Point& p,
        					const Real /* time */)
        {
          Number weight = 0.0;
        
          switch(var)
            {
            case 1:
              {
        	Gradient grad_C = c.point_gradient(3, p);
        
        	weight = grad_C(1);
              }
              break;
        
            case 3:
              {
        	Number v = c.point_value(1, p);
        
        	weight = v;
              }
              break;
        
            default:
              {
        	std::cout&lt;&lt;"Wrong variable number "&lt;&lt;var&lt;&lt;" passed to CoupledFEMFunctionsy object ! Quitting !"&lt;&lt;std::endl;
        	libmesh_error();
              }
            }
        
          return weight;
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file domain.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh.h"
        #include "libmesh/serial_mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/mesh_modification.h"
        #include "libmesh/mesh_refinement.h"
        
        #include "domain.h"
        
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
Mesh construction
</div>

<div class ="fragment">
<pre>
        void build_domain (Mesh &mesh, FEMParameters &param)
        {
          mesh.read(param.domainfile);
        
          std::cout&lt;&lt;"Making elements 2nd order"&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
Right now we are setting approximation orders in the code, rather than reading them in
That needs to be fixed and the second ordering should be done only if one of the
approximation orders is greater than 1
</div>

<div class ="fragment">
<pre>
          mesh.all_second_order();
        
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
        
        using namespace libMesh;
        
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
        
            GETPOT_INT_INPUT(initial_timestep);
            GETPOT_INT_INPUT(n_timesteps);
            GETPOT_INPUT(transient);
            GETPOT_INT_INPUT(deltat_reductions);
            GETPOT_INPUT(timesolver_core);
            GETPOT_INPUT(end_time);
            GETPOT_INPUT(deltat);
            GETPOT_INPUT(timesolver_theta);
            GETPOT_INPUT(timesolver_maxgrowth);
            GETPOT_INPUT(timesolver_upper_tolerance);
            GETPOT_INPUT(timesolver_tolerance);
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
            GETPOT_INPUT(fine_mesh_file_primal);
            GETPOT_INPUT(fine_mesh_soln_primal);
            GETPOT_INPUT(fine_mesh_file_adjoint);
            GETPOT_INPUT(fine_mesh_soln_adjoint);
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
            GETPOT_INT_INPUT(coarsecoarsenings);
            GETPOT_INT_INPUT(extrarefinements);
            GETPOT_INPUT(use_petsc_snes);
            GETPOT_INPUT(time_solver_quiet);
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
            GETPOT_INT_INPUT(initial_adaptivesteps);
            GETPOT_INT_INPUT(initial_sobolev_order);
            GETPOT_INT_INPUT(initial_extra_quadrature);
            GETPOT_INPUT(refine_uniformly);
            GETPOT_INPUT(indicator_type);
            GETPOT_INPUT(adjoint_residual_type);
            GETPOT_INPUT(patch_reuse);
            GETPOT_INT_INPUT(sobolev_order);
            GETPOT_INPUT(alternate_with_uniform_steps);
            GETPOT_INPUT(alternate_step_number);
            GETPOT_INPUT(component_wise_error);
            GETPOT_INPUT(compute_sensitivities);
            GETPOT_INPUT(compare_to_fine_solution_primal);
            GETPOT_INPUT(compare_to_fine_solution_adjoint);
            GETPOT_INPUT(do_forward_sensitivity);
            GETPOT_INPUT(do_adjoint_sensitivity);
            GETPOT_INPUT(postprocess_adjoint);
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
        
            GETPOT_INPUT(run_simulation);
            GETPOT_INPUT(run_postprocess);
        
            run_simulation              = input("run_simulation", run_simulation);
            run_postprocess             = input("run_postprocess", run_postprocess);
        
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
            GETPOT_INPUT(print_solution_norms);
            GETPOT_INPUT(print_solutions);
            GETPOT_INPUT(print_residual_norms);
            GETPOT_INPUT(print_residuals);
            GETPOT_INPUT(print_jacobian_norms);
            GETPOT_INPUT(print_jacobians);
        
        
</pre>
</div>
<div class = "comment">
std::vector<std::string> bad_variables =
input.unidentified_arguments(variable_names);


<br><br>if (libMesh::processor_id() == 0 && !bad_variables.empty())
{
std::cerr << "ERROR: Unrecognized variables:" << std::endl;
for (unsigned int i = 0; i != bad_variables.size(); ++i)
std::cerr << bad_variables[i] << std::endl;
std::cerr << "not found among recognized variables." << std::endl;
for (unsigned int i = 0; i != variable_names.size(); ++i)
std::cerr << variable_names[i] << std::endl;
libmesh_error();
}
</div>

<div class ="fragment">
<pre>
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file H-qoi.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "H-qoi.h"
        
</pre>
</div>
<div class = "comment">
Here we define the functions to compute the QoI (side_qoi)
and supply the right hand side for the associated adjoint problem (side_qoi_derivative)


<br><br></div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        void CoupledSystemQoI::init_qoi( std::vector&lt;Number&gt;& sys_qoi)
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
This function supplies the right hand side for the adjoint problem
We only have one QoI, so we don't bother checking the qois argument
to see if it was requested from us
</div>

<div class ="fragment">
<pre>
        void CoupledSystemQoI::side_qoi_derivative (DiffContext &context,
        					 const QoISet & /* qois */)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.side_fe_var[0]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
Get velocity basis functions phi
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &phi = c.side_fe_var[0]-&gt;get_phi();
        
          const std::vector&lt;Point &gt; &q_point = c.side_fe_var[0]-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[1].size();
        
          DenseSubVector&lt;Number&gt; &Qu = *c.elem_qoi_subderivatives[0][0];
          DenseSubVector&lt;Number&gt; &QC = *c.elem_qoi_subderivatives[0][3];
        
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
          unsigned int n_qpoints = c.side_qrule-&gt;n_points();
        
          Number u = 0. ;
          Number C = 0. ;
        
</pre>
</div>
<div class = "comment">
// If side is on outlet
</div>

<div class ="fragment">
<pre>
          if(c.has_side_boundary_id(2))
            {
</pre>
</div>
<div class = "comment">
Loop over all the qps on this side
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp != n_qpoints; qp++)
          	{
        	  Real x = q_point[qp](0);
        
</pre>
</div>
<div class = "comment">
If side is on left outlet
</div>

<div class ="fragment">
<pre>
                  if(x &lt; 0.)
        	    {
</pre>
</div>
<div class = "comment">
Get u at the qp
</div>

<div class ="fragment">
<pre>
                      u = c.side_value(0,qp);
        	      C = c.side_value(3,qp);
        
</pre>
</div>
<div class = "comment">
Add the contribution from each basis function
</div>

<div class ="fragment">
<pre>
                      for (unsigned int i=0; i != n_u_dofs; i++)
        		{
        		  Qu(i) += JxW[qp] * -phi[i][qp] * C;
        		  QC(i) += JxW[qp] * phi[i][qp] * -u;
        		}
        	    } // end if
        
          	} // end quadrature loop
        
            } // end if on outlet
        }
        
</pre>
</div>
<div class = "comment">
This function computes the actual QoI
</div>

<div class ="fragment">
<pre>
        void CoupledSystemQoI::side_qoi(DiffContext &context, const QoISet & /* qois */)
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
        
</pre>
</div>
<div class = "comment">
Loop over qp's, compute the function at each qp and add
to get the QoI


<br><br></div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = c.side_qrule-&gt;n_points();
        
          Number dQoI_0 = 0. ;
          Number u = 0. ;
          Number C = 0. ;
        
</pre>
</div>
<div class = "comment">
If side is on the left outlet
</div>

<div class ="fragment">
<pre>
          if(c.has_side_boundary_id(2))
            {
</pre>
</div>
<div class = "comment">
Loop over all the qps on this side
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp != n_qpoints; qp++)
        	{
        	  Real x = q_point[qp](0);
        
        	  if(x &lt; 0.)
        	    {
</pre>
</div>
<div class = "comment">
Get u and C at the qp
</div>

<div class ="fragment">
<pre>
                      u = c.side_value(0,qp);
        	      C = c.side_value(3,qp);
        
        	      dQoI_0 += JxW[qp] * -u * C;
        	    } // end if
        
          	} // end quadrature loop
        
            } // end if on bdry
        
          c.elem_qoi[0] += dQoI_0;
        
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
        #include "femparameters.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
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
        Number initial_value(const Point& /* p */,
                             const Parameters& /* param */,
                             const std::string&,
                             const std::string&)
        {
        
          return 1.;
        
        }
        
        
        
        Gradient initial_grad(const Point& /* p */,
                              const Parameters& /* param */,
                              const std::string&,
                              const std::string&)
        {
          return 0.;
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file output.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include &lt;fstream&gt;
        #include &lt;unistd.h&gt;
        
        #include "libmesh/libmesh_common.h"
        
        #include "output.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        void start_output(unsigned int timesteps,
                          std::string filename,
                          std::string varname)
        {
</pre>
</div>
<div class = "comment">
Truncate then append if we're doing a restart
</div>

<div class ="fragment">
<pre>
          if (timesteps)
            {
              std::ifstream infile (filename.c_str());
              char buf[1025];
              if (!infile.is_open())
                libmesh_error();
        
</pre>
</div>
<div class = "comment">
Count off the lines we want to save, including
the original header
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != timesteps+1; ++i)
                infile.getline(buf, 1024);
              if (!infile.good())
                libmesh_error();
              unsigned int length = infile.tellg();
              infile.close();
        
</pre>
</div>
<div class = "comment">
Then throw away the rest
</div>

<div class ="fragment">
<pre>
              int err = truncate(filename.c_str(), length);
              if (err != 0)
                libmesh_error();
            }
</pre>
</div>
<div class = "comment">
Write new headers if we're starting a new simulation
</div>

<div class ="fragment">
<pre>
          else
            {
              std::ofstream outfile (filename.c_str(), std::ios_base::out);
              outfile &lt;&lt; varname &lt;&lt; " = [" &lt;&lt; std::endl;
            }
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file coupled_system.h without comments: </h1> 
<pre> 
  #ifndef __coupled_system_h__
  #define __coupled_system_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_function_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameter_vector.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> CoupledSystem : <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    CoupledSystem(EquationSystems&amp; es,
                 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                 <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
      : FEMSystem(es, name_in, number_in), Peclet(1.) {qoi.resize(1);}
  
  
    Number &amp;get_QoI_value()
      {
        <B><FONT COLOR="#A020F0">return</FONT></B> computed_QoI;
      }
  
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
  
    Real &amp;get_Pe()
      {
        <B><FONT COLOR="#A020F0">return</FONT></B> Peclet;
      }
  
   <B><FONT COLOR="#228B22">protected</FONT></B>:
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                          DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                     DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> postprocess ();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; parameters;
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var, u_var, v_var, C_var;
  
    ParameterVector parameter_vector;
  
    Real Peclet;
  
    Number computed_QoI;
  
  };
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> CoupledFEMFunctionsx : <B><FONT COLOR="#228B22">public</FONT></B> FEMFunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    CoupledFEMFunctionsx(System&amp; <I><FONT COLOR="#B22222">/* sys */</FONT></I>, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var_number) {var = var_number;}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~CoupledFEMFunctionsx () {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt; clone () <B><FONT COLOR="#228B22">const</FONT></B>
    {<B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt;( <B><FONT COLOR="#A020F0">new</FONT></B> CoupledFEMFunctionsx(*<B><FONT COLOR="#A020F0">this</FONT></B>) ); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp;,
  			   <B><FONT COLOR="#228B22">const</FONT></B> Real, DenseVector&lt;Number&gt;&amp;)
    {libmesh_error();}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			     <B><FONT COLOR="#228B22">const</FONT></B> Real time = 0.);
  
   <B><FONT COLOR="#228B22">private</FONT></B>:
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var;
  
  };
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> CoupledFEMFunctionsy : <B><FONT COLOR="#228B22">public</FONT></B> FEMFunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    CoupledFEMFunctionsy(System&amp; <I><FONT COLOR="#B22222">/* sys */</FONT></I>, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var_number) {var = var_number;}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~CoupledFEMFunctionsy () {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt; clone () <B><FONT COLOR="#228B22">const</FONT></B>
    {<B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FEMFunctionBase&lt;Number&gt; &gt;( <B><FONT COLOR="#A020F0">new</FONT></B> CoupledFEMFunctionsy(*<B><FONT COLOR="#A020F0">this</FONT></B>) ); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp;,
  			   <B><FONT COLOR="#228B22">const</FONT></B> Real,
  			   DenseVector&lt;Number&gt;&amp;)
    {libmesh_error();}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			     <B><FONT COLOR="#228B22">const</FONT></B> Real time = 0.);
  
   <B><FONT COLOR="#228B22">private</FONT></B>:
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var;
  
  };
  
  #endif <I><FONT COLOR="#B22222">//__coupled_system_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file domain.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">class</FONT></B> FEMParameters;
  
  <B><FONT COLOR="#228B22">void</FONT></B> build_domain (Mesh &amp;mesh, FEMParameters &amp;param);
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.h without comments: </h1> 
<pre> 
  #ifndef __fem_parameters_h__
  #define __fem_parameters_h__
  
  #include &lt;limits&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/enum_norm_type.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> FEMParameters
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
      FEMParameters() :
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
  	fine_mesh_file_primal(<B><FONT COLOR="#BC8F8F">&quot;fine_mesh.xda&quot;</FONT></B>), fine_mesh_soln_primal(<B><FONT COLOR="#BC8F8F">&quot;fine_mesh_soln.xda&quot;</FONT></B>),
  	fine_mesh_file_adjoint(<B><FONT COLOR="#BC8F8F">&quot;fine_mesh.xda&quot;</FONT></B>), fine_mesh_soln_adjoint(<B><FONT COLOR="#BC8F8F">&quot;fine_mesh_soln.xda&quot;</FONT></B>),
        elementorder(2),
        domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
        domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
        coarsegridx(1), coarsegridy(1), coarsegridz(1),
  	coarserefinements(0), coarsecoarsenings(0), extrarefinements(0),
        use_petsc_snes(false),
        time_solver_quiet(true), solver_quiet(true),
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
        initial_adaptivesteps(0), initial_sobolev_order(1),
  	initial_extra_quadrature(0),
  	indicator_type(<B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>), patch_reuse(true), sobolev_order(1), adjoint_residual_type(<B><FONT COLOR="#BC8F8F">&quot;patchpatch&quot;</FONT></B>), alternate_with_uniform_steps(<B><FONT COLOR="#BC8F8F">&quot;false&quot;</FONT></B>), alternate_step_number(2),
  	component_wise_error(<B><FONT COLOR="#BC8F8F">&quot;false&quot;</FONT></B>), compare_to_fine_solution_primal(false), compare_to_fine_solution_adjoint(false), compute_sensitivities(false), do_forward_sensitivity(true), do_adjoint_sensitivity(false),
  	postprocess_adjoint(false),
        write_interval(10),
        write_gmv_error(false), write_tecplot_error(false),
        output_xda(false), output_xdr(false),
        output_bz2(true), output_gz(true),
        output_gmv(false), output_tecplot(false),
        run_simulation(true), run_postprocess(false),
        fe_family(1, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>), fe_order(1, 1),
        extra_quadrature_order(0),
        analytic_jacobians(true), verify_analytic_jacobians(0.0),
        print_solution_norms(false), print_solutions(false),
        print_residual_norms(false), print_residuals(false),
        print_jacobian_norms(false), print_jacobians(false) {}
  
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
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string domaintype, domainfile, elementtype, fine_mesh_file_primal, fine_mesh_soln_primal, fine_mesh_file_adjoint, fine_mesh_soln_adjoint;
      Real elementorder;
      Real domain_xmin, domain_ymin, domain_zmin;
      Real domain_edge_width, domain_edge_length, domain_edge_height;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarsegridx, coarsegridy, coarsegridz;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarserefinements, coarsecoarsenings, extrarefinements;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> use_petsc_snes;
      <B><FONT COLOR="#228B22">bool</FONT></B> time_solver_quiet, solver_quiet, reuse_preconditioner, require_residual_reduction;
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
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_adaptivesteps;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_sobolev_order;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> initial_extra_quadrature;
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string indicator_type;
      <B><FONT COLOR="#228B22">bool</FONT></B> patch_reuse;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> sobolev_order;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string adjoint_residual_type;
      <B><FONT COLOR="#228B22">bool</FONT></B> alternate_with_uniform_steps;
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> alternate_step_number;
      <B><FONT COLOR="#228B22">bool</FONT></B> component_wise_error;
      <B><FONT COLOR="#228B22">bool</FONT></B> compare_to_fine_solution_primal, compare_to_fine_solution_adjoint;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> compute_sensitivities;
      <B><FONT COLOR="#228B22">bool</FONT></B> do_forward_sensitivity;
      <B><FONT COLOR="#228B22">bool</FONT></B> do_adjoint_sensitivity;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> postprocess_adjoint;
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> write_interval;
      <B><FONT COLOR="#228B22">bool</FONT></B> write_gmv_error, write_tecplot_error, output_xda, output_xdr,
  	 output_bz2, output_gz, output_gmv, output_tecplot;
      <B><FONT COLOR="#228B22">bool</FONT></B> run_simulation, run_postprocess;
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; fe_family;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; fe_order;
      <B><FONT COLOR="#228B22">int</FONT></B> extra_quadrature_order;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> analytic_jacobians;
      Real verify_analytic_jacobians;
  
      <B><FONT COLOR="#228B22">bool</FONT></B> print_solution_norms, print_solutions,
           print_residual_norms, print_residuals,
           print_jacobian_norms, print_jacobians;
  };
  
  #endif <I><FONT COLOR="#B22222">// __fem_parameters_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file H-qoi.h without comments: </h1> 
<pre> 
  #ifndef H_QOI_H
  #define H_QOI_H
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/point.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/diff_qoi.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> CoupledSystemQoI : <B><FONT COLOR="#228B22">public</FONT></B> DifferentiableQoI
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    CoupledSystemQoI(){}
    <B><FONT COLOR="#228B22">virtual</FONT></B> ~CoupledSystemQoI(){}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_qoi( std::vector&lt;Number&gt;&amp; sys_qoi);
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> postprocess( ){}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> side_qoi_derivative(DiffContext &amp;context, <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> side_qoi(DiffContext &amp;context, <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; qois);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;DifferentiableQoI&gt; clone( )
    { AutoPtr&lt;DifferentiableQoI&gt; my_clone( <B><FONT COLOR="#A020F0">new</FONT></B> CoupledSystemQoI );
      *my_clone = *<B><FONT COLOR="#A020F0">this</FONT></B>;
      <B><FONT COLOR="#A020F0">return</FONT></B> my_clone;
    }
  
  };
  #endif <I><FONT COLOR="#B22222">// H_QOI_H
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
<br><br><br> <h1> The source file mysystems.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;coupled_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  FEMSystem &amp;build_system(EquationSystems &amp;es, GetPot &amp;, FEMParameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>)
  {
    CoupledSystem &amp;system = es.add_system&lt;CoupledSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;CoupledSystem&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> system;
  }
  
  Number exact_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; <I><FONT COLOR="#B22222">/* p */</FONT></I>,       <I><FONT COLOR="#B22222">// xyz location
</FONT></I>                     <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>,  <I><FONT COLOR="#B22222">// EquationSystem parameters
</FONT></I>                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name
</FONT></I>                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unknown_name
</FONT></I>  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Warning ! No exact value function specified ! Returning 0.&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  Gradient exact_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; <I><FONT COLOR="#B22222">/* p */</FONT></I>,       <I><FONT COLOR="#B22222">// xyz location
</FONT></I>                      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>,  <I><FONT COLOR="#B22222">// EquationSystems parameters
</FONT></I>                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name
</FONT></I>                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unknown_name
</FONT></I>  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Warning ! No exact value function specified ! Returning 0.&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file output.h without comments: </h1> 
<pre> 
  #include &lt;string&gt;
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> start_output(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> timesteps,
                    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string filename,
                    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string varname);
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file adjoints_ex3.C without comments: </h1> 
<pre> 
  #include &lt;iostream&gt;
  #include &lt;sys/time.h&gt;
  #include &lt;iomanip&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/twostep_time_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/euler_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/euler2_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/petsc_diff_solver.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_tools.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system_norm.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/adjoint_residual_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/const_fem_function.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_function_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parameter_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/patch_recovery_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/petsc_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sensitivity_data.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/uniform_refinement_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/qoi_set.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/weighted_patch_recovery_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;coupled_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;domain.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;initial.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mysystems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;output.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;H-qoi.h&quot;</FONT></B>
  
  
  
  
  
  
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string numbered_filename(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step, <I><FONT COLOR="#B22222">// The timestep count
</FONT></I>                                <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step, <I><FONT COLOR="#B22222">// The adaptive step count
</FONT></I>                                <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type, <I><FONT COLOR="#B22222">// primal or adjoint solve
</FONT></I>                                <B><FONT COLOR="#5F9EA0">std</FONT></B>::string type,
                                <B><FONT COLOR="#5F9EA0">std</FONT></B>::string extension,
                                FEMParameters &amp;param)
  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name;
    file_name &lt;&lt; solution_type
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.out.&quot;</FONT></B>
              &lt;&lt; type
              &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B>
              &lt;&lt; std::setw(3)
              &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
              &lt;&lt; std::right
              &lt;&lt; t_step
              &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B>
              &lt;&lt; std::setw(2)
              &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
              &lt;&lt; std::right
              &lt;&lt; a_step
              &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B>
              &lt;&lt; extension;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_bz2)
      file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.bz2&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_gz)
      file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.gz&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">return</FONT></B> file_name.str();
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output(EquationSystems &amp;es,
                    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step, <I><FONT COLOR="#B22222">// The timestep count
</FONT></I>                    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step, <I><FONT COLOR="#B22222">// The adaptive step count
</FONT></I>                    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type, <I><FONT COLOR="#B22222">// primal or adjoint solve
</FONT></I>                    FEMParameters &amp;param)
  {
    MeshBase &amp;mesh = es.get_mesh();
  
  #ifdef LIBMESH_HAVE_GMV
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_gmv)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name_gmv;
        file_name_gmv &lt;&lt; solution_type
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.out.gmv.&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_step
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B>
                      &lt;&lt; std::setw(2)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; a_step;
  
        GMVIO(mesh).write_equation_systems
          (file_name_gmv.str(), es);
      }
  #endif
  
  #ifdef LIBMESH_HAVE_TECPLOT_API
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_tecplot)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name_tecplot;
        file_name_tecplot &lt;&lt; solution_type
                          &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.out.&quot;</FONT></B>
                          &lt;&lt; std::setw(3)
                          &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                          &lt;&lt; std::right
                          &lt;&lt; t_step
                          &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B>
                          &lt;&lt; std::setw(2)
                          &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                          &lt;&lt; std::right
                          &lt;&lt; a_step
                          &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.plt&quot;</FONT></B>;
  
        TecplotIO(mesh).write_equation_systems
          (file_name_tecplot.str(), es);
      }
  #endif
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_xda || param.output_xdr)
      mesh.renumber_nodes_and_elements();
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_xda)
      {
        mesh.write(numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;mesh&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xda&quot;</FONT></B>, param));
        es.write(numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;soln&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xda&quot;</FONT></B>, param),
                 <B><FONT COLOR="#5F9EA0">libMeshEnums</FONT></B>::WRITE, EquationSystems::WRITE_DATA |
                 <B><FONT COLOR="#5F9EA0">EquationSystems</FONT></B>::WRITE_ADDITIONAL_DATA);
      }
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_xdr)
      {
        mesh.write(numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;mesh&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xdr&quot;</FONT></B>, param));
        es.write(numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;soln&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xdr&quot;</FONT></B>, param),
                 <B><FONT COLOR="#5F9EA0">libMeshEnums</FONT></B>::WRITE, EquationSystems::WRITE_DATA |
                 <B><FONT COLOR="#5F9EA0">EquationSystems</FONT></B>::WRITE_ADDITIONAL_DATA);
      }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output_headers(FEMParameters &amp;param)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() != 0)
      <B><FONT COLOR="#A020F0">return</FONT></B>;
  
    start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_clocktime.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;vector_clocktime&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.run_simulation)
      {
        start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_activemesh.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;table_activemesh&quot;</FONT></B>);
        start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_solvesteps.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;table_solvesteps&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_tolerance)
          {
            start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_time.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;vector_time&quot;</FONT></B>);
            start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_timesteps.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;vector_timesteps&quot;</FONT></B>);
          }
        <B><FONT COLOR="#A020F0">if</FONT></B> (param.steadystate_tolerance)
          start_output(param.initial_timestep, <B><FONT COLOR="#BC8F8F">&quot;out_changerate.m&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;vector_changerate&quot;</FONT></B>);
      }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output_solvedata(EquationSystems &amp;es,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> newton_steps,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> krylov_steps,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> tv_sec,
                              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> tv_usec)
  {
    MeshBase &amp;mesh = es.get_mesh();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_active_elem = mesh.n_active_elem();
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_active_dofs = es.n_active_dofs();
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream activemesh (<B><FONT COLOR="#BC8F8F">&quot;out_activemesh.m&quot;</FONT></B>,
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
        activemesh.precision(17);
        activemesh &lt;&lt; (a_step + 1) &lt;&lt; <B><FONT COLOR="#BC8F8F">' '</FONT></B>
                   &lt;&lt; n_active_elem &lt;&lt; <B><FONT COLOR="#BC8F8F">' '</FONT></B>
                   &lt;&lt; n_active_dofs &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream solvesteps (<B><FONT COLOR="#BC8F8F">&quot;out_solvesteps.m&quot;</FONT></B>,
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
        solvesteps.precision(17);
        solvesteps &lt;&lt; newton_steps &lt;&lt; <B><FONT COLOR="#BC8F8F">' '</FONT></B>
                   &lt;&lt; krylov_steps &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream clocktime (<B><FONT COLOR="#BC8F8F">&quot;out_clocktime.m&quot;</FONT></B>,
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
        clocktime.precision(17);
        clocktime &lt;&lt; tv_sec &lt;&lt; <B><FONT COLOR="#BC8F8F">'.'</FONT></B> &lt;&lt; tv_usec &lt;&lt; std::endl;
      }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> write_output_footers(FEMParameters &amp;param)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream clocktime (<B><FONT COLOR="#BC8F8F">&quot;out_clocktime.m&quot;</FONT></B>,
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
        clocktime &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (param.run_simulation)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream activemesh (<B><FONT COLOR="#BC8F8F">&quot;out_activemesh.m&quot;</FONT></B>,
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
            activemesh &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream solvesteps (<B><FONT COLOR="#BC8F8F">&quot;out_solvesteps.m&quot;</FONT></B>,
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
            solvesteps &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_tolerance)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream times (<B><FONT COLOR="#BC8F8F">&quot;out_time.m&quot;</FONT></B>,
                  <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
                times &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream timesteps (<B><FONT COLOR="#BC8F8F">&quot;out_timesteps.m&quot;</FONT></B>,
                  <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
                timesteps &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
              }
            <B><FONT COLOR="#A020F0">if</FONT></B> (param.steadystate_tolerance)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream changerate (<B><FONT COLOR="#BC8F8F">&quot;out_changerate.m&quot;</FONT></B>,
                  <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::app | std::ios_base::out);
                changerate &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
              }
          }
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream complete (<B><FONT COLOR="#BC8F8F">&quot;complete&quot;</FONT></B>);
        complete &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;complete&quot;</FONT></B> &lt;&lt; std::endl;
      }
  }
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> defined(LIBMESH_HAVE_GMV) || defined(LIBMESH_HAVE_TECPLOT_API)
  <B><FONT COLOR="#228B22">void</FONT></B> write_error(EquationSystems &amp;es,
                   ErrorVector &amp;error,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_number,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_number,
                   FEMParameters &amp;param,
                   <B><FONT COLOR="#5F9EA0">std</FONT></B>::string error_type)
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  <B><FONT COLOR="#228B22">void</FONT></B> write_error(EquationSystems &amp;,
                   ErrorVector &amp;,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>,
                   FEMParameters &amp;,
                   <B><FONT COLOR="#5F9EA0">std</FONT></B>::string)
  #endif
  {
  #ifdef LIBMESH_HAVE_GMV
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.write_gmv_error)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream error_gmv;
        error_gmv &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;error.gmv.&quot;</FONT></B>
                  &lt;&lt; std::setw(3)
                  &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                  &lt;&lt; std::right
                  &lt;&lt; a_number
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.&quot;</FONT></B>
                  &lt;&lt; std::setw(2)
                  &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                  &lt;&lt; std::right
                  &lt;&lt; t_number
                  &lt;&lt; error_type;
  
        error.plot_error(error_gmv.str(), es.get_mesh());
      }
  #endif
  
  #ifdef LIBMESH_HAVE_TECPLOT_API
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.write_tecplot_error)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream error_tecplot;
        error_tecplot &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;error.plt.&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; a_number
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.&quot;</FONT></B>
                      &lt;&lt; std::setw(2)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_number
                      &lt;&lt; error_type;
  
        error.plot_error(error_tecplot.str(), es.get_mesh());
      }
  #endif
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> read_output(EquationSystems &amp;es,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step,
                   <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step,
                   <B><FONT COLOR="#5F9EA0">std</FONT></B>::string solution_type,
                   FEMParameters &amp;param)
  {
    MeshBase &amp;mesh = es.get_mesh();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string file_name_mesh, file_name_soln;
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_xda)
      {
        file_name_mesh = numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;mesh&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xda&quot;</FONT></B>, param);
        file_name_soln = numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;soln&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xda&quot;</FONT></B>, param);
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (param.output_xdr)
      {
        file_name_mesh = numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;mesh&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xdr&quot;</FONT></B>, param);
        file_name_soln = numbered_filename(t_step, a_step, solution_type, <B><FONT COLOR="#BC8F8F">&quot;soln&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;xdr&quot;</FONT></B>, param);
      }
  
    mesh.read(file_name_mesh);
  
    es.read(file_name_soln, libMeshEnums::READ,
            <B><FONT COLOR="#5F9EA0">EquationSystems</FONT></B>::READ_HEADER |
            <B><FONT COLOR="#5F9EA0">EquationSystems</FONT></B>::READ_DATA |
            <B><FONT COLOR="#5F9EA0">EquationSystems</FONT></B>::READ_ADDITIONAL_DATA);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != es.n_systems(); ++i)
      es.get_system&lt;FEMSystem&gt;(i).update();
  
    Real current_time = 0., current_timestep = 0.;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_tolerance)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ifstream times (<B><FONT COLOR="#BC8F8F">&quot;out_time.m&quot;</FONT></B>);
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ifstream timesteps (<B><FONT COLOR="#BC8F8F">&quot;out_timesteps.m&quot;</FONT></B>);
        <B><FONT COLOR="#A020F0">if</FONT></B> (times.is_open() &amp;&amp; timesteps.is_open())
          {
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> headersize = 25;
            <B><FONT COLOR="#228B22">char</FONT></B> header[headersize];
            timesteps.getline (header, headersize);
            <B><FONT COLOR="#A020F0">if</FONT></B> (strcmp(header, <B><FONT COLOR="#BC8F8F">&quot;vector_timesteps = [&quot;</FONT></B>) != 0)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Bad header in out_timesteps.m:&quot;</FONT></B> &lt;&lt; std::endl
                          &lt;&lt; header
                          &lt;&lt; std::endl;
                libmesh_error();
              }
  
            times.getline (header, headersize);
            <B><FONT COLOR="#A020F0">if</FONT></B> (strcmp(header, <B><FONT COLOR="#BC8F8F">&quot;vector_time = [&quot;</FONT></B>) != 0)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Bad header in out_time.m:&quot;</FONT></B> &lt;&lt; std::endl
                          &lt;&lt; header
                          &lt;&lt; std::endl;
                libmesh_error();
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != t_step; ++i)
              {
                <B><FONT COLOR="#A020F0">if</FONT></B> (!times.good())
                  libmesh_error();
                times &gt;&gt; current_time;
                timesteps &gt;&gt; current_timestep;
              }
            current_time += current_timestep;
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          libmesh_error();
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      current_time = t_step * param.deltat;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != es.n_systems(); ++i)
      es.get_system&lt;FEMSystem&gt;(i).time = current_time;
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
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (param.transient)
      {
        UnsteadySolver *innersolver;
        <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_core == <B><FONT COLOR="#BC8F8F">&quot;euler2&quot;</FONT></B>)
          {
            Euler2Solver *euler2solver =
              <B><FONT COLOR="#A020F0">new</FONT></B> Euler2Solver(system);
  
            euler2solver-&gt;theta = param.timesolver_theta;
            innersolver = euler2solver;
          }
        <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_core == <B><FONT COLOR="#BC8F8F">&quot;euler&quot;</FONT></B>)
          {
            EulerSolver *eulersolver =
              <B><FONT COLOR="#A020F0">new</FONT></B> EulerSolver(system);
  
            eulersolver-&gt;theta = param.timesolver_theta;
            innersolver = eulersolver;
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Don't recognize core TimeSolver type: &quot;</FONT></B>
                      &lt;&lt; param.timesolver_core &lt;&lt; std::endl;
            libmesh_error();
          }
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (param.timesolver_tolerance)
          {
            TwostepTimeSolver *timesolver =
              <B><FONT COLOR="#A020F0">new</FONT></B> TwostepTimeSolver(system);
  
            timesolver-&gt;max_growth       = param.timesolver_maxgrowth;
            timesolver-&gt;target_tolerance = param.timesolver_tolerance;
            timesolver-&gt;upper_tolerance  = param.timesolver_upper_tolerance;
            timesolver-&gt;component_norm   = SystemNorm(param.timesolver_norm);
  
            timesolver-&gt;core_time_solver =
              AutoPtr&lt;UnsteadySolver&gt;(innersolver);
            system.time_solver =
              AutoPtr&lt;UnsteadySolver&gt;(timesolver);
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          system.time_solver =
            AutoPtr&lt;TimeSolver&gt;(innersolver);
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      system.time_solver =
        AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
  
    system.time_solver-&gt;reduce_deltat_on_diffsolver_failure =
                                          param.deltat_reductions;
    system.time_solver-&gt;quiet           = param.time_solver_quiet;
  
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
  AutoPtr&lt;ErrorEstimator&gt; build_error_estimator(FEMParameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>)
  {
    AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
    error_estimator.reset(<B><FONT COLOR="#A020F0">new</FONT></B> KellyErrorEstimator);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> error_estimator;
  }
  
  AutoPtr&lt;ErrorEstimator&gt;
  build_error_estimator_component_wise
    (FEMParameters &amp;param,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt; &amp;term_weights,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt; &amp;primal_error_norm_type,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt; &amp;dual_error_norm_type)
  {
    AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
    AdjointResidualErrorEstimator *adjoint_residual_estimator = <B><FONT COLOR="#A020F0">new</FONT></B> AdjointResidualErrorEstimator;
  
    error_estimator.reset (adjoint_residual_estimator);
  
    PatchRecoveryErrorEstimator *p1 =
      <B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
    adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
  
    PatchRecoveryErrorEstimator *p2 =
      <B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
    adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);
  
    p1-&gt;set_patch_reuse(param.patch_reuse);
    p2-&gt;set_patch_reuse(param.patch_reuse);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> size = primal_error_norm_type.size();
  
    libmesh_assert_equal_to (size, dual_error_norm_type.size());
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != size; ++i)
      {
        adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(i, primal_error_norm_type[i]);
        adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(i, dual_error_norm_type[i]);
      }
  
    libmesh_assert_equal_to (size, term_weights.size());
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != term_weights.size(); ++i)
      {
        libmesh_assert_equal_to (size, term_weights[i].size());
        adjoint_residual_estimator-&gt;error_norm.set_weight(i, term_weights[i][i]);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j != size; ++j)
          <B><FONT COLOR="#A020F0">if</FONT></B> (i != j)
            adjoint_residual_estimator-&gt;error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> error_estimator;
  }
  
  AutoPtr&lt;ErrorEstimator&gt;
  build_weighted_error_estimator_component_wise
    (FEMParameters &amp;param,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt; &amp;term_weights,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt; &amp;primal_error_norm_type,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt; &amp;dual_error_norm_type,
     <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FEMFunctionBase&lt;Number&gt;*&gt; coupled_system_weight_functions)
  {
    AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
    AdjointResidualErrorEstimator *adjoint_residual_estimator = <B><FONT COLOR="#A020F0">new</FONT></B> AdjointResidualErrorEstimator;
  
    error_estimator.reset (adjoint_residual_estimator);
  
  
    WeightedPatchRecoveryErrorEstimator *p1 =
      <B><FONT COLOR="#A020F0">new</FONT></B> WeightedPatchRecoveryErrorEstimator;
    adjoint_residual_estimator-&gt;primal_error_estimator().reset(p1);
  
    PatchRecoveryErrorEstimator *p2 =
      <B><FONT COLOR="#A020F0">new</FONT></B> PatchRecoveryErrorEstimator;
    adjoint_residual_estimator-&gt;dual_error_estimator().reset(p2);
  
    p1-&gt;set_patch_reuse(param.patch_reuse);
    p2-&gt;set_patch_reuse(param.patch_reuse);
  
    p1-&gt;weight_functions.clear();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> size = primal_error_norm_type.size();
  
    libmesh_assert(coupled_system_weight_functions.size() == size);
    libmesh_assert(dual_error_norm_type.size() == size);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != size; ++i)
      {
        p1-&gt;weight_functions.push_back(coupled_system_weight_functions[i]);
        adjoint_residual_estimator-&gt;primal_error_estimator()-&gt;error_norm.set_type(i, primal_error_norm_type[i]);
        adjoint_residual_estimator-&gt;dual_error_estimator()-&gt;error_norm.set_type(i, dual_error_norm_type[i]);
      }
  
    libmesh_assert_equal_to (size, term_weights.size());
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != term_weights.size(); ++i)
      {
        libmesh_assert_equal_to (size, term_weights[i].size());
        adjoint_residual_estimator-&gt;error_norm.set_weight(i, term_weights[i][i]);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j != size; ++j)
          <B><FONT COLOR="#A020F0">if</FONT></B> (i != j)
            adjoint_residual_estimator-&gt;error_norm.set_off_diagonal_weight(i, j, term_weights[i][j]);
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> error_estimator;
  }
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_example_assert(libMesh::default_solver_package() != EIGEN_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
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
  
    Mesh mesh(init.comm(), param.dimension);
  
    AutoPtr&lt;MeshRefinement&gt; mesh_refinement =
      build_mesh_refinement(mesh, param);
  
    EquationSystems equation_systems (mesh);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Building mesh&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (!param.initial_timestep &amp;&amp; param.run_simulation)
      build_domain(mesh, param);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Building system&quot;</FONT></B> &lt;&lt; std::endl;
  
    FEMSystem &amp;system = build_system(equation_systems, infile, param);
  
    set_system_parameters(system, param);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initializing systems&quot;</FONT></B> &lt;&lt; std::endl;
  
    {
      CoupledSystemQoI qoi;
      qoi.assemble_qoi_sides = true;
      system.attach_qoi(&amp;qoi);
    }
  
    equation_systems.init ();
  
    mesh.print_info();
    equation_systems.print_info();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Starting adaptive loop&quot;</FONT></B>&lt;&lt;std::endl&lt;&lt;std::endl;
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step = 0;
  
    LinearSolver&lt;Number&gt; *linear_solver = system.get_linear_solver();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (; a_step &lt;= param.max_adaptivesteps; ++a_step)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B>(param.global_tolerance &gt; 0 &amp;&amp; param.nelem_target &gt; 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;We cant refine to both a non-zero tolerance and a target number of elements, EXITING adaptive loop. &quot;</FONT></B>&lt;&lt;std::endl&lt;&lt;std::endl;
            <B><FONT COLOR="#A020F0">break</FONT></B>;
          }
  
        linear_solver-&gt;reuse_preconditioner(false);
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; std::endl;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solving the forward problem&quot;</FONT></B> &lt;&lt;std::endl;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;We have &quot;</FONT></B> &lt;&lt; mesh.n_active_elem()
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B> &lt;&lt; equation_systems.n_active_dofs()
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl;
  
        system.solve();
  
        write_output(equation_systems, 0, a_step, <B><FONT COLOR="#BC8F8F">&quot;primal&quot;</FONT></B>, param);
  
        system.assemble_qoi();
  
        system.postprocess();
  
        Number QoI_0_computed = (dynamic_cast&lt;CoupledSystem&amp;&gt;(system)).get_QoI_value();
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The boundary QoI is &quot;</FONT></B> &lt;&lt; std::setprecision(17) &lt;&lt; QoI_0_computed &lt;&lt; std::endl &lt;&lt; std::endl;
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(a_step!=param.max_adaptivesteps)
          {
            NumericVector&lt;Number&gt; &amp;primal_solution = (*system.solution);
  
            system.assemble_qoi_sides = true;
  
            linear_solver-&gt;reuse_preconditioner(param.reuse_preconditioner);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solving the adjoint problem&quot;</FONT></B> &lt;&lt;std::endl;
            system.adjoint_solve();
  
  	system.set_adjoint_already_solved(true);
  
            NumericVector&lt;Number&gt; &amp;dual_solution = system.get_adjoint_solution();
            primal_solution.swap(dual_solution);
  
            write_output(equation_systems, 0, a_step, <B><FONT COLOR="#BC8F8F">&quot;adjoint&quot;</FONT></B>, param);
  
            primal_solution.swap(dual_solution);
  
            Real Pe = (dynamic_cast&lt;CoupledSystem&amp;&gt;(system)).get_Pe();
  
            ErrorVector error;
  
            ErrorVector error_non_pressure;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              primal_norm_type_vector_non_pressure;
            primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
            primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
            primal_norm_type_vector_non_pressure.push_back(L2);
            primal_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              dual_norm_type_vector_non_pressure;
            dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
            dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
            dual_norm_type_vector_non_pressure.push_back(L2);
            dual_norm_type_vector_non_pressure.push_back(H1_SEMINORM);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt;
              weights_matrix_non_pressure(system.n_vars(),
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt;(system.n_vars(), 0.0));
            weights_matrix_non_pressure[0][0] = 1.;
            weights_matrix_non_pressure[1][1] = 1.;
            weights_matrix_non_pressure[3][3] = 1./Pe;
  
            AutoPtr&lt;ErrorEstimator&gt; error_estimator_non_pressure =
              build_error_estimator_component_wise
                (param, weights_matrix_non_pressure,
                 primal_norm_type_vector_non_pressure,
                 dual_norm_type_vector_non_pressure);
  
            error_estimator_non_pressure-&gt;estimate_error(system, error_non_pressure);
  
            write_error(equation_systems, error_non_pressure, 0, a_step, param, <B><FONT COLOR="#BC8F8F">&quot;_non_pressure&quot;</FONT></B>);
  
            ErrorVector error_with_pressure;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              primal_norm_type_vector_with_pressure;
            primal_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
            primal_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
            primal_norm_type_vector_with_pressure.push_back(L2);
            primal_norm_type_vector_with_pressure.push_back(L2);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              dual_norm_type_vector_with_pressure;
            dual_norm_type_vector_with_pressure.push_back(H1_X_SEMINORM);
            dual_norm_type_vector_with_pressure.push_back(H1_Y_SEMINORM);
            dual_norm_type_vector_with_pressure.push_back(L2);
            dual_norm_type_vector_with_pressure.push_back(L2);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt;
              weights_matrix_with_pressure
                (system.n_vars(),
                 <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt;(system.n_vars(), 0.0));
            weights_matrix_with_pressure[0][2] = 1.;
  
            weights_matrix_with_pressure[1][2] = 1.;
  
            weights_matrix_with_pressure[2][0] = 1.;
            weights_matrix_with_pressure[2][1] = 1.;
  
            AutoPtr&lt;ErrorEstimator&gt; error_estimator_with_pressure =
              build_error_estimator_component_wise
                (param, weights_matrix_with_pressure,
                 primal_norm_type_vector_with_pressure,
                 dual_norm_type_vector_with_pressure);
  
            error_estimator_with_pressure-&gt;estimate_error(system, error_with_pressure);
  
            write_error(equation_systems, error_with_pressure, 0, a_step, param, <B><FONT COLOR="#BC8F8F">&quot;_with_pressure&quot;</FONT></B>);
  
  
            ErrorVector error_convection_diffusion_x;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              primal_norm_type_vector_convection_diffusion_x;
            primal_norm_type_vector_convection_diffusion_x.push_back(L2);
            primal_norm_type_vector_convection_diffusion_x.push_back(L2);
            primal_norm_type_vector_convection_diffusion_x.push_back(L2);
            primal_norm_type_vector_convection_diffusion_x.push_back(H1_X_SEMINORM);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              dual_norm_type_vector_convection_diffusion_x;
            dual_norm_type_vector_convection_diffusion_x.push_back(L2);
            dual_norm_type_vector_convection_diffusion_x.push_back(L2);
            dual_norm_type_vector_convection_diffusion_x.push_back(L2);
            dual_norm_type_vector_convection_diffusion_x.push_back(L2);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt;
              weights_matrix_convection_diffusion_x
                (system.n_vars(),
                 <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt;(system.n_vars(), 0.0));
  	  weights_matrix_convection_diffusion_x[0][3] = 1.;
            weights_matrix_convection_diffusion_x[3][3] = 1.;
  
  
            ConstFEMFunction&lt;Number&gt; identity(1);
  
            CoupledFEMFunctionsx convdiffx0(system, 0);
  	  CoupledFEMFunctionsx convdiffx3(system, 3);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FEMFunctionBase&lt;Number&gt; *&gt; coupled_system_weight_functions_x;
            coupled_system_weight_functions_x.push_back(&amp;convdiffx0);
            coupled_system_weight_functions_x.push_back(&amp;identity);
            coupled_system_weight_functions_x.push_back(&amp;identity);
            coupled_system_weight_functions_x.push_back(&amp;convdiffx3);
  
            AutoPtr&lt;ErrorEstimator&gt; error_estimator_convection_diffusion_x =
            build_weighted_error_estimator_component_wise
              (param, weights_matrix_convection_diffusion_x,
               primal_norm_type_vector_convection_diffusion_x,
               dual_norm_type_vector_convection_diffusion_x,
               coupled_system_weight_functions_x);
  
            error_estimator_convection_diffusion_x-&gt;estimate_error
              (system, error_convection_diffusion_x);
  
            write_error(equation_systems, error_convection_diffusion_x,
                        0, a_step, param, <B><FONT COLOR="#BC8F8F">&quot;_convection_diffusion_x&quot;</FONT></B>);
  
            ErrorVector error_convection_diffusion_y;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              primal_norm_type_vector_convection_diffusion_y;
            primal_norm_type_vector_convection_diffusion_y.push_back(L2);
            primal_norm_type_vector_convection_diffusion_y.push_back(L2);
            primal_norm_type_vector_convection_diffusion_y.push_back(L2);
            primal_norm_type_vector_convection_diffusion_y.push_back(H1_Y_SEMINORM);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;libMeshEnums::FEMNormType&gt;
              dual_norm_type_vector_convection_diffusion_y;
            dual_norm_type_vector_convection_diffusion_y.push_back(L2);
            dual_norm_type_vector_convection_diffusion_y.push_back(L2);
            dual_norm_type_vector_convection_diffusion_y.push_back(L2);
            dual_norm_type_vector_convection_diffusion_y.push_back(L2);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::vector&lt;Real&gt; &gt;
              weights_matrix_convection_diffusion_y
                (system.n_vars(), std::vector&lt;Real&gt;(system.n_vars(), 0.0));
  	  weights_matrix_convection_diffusion_y[1][3] = 1.;
            weights_matrix_convection_diffusion_y[3][3] = 1.;
  
            CoupledFEMFunctionsy convdiffy1(system, 1);
  	  CoupledFEMFunctionsy convdiffy3(system, 3);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FEMFunctionBase&lt;Number&gt; *&gt; coupled_system_weight_functions_y;
            coupled_system_weight_functions_y.push_back(&amp;identity);
            coupled_system_weight_functions_y.push_back(&amp;convdiffy1);
            coupled_system_weight_functions_y.push_back(&amp;identity);
            coupled_system_weight_functions_y.push_back(&amp;convdiffy3);
  
            AutoPtr&lt;ErrorEstimator&gt; error_estimator_convection_diffusion_y =
              build_weighted_error_estimator_component_wise
                (param, weights_matrix_convection_diffusion_y,
                 primal_norm_type_vector_convection_diffusion_y,
                 dual_norm_type_vector_convection_diffusion_y,
                 coupled_system_weight_functions_y);
  
            error_estimator_convection_diffusion_y-&gt;estimate_error(system, error_convection_diffusion_y);
  
            write_error(equation_systems, error_convection_diffusion_y, 0, a_step, param, <B><FONT COLOR="#BC8F8F">&quot;_convection_diffusion_y&quot;</FONT></B>);
  
            <B><FONT COLOR="#A020F0">if</FONT></B>(param.indicator_type == <B><FONT COLOR="#BC8F8F">&quot;adjoint_residual&quot;</FONT></B>)
              {
                error.resize(error_non_pressure.size());
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; error.size(); i++)
                  {
                    error[i] = error_non_pressure[i] + error_with_pressure[i] + error_convection_diffusion_x[i] + error_convection_diffusion_y[i];
                  }
              }
            <B><FONT COLOR="#A020F0">else</FONT></B>
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Using Kelly Estimator&quot;</FONT></B>&lt;&lt;std::endl;
                AutoPtr&lt;ErrorEstimator&gt; error_estimator =  build_error_estimator(param);
  
                error_estimator-&gt;estimate_error(system, error);
              }
  
            write_error(equation_systems, error, 0, a_step, param, <B><FONT COLOR="#BC8F8F">&quot;_total&quot;</FONT></B>);
  
  
  
            <B><FONT COLOR="#A020F0">if</FONT></B>(param.refine_uniformly)
              {
                mesh_refinement-&gt;uniformly_refine(1);
              }
            <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B>(param.global_tolerance &gt;= 0. &amp;&amp; param.nelem_target == 0.) <I><FONT COLOR="#B22222">// Refine based on reaching an error tolerance
</FONT></I>              {
                mesh_refinement-&gt;flag_elements_by_error_tolerance (error);
  
                mesh_refinement-&gt;refine_and_coarsen_elements();
              }
            <B><FONT COLOR="#A020F0">else</FONT></B> <I><FONT COLOR="#B22222">// Refine based on reaching a target number of elements
</FONT></I>              {
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
  
      }
  
  
    write_output_footers(param);
  
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file coupled_system.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parallel.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/zero_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;coupled_system.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> BdyFunction : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    BdyFunction (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var, <B><FONT COLOR="#228B22">int</FONT></B> sign)
      : _u_var(u_var), _v_var(v_var), _sign(sign)
      { <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;_initialized = true; }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Real = 0)
      { libmesh_not_implemented(); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                             <B><FONT COLOR="#228B22">const</FONT></B> Real,
                             DenseVector&lt;Number&gt;&amp; output)
      {
        output.resize(2);
        output.zero();
        <B><FONT COLOR="#228B22">const</FONT></B> Real y=p(1);
        output(_u_var) = (_sign)*((y-2) * (y-3));
        output(_v_var) = 0;
      }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() <B><FONT COLOR="#228B22">const</FONT></B>
      { <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (<B><FONT COLOR="#A020F0">new</FONT></B> BdyFunction(_u_var, _v_var, _sign)); }
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _u_var, _v_var;
    <B><FONT COLOR="#228B22">const</FONT></B> Real _sign;
  };
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystem::init_data ()
  {
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;coupled_system.in&quot;</FONT></B>);
    Peclet = infile(<B><FONT COLOR="#BC8F8F">&quot;Peclet&quot;</FONT></B>, 1.);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> pressure_p = infile(<B><FONT COLOR="#BC8F8F">&quot;pressure_p&quot;</FONT></B>, 1);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string fe_family = infile(<B><FONT COLOR="#BC8F8F">&quot;fe_family&quot;</FONT></B>, std::string(<B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>));
  
    libmesh_assert((pressure_p == 1) || (fe_family != <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>));
  
    FEFamily fefamily = Utility::string_to_enum&lt;FEFamily&gt;(fe_family);
  
    u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, static_cast&lt;Order&gt;(pressure_p+1),
  			      fefamily);
    v_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, static_cast&lt;Order&gt;(pressure_p+1),
  			      fefamily);
  
    p_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, static_cast&lt;Order&gt;(pressure_p),
  			      fefamily);
  
    C_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;C&quot;</FONT></B>, static_cast&lt;Order&gt;(pressure_p+1),
  			      fefamily);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(u_var);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(v_var);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(C_var);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;verify_analytic_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;verify_analytic_jacobians&quot;</FONT></B>, 0.);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_jacobians&quot;</FONT></B>, false);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_element_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_element_jacobians&quot;</FONT></B>, false);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type left_inlet_id = 0;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; left_inlet_bdy;
    left_inlet_bdy.insert(left_inlet_id);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type right_inlet_id = 1;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; right_inlet_bdy;
    right_inlet_bdy.insert(right_inlet_id);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type outlets_id = 2;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; outlets_bdy;
    outlets_bdy.insert(outlets_id);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type wall_id = 3;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; wall_bdy;
    wall_bdy.insert(wall_id);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; uv(1, u_var);
    uv.push_back(v_var);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; C_only(1, C_var);
  
    ZeroFunction&lt;Number&gt; zero;
    ConstFunction&lt;Number&gt; one(1);
  
    <B><FONT COLOR="#228B22">int</FONT></B> velocity_sign = 1;
    BdyFunction inflow_left(u_var, v_var, -velocity_sign);
    BdyFunction inflow_right(u_var, v_var, velocity_sign);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (wall_bdy, uv, &amp;zero));
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (left_inlet_bdy, uv, &amp;inflow_left));
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (left_inlet_bdy, C_only, &amp;one));
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (right_inlet_bdy, uv, &amp;inflow_right));
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
      (DirichletBoundary (right_inlet_bdy, C_only, &amp;zero));
  
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystem::init_context(DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    c.element_fe_var[u_var]-&gt;get_JxW();
    c.element_fe_var[u_var]-&gt;get_phi();
    c.element_fe_var[u_var]-&gt;get_dphi();
    c.element_fe_var[u_var]-&gt;get_xyz();
  
    c.element_fe_var[p_var]-&gt;get_phi();
  
    c.side_fe_var[u_var]-&gt;get_JxW();
    c.side_fe_var[u_var]-&gt;get_phi();
    c.side_fe_var[u_var]-&gt;get_xyz();
  }
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> CoupledSystem::element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                              DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW =
      c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi =
      c.element_fe_var[u_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi =
      c.element_fe_var[p_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
    libmesh_assert_equal_to (n_u_dofs, c.dof_indices_var[v_var].size());
  
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = *c.elem_subjacobians[u_var][u_var];
    DenseSubMatrix&lt;Number&gt; &amp;Kup = *c.elem_subjacobians[u_var][p_var];
    DenseSubVector&lt;Number&gt; &amp;Fu = *c.elem_subresiduals[u_var];
  
    DenseSubMatrix&lt;Number&gt; &amp;Kvv = *c.elem_subjacobians[v_var][v_var];
    DenseSubMatrix&lt;Number&gt; &amp;Kvp = *c.elem_subjacobians[v_var][p_var];
    DenseSubVector&lt;Number&gt; &amp;Fv = *c.elem_subresiduals[v_var];
  
    DenseSubMatrix&lt;Number&gt; &amp;KCu = *c.elem_subjacobians[C_var][u_var];
    DenseSubMatrix&lt;Number&gt; &amp;KCv = *c.elem_subjacobians[C_var][v_var];
    DenseSubMatrix&lt;Number&gt; &amp;KCC = *c.elem_subjacobians[C_var][C_var];
    DenseSubVector&lt;Number&gt; &amp;FC = *c.elem_subresiduals[C_var];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Number p = c.interior_value(p_var, qp),
  	u = c.interior_value(u_var, qp),
  	v = c.interior_value(v_var, qp);
        Gradient grad_u = c.interior_gradient(u_var, qp),
  	grad_v = c.interior_gradient(v_var, qp),
  	grad_C = c.interior_gradient(C_var, qp);
  
        NumberVectorValue U     (u,     v);
        <B><FONT COLOR="#228B22">const</FONT></B> Number C_x = grad_C(0);
        <B><FONT COLOR="#228B22">const</FONT></B> Number C_y = grad_C(1);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          {
            Fu(i) += JxW[qp] *
                     (p*dphi[i][qp](0) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>  		    (grad_u*dphi[i][qp]));            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  
            Fv(i) += JxW[qp] *
                     (p*dphi[i][qp](1) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>  		    (grad_v*dphi[i][qp]));            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  
  	  FC(i) += JxW[qp] *
  	           ( (U*grad_C)*phi[i][qp] +                <I><FONT COLOR="#B22222">// convection term
</FONT></I>  	            (1./Peclet)*(grad_C*dphi[i][qp]) );     <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.elem_solution_derivative)
              {
                libmesh_assert (c.elem_solution_derivative == 1.0);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>
  
                    Kvv(i,j) += JxW[qp] * (-(dphi[i][qp]*dphi[j][qp])); <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>
  
  		  KCu(i,j) += JxW[qp]* ( (phi[j][qp]*C_x)*phi[i][qp] ); <I><FONT COLOR="#B22222">/* convection term */</FONT></I>
  
  		  KCv(i,j) += JxW[qp]*( (phi[j][qp]*C_y)*phi[i][qp] );  <I><FONT COLOR="#B22222">/* convection term */</FONT></I>
  
  		  KCC(i,j) += JxW[qp]*
  		              ( (U*dphi[j][qp])*phi[i][qp] +      <I><FONT COLOR="#B22222">/* nonlinear term (convection) */</FONT></I>
  		              (1./Peclet)*(dphi[j][qp]*dphi[i][qp]) ); <I><FONT COLOR="#B22222">/* diffusion term */</FONT></I>
  		}
  
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
  		{
  		  Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
  		  Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
  		}
  	    }
  	}
  
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
        <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> CoupledSystem::element_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                         DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[u_var]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
      c.element_fe_var[u_var]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi =
      c.element_fe_var[p_var]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.dof_indices_var[p_var].size();
  
    DenseSubMatrix&lt;Number&gt; &amp;Kpu = *c.elem_subjacobians[p_var][u_var];
    DenseSubMatrix&lt;Number&gt; &amp;Kpv = *c.elem_subjacobians[p_var][v_var];
    DenseSubVector&lt;Number&gt; &amp;Fp = *c.elem_subresiduals[p_var];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Gradient grad_u = c.interior_gradient(u_var, qp),
  	grad_v = c.interior_gradient(v_var, qp);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          {
            Fp(i) += JxW[qp] * psi[i][qp] *
                     (grad_u(0) + grad_v(1));
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.elem_solution_derivative)
              {
                libmesh_assert (c.elem_solution_derivative == 1.0);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                    Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                  }
              }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystem::postprocess()
  {
  
    computed_QoI = 0.0;
  
    computed_QoI = System::qoi[0];
  }
  
  Number CoupledFEMFunctionsx::<B><FONT COLOR="#A020F0">operator</FONT></B>()(<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp; c, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  					<B><FONT COLOR="#228B22">const</FONT></B> Real <I><FONT COLOR="#B22222">/* time */</FONT></I>)
  {
    Number weight = 0.0;
  
    <B><FONT COLOR="#A020F0">switch</FONT></B>(var)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
        {
  	Gradient grad_C = c.point_gradient(3, p);
  
  	weight = grad_C(0);
        }
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
        {
  	Number u = c.point_value(0, p);
  
  	weight = u;
        }
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Wrong variable number&quot;</FONT></B>&lt;&lt;var&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot; passed to CoupledFEMFunctionsx object ! Quitting !&quot;</FONT></B>&lt;&lt;std::endl;
  	libmesh_error();
        }
  
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> weight;
  }
  
  Number CoupledFEMFunctionsy::<B><FONT COLOR="#A020F0">operator</FONT></B>()(<B><FONT COLOR="#228B22">const</FONT></B> FEMContext&amp; c, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  					<B><FONT COLOR="#228B22">const</FONT></B> Real <I><FONT COLOR="#B22222">/* time */</FONT></I>)
  {
    Number weight = 0.0;
  
    <B><FONT COLOR="#A020F0">switch</FONT></B>(var)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
        {
  	Gradient grad_C = c.point_gradient(3, p);
  
  	weight = grad_C(1);
        }
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
        {
  	Number v = c.point_value(1, p);
  
  	weight = v;
        }
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Wrong variable number &quot;</FONT></B>&lt;&lt;var&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot; passed to CoupledFEMFunctionsy object ! Quitting !&quot;</FONT></B>&lt;&lt;std::endl;
  	libmesh_error();
        }
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> weight;
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file domain.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;domain.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> build_domain (Mesh &amp;mesh, FEMParameters &amp;param)
  {
    mesh.read(param.domainfile);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Making elements 2nd order&quot;</FONT></B>&lt;&lt;std::endl;
  
    mesh.all_second_order();
  
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file femparameters.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  
  using namespace libMesh;
  
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
  
      GETPOT_INT_INPUT(initial_timestep);
      GETPOT_INT_INPUT(n_timesteps);
      GETPOT_INPUT(transient);
      GETPOT_INT_INPUT(deltat_reductions);
      GETPOT_INPUT(timesolver_core);
      GETPOT_INPUT(end_time);
      GETPOT_INPUT(deltat);
      GETPOT_INPUT(timesolver_theta);
      GETPOT_INPUT(timesolver_maxgrowth);
      GETPOT_INPUT(timesolver_upper_tolerance);
      GETPOT_INPUT(timesolver_tolerance);
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
      GETPOT_INPUT(fine_mesh_file_primal);
      GETPOT_INPUT(fine_mesh_soln_primal);
      GETPOT_INPUT(fine_mesh_file_adjoint);
      GETPOT_INPUT(fine_mesh_soln_adjoint);
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
      GETPOT_INT_INPUT(coarsecoarsenings);
      GETPOT_INT_INPUT(extrarefinements);
      GETPOT_INPUT(use_petsc_snes);
      GETPOT_INPUT(time_solver_quiet);
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
      GETPOT_INT_INPUT(initial_adaptivesteps);
      GETPOT_INT_INPUT(initial_sobolev_order);
      GETPOT_INT_INPUT(initial_extra_quadrature);
      GETPOT_INPUT(refine_uniformly);
      GETPOT_INPUT(indicator_type);
      GETPOT_INPUT(adjoint_residual_type);
      GETPOT_INPUT(patch_reuse);
      GETPOT_INT_INPUT(sobolev_order);
      GETPOT_INPUT(alternate_with_uniform_steps);
      GETPOT_INPUT(alternate_step_number);
      GETPOT_INPUT(component_wise_error);
      GETPOT_INPUT(compute_sensitivities);
      GETPOT_INPUT(compare_to_fine_solution_primal);
      GETPOT_INPUT(compare_to_fine_solution_adjoint);
      GETPOT_INPUT(do_forward_sensitivity);
      GETPOT_INPUT(do_adjoint_sensitivity);
      GETPOT_INPUT(postprocess_adjoint);
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
  
      GETPOT_INPUT(run_simulation);
      GETPOT_INPUT(run_postprocess);
  
      run_simulation              = input(<B><FONT COLOR="#BC8F8F">&quot;run_simulation&quot;</FONT></B>, run_simulation);
      run_postprocess             = input(<B><FONT COLOR="#BC8F8F">&quot;run_postprocess&quot;</FONT></B>, run_postprocess);
  
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
      GETPOT_INPUT(print_solution_norms);
      GETPOT_INPUT(print_solutions);
      GETPOT_INPUT(print_residual_norms);
      GETPOT_INPUT(print_residuals);
      GETPOT_INPUT(print_jacobian_norms);
      GETPOT_INPUT(print_jacobians);
  
  
  
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file H-qoi.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;H-qoi.h&quot;</FONT></B>
  
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystemQoI::init_qoi( std::vector&lt;Number&gt;&amp; sys_qoi)
  {
    sys_qoi.resize(1);
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystemQoI::side_qoi_derivative (DiffContext &amp;context,
  					 <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.side_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;phi = c.side_fe_var[0]-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;q_point = c.side_fe_var[0]-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[1].size();
  
    DenseSubVector&lt;Number&gt; &amp;Qu = *c.elem_qoi_subderivatives[0][0];
    DenseSubVector&lt;Number&gt; &amp;QC = *c.elem_qoi_subderivatives[0][3];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.side_qrule-&gt;n_points();
  
    Number u = 0. ;
    Number C = 0. ;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(c.has_side_boundary_id(2))
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
    	{
  	  Real x = q_point[qp](0);
  
  	  <B><FONT COLOR="#A020F0">if</FONT></B>(x &lt; 0.)
  	    {
  	      u = c.side_value(0,qp);
  	      C = c.side_value(3,qp);
  
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
  		{
  		  Qu(i) += JxW[qp] * -phi[i][qp] * C;
  		  QC(i) += JxW[qp] * phi[i][qp] * -u;
  		}
  	    } <I><FONT COLOR="#B22222">// end if
</FONT></I>  
    	} <I><FONT COLOR="#B22222">// end quadrature loop
</FONT></I>  
      } <I><FONT COLOR="#B22222">// end if on outlet
</FONT></I>  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> CoupledSystemQoI::side_qoi(DiffContext &amp;context, <B><FONT COLOR="#228B22">const</FONT></B> QoISet &amp; <I><FONT COLOR="#B22222">/* qois */</FONT></I>)
  {
  
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.side_fe_var[0]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt; &amp;q_point = c.side_fe_var[0]-&gt;get_xyz();
  
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.side_qrule-&gt;n_points();
  
    Number dQoI_0 = 0. ;
    Number u = 0. ;
    Number C = 0. ;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(c.has_side_boundary_id(2))
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
  	{
  	  Real x = q_point[qp](0);
  
  	  <B><FONT COLOR="#A020F0">if</FONT></B>(x &lt; 0.)
  	    {
  	      u = c.side_value(0,qp);
  	      C = c.side_value(3,qp);
  
  	      dQoI_0 += JxW[qp] * -u * C;
  	    } <I><FONT COLOR="#B22222">// end if
</FONT></I>  
    	} <I><FONT COLOR="#B22222">// end quadrature loop
</FONT></I>  
      } <I><FONT COLOR="#B22222">// end if on bdry
</FONT></I>  
    c.elem_qoi[0] += dQoI_0;
  
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file initial.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;initial.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;femparameters.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> read_initial_parameters()
  {
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> finish_initialization()
  {
  }
  
  
  
  Number initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; <I><FONT COLOR="#B22222">/* p */</FONT></I>,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 1.;
  
  }
  
  
  
  Gradient initial_grad(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; <I><FONT COLOR="#B22222">/* p */</FONT></I>,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; <I><FONT COLOR="#B22222">/* param */</FONT></I>,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> 0.;
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file output.C without comments: </h1> 
<pre> 
  #include &lt;fstream&gt;
  #include &lt;unistd.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;output.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> start_output(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> timesteps,
                    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string filename,
                    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string varname)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (timesteps)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ifstream infile (filename.c_str());
        <B><FONT COLOR="#228B22">char</FONT></B> buf[1025];
        <B><FONT COLOR="#A020F0">if</FONT></B> (!infile.is_open())
          libmesh_error();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != timesteps+1; ++i)
          infile.getline(buf, 1024);
        <B><FONT COLOR="#A020F0">if</FONT></B> (!infile.good())
          libmesh_error();
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> length = infile.tellg();
        infile.close();
  
        <B><FONT COLOR="#228B22">int</FONT></B> err = truncate(filename.c_str(), length);
        <B><FONT COLOR="#A020F0">if</FONT></B> (err != 0)
          libmesh_error();
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream outfile (filename.c_str(), std::ios_base::out);
        outfile &lt;&lt; varname &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; = [&quot;</FONT></B> &lt;&lt; std::endl;
      }
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex3'
***************************************************************
* Running Example adjoints_ex3:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Started /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex3/.libs/lt-example-devel
Building mesh
Making elements 2nd order
Building system
Initializing systems
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=99
    n_local_nodes()=27
  n_elem()=16
    n_local_elem()=4
    n_active_elem()=16
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "CoupledSystem"
    Type "Implicit"
    Variables={ "u" "v" } "p" "C" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=331
    n_local_dofs()=91
    n_constrained_dofs()=138
    n_local_constrained_dofs()=34
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 36.7855
      Average Off-Processor Bandwidth <= 5.91541
      Maximum  On-Processor Bandwidth <= 62
      Maximum Off-Processor Bandwidth <= 40
    DofMap Constraints
      Number of DoF Constraints = 138
      Number of Heterogenous Constraints= 5
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

Starting adaptive loop

Adaptive step 0
Solving the forward problem
We have 16 active elements and 193 active dofs.

  Nonlinear solver converged, step 1, residual reduction 9.03212e-14 < 1e-05
The boundary QoI is 0.083342124064158307

Solving the adjoint problem
Adaptive step 1
Solving the forward problem
We have 22 active elements and 259 active dofs.

  Nonlinear solver converged, step 0, residual reduction 5.6868816087085118e-07 < 1.0000000000000001e-05
The boundary QoI is 0.082286987672833115

Solving the adjoint problem
Adaptive step 2
Solving the forward problem
We have 31 active elements and 363 active dofs.

  Nonlinear solver converged, step 1, residual reduction 3.8583229881642115e-14 < 1.0000000000000001e-05
  Nonlinear solver converged, step 1, relative step size 7.5013647980960183e-07 < 1.0000000000000001e-05
The boundary QoI is 0.083169669189009962

Solving the adjoint problem
Adaptive step 3
Solving the forward problem
We have 43 active elements and 522 active dofs.

  Nonlinear solver converged, step 1, residual reduction 4.3999517315436118e-14 < 1.0000000000000001e-05
  Nonlinear solver converged, step 1, relative step size 6.565604269298289e-06 < 1.0000000000000001e-05
The boundary QoI is 0.083339483627459909

Solving the adjoint problem
Adaptive step 4
Solving the forward problem
We have 58 active elements and 671 active dofs.

  Nonlinear solver converged, step 0, residual reduction 3.9242900634249945e-06 < 1.0000000000000001e-05
The boundary QoI is 0.083044585890271971

Solving the adjoint problem
Adaptive step 5
Solving the forward problem
We have 76 active elements and 898 active dofs.

  Nonlinear solver converged, step 1, residual reduction 5.9607546528653657e-14 < 1.0000000000000001e-05
  Nonlinear solver converged, step 1, relative step size 2.0501238293058485e-07 < 1.0000000000000001e-05
The boundary QoI is 0.08333916169049338


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:47:50 2013                                                                          |
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
 ---------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.09678, Active time=1.97643                                                        |
 ---------------------------------------------------------------------------------------------------------------------
| Event                                   nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                   w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------------|
|                                                                                                                     |
|                                                                                                                     |
| AdjointResidualErrorEstimator                                                                                       |
|   estimate_error()                      20        0.0014      0.000069    0.9382      0.046909    0.07     47.47    |
|                                                                                                                     |
| DofMap                                                                                                              |
|   add_neighbors_to_send_list()          31        0.0024      0.000079    0.0046      0.000149    0.12     0.23     |
|   build_constraint_matrix()             298       0.0030      0.000010    0.0030      0.000010    0.15     0.15     |
|   build_sparsity()                      6         0.0024      0.000406    0.0075      0.001243    0.12     0.38     |
|   cnstrn_elem_mat_vec()                 105       0.0072      0.000068    0.0072      0.000068    0.36     0.36     |
|   constrain_elem_matrix()               44        0.0013      0.000030    0.0013      0.000030    0.07     0.07     |
|   constrain_elem_vector()               149       0.0007      0.000005    0.0007      0.000005    0.04     0.04     |
|   create_dof_constraints()              31        0.0156      0.000503    0.0607      0.001959    0.79     3.07     |
|   distribute_dofs()                     31        0.0056      0.000180    0.0288      0.000928    0.28     1.46     |
|   dof_indices()                         25520     0.2032      0.000008    0.2032      0.000008    10.28    10.28    |
|   enforce_constraints_exactly()         37        0.0068      0.000184    0.0068      0.000184    0.34     0.34     |
|   old_dof_indices()                     590       0.0067      0.000011    0.0067      0.000011    0.34     0.34     |
|   prepare_send_list()                   31        0.0001      0.000003    0.0001      0.000003    0.01     0.01     |
|   reinit()                              31        0.0060      0.000193    0.0060      0.000193    0.30     0.30     |
|                                                                                                                     |
| EquationSystems                                                                                                     |
|   build_discontinuous_solution_vector() 25        0.0019      0.000078    0.0033      0.000134    0.10     0.17     |
|   build_solution_vector()               11        0.0012      0.000113    0.0099      0.000899    0.06     0.50     |
|                                                                                                                     |
| FE                                                                                                                  |
|   compute_shape_functions()             21372     0.1065      0.000005    0.1065      0.000005    5.39     5.39     |
|   init_shape_functions()                5032      0.0401      0.000008    0.0401      0.000008    2.03     2.03     |
|   inverse_map()                         10827     0.0364      0.000003    0.0364      0.000003    1.84     1.84     |
|                                                                                                                     |
| FEMSystem                                                                                                           |
|   assemble_qoi()                        6         0.0031      0.000518    0.0486      0.008096    0.16     2.46     |
|   assemble_qoi_derivative()             5         0.0266      0.005329    0.0438      0.008754    1.35     2.21     |
|   assembly()                            10        0.0316      0.003164    0.0721      0.007211    1.60     3.65     |
|   assembly(get_jacobian)                5         0.0132      0.002642    0.0295      0.005898    0.67     1.49     |
|   assembly(get_residual)                10        0.0094      0.000940    0.0430      0.004302    0.48     2.18     |
|                                                                                                                     |
| FEMap                                                                                                               |
|   compute_affine_map()                  21372     0.0648      0.000003    0.0648      0.000003    3.28     3.28     |
|   compute_face_map()                    1896      0.0142      0.000007    0.0399      0.000021    0.72     2.02     |
|   init_face_shape_functions()           812       0.0013      0.000002    0.0013      0.000002    0.06     0.06     |
|   init_reference_to_physical_map()      5032      0.0538      0.000011    0.0538      0.000011    2.72     2.72     |
|                                                                                                                     |
| GMVIO                                                                                                               |
|   write_nodal_data()                    11        0.0138      0.001252    0.0145      0.001321    0.70     0.74     |
|                                                                                                                     |
| ImplicitSystem                                                                                                      |
|   adjoint_solve()                       5         0.0027      0.000547    0.1677      0.033542    0.14     8.49     |
|                                                                                                                     |
| LocationMap                                                                                                         |
|   find()                                400       0.0005      0.000001    0.0005      0.000001    0.02     0.02     |
|   init()                                10        0.0006      0.000062    0.0006      0.000062    0.03     0.03     |
|                                                                                                                     |
| Mesh                                                                                                                |
|   all_first_order()                     25        0.0016      0.000064    0.0016      0.000064    0.08     0.08     |
|   all_second_order()                    1         0.0001      0.000142    0.0001      0.000142    0.01     0.01     |
|   contract()                            5         0.0000      0.000009    0.0001      0.000024    0.00     0.01     |
|   find_neighbors()                      57        0.0050      0.000088    0.0121      0.000212    0.25     0.61     |
|   renumber_nodes_and_elem()             69        0.0005      0.000007    0.0005      0.000007    0.03     0.03     |
|                                                                                                                     |
| MeshCommunication                                                                                                   |
|   compute_hilbert_indices()             32        0.0030      0.000093    0.0030      0.000093    0.15     0.15     |
|   find_global_indices()                 32        0.0011      0.000034    0.0125      0.000391    0.06     0.63     |
|   parallel_sort()                       32        0.0021      0.000065    0.0054      0.000168    0.10     0.27     |
|                                                                                                                     |
| MeshOutput                                                                                                          |
|   write_equation_systems()              11        0.0003      0.000023    0.0255      0.002315    0.01     1.29     |
|                                                                                                                     |
| MeshRefinement                                                                                                      |
|   _coarsen_elements()                   10        0.0001      0.000005    0.0002      0.000015    0.00     0.01     |
|   _refine_elements()                    10        0.0007      0.000070    0.0033      0.000330    0.04     0.17     |
|   add_point()                           400       0.0006      0.000002    0.0012      0.000003    0.03     0.06     |
|   make_coarsening_compatible()          20        0.0005      0.000027    0.0005      0.000027    0.03     0.03     |
|   make_flags_parallel_consistent()      15        0.0008      0.000055    0.0055      0.000367    0.04     0.28     |
|   make_refinement_compatible()          20        0.0001      0.000004    0.0002      0.000011    0.00     0.01     |
|                                                                                                                     |
| MetisPartitioner                                                                                                    |
|   partition()                           32        0.0107      0.000333    0.0273      0.000854    0.54     1.38     |
|                                                                                                                     |
| NewtonSolver                                                                                                        |
|   solve()                               6         0.0940      0.015661    0.5123      0.085382    4.75     25.92    |
|                                                                                                                     |
| Parallel                                                                                                            |
|   allgather()                           169       0.0067      0.000039    0.0071      0.000042    0.34     0.36     |
|   broadcast()                           22        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   max(bool)                             52        0.0018      0.000036    0.0018      0.000036    0.09     0.09     |
|   max(scalar)                           3484      0.0150      0.000004    0.0150      0.000004    0.76     0.76     |
|   max(vector)                           787       0.0047      0.000006    0.0143      0.000018    0.24     0.72     |
|   max(vector<bool>)                     10        0.0001      0.000007    0.0003      0.000029    0.00     0.01     |
|   min(bool)                             4177      0.0168      0.000004    0.0168      0.000004    0.85     0.85     |
|   min(scalar)                           3367      0.4617      0.000137    0.4617      0.000137    23.36    23.36    |
|   min(vector)                           792       0.0052      0.000007    0.0154      0.000019    0.26     0.78     |
|   probe()                               1128      0.0072      0.000006    0.0072      0.000006    0.37     0.37     |
|   receive()                             1128      0.0023      0.000002    0.0096      0.000009    0.11     0.49     |
|   send()                                1128      0.0011      0.000001    0.0011      0.000001    0.06     0.06     |
|   send_receive()                        1192      0.0033      0.000003    0.0151      0.000013    0.17     0.76     |
|   sum()                                 354       0.0050      0.000014    0.3842      0.001085    0.25     19.44    |
|                                                                                                                     |
| Parallel::Request                                                                                                   |
|   wait()                                1128      0.0008      0.000001    0.0008      0.000001    0.04     0.04     |
|                                                                                                                     |
| Partitioner                                                                                                         |
|   set_node_processor_ids()              58        0.0035      0.000060    0.0195      0.000337    0.18     0.99     |
|   set_parent_processor_ids()            32        0.0003      0.000010    0.0003      0.000010    0.02     0.02     |
|                                                                                                                     |
| PatchRecoveryErrorEstimator                                                                                         |
|   estimate_error()                      120       0.1552      0.001293    0.5901      0.004917    7.85     29.86    |
|                                                                                                                     |
| PetscLinearSolver                                                                                                   |
|   solve()                               15        0.3861      0.025737    0.3861      0.025737    19.53    19.53    |
|                                                                                                                     |
| ProjectVector                                                                                                       |
|   operator()                            10        0.0028      0.000280    0.0135      0.001348    0.14     0.68     |
|                                                                                                                     |
| System                                                                                                              |
|   project_vector()                      10        0.0122      0.001222    0.0321      0.003208    0.62     1.62     |
|                                                                                                                     |
| WeightedPatchRecoveryErrorEstimator                                                                                 |
|   estimate_error()                      40        0.0787      0.001968    0.3467      0.008667    3.98     17.54    |
|                                                                                                                     |
| XdrIO                                                                                                               |
|   read()                                1         0.0003      0.000273    0.0003      0.000305    0.01     0.02     |
 ---------------------------------------------------------------------------------------------------------------------
| Totals:                                 113746    1.9764                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adjoints_ex3:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adjoints/adjoints_ex3'
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
