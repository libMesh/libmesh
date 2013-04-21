<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("fem_system_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file naviersystem.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fem_system.h"
        
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
        class NavierSystem : public FEMSystem
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          NavierSystem(EquationSystems& es,
                       const std::string& name_in,
                       const unsigned int number_in)
            : FEMSystem(es, name_in, number_in), Reynolds(1.), application(0) {}
        
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
          virtual bool side_constraint (bool request_jacobian,
                                        DiffContext& context);
        
</pre>
</div>
<div class = "comment">
Mass matrix part
</div>

<div class ="fragment">
<pre>
          virtual bool mass_residual (bool request_jacobian,
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
Indices for each variable;
</div>

<div class ="fragment">
<pre>
          unsigned int p_var, u_var, v_var, w_var;
        
</pre>
</div>
<div class = "comment">
The Reynolds number to solve for
</div>

<div class ="fragment">
<pre>
          Real Reynolds;
        
</pre>
</div>
<div class = "comment">
The application number controls what boundary conditions and/or
forcing functions are applied.  Current options are:
0 - discontinuous lid velociy driven cavity
1 - homogeneous Dirichlet BC with smooth forcing
</div>

<div class ="fragment">
<pre>
          unsigned int application;
        
</pre>
</div>
<div class = "comment">
Returns the value of a forcing function at point p.  This value
depends on which application is being used.
</div>

<div class ="fragment">
<pre>
          Point forcing(const Point& p);
        };
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file fem_system_ex1.C with comments: </h1> 
<div class = "comment">
<h1>FEMSystem Example 1 - Unsteady Navier-Stokes Equations with
FEMSystem</h1>

<br><br>This example shows how the transient nonlinear problem from
example 13 can be solved using the
DifferentiableSystem class framework


<br><br>C++ includes
</div>

<div class ="fragment">
<pre>
        #include &lt;iomanip&gt;
        
</pre>
</div>
<div class = "comment">
Basic include files
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/getpot.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/uniform_refinement_estimator.h"
        
</pre>
</div>
<div class = "comment">
The systems and solvers we may use
</div>

<div class ="fragment">
<pre>
        #include "naviersystem.h"
        #include "libmesh/diff_solver.h"
        #include "libmesh/euler_solver.h"
        #include "libmesh/steady_solver.h"
        
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
This example fails without at least double precision FP
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
          libmesh_example_assert(false, "--disable-singleprecision");
        #endif
        
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
        
</pre>
</div>
<div class = "comment">
Trilinos and Eigen solvers NaN by default here.
We'll skip this example for now.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc");
          libmesh_example_assert(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc");
        
</pre>
</div>
<div class = "comment">
Parse the input file
</div>

<div class ="fragment">
<pre>
          GetPot infile("fem_system_ex1.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
          const Real global_tolerance          = infile("global_tolerance", 0.);
          const unsigned int nelem_target      = infile("n_elements", 400);
          const bool transient                 = infile("transient", true);
          const Real deltat                    = infile("deltat", 0.005);
          unsigned int n_timesteps             = infile("n_timesteps", 20);
          const unsigned int write_interval    = infile("write_interval", 5);
          const unsigned int coarsegridsize    = infile("coarsegridsize", 1);
          const unsigned int coarserefinements = infile("coarserefinements", 0);
          const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 10);
          const unsigned int dim               = infile("dimension", 2);
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D/3D support");
        
</pre>
</div>
<div class = "comment">
We have only defined 2 and 3 dimensional problems
</div>

<div class ="fragment">
<pre>
          libmesh_assert (dim == 2 || dim == 3);
        
</pre>
</div>
<div class = "comment">
Create a mesh, with dimension to be overridden later, distributed
across the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
And an object to refine it
</div>

<div class ="fragment">
<pre>
          MeshRefinement mesh_refinement(mesh);
          mesh_refinement.coarsen_by_parents() = true;
          mesh_refinement.absolute_global_tolerance() = global_tolerance;
          mesh_refinement.nelem_target() = nelem_target;
          mesh_refinement.refine_fraction() = 0.3;
          mesh_refinement.coarsen_fraction() = 0.3;
          mesh_refinement.coarsen_threshold() = 0.1;
        
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
          if (dim == 2)
            MeshTools::Generation::build_square (mesh,
                                                 coarsegridsize,
                                                 coarsegridsize,
                                                 0., 1.,
                                                 0., 1.,
                                                 QUAD9);
          else if (dim == 3)
            MeshTools::Generation::build_cube (mesh,
                                               coarsegridsize,
                                               coarsegridsize,
                                               coarsegridsize,
                                               0., 1.,
                                               0., 1.,
                                               0., 1.,
                                               HEX27);
        
          mesh_refinement.uniformly_refine(coarserefinements);
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
          mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
</pre>
</div>
<div class = "comment">
Declare the system "Navier-Stokes" and its variables.
</div>

<div class ="fragment">
<pre>
          NavierSystem & system =
            equation_systems.add_system&lt;NavierSystem&gt; ("Navier-Stokes");
        
</pre>
</div>
<div class = "comment">
Solve this as a time-dependent or steady system
</div>

<div class ="fragment">
<pre>
          if (transient)
            system.time_solver =
              AutoPtr&lt;TimeSolver&gt;(new EulerSolver(system));
          else
            {
              system.time_solver =
                AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
              libmesh_assert_equal_to (n_timesteps, 1);
            }
        
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
Set the time stepping options
</div>

<div class ="fragment">
<pre>
          system.deltat = deltat;
        
</pre>
</div>
<div class = "comment">
And the nonlinear solver options
</div>

<div class ="fragment">
<pre>
          DiffSolver &solver = *(system.time_solver-&gt;diff_solver().get());
          solver.quiet = infile("solver_quiet", true);
          solver.verbose = !solver.quiet;
          solver.max_nonlinear_iterations =
            infile("max_nonlinear_iterations", 15);
          solver.relative_step_tolerance =
            infile("relative_step_tolerance", 1.e-3);
          solver.relative_residual_tolerance =
            infile("relative_residual_tolerance", 0.0);
          solver.absolute_residual_tolerance =
            infile("absolute_residual_tolerance", 0.0);
        
</pre>
</div>
<div class = "comment">
And the linear solver options
</div>

<div class ="fragment">
<pre>
          solver.max_linear_iterations =
            infile("max_linear_iterations", 50000);
          solver.initial_linear_tolerance =
            infile("initial_linear_tolerance", 1.e-3);
        
</pre>
</div>
<div class = "comment">
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Now we begin the timestep loop to compute the time-accurate
solution of the equations.
</div>

<div class ="fragment">
<pre>
          for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
            {
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; "\n\nSolving time step " &lt;&lt; t_step &lt;&lt; ", time = "
                        &lt;&lt; system.time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Adaptively solve the timestep
</div>

<div class ="fragment">
<pre>
              unsigned int a_step = 0;
              for (; a_step != max_adaptivesteps; ++a_step)
                {
                  system.solve();
        
                  system.postprocess();
        
                  ErrorVector error;
        
                  AutoPtr&lt;ErrorEstimator&gt; error_estimator;
        
</pre>
</div>
<div class = "comment">
To solve to a tolerance in this problem we
need a better estimator than Kelly
</div>

<div class ="fragment">
<pre>
                  if (global_tolerance != 0.)
                    {
</pre>
</div>
<div class = "comment">
We can't adapt to both a tolerance and a mesh
size at once
</div>

<div class ="fragment">
<pre>
                      libmesh_assert_equal_to (nelem_target, 0);
        
                      UniformRefinementEstimator *u =
                        new UniformRefinementEstimator;
        
</pre>
</div>
<div class = "comment">
The lid-driven cavity problem isn't in H1, so
lets estimate L2 error
</div>

<div class ="fragment">
<pre>
                      u-&gt;error_norm = L2;
        
                      error_estimator.reset(u);
                    }
                  else
                    {
</pre>
</div>
<div class = "comment">
If we aren't adapting to a tolerance we need a
target mesh size
</div>

<div class ="fragment">
<pre>
                      libmesh_assert_greater (nelem_target, 0);
        
</pre>
</div>
<div class = "comment">
Kelly is a lousy estimator to use for a problem
not in H1 - if we were doing more than a few
timesteps we'd need to turn off or limit the
maximum level of our adaptivity eventually
</div>

<div class ="fragment">
<pre>
                      error_estimator.reset(new KellyErrorEstimator);
                    }
        
</pre>
</div>
<div class = "comment">
Calculate error based on u and v (and w?) but not p
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;Real&gt; weights(2,1.0);  // u, v
                  if (dim == 3)
                    weights.push_back(1.0);          // w
                  weights.push_back(0.0);            // p
</pre>
</div>
<div class = "comment">
Keep the same default norm type.
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;FEMNormType&gt;
        	    norms(1, error_estimator-&gt;error_norm.type(0));
        	  error_estimator-&gt;error_norm = SystemNorm(norms, weights);
        
                  error_estimator-&gt;estimate_error(system, error);
        
</pre>
</div>
<div class = "comment">
Print out status at each adaptive step.
</div>

<div class ="fragment">
<pre>
                  Real global_error = error.l2_norm();
                  std::cout &lt;&lt; "Adaptive step " &lt;&lt; a_step &lt;&lt; ": " &lt;&lt; std::endl;
                  if (global_tolerance != 0.)
                    std::cout &lt;&lt; "Global_error = " &lt;&lt; global_error
                              &lt;&lt; std::endl;
                  if (global_tolerance != 0.)
                    std::cout &lt;&lt; "Worst element error = " &lt;&lt; error.maximum()
                              &lt;&lt; ", mean = " &lt;&lt; error.mean() &lt;&lt; std::endl;
        
                  if (global_tolerance != 0.)
                    {
</pre>
</div>
<div class = "comment">
If we've reached our desired tolerance, we
don't need any more adaptive steps
</div>

<div class ="fragment">
<pre>
                      if (global_error &lt; global_tolerance)
                        break;
                      mesh_refinement.flag_elements_by_error_tolerance(error);
                    }
                  else
                    {
</pre>
</div>
<div class = "comment">
If flag_elements_by_nelem_target returns true, this
should be our last adaptive step.
</div>

<div class ="fragment">
<pre>
                      if (mesh_refinement.flag_elements_by_nelem_target(error))
                        {
                          mesh_refinement.refine_and_coarsen_elements();
                          equation_systems.reinit();
                          a_step = max_adaptivesteps;
                          break;
                        }
                    }
        
</pre>
</div>
<div class = "comment">
Carry out the adaptive mesh refinement/coarsening
</div>

<div class ="fragment">
<pre>
                  mesh_refinement.refine_and_coarsen_elements();
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
              if (a_step == max_adaptivesteps)
                {
                  system.solve();
        
                  system.postprocess();
                }
        
</pre>
</div>
<div class = "comment">
Advance to the next timestep in a transient problem
</div>

<div class ="fragment">
<pre>
              system.time_solver-&gt;advance_timestep();
        
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
Write out this timestep if we're requested to
</div>

<div class ="fragment">
<pre>
              if ((t_step+1)%write_interval == 0)
                {
                  std::ostringstream file_name;
        
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
                  file_name &lt;&lt; "out_"
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_step+1
                            &lt;&lt; ".e";
        
                  ExodusII_IO(mesh).write_timestep(file_name.str(),
        					   equation_systems,
        					   1, /* This number indicates how many time steps
        						 are being written to the file */
        					   system.time);
                }
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
            }
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
<br><br><br> <h1> The source file naviersystem.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/getpot.h"
        
        #include "naviersystem.h"
        
        #include "libmesh/boundary_info.h"
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/mesh.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/zero_function.h"
        
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
Boundary conditions for the 3D test case
</div>

<div class ="fragment">
<pre>
        class BdyFunction : public FunctionBase&lt;Number&gt;
        {
        public:
          BdyFunction (unsigned int u_var,
                       unsigned int v_var,
                       unsigned int w_var,
                       Real Reynolds)
            : _u_var(u_var), _v_var(v_var), _w_var(w_var), _Re(Reynolds)
            { this-&gt;_initialized = true; }
        
          virtual Number operator() (const Point&, const Real = 0)
            { libmesh_not_implemented(); }
        
          virtual void operator() (const Point& p,
                                   const Real,
                                   DenseVector&lt;Number&gt;& output)
            {
              output.zero();
              const Real x=p(0), y=p(1), z=p(2);
              output(_u_var) = (_Re+1)*(y*y + z*z);
              output(_v_var) = (_Re+1)*(x*x + z*z);
              output(_w_var) = (_Re+1)*(x*x + y*y);
            }
        
          virtual AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() const
            { return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (new BdyFunction(_u_var, _v_var, _w_var, _Re)); }
        
        private:
          const unsigned int _u_var, _v_var, _w_var;
          const Real _Re;
        };
        
        
        void NavierSystem::init_data ()
        {
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Check the input file for Reynolds number, application type,
approximation type
</div>

<div class ="fragment">
<pre>
          GetPot infile("navier.in");
          Reynolds = infile("Reynolds", 1.);
          application = infile("application", 0);
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
          v_var = this-&gt;add_variable ("v",
        			      static_cast&lt;Order&gt;(pressure_p+1),
        			      fefamily);
        
          if (dim == 3)
            w_var = this-&gt;add_variable ("w",
        			        static_cast&lt;Order&gt;(pressure_p+1),
        			        fefamily);
          else
            w_var = u_var;
        
</pre>
</div>
<div class = "comment">
Add the pressure variable "p". This will
be approximated with a first-order basis,
providing an LBB-stable pressure-velocity pair.
</div>

<div class ="fragment">
<pre>
          p_var = this-&gt;add_variable ("p",
        			      static_cast&lt;Order&gt;(pressure_p),
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
          if (dim == 3)
            this-&gt;time_evolving(w_var);
        
</pre>
</div>
<div class = "comment">
Useful debugging options
Set verify_analytic_jacobians to 1e-6 to use
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
          const boundary_id_type top_id = (dim==3) ? 5 : 2;
        
          std::set&lt;boundary_id_type&gt; top_bdys;
          top_bdys.insert(top_id);
        
          const boundary_id_type all_ids[6] = {0, 1, 2, 3, 4, 5};
          std::set&lt;boundary_id_type&gt; all_bdys(all_ids, all_ids+(dim*2));
        
          std::set&lt;boundary_id_type&gt; nontop_bdys = all_bdys;
          nontop_bdys.erase(top_id);
        
          std::vector&lt;unsigned int&gt; u_only(1, u_var);
          std::vector&lt;unsigned int&gt; vw(1, v_var), uvw(1, u_var);
          uvw.push_back(v_var);
          if (dim == 3)
            {
              vw.push_back(w_var);
              uvw.push_back(w_var);
            }
        
          ZeroFunction&lt;Number&gt; zero;
          ConstFunction&lt;Number&gt; one(1);
</pre>
</div>
<div class = "comment">
For lid-driven cavity, set u=1,v=w=0 on the lid and u=v=w=0 elsewhere
</div>

<div class ="fragment">
<pre>
          if (application == 0)
            {
              this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (top_bdys, u_only, &one));
              this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (top_bdys, vw, &zero));
              this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (nontop_bdys, uvw, &zero));
            }
</pre>
</div>
<div class = "comment">
For forcing with zero wall velocity, set homogeneous Dirichlet BCs
</div>

<div class ="fragment">
<pre>
          else if (application == 1)
            {
              this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (all_bdys, uvw, &zero));
            }
</pre>
</div>
<div class = "comment">
For 3D test case with quadratic velocity field, set that field on walls
</div>

<div class ="fragment">
<pre>
          else if (application == 2)
            {
              BdyFunction bdy(u_var,v_var,w_var,Reynolds);
              this-&gt;get_dof_map().add_dirichlet_boundary
                (DirichletBoundary (all_bdys, uvw, &bdy));
            }
        
</pre>
</div>
<div class = "comment">
Do the parent's initialization after variables and boundary constraints are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        }
        
        
        
        void NavierSystem::init_context(DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
          FEBase* u_elem_fe;
          FEBase* p_elem_fe;
          FEBase* u_side_fe;
        
          c.get_element_fe( u_var, u_elem_fe );
          c.get_element_fe( p_var, p_elem_fe );
          c.get_side_fe( u_var, u_side_fe );
        
</pre>
</div>
<div class = "comment">
We should prerequest all the data
we will need to build the linear system.
</div>

<div class ="fragment">
<pre>
          u_elem_fe-&gt;get_JxW();
          u_elem_fe-&gt;get_phi();
          u_elem_fe-&gt;get_dphi();
          u_elem_fe-&gt;get_xyz();
        
          p_elem_fe-&gt;get_phi();
        
          u_side_fe-&gt;get_JxW();
          u_side_fe-&gt;get_phi();
          u_side_fe-&gt;get_xyz();
        }
        
        
        bool NavierSystem::element_time_derivative (bool request_jacobian,
                                                    DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
          FEBase* u_elem_fe;
          FEBase* p_elem_fe;
        
          c.get_element_fe( u_var, u_elem_fe );
          c.get_element_fe( p_var, p_elem_fe );
        
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = u_elem_fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = u_elem_fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = u_elem_fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi = p_elem_fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
Physical location of the quadrature points
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt;& qpoint = u_elem_fe-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_p_dofs = c.get_dof_indices( p_var ).size();
          const unsigned int n_u_dofs = c.get_dof_indices( u_var ).size();
          libmesh_assert_equal_to (n_u_dofs, c.get_dof_indices( v_var ).size());
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
          DenseSubMatrix&lt;Number&gt; &Kuu = c.get_elem_jacobian( u_var, u_var );
          DenseSubMatrix&lt;Number&gt; &Kvv = c.get_elem_jacobian( v_var, v_var );
          DenseSubMatrix&lt;Number&gt; &Kww = c.get_elem_jacobian( w_var, w_var );
          DenseSubMatrix&lt;Number&gt; &Kuv = c.get_elem_jacobian( u_var, v_var );
          DenseSubMatrix&lt;Number&gt; &Kuw = c.get_elem_jacobian( u_var, w_var );
          DenseSubMatrix&lt;Number&gt; &Kvu = c.get_elem_jacobian( v_var, u_var );
          DenseSubMatrix&lt;Number&gt; &Kvw = c.get_elem_jacobian( v_var, w_var );
          DenseSubMatrix&lt;Number&gt; &Kwu = c.get_elem_jacobian( w_var, u_var );
          DenseSubMatrix&lt;Number&gt; &Kwv = c.get_elem_jacobian( w_var, v_var );
          DenseSubMatrix&lt;Number&gt; &Kup = c.get_elem_jacobian( u_var, p_var );
          DenseSubMatrix&lt;Number&gt; &Kvp = c.get_elem_jacobian( v_var, p_var );
          DenseSubMatrix&lt;Number&gt; &Kwp = c.get_elem_jacobian( w_var, p_var );
          DenseSubVector&lt;Number&gt; &Fu = c.get_elem_residual( u_var );
          DenseSubVector&lt;Number&gt; &Fv = c.get_elem_residual( v_var );
          DenseSubVector&lt;Number&gt; &Fw = c.get_elem_residual( w_var );
        
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
Compute the solution & its gradient at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Number p = c.interior_value(p_var, qp),
                     u = c.interior_value(u_var, qp),
                     v = c.interior_value(v_var, qp),
                     w = c.interior_value(w_var, qp);
              Gradient grad_u = c.interior_gradient(u_var, qp),
                       grad_v = c.interior_gradient(v_var, qp),
                       grad_w = c.interior_gradient(w_var, qp);
        
</pre>
</div>
<div class = "comment">
Definitions for convenience.  It is sometimes simpler to do a
dot product if you have the full vector at your disposal.
</div>

<div class ="fragment">
<pre>
              NumberVectorValue U     (u,     v);
              if (dim == 3)
                U(2) = w;
              const Number  u_x = grad_u(0);
              const Number  u_y = grad_u(1);
              const Number  u_z = (dim == 3)?grad_u(2):0;
              const Number  v_x = grad_v(0);
              const Number  v_y = grad_v(1);
              const Number  v_z = (dim == 3)?grad_v(2):0;
              const Number  w_x = (dim == 3)?grad_w(0):0;
              const Number  w_y = (dim == 3)?grad_w(1):0;
              const Number  w_z = (dim == 3)?grad_w(2):0;
        
</pre>
</div>
<div class = "comment">
Value of the forcing function at this quadrature point
</div>

<div class ="fragment">
<pre>
              Point f = this-&gt;forcing(qpoint[qp]);
        
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
                  Fu(i) += JxW[qp] *
                           (-Reynolds*(U*grad_u)*phi[i][qp] + // convection term
                            p*dphi[i][qp](0) -                // pressure term
        		    (grad_u*dphi[i][qp]) +            // diffusion term
        		    f(0)*phi[i][qp]                   // forcing function
        		    );
        
        
                  Fv(i) += JxW[qp] *
                           (-Reynolds*(U*grad_v)*phi[i][qp] + // convection term
                            p*dphi[i][qp](1) -                // pressure term
        		    (grad_v*dphi[i][qp]) +            // diffusion term
        		    f(1)*phi[i][qp]                   // forcing function
        		    );
        
        
                  if (dim == 3)
                  Fw(i) += JxW[qp] *
                           (-Reynolds*(U*grad_w)*phi[i][qp] + // convection term
                            p*dphi[i][qp](2) -                // pressure term
        		    (grad_w*dphi[i][qp]) +            // diffusion term
        		    f(2)*phi[i][qp]                   // forcing function
        		    );
        
        
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
                      libmesh_assert_equal_to (c.elem_solution_derivative, 1.0);
        
</pre>
</div>
<div class = "comment">
Matrix contributions for the uu and vv couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        {
                          Kuu(i,j) += JxW[qp] *
           /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
           /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
           /* Newton term     */       Reynolds*u_x*phi[i][qp]*phi[j][qp]);
        
                          Kuv(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*u_y*phi[i][qp]*phi[j][qp];
        
                          Kvv(i,j) += JxW[qp] *
           /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
           /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
           /* Newton term     */       Reynolds*v_y*phi[i][qp]*phi[j][qp]);
        
                          Kvu(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*v_x*phi[i][qp]*phi[j][qp];
                          if (dim == 3)
                            {
                              Kww(i,j) += JxW[qp] *
           /* convection term */          (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
           /* diffusion term  */           (dphi[i][qp]*dphi[j][qp]) -
           /* Newton term     */           Reynolds*w_z*phi[i][qp]*phi[j][qp]);
                              Kuw(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*u_z*phi[i][qp]*phi[j][qp];
                              Kvw(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*v_z*phi[i][qp]*phi[j][qp];
                              Kwu(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*w_x*phi[i][qp]*phi[j][qp];
                              Kwv(i,j) += JxW[qp] *
           /* Newton term     */      -Reynolds*w_y*phi[i][qp]*phi[j][qp];
                            }
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
                          if (dim == 3)
                            Kwp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);
                        }
                    }
                }
            } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
        
        
        bool NavierSystem::element_constraint (bool request_jacobian,
                                               DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
          FEBase* u_elem_fe;
          FEBase* p_elem_fe;
        
          c.get_element_fe( u_var, u_elem_fe );
          c.get_element_fe( p_var, p_elem_fe );
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weight for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = u_elem_fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = u_elem_fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi = p_elem_fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.get_dof_indices( u_var ).size();
          const unsigned int n_p_dofs = c.get_dof_indices( p_var ).size();
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
          DenseSubMatrix&lt;Number&gt; &Kpu = c.get_elem_jacobian( p_var, u_var );
          DenseSubMatrix&lt;Number&gt; &Kpv = c.get_elem_jacobian( p_var, v_var );
          DenseSubMatrix&lt;Number&gt; &Kpw = c.get_elem_jacobian( p_var, w_var );
          DenseSubVector&lt;Number&gt; &Fp = c.get_elem_residual( p_var );
        
</pre>
</div>
<div class = "comment">
Add the constraint given by the continuity equation
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
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
                       grad_v = c.interior_gradient(v_var, qp),
                       grad_w = c.interior_gradient(w_var, qp);
        
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
                  if (dim == 3)
                    Fp(i) += JxW[qp] * psi[i][qp] *
                             (grad_w(2));
        
                  if (request_jacobian && c.get_elem_solution_derivative())
                    {
                      libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
        
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        {
                          Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                          Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                          if (dim == 3)
                            Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2);
                        }
                    }
                }
            } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
        
        
        bool NavierSystem::side_constraint (bool request_jacobian,
                                            DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
          FEBase* p_elem_fe;
        
          c.get_element_fe( p_var, p_elem_fe );
        
</pre>
</div>
<div class = "comment">
Pin p = 0 at the origin
</div>

<div class ="fragment">
<pre>
          const Point zero(0.,0.);
        
          if (c.elem-&gt;contains_point(zero))
            {
</pre>
</div>
<div class = "comment">
The pressure penalty value.  \f$ \frac{1}{\epsilon} \f$
</div>

<div class ="fragment">
<pre>
              const Real penalty = 1.e9;
        
              DenseSubMatrix&lt;Number&gt; &Kpp = c.get_elem_jacobian( p_var, p_var );
              DenseSubVector&lt;Number&gt; &Fp = c.get_elem_residual( p_var );
              const unsigned int n_p_dofs = c.get_dof_indices( p_var ).size();
        
              Number p = c.point_value(p_var, zero);
              Number p_value = 0.;
        
              unsigned int dim = get_mesh().mesh_dimension();
              FEType fe_type = p_elem_fe-&gt;get_fe_type();
              Point p_master = FEInterface::inverse_map(dim, fe_type, c.elem, zero);
        
              std::vector&lt;Real&gt; point_phi(n_p_dofs);
              for (unsigned int i=0; i != n_p_dofs; i++)
                {
                  point_phi[i] = FEInterface::shape(dim, fe_type, c.elem, i, p_master);
                }
        
              for (unsigned int i=0; i != n_p_dofs; i++)
                {
                  Fp(i) += penalty * (p - p_value) * point_phi[i];
                  if (request_jacobian && c.get_elem_solution_derivative())
                    {
                      libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
        
                      for (unsigned int j=0; j != n_p_dofs; j++)
        		Kpp(i,j) += penalty * point_phi[i] * point_phi[j];
                    }
                }
            }
        
          return request_jacobian;
        }
        
        
</pre>
</div>
<div class = "comment">
We override the default mass_residual function,
because in the non-dimensionalized Navier-Stokes equations
the time derivative of velocity has a Reynolds number coefficient.
Alternatively we could divide the whole equation by
Reynolds number (or choose a more complicated non-dimensionalization
of time), but this gives us an opportunity to demonstrate overriding
FEMSystem::mass_residual()
</div>

<div class ="fragment">
<pre>
        bool NavierSystem::mass_residual (bool request_jacobian,
                                          DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
          FEBase* u_elem_fe;
        
          c.get_element_fe( u_var, u_elem_fe );
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Element Jacobian * quadrature weight for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = u_elem_fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = u_elem_fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubVector&lt;Number&gt; &Fu = c.get_elem_residual( u_var );
          DenseSubVector&lt;Number&gt; &Fv = c.get_elem_residual( v_var );
          DenseSubVector&lt;Number&gt; &Fw = c.get_elem_residual( w_var );
          DenseSubMatrix&lt;Number&gt; &Kuu = c.get_elem_jacobian( u_var, u_var );
          DenseSubMatrix&lt;Number&gt; &Kvv = c.get_elem_jacobian( v_var, v_var );
          DenseSubMatrix&lt;Number&gt; &Kww = c.get_elem_jacobian( w_var, w_var );
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in velocity
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.get_dof_indices( u_var ).size();
        
          unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
          for (unsigned int qp = 0; qp != n_qpoints; ++qp)
            {
              Number u = c.interior_value(u_var, qp),
                     v = c.interior_value(v_var, qp),
                     w = c.interior_value(w_var, qp);
        
</pre>
</div>
<div class = "comment">
We pull as many calculations as possible outside of loops
</div>

<div class ="fragment">
<pre>
              Number JxWxRe   = JxW[qp] * Reynolds;
              Number JxWxRexU = JxWxRe * u;
              Number JxWxRexV = JxWxRe * v;
              Number JxWxRexW = JxWxRe * w;
        
              for (unsigned int i = 0; i != n_u_dofs; ++i)
                {
                  Fu(i) += JxWxRexU * phi[i][qp];
                  Fv(i) += JxWxRexV * phi[i][qp];
                  if (dim == 3)
                    Fw(i) += JxWxRexW * phi[i][qp];
        
                  if (request_jacobian && c.get_elem_solution_derivative())
                    {
                      libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
        
                      Number JxWxRexPhiI = JxWxRe * phi[i][qp];
                      Number JxWxRexPhiII = JxWxRexPhiI * phi[i][qp];
                      Kuu(i,i) += JxWxRexPhiII;
                      Kvv(i,i) += JxWxRexPhiII;
                      if (dim == 3)
                        Kww(i,i) += JxWxRexPhiII;
        
</pre>
</div>
<div class = "comment">
The mass matrix is symmetric, so we calculate
one triangle and add it to both upper and lower
triangles
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j = i+1; j != n_u_dofs; ++j)
                        {
                          Number Kij = JxWxRexPhiI * phi[j][qp];
                          Kuu(i,j) += Kij;
                          Kuu(j,i) += Kij;
                          Kvv(i,j) += Kij;
                          Kvv(j,i) += Kij;
                          if (dim == 3)
                            {
                              Kww(i,j) += Kij;
                              Kww(j,i) += Kij;
                            }
                        }
                    }
                }
            }
        
          return request_jacobian;
        }
        
        
        
        void NavierSystem::postprocess()
        {
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
          Point p(1./3., 1./3.);
          Number u = point_value(u_var, p),
                 v = point_value(v_var, p),
                 w = (dim == 3) ? point_value(w_var, p) : 0;
        
          std::cout &lt;&lt; "u(1/3,1/3) = ("
                    &lt;&lt; u &lt;&lt; ", "
                    &lt;&lt; v &lt;&lt; ", "
                    &lt;&lt; w &lt;&lt; ")" &lt;&lt; std::endl;
        }
        
        
        
        
        Point NavierSystem::forcing(const Point& p)
        {
          switch (application)
            {
</pre>
</div>
<div class = "comment">
lid driven cavity
</div>

<div class ="fragment">
<pre>
            case 0:
              {
</pre>
</div>
<div class = "comment">
No forcing
</div>

<div class ="fragment">
<pre>
                return Point(0.,0.,0.);
              }
        
</pre>
</div>
<div class = "comment">
Homogeneous Dirichlet BCs + sinusoidal forcing
</div>

<div class ="fragment">
<pre>
            case 1:
              {
        	const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
This assumes your domain is defined on [0,1]^dim.
</div>

<div class ="fragment">
<pre>
                Point f;
        
</pre>
</div>
<div class = "comment">
Counter-Clockwise vortex in x-y plane
</div>

<div class ="fragment">
<pre>
                if (dim==2)
        	  {
        	    f(0) =  std::sin(2.*libMesh::pi*p(1));
        	    f(1) = -std::sin(2.*libMesh::pi*p(0));
        	  }
        
</pre>
</div>
<div class = "comment">
Counter-Clockwise vortex in x-z plane
</div>

<div class ="fragment">
<pre>
                else if (dim==3)
        	  {
        	    f(0) =  std::sin(2.*libMesh::pi*p(1));
        	    f(1) = 0.;
        	    f(2) = -std::sin(2.*libMesh::pi*p(0));
        	  }
        
        	return f;
              }
        
</pre>
</div>
<div class = "comment">
3D test case with quadratic velocity and linear pressure field
</div>

<div class ="fragment">
<pre>
            case 2:
              {
</pre>
</div>
<div class = "comment">
This problem doesn't make sense in 1D...
</div>

<div class ="fragment">
<pre>
                libmesh_assert_not_equal_to (1, this-&gt;get_mesh().mesh_dimension());
        	const Real x=p(0), y=p(1), z=p(2);
        	const Real
        	  u=(Reynolds+1)*(y*y + z*z),
        	  v=(Reynolds+1)*(x*x + z*z),
        	  w=(Reynolds+1)*(x*x + y*y);
        
        	if (this-&gt;get_mesh().mesh_dimension() == 2)
        	  return 2*Reynolds*(Reynolds+1)*Point(v*y,
        					       u*x);
        	else
        	  return 2*Reynolds*(Reynolds+1)*Point(v*y + w*z,
        					       u*x + w*z,
        					       u*x + v*y);
              }
        
            default:
              {
        	libmesh_error();
              }
            }
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file naviersystem.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> NavierSystem : <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    NavierSystem(EquationSystems&amp; es,
                 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
                 <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
      : FEMSystem(es, name_in, number_in), Reynolds(1.), application(0) {}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                          DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                     DiffContext&amp; context);
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> side_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                  DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> mass_residual (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> postprocess ();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var, u_var, v_var, w_var;
  
    Real Reynolds;
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> application;
  
    Point forcing(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p);
  };
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file fem_system_ex1.C without comments: </h1> 
<pre> 
  
  #include &lt;iomanip&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/uniform_refinement_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;naviersystem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/euler_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_example_assert(libMesh::default_solver_package() != EIGEN_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
    libmesh_example_assert(libMesh::default_solver_package() != TRILINOS_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;fem_system_ex1.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real global_tolerance          = infile(<B><FONT COLOR="#BC8F8F">&quot;global_tolerance&quot;</FONT></B>, 0.);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nelem_target      = infile(<B><FONT COLOR="#BC8F8F">&quot;n_elements&quot;</FONT></B>, 400);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> transient                 = infile(<B><FONT COLOR="#BC8F8F">&quot;transient&quot;</FONT></B>, true);
    <B><FONT COLOR="#228B22">const</FONT></B> Real deltat                    = infile(<B><FONT COLOR="#BC8F8F">&quot;deltat&quot;</FONT></B>, 0.005);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesteps             = infile(<B><FONT COLOR="#BC8F8F">&quot;n_timesteps&quot;</FONT></B>, 20);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> write_interval    = infile(<B><FONT COLOR="#BC8F8F">&quot;write_interval&quot;</FONT></B>, 5);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarsegridsize    = infile(<B><FONT COLOR="#BC8F8F">&quot;coarsegridsize&quot;</FONT></B>, 1);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> coarserefinements = infile(<B><FONT COLOR="#BC8F8F">&quot;coarserefinements&quot;</FONT></B>, 0);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_adaptivesteps = infile(<B><FONT COLOR="#BC8F8F">&quot;max_adaptivesteps&quot;</FONT></B>, 10);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim               = infile(<B><FONT COLOR="#BC8F8F">&quot;dimension&quot;</FONT></B>, 2);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
    libmesh_assert (dim == 2 || dim == 3);
  
    Mesh mesh(init.comm());
  
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.coarsen_by_parents() = true;
    mesh_refinement.absolute_global_tolerance() = global_tolerance;
    mesh_refinement.nelem_target() = nelem_target;
    mesh_refinement.refine_fraction() = 0.3;
    mesh_refinement.coarsen_fraction() = 0.3;
    mesh_refinement.coarsen_threshold() = 0.1;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 2)
      <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                           coarsegridsize,
                                           coarsegridsize,
                                           0., 1.,
                                           0., 1.,
                                           QUAD9);
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                         coarsegridsize,
                                         coarsegridsize,
                                         coarsegridsize,
                                         0., 1.,
                                         0., 1.,
                                         0., 1.,
                                         HEX27);
  
    mesh_refinement.uniformly_refine(coarserefinements);
  
    mesh.print_info();
  
    EquationSystems equation_systems (mesh);
  
    NavierSystem &amp; system =
      equation_systems.add_system&lt;NavierSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (transient)
      system.time_solver =
        AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> EulerSolver(system));
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        system.time_solver =
          AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
        libmesh_assert_equal_to (n_timesteps, 1);
      }
  
    equation_systems.init ();
  
    system.deltat = deltat;
  
    DiffSolver &amp;solver = *(system.time_solver-&gt;diff_solver().get());
    solver.quiet = infile(<B><FONT COLOR="#BC8F8F">&quot;solver_quiet&quot;</FONT></B>, true);
    solver.verbose = !solver.quiet;
    solver.max_nonlinear_iterations =
      infile(<B><FONT COLOR="#BC8F8F">&quot;max_nonlinear_iterations&quot;</FONT></B>, 15);
    solver.relative_step_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;relative_step_tolerance&quot;</FONT></B>, 1.e-3);
    solver.relative_residual_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;relative_residual_tolerance&quot;</FONT></B>, 0.0);
    solver.absolute_residual_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;absolute_residual_tolerance&quot;</FONT></B>, 0.0);
  
    solver.max_linear_iterations =
      infile(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 50000);
    solver.initial_linear_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;initial_linear_tolerance&quot;</FONT></B>, 1.e-3);
  
    equation_systems.print_info();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=0; t_step != n_timesteps; ++t_step)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n\nSolving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B>
                  &lt;&lt; system.time &lt;&lt; std::endl;
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step = 0;
        <B><FONT COLOR="#A020F0">for</FONT></B> (; a_step != max_adaptivesteps; ++a_step)
          {
            system.solve();
  
            system.postprocess();
  
            ErrorVector error;
  
            AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              {
                libmesh_assert_equal_to (nelem_target, 0);
  
                UniformRefinementEstimator *u =
                  <B><FONT COLOR="#A020F0">new</FONT></B> UniformRefinementEstimator;
  
                u-&gt;error_norm = L2;
  
                error_estimator.reset(u);
              }
            <B><FONT COLOR="#A020F0">else</FONT></B>
              {
                libmesh_assert_greater (nelem_target, 0);
  
                error_estimator.reset(<B><FONT COLOR="#A020F0">new</FONT></B> KellyErrorEstimator);
              }
  
  	  <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; weights(2,1.0);  <I><FONT COLOR="#B22222">// u, v
</FONT></I>            <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
              weights.push_back(1.0);          <I><FONT COLOR="#B22222">// w
</FONT></I>            weights.push_back(0.0);            <I><FONT COLOR="#B22222">// p
</FONT></I>  	  <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FEMNormType&gt;
  	    norms(1, error_estimator-&gt;error_norm.type(0));
  	  error_estimator-&gt;error_norm = SystemNorm(norms, weights);
  
            error_estimator-&gt;estimate_error(system, error);
  
            Real global_error = error.l2_norm();
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;: &quot;</FONT></B> &lt;&lt; std::endl;
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Global_error = &quot;</FONT></B> &lt;&lt; global_error
                        &lt;&lt; std::endl;
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Worst element error = &quot;</FONT></B> &lt;&lt; error.maximum()
                        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, mean = &quot;</FONT></B> &lt;&lt; error.mean() &lt;&lt; std::endl;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              {
                <B><FONT COLOR="#A020F0">if</FONT></B> (global_error &lt; global_tolerance)
                  <B><FONT COLOR="#A020F0">break</FONT></B>;
                mesh_refinement.flag_elements_by_error_tolerance(error);
              }
            <B><FONT COLOR="#A020F0">else</FONT></B>
              {
                <B><FONT COLOR="#A020F0">if</FONT></B> (mesh_refinement.flag_elements_by_nelem_target(error))
                  {
                    mesh_refinement.refine_and_coarsen_elements();
                    equation_systems.reinit();
                    a_step = max_adaptivesteps;
                    <B><FONT COLOR="#A020F0">break</FONT></B>;
                  }
              }
  
            mesh_refinement.refine_and_coarsen_elements();
            equation_systems.reinit();
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Refined mesh to &quot;</FONT></B>
                      &lt;&lt; mesh.n_active_elem()
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                      &lt;&lt; equation_systems.n_active_dofs()
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl;
          }
        <B><FONT COLOR="#A020F0">if</FONT></B> (a_step == max_adaptivesteps)
          {
            system.solve();
  
            system.postprocess();
          }
  
        system.time_solver-&gt;advance_timestep();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
        <B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_step+1
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
  
            ExodusII_IO(mesh).write_timestep(file_name.str(),
  					   equation_systems,
  					   1, <I><FONT COLOR="#B22222">/* This number indicates how many time steps
  						 are being written to the file */</FONT></I>
  					   system.time);
          }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>      }
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file naviersystem.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;naviersystem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/zero_function.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> BdyFunction : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    BdyFunction (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var,
                 <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var,
                 <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> w_var,
                 Real Reynolds)
      : _u_var(u_var), _v_var(v_var), _w_var(w_var), _Re(Reynolds)
      { <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;_initialized = true; }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Real = 0)
      { libmesh_not_implemented(); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                             <B><FONT COLOR="#228B22">const</FONT></B> Real,
                             DenseVector&lt;Number&gt;&amp; output)
      {
        output.zero();
        <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
        output(_u_var) = (_Re+1)*(y*y + z*z);
        output(_v_var) = (_Re+1)*(x*x + z*z);
        output(_w_var) = (_Re+1)*(x*x + y*y);
      }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() <B><FONT COLOR="#228B22">const</FONT></B>
      { <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (<B><FONT COLOR="#A020F0">new</FONT></B> BdyFunction(_u_var, _v_var, _w_var, _Re)); }
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _u_var, _v_var, _w_var;
    <B><FONT COLOR="#228B22">const</FONT></B> Real _Re;
  };
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> NavierSystem::init_data ()
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;navier.in&quot;</FONT></B>);
    Reynolds = infile(<B><FONT COLOR="#BC8F8F">&quot;Reynolds&quot;</FONT></B>, 1.);
    application = infile(<B><FONT COLOR="#BC8F8F">&quot;application&quot;</FONT></B>, 0);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> pressure_p = infile(<B><FONT COLOR="#BC8F8F">&quot;pressure_p&quot;</FONT></B>, 1);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string fe_family = infile(<B><FONT COLOR="#BC8F8F">&quot;fe_family&quot;</FONT></B>, std::string(<B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>));
  
    libmesh_assert((pressure_p == 1) || (fe_family != <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>));
  
    FEFamily fefamily = Utility::string_to_enum&lt;FEFamily&gt;(fe_family);
  
    u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, static_cast&lt;Order&gt;(pressure_p+1),
  			      fefamily);
    v_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>,
  			      static_cast&lt;Order&gt;(pressure_p+1),
  			      fefamily);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      w_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>,
  			        static_cast&lt;Order&gt;(pressure_p+1),
  			        fefamily);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      w_var = u_var;
  
    p_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>,
  			      static_cast&lt;Order&gt;(pressure_p),
  			      fefamily);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(u_var);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(v_var);
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(w_var);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;verify_analytic_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;verify_analytic_jacobians&quot;</FONT></B>, 0.);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_jacobians&quot;</FONT></B>, false);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_element_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_element_jacobians&quot;</FONT></B>, false);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type top_id = (dim==3) ? 5 : 2;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; top_bdys;
    top_bdys.insert(top_id);
  
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type all_ids[6] = {0, 1, 2, 3, 4, 5};
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; all_bdys(all_ids, all_ids+(dim*2));
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; nontop_bdys = all_bdys;
    nontop_bdys.erase(top_id);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; u_only(1, u_var);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; vw(1, v_var), uvw(1, u_var);
    uvw.push_back(v_var);
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      {
        vw.push_back(w_var);
        uvw.push_back(w_var);
      }
  
    ZeroFunction&lt;Number&gt; zero;
    ConstFunction&lt;Number&gt; one(1);
    <B><FONT COLOR="#A020F0">if</FONT></B> (application == 0)
      {
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (top_bdys, u_only, &amp;one));
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (top_bdys, vw, &amp;zero));
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (nontop_bdys, uvw, &amp;zero));
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (application == 1)
      {
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (all_bdys, uvw, &amp;zero));
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (application == 2)
      {
        BdyFunction bdy(u_var,v_var,w_var,Reynolds);
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary
          (DirichletBoundary (all_bdys, uvw, &amp;bdy));
      }
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> NavierSystem::init_context(DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEBase* u_elem_fe;
    FEBase* p_elem_fe;
    FEBase* u_side_fe;
  
    c.get_element_fe( u_var, u_elem_fe );
    c.get_element_fe( p_var, p_elem_fe );
    c.get_side_fe( u_var, u_side_fe );
  
    u_elem_fe-&gt;get_JxW();
    u_elem_fe-&gt;get_phi();
    u_elem_fe-&gt;get_dphi();
    u_elem_fe-&gt;get_xyz();
  
    p_elem_fe-&gt;get_phi();
  
    u_side_fe-&gt;get_JxW();
    u_side_fe-&gt;get_phi();
    u_side_fe-&gt;get_xyz();
  }
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> NavierSystem::element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                              DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEBase* u_elem_fe;
    FEBase* p_elem_fe;
  
    c.get_element_fe( u_var, u_elem_fe );
    c.get_element_fe( p_var, p_elem_fe );
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = u_elem_fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = u_elem_fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = u_elem_fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = p_elem_fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qpoint = u_elem_fe-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.get_dof_indices( p_var ).size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.get_dof_indices( u_var ).size();
    libmesh_assert_equal_to (n_u_dofs, c.get_dof_indices( v_var ).size());
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = c.get_elem_jacobian( u_var, u_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kvv = c.get_elem_jacobian( v_var, v_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kww = c.get_elem_jacobian( w_var, w_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kuv = c.get_elem_jacobian( u_var, v_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kuw = c.get_elem_jacobian( u_var, w_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kvu = c.get_elem_jacobian( v_var, u_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kvw = c.get_elem_jacobian( v_var, w_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kwu = c.get_elem_jacobian( w_var, u_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kwv = c.get_elem_jacobian( w_var, v_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kup = c.get_elem_jacobian( u_var, p_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kvp = c.get_elem_jacobian( v_var, p_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kwp = c.get_elem_jacobian( w_var, p_var );
    DenseSubVector&lt;Number&gt; &amp;Fu = c.get_elem_residual( u_var );
    DenseSubVector&lt;Number&gt; &amp;Fv = c.get_elem_residual( v_var );
    DenseSubVector&lt;Number&gt; &amp;Fw = c.get_elem_residual( w_var );
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Number p = c.interior_value(p_var, qp),
               u = c.interior_value(u_var, qp),
               v = c.interior_value(v_var, qp),
               w = c.interior_value(w_var, qp);
        Gradient grad_u = c.interior_gradient(u_var, qp),
                 grad_v = c.interior_gradient(v_var, qp),
                 grad_w = c.interior_gradient(w_var, qp);
  
        NumberVectorValue U     (u,     v);
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
          U(2) = w;
        <B><FONT COLOR="#228B22">const</FONT></B> Number  u_x = grad_u(0);
        <B><FONT COLOR="#228B22">const</FONT></B> Number  u_y = grad_u(1);
        <B><FONT COLOR="#228B22">const</FONT></B> Number  u_z = (dim == 3)?grad_u(2):0;
        <B><FONT COLOR="#228B22">const</FONT></B> Number  v_x = grad_v(0);
        <B><FONT COLOR="#228B22">const</FONT></B> Number  v_y = grad_v(1);
        <B><FONT COLOR="#228B22">const</FONT></B> Number  v_z = (dim == 3)?grad_v(2):0;
        <B><FONT COLOR="#228B22">const</FONT></B> Number  w_x = (dim == 3)?grad_w(0):0;
        <B><FONT COLOR="#228B22">const</FONT></B> Number  w_y = (dim == 3)?grad_w(1):0;
        <B><FONT COLOR="#228B22">const</FONT></B> Number  w_z = (dim == 3)?grad_w(2):0;
  
        Point f = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;forcing(qpoint[qp]);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          {
            Fu(i) += JxW[qp] *
                     (-Reynolds*(U*grad_u)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                      p*dphi[i][qp](0) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>  		    (grad_u*dphi[i][qp]) +            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  		    f(0)*phi[i][qp]                   <I><FONT COLOR="#B22222">// forcing function
</FONT></I>  		    );
  
  
            Fv(i) += JxW[qp] *
                     (-Reynolds*(U*grad_v)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                      p*dphi[i][qp](1) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>  		    (grad_v*dphi[i][qp]) +            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  		    f(1)*phi[i][qp]                   <I><FONT COLOR="#B22222">// forcing function
</FONT></I>  		    );
  
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
            Fw(i) += JxW[qp] *
                     (-Reynolds*(U*grad_w)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                      p*dphi[i][qp](2) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>  		    (grad_w*dphi[i][qp]) +            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  		    f(2)*phi[i][qp]                   <I><FONT COLOR="#B22222">// forcing function
</FONT></I>  		    );
  
  
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.elem_solution_derivative)
              {
                libmesh_assert_equal_to (c.elem_solution_derivative, 1.0);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* convection term */</FONT></I>      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
     <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>       (dphi[i][qp]*dphi[j][qp]) -
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>       Reynolds*u_x*phi[i][qp]*phi[j][qp]);
  
                    Kuv(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*u_y*phi[i][qp]*phi[j][qp];
  
                    Kvv(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* convection term */</FONT></I>      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
     <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>       (dphi[i][qp]*dphi[j][qp]) -
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>       Reynolds*v_y*phi[i][qp]*phi[j][qp]);
  
                    Kvu(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*v_x*phi[i][qp]*phi[j][qp];
                    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                      {
                        Kww(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* convection term */</FONT></I>          (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
     <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>           (dphi[i][qp]*dphi[j][qp]) -
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>           Reynolds*w_z*phi[i][qp]*phi[j][qp]);
                        Kuw(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*u_z*phi[i][qp]*phi[j][qp];
                        Kvw(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*v_z*phi[i][qp]*phi[j][qp];
                        Kwu(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*w_x*phi[i][qp]*phi[j][qp];
                        Kwv(i,j) += JxW[qp] *
     <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*w_y*phi[i][qp]*phi[j][qp];
                      }
                  }
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
                  {
                    Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                    Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
                    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                      Kwp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);
                  }
              }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> NavierSystem::element_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                         DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEBase* u_elem_fe;
    FEBase* p_elem_fe;
  
    c.get_element_fe( u_var, u_elem_fe );
    c.get_element_fe( p_var, p_elem_fe );
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = u_elem_fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = u_elem_fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = p_elem_fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.get_dof_indices( u_var ).size();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.get_dof_indices( p_var ).size();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
    DenseSubMatrix&lt;Number&gt; &amp;Kpu = c.get_elem_jacobian( p_var, u_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kpv = c.get_elem_jacobian( p_var, v_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kpw = c.get_elem_jacobian( p_var, w_var );
    DenseSubVector&lt;Number&gt; &amp;Fp = c.get_elem_residual( p_var );
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Gradient grad_u = c.interior_gradient(u_var, qp),
                 grad_v = c.interior_gradient(v_var, qp),
                 grad_w = c.interior_gradient(w_var, qp);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          {
            Fp(i) += JxW[qp] * psi[i][qp] *
                     (grad_u(0) + grad_v(1));
            <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
              Fp(i) += JxW[qp] * psi[i][qp] *
                       (grad_w(2));
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.get_elem_solution_derivative())
              {
                libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  {
                    Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                    Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                      Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2);
                  }
              }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> NavierSystem::side_constraint (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                      DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEBase* p_elem_fe;
  
    c.get_element_fe( p_var, p_elem_fe );
  
    <B><FONT COLOR="#228B22">const</FONT></B> Point zero(0.,0.);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (c.elem-&gt;contains_point(zero))
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e9;
  
        DenseSubMatrix&lt;Number&gt; &amp;Kpp = c.get_elem_jacobian( p_var, p_var );
        DenseSubVector&lt;Number&gt; &amp;Fp = c.get_elem_residual( p_var );
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = c.get_dof_indices( p_var ).size();
  
        Number p = c.point_value(p_var, zero);
        Number p_value = 0.;
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = get_mesh().mesh_dimension();
        FEType fe_type = p_elem_fe-&gt;get_fe_type();
        Point p_master = FEInterface::inverse_map(dim, fe_type, c.elem, zero);
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; point_phi(n_p_dofs);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          {
            point_phi[i] = FEInterface::shape(dim, fe_type, c.elem, i, p_master);
          }
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_p_dofs; i++)
          {
            Fp(i) += penalty * (p - p_value) * point_phi[i];
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.get_elem_solution_derivative())
              {
                libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_p_dofs; j++)
  		Kpp(i,j) += penalty * point_phi[i] * point_phi[j];
              }
          }
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> NavierSystem::mass_residual (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                    DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEBase* u_elem_fe;
  
    c.get_element_fe( u_var, u_elem_fe );
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = u_elem_fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = u_elem_fe-&gt;get_phi();
  
    DenseSubVector&lt;Number&gt; &amp;Fu = c.get_elem_residual( u_var );
    DenseSubVector&lt;Number&gt; &amp;Fv = c.get_elem_residual( v_var );
    DenseSubVector&lt;Number&gt; &amp;Fw = c.get_elem_residual( w_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = c.get_elem_jacobian( u_var, u_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kvv = c.get_elem_jacobian( v_var, v_var );
    DenseSubMatrix&lt;Number&gt; &amp;Kww = c.get_elem_jacobian( w_var, w_var );
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.get_dof_indices( u_var ).size();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp = 0; qp != n_qpoints; ++qp)
      {
        Number u = c.interior_value(u_var, qp),
               v = c.interior_value(v_var, qp),
               w = c.interior_value(w_var, qp);
  
        Number JxWxRe   = JxW[qp] * Reynolds;
        Number JxWxRexU = JxWxRe * u;
        Number JxWxRexV = JxWxRe * v;
        Number JxWxRexW = JxWxRe * w;
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i != n_u_dofs; ++i)
          {
            Fu(i) += JxWxRexU * phi[i][qp];
            Fv(i) += JxWxRexV * phi[i][qp];
            <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
              Fw(i) += JxWxRexW * phi[i][qp];
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.get_elem_solution_derivative())
              {
                libmesh_assert_equal_to (c.get_elem_solution_derivative(), 1.0);
  
                Number JxWxRexPhiI = JxWxRe * phi[i][qp];
                Number JxWxRexPhiII = JxWxRexPhiI * phi[i][qp];
                Kuu(i,i) += JxWxRexPhiII;
                Kvv(i,i) += JxWxRexPhiII;
                <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                  Kww(i,i) += JxWxRexPhiII;
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = i+1; j != n_u_dofs; ++j)
                  {
                    Number Kij = JxWxRexPhiI * phi[j][qp];
                    Kuu(i,j) += Kij;
                    Kuu(j,i) += Kij;
                    Kvv(i,j) += Kij;
                    Kvv(j,i) += Kij;
                    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                      {
                        Kww(i,j) += Kij;
                        Kww(j,i) += Kij;
                      }
                  }
              }
          }
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> NavierSystem::postprocess()
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    Point p(1./3., 1./3.);
    Number u = point_value(u_var, p),
           v = point_value(v_var, p),
           w = (dim == 3) ? point_value(w_var, p) : 0;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;u(1/3,1/3) = (&quot;</FONT></B>
              &lt;&lt; u &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, &quot;</FONT></B>
              &lt;&lt; v &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, &quot;</FONT></B>
              &lt;&lt; w &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;)&quot;</FONT></B> &lt;&lt; std::endl;
  }
  
  
  
  
  Point NavierSystem::forcing(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
  {
    <B><FONT COLOR="#A020F0">switch</FONT></B> (application)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
        {
  	<B><FONT COLOR="#A020F0">return</FONT></B> Point(0.,0.,0.);
        }
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
        {
  	<B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
  	Point f;
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (dim==2)
  	  {
  	    f(0) =  std::sin(2.*libMesh::pi*p(1));
  	    f(1) = -std::sin(2.*libMesh::pi*p(0));
  	  }
  
  	<B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim==3)
  	  {
  	    f(0) =  std::sin(2.*libMesh::pi*p(1));
  	    f(1) = 0.;
  	    f(2) = -std::sin(2.*libMesh::pi*p(0));
  	  }
  
  	<B><FONT COLOR="#A020F0">return</FONT></B> f;
        }
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
        {
  	libmesh_assert_not_equal_to (1, <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension());
  	<B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
  	<B><FONT COLOR="#228B22">const</FONT></B> Real
  	  u=(Reynolds+1)*(y*y + z*z),
  	  v=(Reynolds+1)*(x*x + z*z),
  	  w=(Reynolds+1)*(x*x + y*y);
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension() == 2)
  	  <B><FONT COLOR="#A020F0">return</FONT></B> 2*Reynolds*(Reynolds+1)*Point(v*y,
  					       u*x);
  	<B><FONT COLOR="#A020F0">else</FONT></B>
  	  <B><FONT COLOR="#A020F0">return</FONT></B> 2*Reynolds*(Reynolds+1)*Point(v*y + w*z,
  					       u*x + w*z,
  					       u*x + v*y);
        }
  
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        {
  	libmesh_error();
        }
      }
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/fem_system/fem_system_ex1'
***************************************************************
* Running Example fem_system_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=445
  n_elem()=400
    n_local_elem()=100
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Navier-Stokes"
    Type "Implicit"
    Variables={ "u" "v" } "p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "FIRST", "THIRD" 
    n_dofs()=3803
    n_local_dofs()=1013
    n_constrained_dofs()=320
    n_local_constrained_dofs()=78
    n_vectors()=2
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 37.4744
      Average Off-Processor Bandwidth <= 2.15514
      Maximum  On-Processor Bandwidth <= 63
      Maximum Off-Processor Bandwidth <= 37
    DofMap Constraints
      Number of DoF Constraints = 320
      Number of Heterogenous Constraints= 41
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0



Solving time step 0, time = 0
u(1/3,1/3) = (-0.0312821, 0.0156421, 0)


Solving time step 1, time = 0.005
u(1/3,1/3) = (-0.0494434, 0.0255575, 0)


Solving time step 2, time = 0.01
u(1/3,1/3) = (-0.0629874, 0.0332361, 0)


Solving time step 3, time = 0.015
u(1/3,1/3) = (-0.0738558, 0.0395136, 0)


Solving time step 4, time = 0.02
u(1/3,1/3) = (-0.0827407, 0.0447345, 0)


Solving time step 5, time = 0.025
u(1/3,1/3) = (-0.0900054, 0.0491005, 0)


Solving time step 6, time = 0.03
u(1/3,1/3) = (-0.0959092, 0.0527519, 0)


Solving time step 7, time = 0.035
u(1/3,1/3) = (-0.100671, 0.0557977, 0)


Solving time step 8, time = 0.04
u(1/3,1/3) = (-0.104486, 0.0583278, 0)


Solving time step 9, time = 0.045
u(1/3,1/3) = (-0.107524, 0.0604198, 0)


Solving time step 10, time = 0.05
u(1/3,1/3) = (-0.109933, 0.0621413, 0)


Solving time step 11, time = 0.055
u(1/3,1/3) = (-0.111839, 0.0635516, 0)


Solving time step 12, time = 0.06
u(1/3,1/3) = (-0.113342, 0.0647021, 0)


Solving time step 13, time = 0.065
u(1/3,1/3) = (-0.114527, 0.0656372, 0)


Solving time step 14, time = 0.07
u(1/3,1/3) = (-0.115461, 0.0663945, 0)

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:50:04 2013                                                                          |
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
| libMesh Performance: Alive time=9.23078, Active time=8.74841                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0027      0.002685    0.0042      0.004199    0.03     0.05     |
|   build_constraint_matrix()        5400      0.0144      0.000003    0.0144      0.000003    0.17     0.17     |
|   build_sparsity()                 1         0.0019      0.001852    0.0042      0.004196    0.02     0.05     |
|   cnstrn_elem_mat_vec()            2700      0.0177      0.000007    0.0177      0.000007    0.20     0.20     |
|   constrain_elem_vector()          2700      0.0043      0.000002    0.0043      0.000002    0.05     0.05     |
|   create_dof_constraints()         1         0.0072      0.007210    0.0354      0.035405    0.08     0.40     |
|   distribute_dofs()                1         0.0043      0.004284    0.0102      0.010175    0.05     0.12     |
|   dof_indices()                    27851     0.2724      0.000010    0.2724      0.000010    3.11     3.11     |
|   enforce_constraints_exactly()    57        0.0167      0.000292    0.0167      0.000292    0.19     0.19     |
|   prepare_send_list()              1         0.0000      0.000035    0.0000      0.000035    0.00     0.00     |
|   reinit()                         1         0.0055      0.005533    0.0055      0.005533    0.06     0.06     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          15        0.0093      0.000621    0.0597      0.003980    0.11     0.68     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               15        1.5110      0.100730    1.5110      0.100730    17.27    17.27    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        69064     0.6402      0.000009    0.6402      0.000009    7.32     7.32     |
|   init_shape_functions()           4480      0.0368      0.000008    0.0368      0.000008    0.42     0.42     |
|   inverse_map()                    14328     0.0617      0.000004    0.0617      0.000004    0.70     0.70     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       27        0.6740      0.024962    1.3891      0.051450    7.70     15.88    |
|   assembly(get_residual)           27        0.2234      0.008274    0.8586      0.031801    2.55     9.81     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             69064     0.2725      0.000004    0.2725      0.000004    3.12     3.12     |
|   compute_face_map()               4264      0.0329      0.000008    0.0888      0.000021    0.38     1.01     |
|   init_face_shape_functions()      376       0.0009      0.000002    0.0009      0.000002    0.01     0.01     |
|   init_reference_to_physical_map() 4480      0.0417      0.000009    0.0417      0.000009    0.48     0.48     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0013      0.001306    0.0015      0.001497    0.01     0.02     |
|   renumber_nodes_and_elem()        2         0.0003      0.000132    0.0003      0.000132    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0037      0.001873    0.0037      0.001873    0.04     0.04     |
|   find_global_indices()            2         0.0005      0.000267    0.0050      0.002498    0.01     0.06     |
|   parallel_sort()                  2         0.0004      0.000180    0.0004      0.000215    0.00     0.00     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         15        0.0004      0.000025    1.5720      0.104803    0.00     17.97    |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0008      0.000829    0.0008      0.000829    0.01     0.01     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0051      0.005124    0.0075      0.007518    0.06     0.09     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          15        0.7238      0.048254    7.1164      0.474428    8.27     81.35    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0001      0.000015    0.0002      0.000023    0.00     0.00     |
|   broadcast()                      30        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   max(bool)                        1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                      903       0.0036      0.000004    0.0036      0.000004    0.04     0.04     |
|   max(vector)                      214       0.0013      0.000006    0.0037      0.000017    0.01     0.04     |
|   min(bool)                        1107      0.0038      0.000003    0.0038      0.000003    0.04     0.04     |
|   min(scalar)                      927       0.0223      0.000024    0.0223      0.000024    0.26     0.26     |
|   min(vector)                      214       0.0015      0.000007    0.0044      0.000020    0.02     0.05     |
|   probe()                          36        0.0002      0.000005    0.0002      0.000005    0.00     0.00     |
|   receive()                        36        0.0001      0.000004    0.0003      0.000009    0.00     0.00     |
|   send()                           36        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   send_receive()                   40        0.0002      0.000006    0.0007      0.000018    0.00     0.01     |
|   sum()                            119       0.0024      0.000020    0.0183      0.000154    0.03     0.21     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           36        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0008      0.000770    0.0009      0.000938    0.01     0.01     |
|   set_parent_processor_ids()       1         0.0001      0.000137    0.0001      0.000137    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          27        4.1211      0.152633    4.1211      0.152633    47.11    47.11    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  1         0.0018      0.001783    0.0023      0.002271    0.02     0.03     |
|   operator()                       30        0.0009      0.000030    0.0031      0.000104    0.01     0.04     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            208663    8.7484                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example fem_system_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/fem_system/fem_system_ex1'
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
