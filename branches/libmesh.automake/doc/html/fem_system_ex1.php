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
<div class = "comment">
<h1>FEMSystem Example 1 - Unsteady Navier-Stokes Equations with
FEMSystem</h1>

<br><br>This example shows how the transient nonlinear problem from
example 13 can be solved using the
DifferentiableSystem class framework


<br><br>Basic include files
</div>

<div class ="fragment">
<pre>
        #include "equation_systems.h"
        #include "error_vector.h"
        #include "getpot.h"
        #include "exodusII_io.h"
        #include "kelly_error_estimator.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "mesh_refinement.h"
        #include "uniform_refinement_estimator.h"
        
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
The systems and solvers we may use
</div>

<div class ="fragment">
<pre>
        #include "naviersystem.h"
        #include "diff_solver.h"
        #include "euler_solver.h"
        #include "steady_solver.h"
        
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
Trilinos solver NaNs by default on the zero pressure block.
We'll skip this example for now.
</div>

<div class ="fragment">
<pre>
          if (libMesh::default_solver_package() == TRILINOS_SOLVERS)
            {
              std::cout &lt;&lt; "We skip example 18 when using the Trilinos solvers.\n"
                        &lt;&lt; std::endl;
              return 0;
            }
        
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
              libmesh_assert(n_timesteps == 1);
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
                      libmesh_assert (nelem_target == 0);
        
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
                      libmesh_assert (nelem_target &gt; 0);
        
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
                  OStringStream file_name;
        
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
                  file_name &lt;&lt; "out_";
                  OSSRealzeroright(file_name,3,0, t_step + 1);
                  file_name &lt;&lt; ".e";
        
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

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;uniform_refinement_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;naviersystem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;euler_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;steady_solver.h&quot;</FONT></B>
  
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
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::default_solver_package() == TRILINOS_SOLVERS)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;We skip example 18 when using the Trilinos solvers.\n&quot;</FONT></B>
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
      }
  
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
  
    Mesh mesh;
    
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
        libmesh_assert(n_timesteps == 1);
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
                libmesh_assert (nelem_target == 0);
  
                UniformRefinementEstimator *u =
                  <B><FONT COLOR="#A020F0">new</FONT></B> UniformRefinementEstimator;
  
                u-&gt;error_norm = L2;
  
                error_estimator.reset(u);
              }
            <B><FONT COLOR="#A020F0">else</FONT></B>
              {
                libmesh_assert (nelem_target &gt; 0);
  
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
            OStringStream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0, t_step + 1);
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
  
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
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking fem_system_ex1-opt...
***************************************************************
* Running Example  mpirun -np 6 ./fem_system_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=305
  n_elem()=400
    n_local_elem()=67
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Navier-Stokes"
    Type "Implicit"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" "FIRST", "THIRD" 
    n_dofs()=3803
    n_local_dofs()=696
    n_constrained_dofs()=321
    n_local_constrained_dofs()=66
    n_vectors()=2
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 36.7155
      Average Off-Processor Bandwidth <= 2.72557
      Maximum  On-Processor Bandwidth <= 59
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
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./fem_system_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:26 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           9.592e+00      1.00013   9.592e+00
Objects:              4.430e+02      1.00000   4.430e+02
Flops:                2.525e+09      1.35781   2.136e+09  1.282e+10
Flops/sec:            2.632e+08      1.35774   2.227e+08  1.336e+09
MPI Messages:         1.244e+04      1.86142   9.541e+03  5.725e+04
MPI Message Lengths:  8.619e+06      1.50071   7.289e+02  4.173e+07
MPI Reductions:       7.208e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.5915e+00 100.0%  1.2818e+10 100.0%  5.725e+04 100.0%  7.289e+02      100.0%  7.052e+03  97.8% 

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

KSPGMRESOrthog      2540 1.0 1.5597e+00 3.2 1.04e+08 1.2 0.0e+00 0.0e+00 2.5e+03  9  4  0  0 35   9  4  0  0 36   364
KSPSetup              54 1.0 5.9843e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              27 1.0 8.1034e+00 1.0 2.52e+09 1.4 4.8e+04 4.5e+02 5.3e+03 84100 83 52 73  84100 83 52 75  1582
PCSetUp               54 1.0 3.5936e+00 1.7 1.08e+09 1.5 0.0e+00 0.0e+00 8.1e+01 29 41  0  0  1  29 41  0  0  1  1474
PCSetUpOnBlocks       27 1.0 3.5932e+00 1.7 1.08e+09 1.5 0.0e+00 0.0e+00 8.1e+01 29 41  0  0  1  29 41  0  0  1  1474
PCApply             2640 1.0 1.2995e+00 1.4 1.18e+09 1.3 0.0e+00 0.0e+00 0.0e+00 11 48  0  0  0  11 48  0  0  0  4715
VecMDot             2540 1.0 1.5364e+00 3.3 5.20e+07 1.2 0.0e+00 0.0e+00 2.5e+03  9  2  0  0 35   9  2  0  0 36   185
VecNorm             2775 1.0 2.2043e+00 3.8 3.86e+06 1.2 0.0e+00 0.0e+00 2.8e+03 13  0  0  0 38  13  0  0  0 39    10
VecScale            2640 1.0 3.7777e-03 1.7 1.84e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2658
VecCopy              427 1.0 5.9485e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              2855 1.0 3.7463e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              227 1.0 5.3239e-04 1.4 3.16e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3243
VecMAXPY            2640 1.0 2.1785e-02 1.1 5.56e+07 1.2 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0 13940
VecAssemblyBegin     421 1.0 7.3882e-01 1.7 0.00e+00 0.0 9.7e+02 4.3e+02 1.1e+03  6  0  2  1 15   6  0  2  1 15     0
VecAssemblyEnd       421 1.0 4.3893e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     2808 1.0 2.2934e-02 1.8 0.00e+00 0.0 5.1e+04 6.2e+02 0.0e+00  0  0 90 76  0   0  0 90 76  0     0
VecScatterEnd       2808 1.0 3.0482e+00 3.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 19  0  0  0  0  19  0  0  0  0     0
VecNormalize        2640 1.0 2.1417e+00 4.1 5.51e+06 1.2 0.0e+00 0.0e+00 2.6e+03 13  0  0  0 37  13  0  0  0 37    14
MatMult             2640 1.0 3.1629e+00 3.1 1.43e+08 1.2 4.8e+04 4.5e+02 0.0e+00 20  6 83 52  0  20  6 83 52  0   244
MatSolve            2640 1.0 1.2728e+00 1.4 1.18e+09 1.3 0.0e+00 0.0e+00 0.0e+00 11 48  0  0  0  11 48  0  0  0  4814
MatLUFactorNum        27 1.0 1.1125e+00 1.7 1.08e+09 1.5 0.0e+00 0.0e+00 0.0e+00  9 41  0  0  0   9 41  0  0  0  4762
MatILUFactorSym       27 1.0 2.4791e+00 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 2.7e+01 20  0  0  0  0  20  0  0  0  0     0
MatAssemblyBegin      54 1.0 1.2764e-02 1.7 0.00e+00 0.0 7.3e+02 7.0e+03 1.1e+02  0  0  1 12  1   0  0  1 12  2     0
MatAssemblyEnd        54 1.0 1.5276e-02 1.3 0.00e+00 0.0 3.6e+01 1.2e+02 6.0e+01  0  0  0  0  1   0  0  0  0  1     0
MatGetRowIJ           27 1.0 1.2398e-05 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        27 1.0 3.6716e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 5.4e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        29 1.0 1.3266e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     3              3        19872     0
      Preconditioner     3              3         2048     0
                 Vec   145            145      2329872     0
         Vec Scatter    74             74        64232     0
           Index Set   172            172       328348     0
   IS L to G Mapping    16             16        63744     0
              Matrix    30             30     73202212     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 0.000487232
Average time for zero size MPI_Send(): 7.93139e-05
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
| Time:           Fri Aug 24 15:21:26 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=9.87113, Active time=9.53144                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000088    0.0001      0.000125    0.00     0.00     |
|   build_constraint_matrix()        3618      0.0118      0.000003    0.0118      0.000003    0.12     0.12     |
|   build_sparsity()                 1         0.0022      0.002173    0.0025      0.002501    0.02     0.03     |
|   cnstrn_elem_mat_vec()            1809      0.0147      0.000008    0.0147      0.000008    0.15     0.15     |
|   constrain_elem_vector()          1809      0.0015      0.000001    0.0015      0.000001    0.02     0.02     |
|   create_dof_constraints()         1         0.0038      0.003773    0.0071      0.007063    0.04     0.07     |
|   distribute_dofs()                1         0.0005      0.000523    0.0025      0.002510    0.01     0.03     |
|   dof_indices()                    19528     0.0106      0.000001    0.0106      0.000001    0.11     0.11     |
|   enforce_constraints_exactly()    57        0.0470      0.000825    0.0470      0.000825    0.49     0.49     |
|   prepare_send_list()              1         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|   reinit()                         1         0.0010      0.000959    0.0010      0.000959    0.01     0.01     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          15        0.0040      0.000263    0.0166      0.001109    0.04     0.17     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               15        0.0137      0.000915    0.0137      0.000915    0.14     0.14     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        47032     0.1914      0.000004    0.1914      0.000004    2.01     2.01     |
|   init_shape_functions()           3832      0.0166      0.000004    0.0166      0.000004    0.17     0.17     |
|   inverse_map()                    12222     0.0195      0.000002    0.0195      0.000002    0.20     0.20     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       27        0.3269      0.012107    0.5188      0.019215    3.43     5.44     |
|   assembly(get_residual)           27        0.1946      0.007209    0.3641      0.013486    2.04     3.82     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             47032     0.0558      0.000001    0.0558      0.000001    0.59     0.59     |
|   compute_face_map()               3616      0.0142      0.000004    0.0320      0.000009    0.15     0.34     |
|   init_face_shape_functions()      376       0.0007      0.000002    0.0007      0.000002    0.01     0.01     |
|   init_reference_to_physical_map() 3832      0.0202      0.000005    0.0202      0.000005    0.21     0.21     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0004      0.000393    0.0015      0.001469    0.00     0.02     |
|   renumber_nodes_and_elem()        2         0.0001      0.000046    0.0001      0.000046    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0015      0.000775    0.0015      0.000775    0.02     0.02     |
|   find_global_indices()            2         0.0002      0.000096    0.0060      0.002980    0.00     0.06     |
|   parallel_sort()                  2         0.0019      0.000956    0.0038      0.001896    0.02     0.04     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         15        0.0003      0.000022    0.0307      0.002047    0.00     0.32     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0004      0.000396    0.0004      0.000396    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0014      0.001356    0.0044      0.004357    0.01     0.05     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          15        0.4028      0.026852    9.4725      0.631500    4.23     99.38    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0007      0.000082    0.0007      0.000082    0.01     0.01     |
|   broadcast()                      31        0.0042      0.000136    0.0042      0.000136    0.04     0.04     |
|   gather()                         1         0.0001      0.000066    0.0001      0.000066    0.00     0.00     |
|   max(scalar)                      2         0.0071      0.003541    0.0071      0.003541    0.07     0.07     |
|   max(vector)                      4         0.0008      0.000190    0.0008      0.000190    0.01     0.01     |
|   min(scalar)                      30        0.0027      0.000090    0.0027      0.000090    0.03     0.03     |
|   min(vector)                      4         0.0003      0.000074    0.0003      0.000074    0.00     0.00     |
|   probe()                          50        0.0005      0.000011    0.0005      0.000011    0.01     0.01     |
|   receive()                        50        0.0001      0.000002    0.0006      0.000013    0.00     0.01     |
|   send()                           50        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|   send_receive()                   54        0.0001      0.000002    0.0009      0.000016    0.00     0.01     |
|   sum()                            109       0.0155      0.000142    0.0155      0.000142    0.16     0.16     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000122    0.0002      0.000170    0.00     0.00     |
|   set_parent_processor_ids()       1         0.0000      0.000027    0.0000      0.000027    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          27        8.1376      0.301391    8.1376      0.301391    85.38    85.38    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  1         0.0014      0.001389    0.0017      0.001718    0.01     0.02     |
|   operator()                       30        0.0005      0.000016    0.0014      0.000046    0.01     0.01     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            145397    9.5314                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./fem_system_ex1-opt
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
