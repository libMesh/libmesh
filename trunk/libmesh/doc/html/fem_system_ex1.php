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
        
                  ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                      equation_systems);
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
  
            ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
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
Compiling C++ (in optimized mode) fem_system_ex1.C...
Compiling C++ (in optimized mode) naviersystem.C...
Linking fem_system_ex1-opt...
***************************************************************
* Running Example  ./fem_system_ex1-opt
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=1681
  n_elem()=400
    n_local_elem()=400
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=3803
    n_constrained_dofs()=320
    n_local_constrained_dofs()=320
    n_vectors()=2
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 38.9716
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 59
      Maximum Off-Processor Bandwidth <= 0
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

-------------------------------------------------------------------
| Time:           Sat Apr  7 15:58:49 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=12.7738, Active time=12.6052                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  1         0.0005      0.000483    0.0005      0.000483    0.00     0.00     |
|   build_constraint_matrix()     21600     0.0784      0.000004    0.0784      0.000004    0.62     0.62     |
|   cnstrn_elem_mat_vec()         10800     0.0730      0.000007    0.0730      0.000007    0.58     0.58     |
|   compute_sparsity()            1         0.0198      0.019806    0.0230      0.022977    0.16     0.18     |
|   constrain_elem_vector()       10800     0.0274      0.000003    0.0274      0.000003    0.22     0.22     |
|   create_dof_constraints()      1         0.0090      0.008978    0.0188      0.018827    0.07     0.15     |
|   distribute_dofs()             1         0.0019      0.001901    0.0054      0.005450    0.02     0.04     |
|   dof_indices()                 106430    0.1740      0.000002    0.1740      0.000002    1.38     1.38     |
|   enforce_constraints_exactly() 57        0.0077      0.000135    0.0077      0.000135    0.06     0.06     |
|   prepare_send_list()           1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                      1         0.0035      0.003545    0.0035      0.003545    0.03     0.03     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       15        0.0220      0.001469    0.0363      0.002419    0.17     0.29     |
|                                                                                                             |
| ExodusII_IO                                                                                                 |
|   write_nodal_data()            15        0.0261      0.001740    0.0261      0.001740    0.21     0.21     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          138430    0.2248      0.000002    0.2248      0.000002    1.78     1.78     |
|   compute_face_map()            8800      0.0501      0.000006    0.1264      0.000014    0.40     1.00     |
|   compute_shape_functions()     138430    0.0938      0.000001    0.0938      0.000001    0.74     0.74     |
|   init_face_shape_functions()   268       0.0006      0.000002    0.0006      0.000002    0.00     0.00     |
|   init_shape_functions()        8938      0.0773      0.000009    0.0773      0.000009    0.61     0.61     |
|   inverse_map()                 31476     0.0849      0.000003    0.0849      0.000003    0.67     0.67     |
|                                                                                                             |
| FEMSystem                                                                                                   |
|   assembly()                    27        1.0053      0.037234    1.4654      0.054273    7.98     11.63    |
|   assembly(get_residual)        27        0.3922      0.014525    0.8667      0.032101    3.11     6.88     |
|                                                                                                             |
| Mesh                                                                                                        |
|   find_neighbors()              1         0.0014      0.001437    0.0014      0.001437    0.01     0.01     |
|   renumber_nodes_and_elem()     2         0.0002      0.000097    0.0002      0.000097    0.00     0.00     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      15        0.0003      0.000017    0.0627      0.004178    0.00     0.50     |
|                                                                                                             |
| MeshTools::Generation                                                                                       |
|   build_cube()                  1         0.0014      0.001417    0.0014      0.001417    0.01     0.01     |
|                                                                                                             |
| NewtonSolver                                                                                                |
|   solve()                       15        0.0071      0.000474    12.5671     0.837805    0.06     99.70    |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   1         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   single_partition()            1         0.0002      0.000158    0.0002      0.000158    0.00     0.00     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       27        10.2201     0.378523    10.2201     0.378523    81.08    81.08    |
|                                                                                                             |
| PointLocatorTree                                                                                            |
|   init(no master)               1         0.0013      0.001333    0.0013      0.001333    0.01     0.01     |
|   operator()                    30        0.0007      0.000023    0.0020      0.000066    0.01     0.02     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         476213    12.6052                                         100.00            |
 -------------------------------------------------------------------------------------------------------------

 
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
