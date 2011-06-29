<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex18",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 18 - Unsteady Navier-Stokes Equations with DiffSystem</h1>

<br><br>This example shows how the transient nonlinear problem from
example 13 can be solved using the new (and experimental)
DiffSystem class framework


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
          GetPot infile("ex18.in");
        
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
                  std::cout &lt;&lt; "adaptive step " &lt;&lt; a_step &lt;&lt; ": ";
                  if (global_tolerance != 0.)
                    std::cout &lt;&lt; "global_error = " &lt;&lt; global_error
                              &lt;&lt; " with ";
                  std::cout &lt;&lt; mesh.n_active_elem()
                            &lt;&lt; " active elements and "
                            &lt;&lt; equation_systems.n_active_dofs()
                            &lt;&lt; " active dofs." &lt;&lt; std::endl;
                  if (global_tolerance != 0.)
                    std::cout &lt;&lt; "worst element error = " &lt;&lt; error.maximum()
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
                  file_name &lt;&lt; ".exd";
        
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
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;ex18.in&quot;</FONT></B>);
  
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
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;adaptive step &quot;</FONT></B> &lt;&lt; a_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;: &quot;</FONT></B>;
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;global_error = &quot;</FONT></B> &lt;&lt; global_error
                        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; with &quot;</FONT></B>;
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; mesh.n_active_elem()
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active elements and &quot;</FONT></B>
                      &lt;&lt; equation_systems.n_active_dofs()
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; active dofs.&quot;</FONT></B> &lt;&lt; std::endl;
            <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;worst element error = &quot;</FONT></B> &lt;&lt; error.maximum()
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
          }
        <B><FONT COLOR="#A020F0">if</FONT></B> (a_step == max_adaptivesteps)
          {
            system.solve();
          }
  
        system.time_solver-&gt;advance_timestep();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
        <B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
          {
            OStringStream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0, t_step + 1);
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.exd&quot;</FONT></B>;
  
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
Compiling C++ (in optimized mode) ex18.C...
Compiling C++ (in optimized mode) naviersystem.C...
Linking ex18-opt...
***************************************************************
* Running Example  mpirun -np 2 ./ex18-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=863
  n_elem()=400
    n_local_elem()=200
    n_active_elem()=400
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Navier-Stokes"
    Type "Implicit"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" 
    Approximation Orders="SECOND" "SECOND" "FIRST" 
    n_dofs()=3803
    n_local_dofs()=1958
    n_constrained_dofs()=0
    n_vectors()=2



Solving time step 0, time = 0


Solving time step 1, time = 0.005


Solving time step 2, time = 0.01


Solving time step 3, time = 0.015


Solving time step 4, time = 0.02


Solving time step 5, time = 0.025


Solving time step 6, time = 0.03


Solving time step 7, time = 0.035


Solving time step 8, time = 0.04


Solving time step 9, time = 0.045


Solving time step 10, time = 0.05


Solving time step 11, time = 0.055


Solving time step 12, time = 0.06


Solving time step 13, time = 0.065


Solving time step 14, time = 0.07
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex18-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Tue Feb 22 12:21:48 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.322e+01      1.00007   1.322e+01
Objects:              4.840e+02      1.00000   4.840e+02
Flops:                8.464e+09      1.14869   7.916e+09  1.583e+10
Flops/sec:            6.403e+08      1.14877   5.988e+08  1.198e+09
MPI Messages:         2.178e+03      1.00000   2.178e+03  4.356e+03
MPI Message Lengths:  5.076e+06      1.01961   2.308e+03  1.005e+07
MPI Reductions:       4.729e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.3218e+01 100.0%  1.5832e+10 100.0%  4.356e+03 100.0%  2.308e+03      100.0%  4.632e+03  97.9% 

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

KSPGMRESOrthog      1454 1.0 1.9099e-01 3.1 1.64e+08 1.1 0.0e+00 0.0e+00 1.5e+03  1  2  0  0 31   1  2  0  0 31  1672
KSPSetup              52 1.0 4.7922e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              26 1.0 1.2732e+01 1.0 8.46e+09 1.1 3.0e+03 1.3e+03 3.1e+03 96100 69 40 65  96100 69 40 66  1243
PCSetUp               52 1.0 1.0410e+01 1.1 5.54e+09 1.2 0.0e+00 0.0e+00 7.8e+01 74 65  0  0  2  74 65  0  0  2   984
PCSetUpOnBlocks       26 1.0 1.0410e+01 1.1 5.54e+09 1.2 0.0e+00 0.0e+00 7.8e+01 74 65  0  0  2  74 65  0  0  2   984
PCApply             1513 1.0 2.0386e+00 1.1 2.52e+09 1.1 0.0e+00 0.0e+00 0.0e+00 15 30  0  0  0  15 30  0  0  0  2351
VecMDot             1454 1.0 1.6936e-01 4.3 8.22e+07 1.1 0.0e+00 0.0e+00 1.5e+03  1  1  0  0 31   1  1  0  0 31   943
VecNorm             1654 1.0 3.5646e-02 1.2 6.48e+06 1.1 0.0e+00 0.0e+00 1.7e+03  0  0  0  0 35   0  0  0  0 36   353
VecScale            1513 1.0 2.4424e-03 1.1 2.96e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2356
VecCopy              193 1.0 4.4966e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              1684 1.0 2.1021e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              144 1.0 4.9233e-04 1.3 5.64e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2225
VecMAXPY            1513 1.0 2.2277e-02 1.1 8.79e+07 1.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  7667
VecAssemblyBegin     283 1.0 8.8170e-03 3.1 0.00e+00 0.0 1.0e+02 1.1e+03 8.5e+02  0  0  2  1 18   0  0  2  1 18     0
VecAssemblyEnd       283 1.0 1.8954e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1729 1.0 1.2525e-02 1.9 0.00e+00 0.0 3.2e+03 1.4e+03 0.0e+00  0  0 74 44  0   0  0 74 44  0     0
VecScatterEnd       1729 1.0 1.3028e+00292.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  5  0  0  0  0   5  0  0  0  0     0
VecNormalize        1513 1.0 3.4015e-02 1.2 8.89e+06 1.1 0.0e+00 0.0e+00 1.5e+03  0  0  0  0 32   0  0  0  0 33   507
MatMult             1513 1.0 1.5012e+00 8.1 2.30e+08 1.1 3.0e+03 1.3e+03 0.0e+00  6  3 69 40  0   6  3 69 40  0   295
MatSolve            1513 1.0 2.0200e+00 1.1 2.52e+09 1.1 0.0e+00 0.0e+00 0.0e+00 15 30  0  0  0  15 30  0  0  0  2373
MatLUFactorNum        26 1.0 4.1903e+00 1.2 5.54e+09 1.2 0.0e+00 0.0e+00 0.0e+00 29 65  0  0  0  29 65  0  0  0  2445
MatILUFactorSym       26 1.0 6.2176e+00 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.6e+01 44  0  0  0  1  44  0  0  0  1     0
MatAssemblyBegin      52 1.0 2.0273e-03 2.9 0.00e+00 0.0 7.8e+01 2.2e+04 1.0e+02  0  0  2 17  2   0  0  2 17  2     0
MatAssemblyEnd        52 1.0 9.5301e-03 1.2 0.00e+00 0.0 4.0e+00 3.3e+02 5.8e+01  0  0  0  0  1   0  0  0  0  1     0
MatGetRowIJ           26 1.0 7.6294e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        26 1.0 3.6812e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 5.2e+01  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        28 1.0 2.0542e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     3              3        19872     0
      Preconditioner     3              3         2048     0
                 Vec    88             88      1211040     0
         Vec Scatter   125            125       108500     0
           Index Set   220            220       739584     0
   IS L to G Mapping    16             16       145344     0
              Matrix    29             29    260814268     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 1.62125e-06
Average time for zero size MPI_Send(): 5.96046e-06
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
 
***************************************************************
* Done Running Example  ./ex18-opt
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
