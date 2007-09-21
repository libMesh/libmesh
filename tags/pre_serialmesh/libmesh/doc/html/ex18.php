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
        #include "gmv_io.h"
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
          libMesh::init (argc, argv);
        
        #ifndef ENABLE_AMR
          std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
                    &lt;&lt; "compiled with AMR support!"
                    &lt;&lt; std::endl;
          return 0;
        #else
        
          {    
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
        
            assert (dim == 2 || dim == 3);
</pre>
</div>
<div class = "comment">
Create a n-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
            
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
                assert(n_timesteps == 1);
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
                std::cout &lt;&lt; " Solving time step " &lt;&lt; t_step &lt;&lt; ", time = "
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
                        assert (nelem_target == 0);
        
                        UniformRefinementEstimator *u =
                          new UniformRefinementEstimator;
        
</pre>
</div>
<div class = "comment">
The lid-driven cavity problem isn't in H1, so
lets estimate H0 (i.e. L2) error
</div>

<div class ="fragment">
<pre>
                        u-&gt;sobolev_order() = 0;
        
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
                        assert (nelem_target &gt; 0);
        
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
Calculate error based on u and v but not p
</div>

<div class ="fragment">
<pre>
                    error_estimator-&gt;component_scale.push_back(1.0); // u
                    error_estimator-&gt;component_scale.push_back(1.0); // v
                    if (dim == 3)
                      error_estimator-&gt;component_scale.push_back(1.0); // w
                    error_estimator-&gt;component_scale.push_back(0.0); // p
        
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
We write the file name in the gmv auto-read format.
</div>

<div class ="fragment">
<pre>
                    file_name &lt;&lt; "out.gmv.";
                    OSSRealzeroright(file_name,3,0, t_step + 1);
        
                    GMVIO(mesh).write_equation_systems (file_name.str(),
                                                        equation_systems);
                  }
              }
          }
        #endif // #ifndef ENABLE_AMR
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
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
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::init (argc, argv);
  
  #ifndef ENABLE_AMR
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with AMR support!&quot;</FONT></B>
              &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    {    
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
  
      assert (dim == 2 || dim == 3);
      Mesh mesh (dim);
      
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
          assert(n_timesteps == 1);
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
  
      solver.max_linear_iterations =
        infile(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 50000);
      solver.initial_linear_tolerance =
        infile(<B><FONT COLOR="#BC8F8F">&quot;initial_linear_tolerance&quot;</FONT></B>, 1.e-3);
  
      equation_systems.print_info();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=0; t_step != n_timesteps; ++t_step)
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B>
                    &lt;&lt; system.time &lt;&lt; std::endl;
  
          <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> a_step = 0;
          <B><FONT COLOR="#A020F0">for</FONT></B> (; a_step != max_adaptivesteps; ++a_step)
            {
              system.solve();
  
              ErrorVector error;
  
              AutoPtr&lt;ErrorEstimator&gt; error_estimator;
  
              <B><FONT COLOR="#A020F0">if</FONT></B> (global_tolerance != 0.)
                {
                  assert (nelem_target == 0);
  
                  UniformRefinementEstimator *u =
                    <B><FONT COLOR="#A020F0">new</FONT></B> UniformRefinementEstimator;
  
                  u-&gt;sobolev_order() = 0;
  
                  error_estimator.reset(u);
                }
              <B><FONT COLOR="#A020F0">else</FONT></B>
                {
                  assert (nelem_target &gt; 0);
  
                  error_estimator.reset(<B><FONT COLOR="#A020F0">new</FONT></B> KellyErrorEstimator);
                }
  
              error_estimator-&gt;component_scale.push_back(1.0); <I><FONT COLOR="#B22222">// u
</FONT></I>              error_estimator-&gt;component_scale.push_back(1.0); <I><FONT COLOR="#B22222">// v
</FONT></I>              <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                error_estimator-&gt;component_scale.push_back(1.0); <I><FONT COLOR="#B22222">// w
</FONT></I>              error_estimator-&gt;component_scale.push_back(0.0); <I><FONT COLOR="#B22222">// p
</FONT></I>  
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
  
          <B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
            {
              OStringStream file_name;
  
              file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.gmv.&quot;</FONT></B>;
              OSSRealzeroright(file_name,3,0, t_step + 1);
  
              GMVIO(mesh).write_equation_systems (file_name.str(),
                                                  equation_systems);
            }
        }
    }
  #endif <I><FONT COLOR="#B22222">// #ifndef ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex18-devel
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
  n_elem()=400
   n_local_elem()=400
   n_active_elem()=400
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

*** Warning, This code is untested, experimental, or likely to see future API changes: src/solvers/diff_system.C, line 29, compiled Jun  1 2007 at 14:30:46 ***
 EquationSystems
  n_systems()=1
   System "Navier-Stokes"
    Type "Implicit"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" 
    Approximation Orders="SECOND" "SECOND" "FIRST" 
    n_dofs()=3803
    n_local_dofs()=3803
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
 
***************************************************************
* Done Running Example  ./ex18-devel
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
