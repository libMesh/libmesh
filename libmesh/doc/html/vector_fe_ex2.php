<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("vector_fe_ex2",$root)?>
 
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
        #include "getpot.h"
        #include "exodusII_io.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exact_solution.h"
        
</pre>
</div>
<div class = "comment">
The systems and solvers we may use
</div>

<div class ="fragment">
<pre>
        #include "laplace_system.h"
        #include "diff_solver.h"
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
Parse the input file
</div>

<div class ="fragment">
<pre>
          GetPot infile("vector_fe_ex2.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
          const unsigned int grid_size = infile( "grid_size", 2 );
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "2D/3D support");
        
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
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_cube (mesh,
        				     grid_size,
        				     grid_size,
        				     grid_size,
        				     -1., 1.,
        				     -1., 1.,
        				     -1., 1.,
        				     HEX27);
        
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
          LaplaceSystem & system = 
            equation_systems.add_system&lt;LaplaceSystem&gt; ("Laplace");
        
</pre>
</div>
<div class = "comment">
This example only implements the steady-state problem
</div>

<div class ="fragment">
<pre>
          system.time_solver =
            AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
</pre>
</div>
<div class = "comment">
Initialize the system
</div>

<div class ="fragment">
<pre>
          equation_systems.init();
        
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
        
          system.solve();
        
          ExactSolution exact_sol( equation_systems );
          
          std::vector&lt;FunctionBase&lt;Number&gt;* &gt; sols;
          std::vector&lt;FunctionBase&lt;Gradient&gt;* &gt; grads;
        
          sols.push_back( new SolutionFunction(system.variable_number("u")) );
          grads.push_back( new SolutionGradient(system.variable_number("u")) );
        
          exact_sol.attach_exact_values(sols);
          exact_sol.attach_exact_derivs(grads);
          
</pre>
</div>
<div class = "comment">
Use higher quadrature order for more accurate error results
</div>

<div class ="fragment">
<pre>
          int extra_error_quadrature = infile("extra_error_quadrature",2);
          exact_sol.extra_quadrature_order(extra_error_quadrature);
        
</pre>
</div>
<div class = "comment">
Compute the error.
</div>

<div class ="fragment">
<pre>
          exact_sol.compute_error("Laplace", "u");
          
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "L2-Error is: "
        	    &lt;&lt; exact_sol.l2_error("Laplace", "u")
        	    &lt;&lt; std::endl;
          std::cout &lt;&lt; "H1-Error is: "
        	    &lt;&lt; exact_sol.h1_error("Laplace", "u")
        	    &lt;&lt; std::endl;
          
        #ifdef LIBMESH_HAVE_EXODUS_API
          
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
          ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
          
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
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
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exact_solution.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;laplace_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;steady_solver.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;vector_fe_ex2.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> grid_size = infile( <B><FONT COLOR="#BC8F8F">&quot;grid_size&quot;</FONT></B>, 2 );
  
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
    Mesh mesh;
    
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
  				     grid_size,
  				     grid_size,
  				     grid_size,
  				     -1., 1.,
  				     -1., 1.,
  				     -1., 1.,
  				     HEX27);
  
    mesh.print_info();
  
    EquationSystems equation_systems (mesh);
  
    LaplaceSystem &amp; system = 
      equation_systems.add_system&lt;LaplaceSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
    system.time_solver =
      AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
  
    equation_systems.init();
  
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
  
    system.solve();
  
    ExactSolution exact_sol( equation_systems );
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FunctionBase&lt;Number&gt;* &gt; sols;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FunctionBase&lt;Gradient&gt;* &gt; grads;
  
    sols.push_back( <B><FONT COLOR="#A020F0">new</FONT></B> SolutionFunction(system.variable_number(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)) );
    grads.push_back( <B><FONT COLOR="#A020F0">new</FONT></B> SolutionGradient(system.variable_number(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)) );
  
    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
    
    <B><FONT COLOR="#228B22">int</FONT></B> extra_error_quadrature = infile(<B><FONT COLOR="#BC8F8F">&quot;extra_error_quadrature&quot;</FONT></B>,2);
    exact_sol.extra_quadrature_order(extra_error_quadrature);
  
    exact_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Error is: &quot;</FONT></B>
  	    &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  	    &lt;&lt; std::endl;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1-Error is: &quot;</FONT></B>
  	    &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  	    &lt;&lt; std::endl;
    
  #ifdef LIBMESH_HAVE_EXODUS_API
    
    ExodusII_IO(mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out.e&quot;</FONT></B>, equation_systems);
    
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) laplace_system.C...
Compiling C++ (in optimized mode) vector_fe_ex2.C...
Linking vector_fe_ex2-opt...
***************************************************************
* Running Example  ./vector_fe_ex2-opt
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1331
    n_local_nodes()=1331
  n_elem()=125
    n_local_elem()=125
    n_active_elem()=125
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Laplace"
    Type "Implicit"
    Variables="u" 
    Finite Element Types="LAGRANGE_VEC" 
    Approximation Orders="SECOND" 
    n_dofs()=3993
    n_local_dofs()=3993
    n_constrained_dofs()=1806
    n_local_constrained_dofs()=1806
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 155.344
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 375
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 1806
      Number of Heterogenous Constraints= 1781
      Average DoF Constraint Length= 0

Assembling the System
Nonlinear Residual: 5.35071
Linear solve starting, tolerance 1e-12
Linear solve finished, step 20, residual 1.72623e-12
Trying full Newton step
  Current Residual: 1.45958e-12
  Nonlinear solver converged, step 0, residual reduction 2.72783e-13 < 1e-12
  Nonlinear solver relative step size 1.18335 > 1e-06
L2-Error is: 0.00428631
H1-Error is: 0.0697201

-------------------------------------------------------------------------------------------------------------------
| Time:           Thu Aug 23 08:08:09 2012                                                                         |
| OS:             Darwin                                                                                           |
| HostName:       www.example.com                                                                                  |
| OS Release:     10.8.0                                                                                           |
| OS Version:     Darwin Kernel Version 10.8.0: Tue Jun  7 16:32:41 PDT 2011; root:xnu-1504.15.3~1/RELEASE_X86_64  |
| Machine:        x86_64                                                                                           |
| Username:       jwpeterson                                                                                       |
| Configuration:  ./configure run on Wed Aug 22 13:54:56 MDT 2012                                                  |
-------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.968076, Active time=0.862715                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000125    0.0001      0.000125    0.01     0.01     |
|   build_constraint_matrix()        250       0.0100      0.000040    0.0100      0.000040    1.16     1.16     |
|   build_sparsity()                 1         0.0222      0.022215    0.0224      0.022417    2.58     2.60     |
|   cnstrn_elem_mat_vec()            125       0.0678      0.000543    0.0678      0.000543    7.86     7.86     |
|   constrain_elem_vector()          125       0.0013      0.000010    0.0013      0.000010    0.15     0.15     |
|   create_dof_constraints()         1         0.0103      0.010280    0.0723      0.072336    1.19     8.38     |
|   distribute_dofs()                1         0.0002      0.000243    0.0010      0.001005    0.03     0.12     |
|   dof_indices()                    1000      0.0017      0.000002    0.0017      0.000002    0.19     0.19     |
|   enforce_constraints_exactly()    3         0.0006      0.000214    0.0006      0.000214    0.07     0.07     |
|   prepare_send_list()              1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                         1         0.0008      0.000761    0.0008      0.000761    0.09     0.09     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0003      0.000314    0.0005      0.000460    0.04     0.05     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0023      0.002317    0.0023      0.002317    0.27     0.27     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        1365      0.3139      0.000230    0.3139      0.000230    36.39    36.39    |
|   init_shape_functions()           993       0.0091      0.000009    0.0091      0.000009    1.06     1.06     |
|   inverse_map()                    1620      0.0158      0.000010    0.0158      0.000010    1.84     1.84     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       1         0.2512      0.251199    0.3631      0.363132    29.12    42.09    |
|   assembly(get_residual)           1         0.0155      0.015502    0.0605      0.060518    1.80     7.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             1365      0.0095      0.000007    0.0095      0.000007    1.10     1.10     |
|   compute_edge_map()               540       0.0006      0.000001    0.0006      0.000001    0.07     0.07     |
|   compute_face_map()               450       0.0019      0.000004    0.0019      0.000004    0.22     0.22     |
|   init_edge_shape_functions()      540       0.0006      0.000001    0.0006      0.000001    0.07     0.07     |
|   init_face_shape_functions()      152       0.0018      0.000012    0.0018      0.000012    0.21     0.21     |
|   init_reference_to_physical_map() 993       0.0515      0.000052    0.0515      0.000052    5.97     5.97     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0006      0.000632    0.0006      0.000632    0.07     0.07     |
|   renumber_nodes_and_elem()        2         0.0001      0.000066    0.0001      0.000066    0.02     0.02     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000043    0.0028      0.002820    0.00     0.33     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0005      0.000525    0.0005      0.000525    0.06     0.06     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          1         0.0005      0.000543    0.4965      0.496459    0.06     57.55    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   single_partition()               1         0.0001      0.000057    0.0001      0.000057    0.01     0.01     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0716      0.071621    0.0716      0.071621    8.30     8.30     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            9540      0.8627                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./vector_fe_ex2-opt
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
