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
Linking vector_fe_ex2-opt...
***************************************************************
* Running Example  mpirun -np 6 ./vector_fe_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=1331
    n_local_nodes()=275
  n_elem()=125
    n_local_elem()=20
    n_active_elem()=125
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Laplace"
    Type "Implicit"
    Variables="u" 
    Finite Element Types="LAGRANGE_VEC", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=3993
    n_local_dofs()=825
    n_constrained_dofs()=1811
    n_local_constrained_dofs()=346
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 130.745
      Average Off-Processor Bandwidth <= 38.5964
      Maximum  On-Processor Bandwidth <= 375
      Maximum Off-Processor Bandwidth <= 294
    DofMap Constraints
      Number of DoF Constraints = 1806
      Number of Heterogenous Constraints= 1783
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

Assembling the System
Nonlinear Residual: 5.35071
Linear solve starting, tolerance 1e-12
Linear solve finished, step 30, residual 1.94384e-12
Trying full Newton step
  Current Residual: 1.94342e-12
  Nonlinear solver converged, step 0, residual reduction 3.63207e-13 < 1e-12
  Nonlinear solver relative step size 1.18335 > 1e-06
L2-Error is: 0.00428631
H1-Error is: 0.0697201
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./vector_fe_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:28:20 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           8.080e-01      1.00245   8.066e-01
Objects:              7.400e+01      1.00000   7.400e+01
Flops:                7.139e+07      2.81323   4.520e+07  2.712e+08
Flops/sec:            8.836e+07      2.80636   5.603e+07  3.362e+08
MPI Messages:         2.575e+02      1.16516   2.442e+02  1.465e+03
MPI Message Lengths:  9.078e+05      1.26672   3.353e+03  4.912e+06
MPI Reductions:       1.690e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 8.0656e-01 100.0%  2.7122e+08 100.0%  1.465e+03 100.0%  3.353e+03      100.0%  1.500e+02  88.8% 

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

KSPGMRESOrthog        30 1.0 9.4635e-03 2.2 1.53e+06 1.8 0.0e+00 0.0e+00 3.0e+01  1  3  0  0 18   1  3  0  0 20   785
KSPSetup               2 1.0 3.1090e-04 7.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 3.2399e-01 1.0 7.14e+07 2.8 8.7e+02 1.3e+03 6.7e+01 40100 59 23 40  40100 59 23 45   837
PCSetUp                2 1.0 2.6970e-01 1.5 3.66e+07 3.6 0.0e+00 0.0e+00 5.0e+00 28 46  0  0  3  28 46  0  0  3   464
PCSetUpOnBlocks        1 1.0 2.6933e-01 1.6 3.66e+07 3.6 0.0e+00 0.0e+00 5.0e+00 28 46  0  0  3  28 46  0  0  3   464
PCApply               31 1.0 2.7084e-02 2.5 2.44e+07 2.4 0.0e+00 0.0e+00 0.0e+00  2 37  0  0  0   2 37  0  0  0  3684
VecMDot               30 1.0 9.2125e-03 2.3 7.67e+05 1.8 0.0e+00 0.0e+00 3.0e+01  1  1  0  0 18   1  1  0  0 20   403
VecNorm               36 1.0 1.6997e-02 2.9 5.94e+04 1.8 0.0e+00 0.0e+00 3.6e+01  1  0  0  0 21   1  0  0  0 24    17
VecScale              31 1.0 7.0810e-05 1.4 2.56e+04 1.8 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1748
VecCopy                7 1.0 2.4319e-05 2.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                39 1.0 4.9353e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                3 1.0 6.1750e-05 1.8 4.95e+03 1.8 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   388
VecMAXPY              31 1.0 3.5286e-04 1.9 8.17e+05 1.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0 11203
VecAssemblyBegin      18 1.0 7.1602e-02 2.5 0.00e+00 0.0 5.2e+01 1.2e+03 4.5e+01  7  0  4  1 27   7  0  4  1 30     0
VecAssemblyEnd        18 1.0 4.8876e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       38 1.0 7.7200e-04 1.9 0.00e+00 0.0 1.1e+03 1.7e+03 0.0e+00  0  0 74 38  0   0  0 74 38  0     0
VecScatterEnd         38 1.0 1.0075e-0137.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  7  0  0  0  0   7  0  0  0  0     0
VecNormalize          31 1.0 1.6758e-02 2.9 7.67e+04 1.8 0.0e+00 0.0e+00 3.1e+01  1  0  0  0 18   1  0  0  0 21    22
MatMult               31 1.0 1.0848e-01 9.4 8.64e+06 2.1 8.7e+02 1.3e+03 0.0e+00  8 14 59 23  0   8 14 59 23  0   353
MatSolve              31 1.0 2.6614e-02 2.6 2.44e+07 2.4 0.0e+00 0.0e+00 0.0e+00  2 37  0  0  0   2 37  0  0  0  3749
MatLUFactorNum         1 1.0 1.8399e-02 2.4 3.66e+07 3.6 0.0e+00 0.0e+00 0.0e+00  2 46  0  0  0   2 46  0  0  0  6794
MatILUFactorSym        1 1.0 2.5595e-01 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00 26  0  0  0  1  26  0  0  0  1     0
MatAssemblyBegin       2 1.0 2.3050e-03 4.2 0.00e+00 0.0 3.9e+01 6.9e+04 4.0e+00  0  0  3 55  2   0  0  3 55  3     0
MatAssemblyEnd         2 1.0 1.3324e-02 1.1 0.00e+00 0.0 5.6e+01 3.3e+02 8.0e+00  2  0  4  0  5   2  0  4  0  5     0
MatGetRowIJ            1 1.0 1.3113e-05 6.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 9.7036e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 4.0e+00  0  0  0  0  2   0  0  0  0  3     0
MatZeroEntries         3 1.0 1.4699e-03 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     3              3        19872     0
      Preconditioner     3              3         2048     0
                 Vec    46             46       425384     0
         Vec Scatter     5              5         4340     0
           Index Set    12             12        26956     0
   IS L to G Mapping     1              1         8048     0
              Matrix     4              4      6452564     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 5.29766e-05
Average time for zero size MPI_Send(): 9.28243e-05
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
| Time:           Fri Aug 24 15:28:20 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.962267, Active time=0.773428                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000109    0.0002      0.000163    0.01     0.02     |
|   build_constraint_matrix()        40        0.0019      0.000046    0.0019      0.000046    0.24     0.24     |
|   build_sparsity()                 1         0.0072      0.007184    0.0074      0.007445    0.93     0.96     |
|   cnstrn_elem_mat_vec()            20        0.0123      0.000613    0.0123      0.000613    1.58     1.58     |
|   constrain_elem_vector()          20        0.0001      0.000005    0.0001      0.000005    0.01     0.01     |
|   create_dof_constraints()         1         0.0093      0.009347    0.0632      0.063217    1.21     8.17     |
|   distribute_dofs()                1         0.0002      0.000237    0.0020      0.001994    0.03     0.26     |
|   dof_indices()                    428       0.0008      0.000002    0.0008      0.000002    0.10     0.10     |
|   enforce_constraints_exactly()    3         0.0028      0.000927    0.0028      0.000927    0.36     0.36     |
|   prepare_send_list()              1         0.0001      0.000119    0.0001      0.000119    0.02     0.02     |
|   reinit()                         1         0.0004      0.000430    0.0004      0.000430    0.06     0.06     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0001      0.000147    0.0047      0.004715    0.02     0.61     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0018      0.001809    0.0018      0.001809    0.23     0.23     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        798       0.0703      0.000088    0.0703      0.000088    9.09     9.09     |
|   init_shape_functions()           741       0.0072      0.000010    0.0072      0.000010    0.94     0.94     |
|   inverse_map()                    1620      0.0128      0.000008    0.0128      0.000008    1.66     1.66     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       1         0.0565      0.056457    0.0785      0.078453    7.30     10.14    |
|   assembly(get_residual)           1         0.0022      0.002175    0.0103      0.010278    0.28     1.33     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             798       0.0034      0.000004    0.0034      0.000004    0.44     0.44     |
|   compute_edge_map()               540       0.0004      0.000001    0.0004      0.000001    0.05     0.05     |
|   compute_face_map()               198       0.0008      0.000004    0.0008      0.000004    0.11     0.11     |
|   init_edge_shape_functions()      540       0.0004      0.000001    0.0004      0.000001    0.05     0.05     |
|   init_face_shape_functions()      152       0.0019      0.000012    0.0019      0.000012    0.24     0.24     |
|   init_reference_to_physical_map() 741       0.0268      0.000036    0.0268      0.000036    3.47     3.47     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0002      0.000200    0.0016      0.001566    0.03     0.20     |
|   renumber_nodes_and_elem()        2         0.0001      0.000039    0.0001      0.000039    0.01     0.01     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0005      0.000268    0.0005      0.000268    0.07     0.07     |
|   find_global_indices()            2         0.0001      0.000055    0.0032      0.001609    0.01     0.42     |
|   parallel_sort()                  2         0.0011      0.000555    0.0020      0.000986    0.14     0.26     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000021    0.0065      0.006545    0.00     0.85     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0003      0.000284    0.0003      0.000284    0.04     0.04     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0007      0.000747    0.0021      0.002149    0.10     0.28     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          1         0.0627      0.062712    0.4926      0.492598    8.11     63.69    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0014      0.000171    0.0014      0.000171    0.18     0.18     |
|   broadcast()                      1         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   gather()                         1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(scalar)                      3         0.1403      0.046751    0.1403      0.046751    18.13    18.13    |
|   max(vector)                      2         0.0002      0.000102    0.0002      0.000102    0.03     0.03     |
|   min(vector)                      2         0.0004      0.000190    0.0004      0.000190    0.05     0.05     |
|   probe()                          50        0.0014      0.000028    0.0014      0.000028    0.18     0.18     |
|   receive()                        50        0.0001      0.000002    0.0015      0.000030    0.01     0.19     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   send_receive()                   54        0.0001      0.000002    0.0017      0.000031    0.01     0.22     |
|   sum()                            14        0.0057      0.000409    0.0057      0.000409    0.74     0.74     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000108    0.0010      0.000982    0.01     0.13     |
|   set_parent_processor_ids()       1         0.0000      0.000014    0.0000      0.000014    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.3381      0.338095    0.3381      0.338095    43.71    43.71    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            6950      0.7734                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./vector_fe_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
