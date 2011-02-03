<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 2 - Defining a Simple System</h1>

<br><br>This is the second example program.  It demonstrates how to
create an equation system for a simple scalar system.  This
example will also introduce some of the issues involved with using Petsc
in your application.

<br><br>This is the first example program that indirectly
uses the Petsc library.  By default equation data is stored
in Petsc vectors, which may span multiple processors.  Before
Petsc is used it must be initialized via libMesh::init().  Note that
by passing argc and argv to Petsc you may specify
command line arguments to Petsc.  For example, you might
try running this example as:

<br><br>./ex2 -log_info

<br><br>to see what Petsc is doing behind the scenes or

<br><br>./ex2 -log_summary

<br><br>to get a summary of what Petsc did.
Among other things, libMesh::init() initializes the MPI
communications library on your system if you haven't already
done so.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
</pre>
</div>
<div class = "comment">
Include file that defines various mesh generation utilities
</div>

<div class ="fragment">
<pre>
        #include "mesh_generation.h"
</pre>
</div>
<div class = "comment">
Include file that defines (possibly multiple) systems of equations.
</div>

<div class ="fragment">
<pre>
        #include "equation_systems.h"
</pre>
</div>
<div class = "comment">
Include files that define a simple steady system
</div>

<div class ="fragment">
<pre>
        #include "linear_implicit_system.h"
        #include "transient_system.h"
        #include "explicit_system.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        
        int main (int argc, char** argv)
        {
          LibMeshInit init (argc, argv);
        
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
A brief message to the user to inform her of the
exact name of the program being run, and its command line.
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Running " &lt;&lt; argv[0];
          for (int i=1; i&lt;argc; i++)
            std::cout &lt;&lt; " " &lt;&lt; argv[i];
          std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
          
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
2D grid on the unit square.  By default a mesh of QUAD4
elements will be created.  We instruct the mesh generator
to build a mesh of 5x5 elements.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_square (mesh, 5, 5);
        
</pre>
</div>
<div class = "comment">
Create an equation systems object. This object can
contain multiple systems of different 
flavors for solving loosely coupled physics.  Each system can 
contain multiple variables of different approximation orders.  
Here we will simply create a single system with one variable.  
Later on, other flavors of systems will be introduced.  For the 
moment, we use the general system.
The EquationSystems object needs a reference to the mesh
object, so the order of construction here is important.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          
</pre>
</div>
<div class = "comment">
Add a flag "test" that is visible for all systems.  This
helps in inter-system communication.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;bool&gt; ("test") = true;
            
</pre>
</div>
<div class = "comment">
Set a simulation-specific parameter visible for all systems.
This helps in inter-system-communication.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt; ("dummy") = 42.;
            
</pre>
</div>
<div class = "comment">
Set another simulation-specific parameter 
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt; ("nobody") = 0.;
          
</pre>
</div>
<div class = "comment">
Now we declare the system and its variables.
We begin by adding a "TransientLinearImplicitSystem" to the
EquationSystems object, and we give it the name
"Simple System".
</div>

<div class ="fragment">
<pre>
          equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; ("Simple System");
            
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Simple System".  "u"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Simple System").add_variable("u", FIRST);
        
</pre>
</div>
<div class = "comment">
Next we'll by add an "ExplicitSystem" to the
EquationSystems object, and we give it the name
"Complex System".
</div>

<div class ="fragment">
<pre>
          equation_systems.add_system&lt;ExplicitSystem&gt; ("Complex System");
        
</pre>
</div>
<div class = "comment">
Give "Complex System" three variables -- each with a different approximation
order.  Variables "c" and "T" will use first-order Lagrange approximation, 
while variable "dv" will use a second-order discontinuous
approximation space.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Complex System").add_variable("c", FIRST);
          equation_systems.get_system("Complex System").add_variable("T", FIRST);
          equation_systems.get_system("Complex System").add_variable("dv", SECOND, MONOMIAL);
            
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          equation_systems.init();
                
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
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Write the equation system if the user specified an
output file name.  Note that there are two possible
formats to write to.  Specifying libMeshEnums::WRITE will
create a formatted ASCII file.  Optionally, you can specify
libMeshEnums::ENCODE and get an XDR-encoded binary file.

<br><br>We will write the data, clear the object, and read the file
we just wrote.  This is simply to demonstrate capability.
Note that you might use this in an application to periodically
dump the state of your simulation.  You can then restart from
this data later.
</div>

<div class ="fragment">
<pre>
          if (argc &gt; 1)
            if (argv[1][0] != '-')
              {
                std::cout &lt;&lt; "&lt;&lt;&lt; Writing system to file " &lt;&lt; argv[1]
                          &lt;&lt; std::endl;
                
</pre>
</div>
<div class = "comment">
Write the system.
</div>

<div class ="fragment">
<pre>
                equation_systems.write (argv[1], libMeshEnums::WRITE);
                
</pre>
</div>
<div class = "comment">
Clear the equation systems data structure.
</div>

<div class ="fragment">
<pre>
                equation_systems.clear ();
        
                std::cout &lt;&lt; "&gt;&gt;&gt; Reading system from file " &lt;&lt; argv[1]
                          &lt;&lt; std::endl &lt;&lt; std::endl;
                
</pre>
</div>
<div class = "comment">
Read the file we just wrote.  This better
work!
</div>

<div class ="fragment">
<pre>
                equation_systems.read (argv[1], libMeshEnums::READ);
        
</pre>
</div>
<div class = "comment">
Print the information again.
</div>

<div class ="fragment">
<pre>
                equation_systems.print_info();
              }
          
</pre>
</div>
<div class = "comment">
All done.  libMesh objects are destroyed here.  Because the
LibMeshInit object was created first, its destruction occurs
last, and it's destructor finalizes any external libraries and
checks for leaked memory.
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
  
  #include &lt;iostream&gt;
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;explicit_system.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
    
    Mesh mesh;
    
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh, 5, 5);
  
    EquationSystems equation_systems (mesh);
    
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt; (<B><FONT COLOR="#BC8F8F">&quot;test&quot;</FONT></B>) = true;
      
    equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dummy&quot;</FONT></B>) = 42.;
      
    equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;nobody&quot;</FONT></B>) = 0.;
    
    equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Simple System&quot;</FONT></B>);
      
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Simple System&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
    equation_systems.add_system&lt;ExplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Complex System&quot;</FONT></B>);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Complex System&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;c&quot;</FONT></B>, FIRST);
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Complex System&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;T&quot;</FONT></B>, FIRST);
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Complex System&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;dv&quot;</FONT></B>, SECOND, MONOMIAL);
      
    equation_systems.init();
          
    mesh.print_info();
    equation_systems.print_info();
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &gt; 1)
      <B><FONT COLOR="#A020F0">if</FONT></B> (argv[1][0] != <B><FONT COLOR="#BC8F8F">'-'</FONT></B>)
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&lt;&lt;&lt; Writing system to file &quot;</FONT></B> &lt;&lt; argv[1]
                    &lt;&lt; std::endl;
          
          equation_systems.write (argv[1], libMeshEnums::WRITE);
          
          equation_systems.clear ();
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&gt;&gt;&gt; Reading system from file &quot;</FONT></B> &lt;&lt; argv[1]
                    &lt;&lt; std::endl &lt;&lt; std::endl;
          
          equation_systems.read (argv[1], libMeshEnums::READ);
  
          equation_systems.print_info();
        }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  mpirun -np 2 ./ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
    n_local_nodes()=22
  n_elem()=25
    n_local_elem()=13
    n_active_elem()=25
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=2
   System "Complex System"
    Type "Explicit"
    Variables="c" "T" "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=122
    n_constrained_dofs()=0
    n_vectors()=1
   System "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=22
    n_constrained_dofs()=0
    n_vectors()=3

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex2-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:11:23 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.876e-03      1.00000   3.876e-03
Objects:              3.300e+01      1.00000   3.300e+01
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         8.000e+00      1.00000   8.000e+00  1.600e+01
MPI Message Lengths:  3.700e+02      1.00000   4.625e+01  7.400e+02
MPI Reductions:       4.400e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.8500e-03  99.3%  0.0000e+00   0.0%  1.600e+01 100.0%  4.625e+01      100.0%  2.400e+01  54.5% 

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

VecSet                10 1.0 5.9605e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 5.0068e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    14             14        23288     0
         Vec Scatter     4              4         3472     0
           Index Set     8              8         4412     0
   IS L to G Mapping     4              4         2652     0
              Matrix     3              3         8996     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.19209e-06
Average time for zero size MPI_Send(): 5.00679e-06
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

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:11:23 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.010154, Active time=0.001482                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 2         0.0000      0.000021    0.0000      0.000025    2.83     3.37     |
|   compute_sparsity()           1         0.0001      0.000076    0.0001      0.000089    5.13     6.01     |
|   create_dof_constraints()     2         0.0000      0.000014    0.0000      0.000014    1.89     1.89     |
|   distribute_dofs()            2         0.0001      0.000055    0.0003      0.000141    7.42     19.03    |
|   dof_indices()                49        0.0000      0.000000    0.0000      0.000000    1.21     1.21     |
|   prepare_send_list()          2         0.0000      0.000003    0.0000      0.000003    0.47     0.47     |
|   reinit()                     2         0.0001      0.000059    0.0001      0.000059    8.03     8.03     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0001      0.000094    0.0001      0.000110    6.34     7.42     |
|   renumber_nodes_and_elem()    2         0.0000      0.000002    0.0000      0.000002    0.20     0.20     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    2         0.0002      0.000088    0.0002      0.000088    11.88    11.88    |
|   find_global_indices()        2         0.0001      0.000032    0.0006      0.000284    4.39     38.39    |
|   parallel_sort()              2         0.0002      0.000101    0.0002      0.000118    13.63    15.92    |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0001      0.000054    0.0001      0.000054    3.64     3.64     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0002      0.000215    0.0003      0.000345    14.51    23.28    |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  9         0.0000      0.000005    0.0000      0.000005    3.04     3.04     |
|   broadcast()                  1         0.0000      0.000005    0.0000      0.000005    0.34     0.34     |
|   gather()                     1         0.0000      0.000002    0.0000      0.000002    0.13     0.13     |
|   max(scalar)                  3         0.0000      0.000009    0.0000      0.000009    1.75     1.75     |
|   max(vector)                  2         0.0000      0.000002    0.0000      0.000002    0.27     0.27     |
|   min(vector)                  2         0.0000      0.000005    0.0000      0.000005    0.67     0.67     |
|   probe()                      14        0.0000      0.000002    0.0000      0.000002    1.48     1.48     |
|   receive()                    14        0.0000      0.000003    0.0001      0.000004    2.63     4.12     |
|   send()                       14        0.0000      0.000001    0.0000      0.000001    1.15     1.15     |
|   send_receive()               18        0.0000      0.000002    0.0001      0.000007    2.23     8.50     |
|   sum()                        6         0.0000      0.000006    0.0000      0.000006    2.23     2.23     |
|   wait()                       14        0.0000      0.000001    0.0000      0.000001    0.67     0.67     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0000      0.000025    0.0000      0.000044    1.69     2.97     |
|   set_parent_processor_ids()   1         0.0000      0.000002    0.0000      0.000002    0.13     0.13     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        171       0.0015                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
Running ./ex2-opt eqn_sys.dat

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
    n_local_nodes()=22
  n_elem()=25
    n_local_elem()=13
    n_active_elem()=25
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=2
   System "Complex System"
    Type "Explicit"
    Variables="c" "T" "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=122
    n_constrained_dofs()=0
    n_vectors()=1
   System "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=22
    n_constrained_dofs()=0
    n_vectors()=3

<<< Writing system to file eqn_sys.dat
>>> Reading system from file eqn_sys.dat

 EquationSystems
  n_systems()=2
   System "Complex System"
    Type "Explicit"
    Variables="c" "T" "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=122
    n_constrained_dofs()=0
    n_vectors()=1
   System "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=22
    n_constrained_dofs()=0
    n_vectors()=3


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:11:23 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.029722, Active time=0.020792                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 4         0.0001      0.000018    0.0001      0.000023    0.35     0.45     |
|   compute_sparsity()           2         0.0001      0.000058    0.0001      0.000069    0.56     0.67     |
|   create_dof_constraints()     4         0.0000      0.000010    0.0000      0.000010    0.20     0.20     |
|   distribute_dofs()            4         0.0002      0.000049    0.0005      0.000133    0.95     2.55     |
|   dof_indices()                98        0.0000      0.000000    0.0000      0.000000    0.19     0.19     |
|   prepare_send_list()          4         0.0000      0.000003    0.0000      0.000003    0.06     0.06     |
|   reinit()                     4         0.0002      0.000058    0.0002      0.000058    1.13     1.13     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   read()                       1         0.0063      0.006268    0.0086      0.008598    30.15    41.35    |
|   update()                     1         0.0002      0.000202    0.0002      0.000202    0.97     0.97     |
|   write()                      1         0.0100      0.009996    0.0101      0.010083    48.08    48.49    |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0001      0.000095    0.0001      0.000111    0.46     0.53     |
|   renumber_nodes_and_elem()    2         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   assign_global_indices()      2         0.0007      0.000370    0.0008      0.000416    3.55     4.01     |
|   compute_hilbert_indices()    2         0.0002      0.000088    0.0002      0.000088    0.85     0.85     |
|   find_global_indices()        2         0.0001      0.000036    0.0006      0.000284    0.35     2.73     |
|   parallel_sort()              2         0.0002      0.000103    0.0002      0.000117    1.00     1.13     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0001      0.000055    0.0001      0.000055    0.26     0.26     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0002      0.000213    0.0004      0.000351    1.02     1.69     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  31        0.0015      0.000050    0.0015      0.000050    7.42     7.42     |
|   broadcast()                  79        0.0001      0.000001    0.0001      0.000001    0.35     0.27     |
|   gather()                     17        0.0001      0.000003    0.0001      0.000003    0.28     0.28     |
|   max(scalar)                  5         0.0000      0.000008    0.0000      0.000008    0.20     0.20     |
|   max(vector)                  4         0.0000      0.000004    0.0000      0.000004    0.08     0.08     |
|   min(vector)                  4         0.0000      0.000005    0.0000      0.000005    0.09     0.09     |
|   probe()                      30        0.0000      0.000001    0.0000      0.000001    0.16     0.16     |
|   receive()                    62        0.0001      0.000001    0.0001      0.000002    0.29     0.47     |
|   send()                       46        0.0000      0.000001    0.0000      0.000001    0.18     0.18     |
|   send_receive()               42        0.0001      0.000002    0.0002      0.000005    0.34     0.95     |
|   sum()                        18        0.0001      0.000003    0.0001      0.000003    0.25     0.25     |
|   wait()                       30        0.0000      0.000000    0.0000      0.000000    0.06     0.06     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0000      0.000029    0.0000      0.000042    0.14     0.20     |
|   set_parent_processor_ids()   1         0.0000      0.000003    0.0000      0.000003    0.01     0.01     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        506       0.0208                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
