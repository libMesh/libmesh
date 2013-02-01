<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("introduction_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file introduction_ex2.C with comments: </h1> 
<div class = "comment">
<h1>Introduction Example 2 - Defining a Simple System</h1>

<br><br>This is the second example program.  It demonstrates how to
create an equation system for a simple scalar system.  This
example will also introduce some of the issues involved with using PETSc
in your application.

<br><br>This is the first example program that indirectly
uses the PETSc library.  By default equation data is stored
in PETSc vectors, which may span multiple processors.  Before
PETSc is used it must be initialized via libMesh::init().  Note that
by passing argc and argv to PETSc you may specify
command line arguments to PETSc.  For example, you might
try running this example as:

<br><br>./introduction_ex2 -log_info

<br><br>to see what PETSc is doing behind the scenes or

<br><br>./introduction_ex2 -log_summary

<br><br>to get a summary of what PETSc did.
Among other things, libMesh::init() initializes the MPI
communications library and PETSc numeric library on your system if
you haven't already done so.


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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
</pre>
</div>
<div class = "comment">
Include file that defines various mesh generation utilities
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh_generation.h"
</pre>
</div>
<div class = "comment">
Include file that defines (possibly multiple) systems of equations.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
</pre>
</div>
<div class = "comment">
Include files that define a simple steady system
</div>

<div class ="fragment">
<pre>
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/transient_system.h"
        #include "libmesh/explicit_system.h"
        
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
<br><br><br> <h1> The source file introduction_ex2.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/explicit_system.h&quot;</FONT></B>
  
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
* Running Example introduction_ex2:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex2/.libs/lt-example-devel -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
    n_local_nodes()=6
  n_elem()=25
    n_local_elem()=2
    n_active_elem()=25
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=2
   System #1, "Complex System"
    Type "Explicit"
    Variables={ "c" "T" } "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=24
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=0
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 0
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 0
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #0, "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=6
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 4.05556
      Average Off-Processor Bandwidth <= 4.22222
      Maximum  On-Processor Bandwidth <= 6
      Maximum Off-Processor Bandwidth <= 10
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:55:48 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.549e-01      1.00000   2.549e-01
Objects:              4.400e+01      1.00000   4.400e+01
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         1.090e+02      3.02778   6.650e+01  7.980e+02
MPI Message Lengths:  1.190e+03      3.52071   9.774e+00  7.800e+03
MPI Reductions:       5.900e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.5480e-01 100.0%  0.0000e+00   0.0%  7.980e+02 100.0%  9.774e+00      100.0%  5.800e+01  98.3% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecSet                10 1.0 1.3590e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         2 1.0 1.5020e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    16             16        25728     0
      Vector Scatter     6              6         6216     0
           Index Set    12             12         9512     0
   IS L to G Mapping     6              6         3384     0
              Matrix     3              3         8504     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 8.39233e-06
Average time for zero size MPI_Send(): 1.35899e-05
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:55:48 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.325193, Active time=0.039528                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 2         0.0008      0.000391    0.0056      0.002812    1.98     14.23    |
|   build_sparsity()             1         0.0008      0.000834    0.0036      0.003611    2.11     9.14     |
|   create_dof_constraints()     2         0.0003      0.000168    0.0003      0.000168    0.85     0.85     |
|   distribute_dofs()            2         0.0025      0.001271    0.0084      0.004208    6.43     21.29    |
|   dof_indices()                34        0.0049      0.000145    0.0049      0.000145    12.51    12.51    |
|   prepare_send_list()          2         0.0001      0.000044    0.0001      0.000044    0.22     0.22     |
|   reinit()                     2         0.0028      0.001400    0.0028      0.001400    7.09     7.09     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0014      0.001438    0.0016      0.001603    3.64     4.06     |
|   renumber_nodes_and_elem()    2         0.0001      0.000071    0.0001      0.000071    0.36     0.36     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    2         0.0005      0.000275    0.0005      0.000275    1.39     1.39     |
|   find_global_indices()        2         0.0008      0.000403    0.0065      0.003260    2.04     16.49    |
|   parallel_sort()              2         0.0024      0.001222    0.0036      0.001811    6.18     9.17     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0007      0.000727    0.0007      0.000727    1.84     1.84     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0058      0.005812    0.0078      0.007820    14.70    19.78    |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  12        0.0004      0.000034    0.0005      0.000038    1.04     1.16     |
|   max(bool)                    1         0.0001      0.000053    0.0001      0.000053    0.13     0.13     |
|   max(scalar)                  117       0.0015      0.000013    0.0015      0.000013    3.87     3.87     |
|   max(vector)                  27        0.0006      0.000021    0.0017      0.000063    1.46     4.27     |
|   min(bool)                    133       0.0017      0.000013    0.0017      0.000013    4.23     4.23     |
|   min(scalar)                  108       0.0035      0.000032    0.0035      0.000032    8.87     8.87     |
|   min(vector)                  27        0.0006      0.000023    0.0021      0.000078    1.59     5.35     |
|   probe()                      176       0.0015      0.000008    0.0015      0.000008    3.68     3.68     |
|   receive()                    176       0.0012      0.000007    0.0027      0.000016    3.09     6.90     |
|   send()                       176       0.0006      0.000004    0.0006      0.000004    1.60     1.60     |
|   send_receive()               180       0.0016      0.000009    0.0056      0.000031    4.12     14.18    |
|   sum()                        24        0.0008      0.000034    0.0013      0.000055    2.07     3.36     |
|                                                                                                            |
| Parallel::Request                                                                                          |
|   wait()                       176       0.0005      0.000003    0.0005      0.000003    1.15     1.15     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0006      0.000561    0.0012      0.001238    1.42     3.13     |
|   set_parent_processor_ids()   1         0.0001      0.000134    0.0001      0.000134    0.34     0.34     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        1391      0.0395                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex2:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
***************************************************************
* Running Example introduction_ex2:
*  mpirun -np 12 example-devel eqn_sys.dat -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex2/.libs/lt-example-devel eqn_sys.dat -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=36
    n_local_nodes()=6
  n_elem()=25
    n_local_elem()=2
    n_active_elem()=25
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=2
   System #1, "Complex System"
    Type "Explicit"
    Variables={ "c" "T" } "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=24
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=0
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 0
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 0
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #0, "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=6
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 4.05556
      Average Off-Processor Bandwidth <= 4.22222
      Maximum  On-Processor Bandwidth <= 6
      Maximum Off-Processor Bandwidth <= 10
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

<<< Writing system to file eqn_sys.dat
>>> Reading system from file eqn_sys.dat

 EquationSystems
  n_systems()=2
   System #0, "Complex System"
    Type "Explicit"
    Variables={ "c" "T" } "dv" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "MONOMIAL", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" "SECOND", "THIRD" 
    n_dofs()=222
    n_local_dofs()=24
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=0
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 0
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 0
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0
   System #1, "Simple System"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=36
    n_local_dofs()=6
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 4.05556
      Average Off-Processor Bandwidth <= 4.22222
      Maximum  On-Processor Bandwidth <= 6
      Maximum Off-Processor Bandwidth <= 10
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:55:49 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.067e-01      1.01965   1.063e-01
Objects:              7.700e+01      1.00000   7.700e+01
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         2.005e+02      3.03788   1.222e+02  1.467e+03
MPI Message Lengths:  3.172e+03      3.43290   1.456e+01  2.136e+04
MPI Reductions:       1.150e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.0630e-01 100.0%  0.0000e+00   0.0%  1.467e+03 100.0%  1.456e+01      100.0%  1.140e+02  99.1% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecCopy                2 1.0 1.0967e-05 5.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                20 1.0 2.1458e-05 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin       4 1.0 1.1959e-03 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  1  0  0  0 10   1  0  0  0 11     0
VecAssemblyEnd         4 1.0 2.8133e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        2 1.0 8.3208e-05 1.2 0.00e+00 0.0 1.4e+02 5.6e+01 0.0e+00  0  0  9 35  0   0  0  9 35  0     0
VecScatterEnd          2 1.0 3.4094e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 1.3828e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    30             30        48368     0
      Vector Scatter    10             10        10360     0
           Index Set    20             20        15952     0
   IS L to G Mapping    10             10         5640     0
              Matrix     6              6        17008     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 8.2016e-06
Average time for zero size MPI_Send(): 1.35104e-05
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:55:49 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.171392, Active time=0.092489                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 4         0.0014      0.000361    0.0112      0.002789    1.56     12.06    |
|   build_sparsity()             2         0.0013      0.000650    0.0040      0.001981    1.40     4.28     |
|   create_dof_constraints()     4         0.0005      0.000129    0.0005      0.000129    0.56     0.56     |
|   distribute_dofs()            4         0.0050      0.001253    0.0166      0.004146    5.42     17.93    |
|   dof_indices()                68        0.0099      0.000146    0.0099      0.000146    10.74    10.74    |
|   prepare_send_list()          4         0.0002      0.000041    0.0002      0.000041    0.18     0.18     |
|   reinit()                     4         0.0056      0.001406    0.0056      0.001406    6.08     6.08     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   read()                       1         0.0055      0.005490    0.0396      0.039576    5.94     42.79    |
|   update()                     1         0.0003      0.000344    0.0003      0.000344    0.37     0.37     |
|   write()                      1         0.0074      0.007398    0.0094      0.009433    8.00     10.20    |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0014      0.001431    0.0016      0.001572    1.55     1.70     |
|   renumber_nodes_and_elem()    2         0.0001      0.000069    0.0001      0.000069    0.15     0.15     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   assign_global_indices()      2         0.0059      0.002972    0.0094      0.004692    6.43     10.15    |
|   compute_hilbert_indices()    2         0.0006      0.000298    0.0006      0.000298    0.65     0.65     |
|   find_global_indices()        2         0.0008      0.000413    0.0056      0.002785    0.89     6.02     |
|   parallel_sort()              2         0.0023      0.001173    0.0026      0.001322    2.54     2.86     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0005      0.000467    0.0005      0.000467    0.50     0.50     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0058      0.005824    0.0078      0.007827    6.30     8.46     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  28        0.0006      0.000021    0.0007      0.000025    0.65     0.76     |
|   broadcast()                  70        0.0006      0.000009    0.0006      0.000008    0.66     0.61     |
|   max(bool)                    2         0.0000      0.000010    0.0000      0.000010    0.02     0.02     |
|   max(scalar)                  235       0.0027      0.000012    0.0027      0.000012    2.96     2.96     |
|   max(vector)                  53        0.0010      0.000019    0.0028      0.000053    1.10     3.05     |
|   min(bool)                    267       0.0030      0.000011    0.0030      0.000011    3.27     3.27     |
|   min(scalar)                  218       0.0090      0.000041    0.0090      0.000041    9.72     9.72     |
|   min(vector)                  53        0.0011      0.000022    0.0031      0.000059    1.24     3.38     |
|   probe()                      662       0.0078      0.000012    0.0078      0.000012    8.48     8.48     |
|   receive()                    522       0.0037      0.000007    0.0058      0.000011    4.00     6.28     |
|   send()                       434       0.0016      0.000004    0.0016      0.000004    1.73     1.73     |
|   send_receive()               386       0.0035      0.000009    0.0105      0.000027    3.81     11.40    |
|   sum()                        60        0.0011      0.000018    0.0021      0.000035    1.18     2.28     |
|                                                                                                            |
| Parallel::Request                                                                                          |
|   wait()                       434       0.0011      0.000003    0.0011      0.000003    1.20     1.20     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0005      0.000545    0.0012      0.001208    0.59     1.31     |
|   set_parent_processor_ids()   1         0.0001      0.000128    0.0001      0.000128    0.14     0.14     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        3532      0.0925                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex2:
*  mpirun -np 12 example-devel eqn_sys.dat -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
