<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("introduction_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file introduction_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Introduction Example 1 - Creation of a Mesh Object</h1>

<br><br>This is the first example program.  It simply demonstrates
how to create a mesh object.  A mesh is read from file,
information is printed to the screen, and the mesh is then
written.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
</pre>
</div>
<div class = "comment">
Functions to initialize the library.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh.h"
        
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
</pre>
</div>
<div class = "comment">
Initialize the library.  This is necessary because the library
may depend on a number of other libraries (i.e. MPI and PETSc)
that require initialization before use.  When the LibMeshInit
object goes out of scope, other libraries and resources are
finalized.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Check for proper usage. The program is designed to be run
as follows:
./ex1 -d DIM input_mesh_name [-o output_mesh_name]
where [output_mesh_name] is an optional parameter giving
a filename to write the mesh into.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 4)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -d 2 in.mesh [-o out.mesh]"
                          &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
This handy function will print the file name, line number,
and then abort.  Currently the library does not use C++
exception handling.
</div>

<div class ="fragment">
<pre>
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
Get the dimensionality of the mesh from argv[2]
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = std::atoi(argv[2]);
        
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
Create a mesh, with dimension to be overridden later, on the
default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
Read the input mesh.
</div>

<div class ="fragment">
<pre>
          mesh.read (argv[3]);
        
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
Write the output mesh if the user specified an
output file name.
</div>

<div class ="fragment">
<pre>
          if (argc &gt;= 6 && std::string("-o") == argv[4])
            mesh.write (argv[5]);
        
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
<br><br><br> <h1> The source file introduction_ex1.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 4)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -d 2 in.mesh [-o out.mesh]&quot;</FONT></B>
                    &lt;&lt; std::endl;
  
        libmesh_error();
      }
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = std::atoi(argv[2]);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
    Mesh mesh(init.comm());
  
    mesh.read (argv[3]);
  
    mesh.print_info();
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &gt;= 6 &amp;&amp; std::string(<B><FONT COLOR="#BC8F8F">&quot;-o&quot;</FONT></B>) == argv[4])
      mesh.write (argv[5]);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/introduction/introduction_ex1'
***************************************************************
* Running Example introduction_ex1:
*  mpirun -np 4 example-devel -d 3 /h2/roystgnr/libmesh/git/devel/../reference_elements/3D/one_hex27.xda -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=27
    n_local_nodes()=27
  n_elem()=1
    n_local_elem()=1
    n_active_elem()=1
  n_subdomains()=1
  n_partitions()=1
  n_processors()=4
  n_threads()=1
  processor_id()=0


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:44:12 2013                                                                          |
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
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.022938, Active time=0.002683                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            1         0.0001      0.000076    0.0001      0.000149    2.83     5.55     |
|   renumber_nodes_and_elem()   2         0.0000      0.000010    0.0000      0.000010    0.78     0.78     |
|                                                                                                           |
| Parallel                                                                                                  |
|   allgather()                 2         0.0001      0.000040    0.0001      0.000056    2.98     4.17     |
|   broadcast()                 22        0.0001      0.000003    0.0001      0.000003    2.61     2.09     |
|   max(scalar)                 37        0.0004      0.000011    0.0004      0.000011    15.21    15.21    |
|   max(vector)                 9         0.0001      0.000010    0.0004      0.000039    3.47     13.12    |
|   min(bool)                   45        0.0004      0.000009    0.0004      0.000009    14.98    14.98    |
|   min(scalar)                 36        0.0008      0.000023    0.0008      0.000023    31.27    31.27    |
|   min(vector)                 9         0.0001      0.000016    0.0004      0.000045    5.40     15.17    |
|   probe()                     6         0.0000      0.000006    0.0000      0.000006    1.45     1.45     |
|   receive()                   6         0.0000      0.000004    0.0001      0.000010    0.82     2.27     |
|   send()                      6         0.0000      0.000002    0.0000      0.000002    0.52     0.52     |
|   send_receive()              6         0.0000      0.000004    0.0001      0.000021    0.97     4.73     |
|                                                                                                           |
| Parallel::Request                                                                                         |
|   wait()                      6         0.0000      0.000004    0.0000      0.000004    0.97     0.97     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0001      0.000071    0.0003      0.000277    2.65     10.32    |
|   single_partition()          1         0.0000      0.000006    0.0000      0.000006    0.22     0.22     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0003      0.000345    0.0004      0.000373    12.86    13.90    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       196       0.0027                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex1:
*  mpirun -np 4 example-devel -d 3 /h2/roystgnr/libmesh/git/devel/../reference_elements/3D/one_hex27.xda -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
***************************************************************
* Running Example introduction_ex1:
*  mpirun -np 4 example-devel -d 3 /h2/roystgnr/libmesh/git/devel/../reference_elements/3D/one_hex27.xda -o output.xda -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=27
    n_local_nodes()=27
  n_elem()=1
    n_local_elem()=1
    n_active_elem()=1
  n_subdomains()=1
  n_partitions()=1
  n_processors()=4
  n_threads()=1
  processor_id()=0


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:44:13 2013                                                                          |
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
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.069081, Active time=0.002596                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            1         0.0001      0.000067    0.0001      0.000140    2.58     5.39     |
|   renumber_nodes_and_elem()   2         0.0000      0.000009    0.0000      0.000009    0.69     0.69     |
|                                                                                                           |
| Parallel                                                                                                  |
|   allgather()                 4         0.0001      0.000020    0.0001      0.000030    3.08     4.70     |
|   barrier()                   1         0.0000      0.000013    0.0000      0.000013    0.50     0.50     |
|   broadcast()                 22        0.0001      0.000003    0.0000      0.000002    2.27     1.81     |
|   gather()                    4         0.0001      0.000014    0.0001      0.000014    2.16     2.16     |
|   max(scalar)                 55        0.0003      0.000006    0.0003      0.000006    12.02    12.02    |
|   max(vector)                 13        0.0001      0.000006    0.0003      0.000019    3.20     9.63     |
|   min(bool)                   65        0.0002      0.000004    0.0002      0.000004    9.44     9.44     |
|   min(scalar)                 52        0.0005      0.000009    0.0005      0.000009    18.07    18.07    |
|   min(vector)                 13        0.0001      0.000009    0.0003      0.000025    4.39     12.48    |
|   probe()                     12        0.0001      0.000006    0.0001      0.000006    2.66     2.66     |
|   receive()                   18        0.0001      0.000003    0.0001      0.000007    1.96     4.70     |
|   send()                      6         0.0000      0.000003    0.0000      0.000003    0.58     0.58     |
|   send_receive()              6         0.0000      0.000004    0.0001      0.000021    0.92     4.97     |
|   sum()                       1         0.0000      0.000008    0.0000      0.000008    0.31     0.31     |
|                                                                                                           |
| Parallel::Request                                                                                         |
|   wait()                      12        0.0000      0.000003    0.0000      0.000003    1.31     1.31     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0001      0.000071    0.0003      0.000348    2.73     13.41    |
|   single_partition()          1         0.0000      0.000008    0.0000      0.000008    0.31     0.31     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0003      0.000318    0.0003      0.000344    12.25    13.25    |
|   write()                     1         0.0005      0.000482    0.0009      0.000917    18.57    35.32    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       291       0.0026                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex1:
*  mpirun -np 4 example-devel -d 3 /h2/roystgnr/libmesh/git/devel/../reference_elements/3D/one_hex27.xda -o output.xda -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/introduction/introduction_ex1'
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
