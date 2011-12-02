<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 1 - Creation of a Mesh Object</h1>

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
        #include "libmesh.h"
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        
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
./ex1 -d DIM input_mesh_name [output_mesh_name]
where [output_mesh_name] is an optional parameter giving
a filename to write the mesh into.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 4)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -d 2 in.mesh [out.mesh]"
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
Create a mesh
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
          
</pre>
</div>
<div class = "comment">
Read the input mesh.
</div>

<div class ="fragment">
<pre>
          mesh.read (argv[3]);
          
          mesh.find_neighbors();
        
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
<br><br><br> <h1> The program without comments: </h1> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex1-dbg -d 3 /home/dknez/software/libmesh/reference_elements/3D/one_hex27.xda [-o output.xda]
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
  n_processors()=1
  n_threads()=1
  processor_id()=0

WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
Option left: name:-d value: 3

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh4ElemE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh4NodeE reference count information:
|  Creations:    27
|  Destructions: 27
| N7libMesh9DofObjectE reference count information:
|  Creations:    29
|  Destructions: 29
 ---------------------------------------------------------------------------- 

-------------------------------------------------------------------
| Time:           Thu Dec  1 21:20:45 2011                         |
| OS:             Linux                                            |
| HostName:       dknez-laptop                                     |
| OS Release:     3.0.0-13-generic                                 |
| OS Version:     #22-Ubuntu SMP Wed Nov 2 13:27:26 UTC 2011       |
| Machine:        x86_64                                           |
| Username:       dknez                                            |
| Configuration:  ./configure run on Tue Nov 22 16:02:41 EST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.003803, Active time=0.000979                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            2         0.0001      0.000053    0.0001      0.000053    10.93    10.93    |
|   renumber_nodes_and_elem()   2         0.0000      0.000007    0.0000      0.000007    1.53     1.53     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0002      0.000151    0.0002      0.000151    15.42    15.42    |
|   single_partition()          1         0.0000      0.000021    0.0000      0.000021    2.15     2.15     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0007      0.000685    0.0007      0.000685    69.97    69.97    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       7         0.0010                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

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
  n_processors()=1
  n_threads()=1
  processor_id()=0

WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
Option left: name:-d value: 3
Option left: name:-o value: output.xda

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh4ElemE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh4NodeE reference count information:
|  Creations:    27
|  Destructions: 27
| N7libMesh9DofObjectE reference count information:
|  Creations:    29
|  Destructions: 29
 ---------------------------------------------------------------------------- 

-------------------------------------------------------------------
| Time:           Thu Dec  1 21:20:45 2011                         |
| OS:             Linux                                            |
| HostName:       dknez-laptop                                     |
| OS Release:     3.0.0-13-generic                                 |
| OS Version:     #22-Ubuntu SMP Wed Nov 2 13:27:26 UTC 2011       |
| Machine:        x86_64                                           |
| Username:       dknez                                            |
| Configuration:  ./configure run on Tue Nov 22 16:02:41 EST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.004027, Active time=0.001408                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            2         0.0001      0.000053    0.0001      0.000053    7.53     7.53     |
|   renumber_nodes_and_elem()   2         0.0000      0.000008    0.0000      0.000008    1.14     1.14     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0001      0.000149    0.0001      0.000149    10.58    10.58    |
|   single_partition()          1         0.0000      0.000021    0.0000      0.000021    1.49     1.49     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0007      0.000713    0.0007      0.000713    50.64    50.64    |
|   write()                     1         0.0004      0.000403    0.0004      0.000403    28.62    28.62    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       8         0.0014                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex1-dbg -d 3 /home/dknez/software/libmesh/reference_elements/3D/one_hex27.xda [-o output.xda]
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
