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
<pre> 
  
  #include &lt;iostream&gt;
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  
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
    
    Mesh mesh;
    
    mesh.read (argv[3]);
    
    mesh.find_neighbors();
  
    mesh.print_info();
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &gt;= 6 &amp;&amp; std::string(<B><FONT COLOR="#BC8F8F">&quot;-o&quot;</FONT></B>) == argv[4])
      mesh.write (argv[5]);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) introduction_ex1.C...
Linking introduction_ex1-opt...
***************************************************************
* Running Example  ./introduction_ex1-opt -d 3 /home/benkirk/codes/libmesh/reference_elements/3D/one_hex27.xda [-o output.xda]
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


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:58:56 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.081479, Active time=0.000461                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            2         0.0001      0.000038    0.0001      0.000038    16.49    16.49    |
|   renumber_nodes_and_elem()   2         0.0000      0.000002    0.0000      0.000002    0.65     0.65     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0000      0.000036    0.0000      0.000036    7.81     7.81     |
|   single_partition()          1         0.0000      0.000004    0.0000      0.000004    0.87     0.87     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0003      0.000342    0.0003      0.000342    74.19    74.19    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       7         0.0005                                          100.00            |
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


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:58:56 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.054565, Active time=0.000855                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Mesh                                                                                                      |
|   find_neighbors()            2         0.0001      0.000041    0.0001      0.000041    9.71     9.71     |
|   renumber_nodes_and_elem()   2         0.0000      0.000001    0.0000      0.000001    0.23     0.23     |
|                                                                                                           |
| Partitioner                                                                                               |
|   set_node_processor_ids()    1         0.0000      0.000039    0.0000      0.000039    4.56     4.56     |
|   single_partition()          1         0.0000      0.000004    0.0000      0.000004    0.47     0.47     |
|                                                                                                           |
| XdrIO                                                                                                     |
|   read()                      1         0.0004      0.000380    0.0004      0.000380    44.44    44.44    |
|   write()                     1         0.0003      0.000347    0.0003      0.000347    40.58    40.58    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       8         0.0009                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./introduction_ex1-opt -d 3 /home/benkirk/codes/libmesh/reference_elements/3D/one_hex27.xda [-o output.xda]
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
