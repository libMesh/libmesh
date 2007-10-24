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
        
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize the library.  This is necessary because the library
may depend on a number of other libraries (i.e. MPI  and Petsc)
that require initialization before use.
</div>

<div class ="fragment">
<pre>
          libMesh::init (argc, argv);
</pre>
</div>
<div class = "comment">
Force all our objects to have local scope.  By declaring
libMesh objects in the next pair of braces we can assert
that they will go out of scope (and should have been deleted)
before we return from main.  This allows the library to do
internal reference counting and assure memory is not leaked.
</div>

<div class ="fragment">
<pre>
          {    
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
                error();
              }
            
</pre>
</div>
<div class = "comment">
Get the dimensionality of the mesh from argv[2]
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = atoi(argv[2]);
            
</pre>
</div>
<div class = "comment">
Create a mesh with the requested dimension.
</div>

<div class ="fragment">
<pre>
            Mesh mesh(dim);
            
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
            if (argc == 5)
              mesh.write (argv[4]);
</pre>
</div>
<div class = "comment">
At this closing brace all of our objects will be forced
out of scope and will get deconstructed.
</div>

<div class ="fragment">
<pre>
          }
</pre>
</div>
<div class = "comment">
All done.  Call the libMesh::close() function to close any
external libraries and check for leaked memory.  To be absolutey
certain this is called last we will return its value.  This
also allows main to return nonzero if memory is leaked, which
can be useful for testing purposes.
</div>

<div class ="fragment">
<pre>
          return libMesh::close();
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    {    
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 4)
        {
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -d 2 in.mesh [out.mesh]&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	error();
        }
      
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = atoi(argv[2]);
      
      Mesh mesh(dim);
      
      mesh.read (argv[3]);
      
      mesh.print_info();
      
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc == 5)
        mesh.write (argv[4]);
    }
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close();
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in debug mode) ex1.C...
Linking ex1...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'

***************************************************************
* Running Example  ./ex1
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=27
  n_elem()=1
   n_local_elem()=1
   n_active_elem()=1
  n_subdomains()=1
  n_processors()=1
  processor_id()=0


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 4Elem reference count information:
|  Creations:    7
|  Destructions: 7
| 4Node reference count information:
|  Creations:    27
|  Destructions: 27
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex1
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
